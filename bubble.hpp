#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <set>
#include <map>
#include <array>
#include "Chaperone/constexpr_math.hpp"
#include "Chaperone/NamedParameter.hpp"
#include "Chaperone/Vector3D.hpp"
#include "Levenberg-Marquardt/LevMar.hpp"
#include "ReadFile.hpp"
#include <stack>

NP_MAKE_NAMED_PARAMETER(theta);//[0:pi]
NP_MAKE_NAMED_PARAMETER(phi);  //[0:2pi)
NP_MAKE_NAMED_PARAMETER(u);//amplitude
NP_MAKE_NAMED_PARAMETER(k);//index of element of P
NP_MAKE_NAMED_PARAMETER(d);//distance between two points

constexpr size_t N_grid_points= 100;
#include "grid.hpp"
constexpr size_t N_SAMPLING_CURVE = 10;
constexpr size_t N_NEIGHBOR = 6;

bool operator<(k_<size_t> a,k_<size_t> b)
{
   return a.value()<b.value();
}

std::vector<std::tuple<theta_<double>,phi_<double> > > grid_points = generate_random_theta_phi();
std::vector<std::array<size_t, N_NEIGHBOR> > network = generate_network(grid_points);

class Quadratic_function
{  //y = a*x*x + b*x + c
   public:
      std::tuple<double,double,double> beta;//a,b,c
      Quadratic_function();
      Quadratic_function(double a_, double b_, double c_);
      void set(size_t index_grid_point, Vector3D base, const Vector3D& neighbor);
      void add(std::tuple<u_<double>,d_<double> > pnt);
      std::tuple<double,double,double> get_parameter()const; 
   private:
      void LM();
      std::vector<double> ys;
      std::vector<double> xs;
      bool f_useable;
};

double solve(const std::tuple<double,double,double>& beta_1, const std::tuple<double,double,double>& beta_2, const double x_init=0.0);

class Voronoi_cell
{
   public:
   k_<size_t> k;
   std::vector<Vector3D> const * P; //pointer of P (ps)
   private:
   std::vector< u_<double> > u; //amplitude(theta,phi)
   std::set<k_<size_t> >        K; //indexies of neighboring subset P in X
   public: 
   Voronoi_cell(const k_<size_t> i, std::vector<Vector3D> const * ps);
   double get_volume()const;
   std::set<k_<size_t> > get_neighbor()const; //return indexes of neighboring voronoi cells (K)
   void change_pointer(std::vector<Vector3D> const * ps);
   private:
   void boundary_fitting();//to determine the cell, fit u to the boundary.
   void solve_u_stack(const size_t pos_origin);
   std::vector<Quadratic_function> qfs;
};

class Voronoi_diagram
{
   private:
   std::vector<Vector3D> const * P; //point to P in X
   std::vector<Voronoi_cell> R; //voronoi cells
   public:
   Voronoi_diagram(std::vector<Vector3D>* const ps);
   void generate(); //for all voronoi cells
   std::vector<Voronoi_cell> get_vertual_cells(k_<size_t> k, const Vector3D& p)const;//for MC
   void change_pointer(std::vector<Vector3D> const * ps);
};

Quadratic_function::Quadratic_function()
{
   beta={0.0,0.0,0.0};
   f_useable=false;
   xs.reserve(N_SAMPLING_CURVE);
   ys.reserve(N_SAMPLING_CURVE);
}

Quadratic_function::Quadratic_function(double a_, double b_, double c_)
{
   beta={a_,b_,c_};
   f_useable=false;
   xs.reserve(N_SAMPLING_CURVE);
   ys.reserve(N_SAMPLING_CURVE);
}

void Quadratic_function::set(size_t index_grid_point, Vector3D base, const Vector3D& neighbor)
{
   const auto& gp = grid_points.at(index_grid_point);
   const double& gp_theta = (std::get<theta_<double> >(gp)).value();
   const double& gp_phi   = (std::get<phi_<double>   >(gp)).value();
   //1. get direction vector
   const double sin_thteta = std::sin(gp_theta);
   const double sin_phi    = std::sin(gp_phi  );
   const double cos_thteta = std::cos(gp_theta);
   const double cos_phi    = std::cos(gp_phi  );
   const Vector3D direction
      (
         sin_thteta*cos_phi,
         sin_thteta*sin_phi,
         cos_thteta
      );
   //2. sampling distance with constant +u (= unit length)
   xs.resize(N_SAMPLING_CURVE);
   ys.resize(N_SAMPLING_CURVE);
   for(size_t s=0;s<N_SAMPLING_CURVE;++s)
   {
      base+=direction;
      xs.at(s)=static_cast<int>(s);
      ys.at(s)=(base-neighbor).norm();
   }
   //3. get beta
   LM();
}

void Quadratic_function::add(std::tuple<u_<double>,d_<double> > pnt)
{
   xs.push_back(std::get<u_<double>>(pnt).value());
   ys.push_back(std::get<d_<double>>(pnt).value());
   f_useable=false;
}

std::tuple<double,double,double> Quadratic_function::get_parameter()const
{
   assert(f_useable);
   return beta;
}

void Quadratic_function::LM()
{
   if(xs.size()<5)return ;
   beta = LevMar(xs,ys, beta); 
   f_useable=true;
}

std::set<k_<size_t> > Voronoi_cell::get_neighbor()const
{//return indexes of neighboring voronoi cells (K)
   return K;
}

void Voronoi_diagram::generate()
{
   R.reserve(P->size());
   for(size_t s=0, size=P->size(); s<size; ++s)
   {
      R.push_back(Voronoi_cell(k_<size_t>(s), P));
   } 
   return ;
}

Voronoi_cell::Voronoi_cell(const k_<size_t> i, std::vector<Vector3D> const * ps)
{
   k=i;
   change_pointer(ps);
   qfs.resize(N_grid_points);
}

void Voronoi_cell::change_pointer(std::vector<Vector3D> const * ps)
{
   P=ps;
}

double Voronoi_cell::get_volume()const
{
   constexpr size_t N = N_grid_points;
   constexpr double V0 = (4.0*M_PI/(3*N))*cexpr_math::sqrt(1.0-4.0/N);
   double sum=0.0;
   for(size_t i=0,size=u.size();i<size;++i)
   {
      sum+=u.at(i).value();
   }
   return sum*V0;
}

void Voronoi_cell::boundary_fitting()
{
   //1. find nearest neighbor
   const auto [min_dist,min_k] = [this]()->std::tuple<double,k_<size_t> >
   {
      const auto& ps = *P;
      double min = DBL_MAX;
      k_<size_t> min_k_ = k;
      const Vector3D& center = ps.at(k.value());
      for(size_t i=0,size=ps.size();i<size;++i)
      {
         if(i!=k.value())
         {
            const Vector3D& pnt = ps.at(i);
            const double dis = (pnt-center).norm();
            if(dis<min){min=dis;min_k_.value()=i;}
         }
      }
      return {min,min_k_};
   }();
   //2. search nearest grid point
   const auto grid_point_init = [min_k=min_k,this]()->size_t
   {
      const auto& ps = *P;
      Vector3D direction = (ps.at(min_k.value())-ps.at(k.value()));
      direction.normalize();
      const auto tp2v = [](const std::tuple<theta_<double>,phi_<double> >& tp)->Vector3D
      {
         const auto theta = (std::get<0>(tp)).value();
         const auto phi   = (std::get<1>(tp)).value();
         return {std::sin(theta)*std::cos(phi),std::sin(theta)*std::sin(phi),std::cos(theta)};
      };
      const auto dist = [&direction,&tp2v]
      (
         const std::tuple<theta_<double>,phi_<double> >& a
      )->double
      {
         return ((tp2v(a)-direction).norm());
      };
      double min = DBL_MAX;
      size_t min_index = 0;
      for(size_t i=0,size=grid_points.size();i<size;++i)
      {
         const double distance = dist(grid_points.at(i));
         if(min>distance)
         {
            min=distance;
            min_index=i;
         }
      }
      return min_index;
   }();
   ////3. fitting using filling algorithm
   //void solve_u_stack(const size_t pos);
   //TODO::u の初期化を行ってそれの値にフラグの意味をもたせてsolve_u_stackで利用する

}

Voronoi_diagram::Voronoi_diagram(std::vector<Vector3D>* const ps)
{
   P=ps;
   generate();
}

void Voronoi_diagram::change_pointer(std::vector<Vector3D> const * ps)
{
   P=ps;
   std::for_each(R.begin(),R.end(),[this](auto& r){r.change_pointer(P);});
}

double solve
(
   const std::tuple<double,double,double>& beta_1, 
   const std::tuple<double,double,double>& beta_2, 
   const double x_init
)
{//Newton's method + damping
   const double A = std::get<0>(beta_1)-std::get<0>(beta_2);
   const double B = std::get<1>(beta_1)-std::get<1>(beta_2);
   const double C = std::get<2>(beta_1)-std::get<2>(beta_2);
   const double A2 = std::pow(A,2);
   const double B2 = std::pow(B,2);
   const double C2 = std::pow(C,2);
   const double AB = A*B;
   const double BC = B*C;
   const double AC = A*C;

   const auto f = [&](const double x1)->std::tuple<double,double>
   {
      const double x2 = x1*x1;
      const double x3 = x1*x2;
      const double x4 = x1*x3;
      return 
      {
         A2*x4 + 2*AB*x3 + 2*AC*x2 + B2*x2 + 2*BC*x1           + C2,  //f(x)
                 4*A2*x3 + 6*AB*x2         + 4*AC*x1 + 2*B2*x1 + 2*BC //f'(x)
      };
   };

   double x   = x_init;
   double xp1 = x_init;
   std::tuple<double,double> ff;
   constexpr double EPS_NEWTON      = 0.000001;
   constexpr double LMD_NEWTON      = 0.1;
   constexpr int NEWTON_COUNT_LIMIT = 100;
   int count=0;

   do
   {
      x   = xp1;
      ff  = f(x);
      xp1 = x - std::get<0>(ff) / (std::get<1>(ff) * (1 + LMD_NEWTON));
      if(++count > NEWTON_COUNT_LIMIT){break;}
   }while( std::abs( (xp1 + x) / x) > EPS_NEWTON);

   return x;
}

[[deprecated]]
void fill_it_stack
(
   std::vector<std::vector<std::vector<Color> > >&array,
   const int& x0, const int& y0, const int& z0
)
{
   const int maxx = array.size();
   const int maxy = array[0].size();
   const int maxz = array[0][0].size();

   typedef std::tuple<int,int,int> ref ;
   std::stack<ref> refs;
   refs.push(ref(x0,y0,z0));
   while(!refs.empty())
   {
      const ref& top = refs.top();
      const int& x = std::get<0>(top);
      const int& y = std::get<1>(top);
      const int& z = std::get<2>(top);
      Color& topc = array[std::get<0>(top)][std::get<1>(top)][std::get<2>(top)];
      if(topc==UNKNOWN)
      {
         topc=PAST;
         if(((x+1)<maxx)&&array[x+1][y][z]==UNKNOWN){ refs.push(ref(x+1,y,z)); }
         if(((y+1)<maxy)&&array[x][y+1][z]==UNKNOWN){ refs.push(ref(x,y+1,z)); }
         if(((z+1)<maxz)&&array[x][y][z+1]==UNKNOWN){ refs.push(ref(x,y,z+1)); }
         if(((x-1)>=0)  &&array[x-1][y][z]==UNKNOWN){ refs.push(ref(x-1,y,z)); }
         if(((y-1)>=0)  &&array[x][y-1][z]==UNKNOWN){ refs.push(ref(x,y-1,z)); }
         if(((z-1)>=0)  &&array[x][y][z-1]==UNKNOWN){ refs.push(ref(x,y,z-1)); }
      }
      refs.pop();
   }
}

void Voronoi_cell::solve_u_stack
(
   const size_t idx_grid_point_init,
   const k_<size_t>& idx_neighbor,
   const d_<double>& min_distance// ||P_k-P_neighbor||
)
{
   #warning //まだだよ
   const u_<double> u0(DBL_MAX);
   std::fill(u.begin(),u.end(),u0);
   const std::vector<Vector3D>& ps = *P;
   const Vector3D& base = ps.at(k.value());
   qfs.at(k.value()   ).set(idx_grid_point_init,base,k.value()   );
   qfs.at(idx_neighbor).set(idx_grid_point_init,base,idx_neighbor);
   u.at(idx_grid_point_init)=solve(qfs.at(k.value()).get_parameter(),qfs.at(idx_neighbor).get_parameter());
   std::stack<std::tuple<size_t,size_t> > stack;//次計算するgpと、そのとき参照するgp
   for(size_t i=0;i<N_NEIGHBOR;++i)
   {
      stack.push({network.at(idx_grid_point_init),idx_grid_point_init});
   }
   while(!stack.empty())
   {
      //1. pop
      const auto [pos,ref] = stack.pop();
      if(DBL_MAX!=u.at(pos).value()){continue;}
      //2. init qfs
      for(size_t i=0;i<N_grid_points;++i)
      {
         qfs.set(pos,base,ps.at(i));
      }
      //3. find min_arg(solve(i,j))
      double min(DBL_MAX);
      size_t min_idx;
      for(size_t i=0;i<N_grid_points;++i)
      {
        if(k.value()!=i) 
        {
           //solve
           u_res = solve(qfs.at(k.value()),qfs.at(i));
           if(u_res<min)
           {
              min     = u_res;
              min_idx = i;
           }
        }
      }
      u.at(pos)=min;
      //TODO::push with min_idx
   }
}
