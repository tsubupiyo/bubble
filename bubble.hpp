#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <set>
#include <map>
#include <array>
#include <boost/dynamic_bitset.hpp>
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

typedef std::tuple<double,double,double> Beta;

constexpr size_t N_GRID_POINTS= 50;
constexpr size_t N_SAMPLING_CURVE = 5;
constexpr double FIT_LIMIT_VALUE = 10*DBL_EPSILON;
constexpr double BETA_EPS = 10*DBL_EPSILON;
#include "grid.hpp"

bool operator<(const k_<size_t>& a, const k_<size_t>& b)
{
   return a.value()<b.value();
}

std::vector<std::tuple<theta_<double>,phi_<double> > > grid_points = generate_random_theta_phi();
std::vector<std::list<size_t>> network = get_network(grid_points);

class Quadratic_function
{  //y = a*x*x + b*x + c
   public:
      Beta beta;//a,b,c
      Quadratic_function();
      Quadratic_function(const double& a_, const double& b_, const double& c_);
      void set(const size_t& index_grid_point, const Vector3D& base, const Vector3D& neighbor, const Beta& ref_beta, bool sign_plus);
      //void add(const std::tuple<u_<double>,d_<double> >& pnt);
      const Beta& get_parameter()const; 
      double solve_coross_point()const;
   private:
      double limit;
      bool f_useable;
   public:
      double mismatch()const;
};

std::tuple<double,double> solve(const Beta& beta);

class Voronoi_cell
{
   public:
   k_<size_t> k;
   std::vector<Vector3D> const * P; //pointer of P (ps)
   private:
   std::vector< u_<double> >     u; //amplitude(theta,phi)
   std::set<k_<size_t> >         K; //indexies of neighboring subset P in X
   public: 
   Voronoi_cell(const k_<size_t>& i, std::vector<Vector3D> const * ps);
   double get_volume()const;
   std::set<k_<size_t> > get_neighbor()const; //return indexes of neighboring voronoi cells (K)
   void change_pointer(std::vector<Vector3D> const * ps);
   private:
   void boundary_fitting();//to determine the cell, fit u to the boundary.
   std::vector<Quadratic_function> qfs;
   std::vector<Beta> ref_beta;
};

class Voronoi_diagram
{
   private:
   std::vector<Vector3D> const * P; //point to P in X
   std::vector<Voronoi_cell>     R; //voronoi cells
   public:
   Voronoi_diagram(std::vector<Vector3D>* const ps);
   void generate(); //for all voronoi cells
   std::vector<Voronoi_cell> get_vertual_cells(k_<size_t> k, const Vector3D& p)const;//for MC
   void change_pointer(std::vector<Vector3D> const * ps);
   std::map<k_<size_t>,std::list<k_<size_t> > > get_Delaunay_diagram()const;
   std::list<k_<size_t> > get_partial_Delaunay_diagram(const k_<size_t>& center)const;
};

Quadratic_function::Quadratic_function()
{
   beta={0.0,0.0,0.0};
   f_useable=false;
   xs.resize(N_SAMPLING_CURVE);
   ys.resize(N_SAMPLING_CURVE);
   limit=DBL_MAX;
}

Quadratic_function::Quadratic_function(const double& a_, const double& b_, const double& c_)
{
   beta={a_,b_,c_};
   f_useable=false;
   xs.resize(N_SAMPLING_CURVE);
   ys.resize(N_SAMPLING_CURVE);
   limit=DBL_MAX;
}

void Quadratic_function::set
(
   const size_t&   index_grid_point,
   const Vector3D& base,
   const Vector3D& neighbor,
   const Beta&     ref_beta,
   bool sign_plus
)
{
   f_useable=false;
   beta=ref_beta;
   limit=DBL_MAX;
   const auto&   gp       = grid_points.at(index_grid_point);
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

   ////2. sampling distance with constant u (= unit length)
   double current_y = (base-neighbor).norm();
   std::vector<double> xs(N_SAMPLING_CURVE);
   std::vector<double> ys(N_SAMPLING_CURVE);

   const auto get_sample = [&](bool sign_plus)->void
   {
      auto base_ = base;
      int  x       = 0;
      int  counter = 0;
      bool   f_sampling = false;
      double current_y  = (base-neighbor).norm();
      do
      {
         if(sign_plus){++x; base += direction;}
         else         {--x; base -= direction;}

         const double y = (base-neighbor).norm();

         if(f_sampling)
         {
         //std::cout<<x<<" "<<y<<std::endl;
            xs.at(counter)=x;
            ys.at(counter)=y;   
            ++counter;
         }
         else
         {
            f_sampling=(sign_plus)?((current_y-y)<=0.0):((current_y-y)>=0.0);
         }
         current_y=y;
      }while(counter!=N_SAMPLING_CURVE);
   };

   get_sample(false);
   const Beta lower_side  = LevMar(xs,ys,ref_beta,FIT_LIMIT_VALUE);
   get_sample(true);
   const Beta higher_side = LevMar(xs,ys,ref_beta,FIT_LIMIT_VALUE);

   //lineality check
   if
   (//type-U
      (std::abs(std::get<0>(lower_side)-std::get<0>(higher_side))<BETA_EPS)&&
      (std::abs(std::get<1>(lower_side)-std::get<1>(higher_side))<BETA_EPS)&&
      (std::abs(std::get<2>(lower_side)-std::get<2>(higher_side))<BETA_EPS)
   )
   {
      beta=lower_side;  
      return ;
   }

   assert(std::abs(std::get<0>(lower_side ))<BETA_EPS);
   assert(std::abs(std::get<0>(higher_side))<BETA_EPS);
   assert(std::abs(1 + std::get<1>(lower_side ))<BETA_EPS);
   assert(std::abs(1 - std::get<1>(higher_side))<BETA_EPS);

   limit=(std::get<2>(higher_side)-std::get<2>(lower_side))/(std::get<1>(lower_side)-std::get<1>(higher_side)); 
   beta=lower_side;

   return ;
}

double Quadratic_function::solve()const
{
   const double& a = std::get<0>(beta);
   const double  b = std::get<1>(beta)-1; 
   const double& c = std::get<2>(beta);
   const double root = std::sqrt(b*b-4*a*c);
   const auto [res_min,res_max] = std::minmax(-(root+b)/(2*a),+(root-b)/(2*a));
   if(0.0<res_min && res_min<limit){return res_min;}
   if(0.0<res_max && res_max<limit){return res_max;}
   return DBL_MAX;
}

//void Quadratic_function::add(const std::tuple<u_<double>,d_<double> >& pnt)
//{
//   xs.push_back(std::get<u_<double>>(pnt).value());
//   ys.push_back(std::get<d_<double>>(pnt).value());
//   f_useable=false;
//}

const Beta& Quadratic_function::get_parameter()const
{
   assert(f_useable);
   return beta;
}

//void Quadratic_function::LM()
//{
//   if(xs.size()<5)return ;
//   beta = LevMar(xs,ys, beta); 
//   f_useable=true;
//}

double Quadratic_function::mismatch()const
{
   double sum_sqr_error=0.0;
   for(size_t i=0,size=xs.size();i<size;++i)
   {
      sum_sqr_error+=
      std::pow
      (
         ys.at(i)-(std::get<0>(beta)*sqr(xs.at(i))+std::get<1>(beta)*xs.at(i)+std::get<2>(beta)),2
      );
   }
   return std::sqrt(sum_sqr_error/(int)(xs.size()));
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

Voronoi_cell::Voronoi_cell(const k_<size_t>& i, std::vector<Vector3D> const * ps)
{
   k=i;
   change_pointer(ps);
   qfs.resize((*P).size());
   u.resize(N_GRID_POINTS);
   ref_beta.resize((*P).size());
   boundary_fitting();
   for(size_t i_=0,size=u.size();i_<size;++i_)
   {
      std::cout<<i_<<" u: "<<u.at(i_).value()<<std::endl;
   }
   exit(0);
}

void Voronoi_cell::change_pointer(std::vector<Vector3D> const * ps)
{
   P=ps;
}

double Voronoi_cell::get_volume()const
{
   constexpr size_t N = N_GRID_POINTS;
   constexpr double V0 = (4.0*M_PI/(3*N))*cexpr_math::sqrt(1.0-4.0/N);
   //here V0 is a volume of unit-cone (bottom area is Area_of_unit_sphere/N, length of cone is unit length(1)).
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
   const auto [min_distance,min_k] = [this]()->std::tuple<double,k_<size_t> >
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
            //std::cout<<pnt<<std::endl;
         }
      }
      return {min,min_k_};
   }();
   //std::cout<<"min_distance: "<<min_distance<<std::endl;
   //std::cout<<"min_k: "<<min_k.value()<<std::endl;
   //2. search nearest grid point
   const auto idx_grid_point_init = [min_k=min_k,this]()->size_t
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
   std::cout<<"idx_grid_point_init: "<<idx_grid_point_init<<std::endl;
   ////3. fitting using filling algorithm
   boost::dynamic_bitset<> filled((*P).size());
   constexpr u_<double> u0(DBL_MAX);
   for(size_t i=0;i<N_GRID_POINTS;++i)
   {
      u.at(i)=u0;
      filled[i]=false;
   }
   u.at(idx_grid_point_init).value()=0.5*min_distance;
   filled[idx_grid_point_init]=true;
   //std::cout<<"u init: "<<u.at(idx_grid_point_init).value()<<std::endl;
   const std::vector<Vector3D>& ps = *P;
   for(size_t i=0,size=ps.size();i<size;++i)
   {
      qfs.at(i).set(idx_grid_point_init,ps.at(k.value()),ps.at(i));
      ref_beta.at(i)=qfs.at(i).get_parameter();
   }
   std::stack<std::tuple<size_t,size_t> > stack;//2nd is ref-index of 1st(k)
   stack.push({idx_grid_point_init,idx_grid_point_init});
   K.clear();
   while(!stack.empty())
   {
      const auto [top,idx_ref_u] = stack.top(); stack.pop();
      filled[top]=true;
      //if(top!=idx_grid_point_init && DBL_MAX!=u.at(top).value()){continue;}//In this case, u(pos) is calculated-grid-point.
      for(size_t i=0,size=ps.size();i<size;++i)//estimation of distance function for each k
      {
         qfs.at(i).set(top,ps.at(k.value()),ps.at(i),ref_beta.at(i));
      }
      for(size_t i=0,size=ps.size();i<size;++i)//To stock & To lighten the next estimation
      {
         ref_beta.at(i)=qfs.at(i).get_parameter();
      }
      //find min(positive u)
      double min(DBL_MAX);//positive, and 0.5*min_distance<=
      size_t k_neighbor = std::numeric_limits<std::size_t>::max();
      for(size_t j=0,size=ps.size();j<size;++j)
      {
         if(k.value()==j){continue;}
         std::cout<<"j:"<<j<<" ";
         const auto [alpha,bravo] = solve(ref_beta.at(j));
         if((alpha>=0.0) && min>alpha)
         {
            min=alpha;
            k_neighbor=j;
         }
         if((bravo>=0.0) && min>bravo)
         {
            min=bravo;
            k_neighbor=j;
         }
      }
      u.at(top).value()=min;
      if(std::numeric_limits<std::size_t>::max()!=k_neighbor){K.insert(k_<size_t>(k_neighbor));};
      //If u is still DBL_MAX, the space in the direction is open.
      for(auto it=network.at(top).begin();
      it!=network.at(top).end();++it)
      {
         if(!filled[*it])
         stack.push({*it,top});
      }
   }
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

std::tuple<double,double> solve
(
   const Beta& beta
)
{
   const double& a = std::get<0>(beta);
   const double  b = std::get<1>(beta)-1; 
   const double& c = std::get<2>(beta);
   const double root = std::sqrt(b*b-4*a*c);
   std::cout<<"solv:"
      << -(root+b)/(2*a)<<" "
      << +(root-b)/(2*a)<<" "<<a<<" "<<b<<" "<<c<<std::endl;
   return 
   {
      -(root+b)/(2*a),
      +(root-b)/(2*a)
   };
}

std::list<k_<size_t> > Voronoi_diagram::get_partial_Delaunay_diagram(const k_<size_t>& center)const
{
   std::list<k_<size_t> > res;
   const Voronoi_cell& r = R.at(center.value());
   const std::set<k_<size_t> >& set = r.get_neighbor();
   for(auto& a: set)
   {
      res.push_back(a); 
   }
   return res;
}

std::map<k_<size_t>,std::list<k_<size_t> > > Voronoi_diagram::get_Delaunay_diagram()const
{
   std::map<k_<size_t>,std::list<k_<size_t> > > res;
   for(size_t r=0,size=R.size();r<size;++r)
   {
      const k_<size_t> k_(r);
      res.insert(std::make_pair(k_,get_partial_Delaunay_diagram(k_)));
   }
   return res;
}
