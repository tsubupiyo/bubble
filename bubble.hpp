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

NP_MAKE_NAMED_PARAMETER(theta);//[0:pi]
NP_MAKE_NAMED_PARAMETER(phi);  //[0:2pi)
NP_MAKE_NAMED_PARAMETER(u);//amplitude
NP_MAKE_NAMED_PARAMETER(k);//index of element of P
NP_MAKE_NAMED_PARAMETER(d);//distance between two points

constexpr int N_grid_points= 100;
#include "grid.hpp"

bool operator<(k_<int> a,k_<int> b)
{
   return a.value()<b.value();
}

std::vector<std::tuple<theta_<double>,phi_<double> > > grid_points = generate_random_theta_phi();
std::vector<std::array<size_t, 6> > network = generate_network(grid_points);

class Quadratic_function
{  //y = a*x*x + b*x + c
   public:
      std::tuple<double,double,double> beta;//a,b,c
      Quadratic_function(double a_, double b_, double c_);
      void add(std::tuple<u_<double>,d_<double> > pnt);
      std::tuple<double,double,double> get_parameter()const; 
   private:
      void LM();
      std::vector<double> ys;
      std::vector<double> xs;
      bool f_useable;
};

class Voronoi_cell
{
   public:
   k_<int> k;
   std::vector<Vector3D> const * P; //pointer of P (ps)
   private:
   std::vector< u_<double> > u; //amplitude(theta,phi)
   std::set<k_<int> >        K; //indexies of neighboring subset P in X
   public: 
   Voronoi_cell(const k_<int> i, std::vector<Vector3D> const * ps);
   double get_volume()const;
   std::set<k_<int> > get_neighbor()const; //return indexes of neighboring voronoi cells (K)
   void change_pointer(std::vector<Vector3D> const * ps);
   private:
   void boundary_fitting();//to determine the cell, fit u to the boundary.
};

class Voronoi_diagram
{
   private:
   std::vector<Vector3D> const * P; //point to P in X
   std::vector<Voronoi_cell> R; //voronoi cells
   public:
   Voronoi_diagram(std::vector<Vector3D>* const ps);
   void generate(); //for all voronoi cells
   std::vector<Voronoi_cell> get_vertual_cells(k_<int> k, const Vector3D& p)const;//for MC
   void change_pointer(std::vector<Vector3D> const * ps);
};

Quadratic_function::Quadratic_function(double a_, double b_, double c_)
{
   beta={a_,b_,c_};
   f_useable=false;
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

std::set<k_<int> > Voronoi_cell::get_neighbor()const
{//return indexes of neighboring voronoi cells (K)
   return K;
}

void Voronoi_diagram::generate()
{
   R.reserve(P->size());
   for(size_t s=0, size=P->size(); s<size; ++s)
   {
      R.push_back(Voronoi_cell(k_<int>(static_cast<int>(s)), P));
   } 
   return ;
}

Voronoi_cell::Voronoi_cell(const k_<int> i, std::vector<Vector3D> const * ps)
{
   k=i;
   change_pointer(ps);
}

void Voronoi_cell::change_pointer(std::vector<Vector3D> const * ps)
{
   P=ps;
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

