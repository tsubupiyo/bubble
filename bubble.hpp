#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <set>
#include <map>
#include <array>
#include "NamedParameter.hpp"
#include "Vector3D.hpp"

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

std::vector<std::tuple<theta_<double>,phi_<double> > > grid_points;
std::map<k_<int>,std::array<k_<int>, 6> > network;

class Quadratic_function
{  //y = a*x*x + b*x + c
   public:
      double a;
      double b;
      double c;
   Quadratic_function(double a_, double b_, double c_);
   d_<double> get(std::tuple<u_<double>,d_<double> > point)const;
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


std::set<k_<int> > Voronoi_cell::get_neighbor()const
{//return indexes of neighboring voronoi cells (K)
   return K;
}

Quadratic_function::Quadratic_function(double a_, double b_, double c_):a(a_),b(b_),c(c_){}

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

