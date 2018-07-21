#pragma once

#include <iostream>
#include <cmath>
#include <vector>
#include <set>
#include "NamedParameter.hpp"
#include "Vector3D.hpp"

NP_MAKE_NAMED_PARAMETER(theta);//[-0.5pi:0.5pi]
NP_MAKE_NAMED_PARAMETER(phi);  //[0:2pi)
NP_MAKE_NAMED_PARAMETER(u);//amplitude
NP_MAKE_NAMED_PARAMETER(k);//index of element of P
NP_MAKE_NAMED_PARAMETER(d);//distance between two points
constexpr int N_grid_point = 100;
std::vector<std::tuple<theta_<double>,phi_<double> > > grid_points;

class Quadratic_function
{  //y = a*x*x + b*x + c
   public:
      double a;
      double b;
      double c;
   Quadratic_function(double a_, double b_, double c_);
   Quadratic_function(const std::vector<std::tuple<u_<double>,d_<double> > >& points);
   d_<double> get(std::tuple<u_<double>,d_<double> > point)const;
};

class Voronoi_cell
{
   public:
   k_<int> k;
   std::vector<Vector3D>* const P; //pointer of P (ps)
   private:
   std::vector< u_<double> > u; //amplitude(theta,phi)
   std::set<k_<int> >        K; //indexies of neighboring subset P in X
   public: 
   Voronoi_cell(const k_<int> i, const std::vector<Vector3D>* const ps);
   double get_volume()const;
   std::set<k_<int> > get_neighbor()const; //return indexes of neighboring voronoi cells (K)
   private:
   std::tuple<u_<double>,k_<int> > get_u_min()const; //return closest distance between Pk and Pother, and index of Pother
   void boundary_fitting();//fitting for u
   
};

class Voronoi_diagram
{
   private:
   std::vector<Vector3D>* const P; //point to P in X
   std::vector<Voronoi_cell> R; //voronoi cells
   public:
   Voronoi_diagram(std::vector<Vector3D>* const ps);
   void generate(); //for all voronoi cells
   std::vector<Voronoi_cell> get_vertual_cells(k_<int> k, const Vector3D& p)const;//for MC
};

