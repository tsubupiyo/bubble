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
#include "ReadFile.hpp"
#include <stack>

NP_MAKE_NAMED_PARAMETER(theta);//[0:pi]
NP_MAKE_NAMED_PARAMETER(phi);  //[0:2pi)
NP_MAKE_NAMED_PARAMETER(u);//amplitude
NP_MAKE_NAMED_PARAMETER(k);//index of element of P
NP_MAKE_NAMED_PARAMETER(d);//distance between two points

typedef std::tuple<double,double,double> Beta;

constexpr size_t N_GRID_POINTS    = 50;
#include "grid.hpp"

bool operator<(const k_<size_t>& a, const k_<size_t>& b)
{
   return a.value()<b.value();
}

std::vector<std::tuple<theta_<double>,phi_<double> > > grid_points = generate_random_theta_phi();

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
   u.resize(N_GRID_POINTS);
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
   #warning //書いてね
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
