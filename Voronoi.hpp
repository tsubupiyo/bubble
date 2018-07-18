#pragma once

class Plane
{
   public:
      Vector3D n;
      Vector3D o; 
  //名前付しておけばベクトル２つとる関数２つかけるのでは？２点から平面を生成することもできるぞ  
};

template <typename SUBSET>//Pooint, Segment, etc norm()あればいい
class Voronoi_cell
{
   private:
   int k;
   SUBSET* const P;   // in the space X
};

template <typename SUBSET>
class Voronoi_diagram
{
   private:
      std::vector<SUBSET>* const Ps;
      std::vector<Voronoi_cell>  Rs;
      std::map<std::tuple<int,int>,Plane>
};

