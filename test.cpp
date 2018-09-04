#include "bubble.hpp"

int main()
{
   constexpr int N=5;
   std::vector<Vector3D> cubic;
   {
      for(int x=0;x<N;++x)
      for(int y=0;y<N;++y)
      for(int z=0;z<N;++z)
      {
         cubic.push_back(Vector3D(x,y,z));   
      }
   }
   //std::vector<Vector3D> hcp;
   //{
   //   
   //}
   //std::vector<Vector3D> random;
   //{
   //   std::mt19937 mt;
   //   std::uniform_real_distribution<double> dist(0,N);    
   //   for(int x=0;x<N;++x)
   //   for(int y=0;y<N;++y)
   //   for(int z=0;z<N;++z)
   //   random.push_back(Vector3D(dist(mt),dist(mt),dist(mt)));   
   //}
   Voronoi_cell vc(k=0,&cubic);
//   std::cout<<vc.get_volume()<<std::endl;
   
   return EXIT_SUCCESS;   
}
