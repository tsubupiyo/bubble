#include "bubble.hpp"

int main()
{
   std::vector<Vector3D> cubic_lattice;
   {
      constexpr int N=5;
      for(int x=0;x<N;++x)
      for(int y=0;y<N;++y)
      for(int z=0;z<N;++z)
      {
         cubic_lattice.push_back(Vector3D(x,y,z));   
      }
   }
   
   return EXIT_SUCCESS;   
}
