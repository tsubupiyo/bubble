#include "bubble.hpp"

int main()
{
   ////preview network
   //std::cout<<"# Bead 0"<<std::endl;
   //for(size_t i=0,size=grid_points.size();i<size;++i)
   //{
   //   std::cout<<i<<" "<<S2R(grid_points.at(i))<<std::endl;
   //}
   //std::cout<<"# Bond 0"<<std::endl;
   //for(size_t i=0,size=grid_points.size();i<size;++i)
   //{
   //   const auto& nbrs = network.at(i);   
   //   for(auto it=nbrs.begin();it!=nbrs.end();++it)
   //   {
   //      std::cout<<i<<" "<<*it<<std::endl;
   //   }
   //}
   //exit(0);

   constexpr int N=5;
   std::vector<Vector3D> cubic;
   {
      for(int x=0;x<N;++x)
      for(int y=0;y<N;++y)
      for(int z=0;z<N;++z)
      {
         cubic.push_back(Vector3D(x,y,z));   
         //if(x==3 && y==3 && z==3)
         //{
         //   std::cout<<"size:"<<cubic.size()<<std::endl;
         //   exit(0);
         //}
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

   //for(size_t i=0,size=cubic.size();i<size;++i)
   //{
      Voronoi_cell vc(k=1,&cubic);
      std::cout<<vc.get_volume()<<std::endl;
   //}
   
   return EXIT_SUCCESS;   
}
