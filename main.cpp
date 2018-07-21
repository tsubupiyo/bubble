#include <random>
#include "bubble.hpp"

constexpr std::uint_fast64_t SEED = static_cast<std::uint_fast64_t>(34609273465027436);
constexpr int N = 3;
constexpr double L=10;

int main()
{
   std::mt19937_64 mt(SEED); 
   const std::vector<Vector3D> ps=[&mt]()
   {
      std::uniform_real_distribution<double>  dist_d(0,L);   
      std::vector<Vector3D> result;
      for(int n=0; n<N; ++n)
      {
         result.push_back({dist_d(mt),dist_d(mt),dist_d(mt)});
      }
      return result;
   }();   

   for(int i=0; i<N; ++i)
   {
      std::cout<<ps.at(i)<<std::endl;
   }

   return EXIT_SUCCESS;
}
