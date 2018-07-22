#pragma once
#include <random>
//std::vector<std::tuple<theta_<double>,phi_<double> > > grid_points;

std::mt19937_64 mt_grid(265278465287);

std::vector<std::tuple<theta_<double>,phi_<double> > > generate_random_theta_phi()
{
   std::vector< std::tuple<theta_<double>,phi_<double> > > random_points(N_grid_points);

   std::uniform_real_distribution<double> dist_theta(0.0,  M_PI);
   std::uniform_real_distribution<double> dist_phi(  0.0,2*M_PI);
   for(int i=0;i<N_grid_points;++i)
   {
      std::get<theta_<double> >(random_points.at(i)).value() = dist_theta(mt_grid);
      std::get<phi_<double>   >(random_points.at(i)).value() = dist_phi(mt_grid);
   }
   const auto S2R = [](const std::tuple<theta_<double>,phi_<double> >& p)->Vector3D
   {
      return Vector3D
         (
            std::sin(std::get<theta_<double>>(p).value())*std::cos(std::get<phi_<double> >(p).value()),
            std::sin(std::get<theta_<double>>(p).value())*std::sin(std::get<phi_<double> >(p).value()),
            std::cos(std::get<theta_<double>>(p).value())
         ); 
   };

   const auto LJ = [](const double& distance)->double
   {
      constexpr double EPS = 1.0;
      constexpr double sigma = cexpr_math::sqrt(1.0/N_grid_points);
      const double six = std::pow(sigma/distance,6);
      const double twl = std::pow(six,2);
      const auto res = 4.0*EPS*(twl-six);
      return res;
   };

   const auto E_system = [&]
   (
      bool is_virtual, size_t pos=0, 
      std::tuple<theta_<double>,phi_<double> > p = std::tuple<theta_<double>,phi_<double> >(0.0,0.0)
   )->double
   {
      double E=0.0; 
      for(size_t i=0,size=random_points.size();i<size;++i)
      for(size_t j=i+1;j<size;++j)
      {
         const std::tuple<theta_<double>,phi_<double> >& p_i = (!is_virtual)?(random_points.at(i)):((pos==i)?p:random_points.at(i)); 
         const std::tuple<theta_<double>,phi_<double> >& p_j = (!is_virtual)?(random_points.at(j)):((pos==j)?p:random_points.at(j)); 
         E+=LJ((S2R(p_i)-S2R(p_j)).norm());
      }
      return E;
   };

   constexpr int N_STEP    = 10000*N_grid_points;
   constexpr int DOWN_STEP = N_STEP/10;
   double kT = 1000.0;
   double kT_=kT/10;
   double E_current = E_system(false);
   constexpr double DELTA_THETA =   M_PI/100;
   constexpr double DELTA_PHI   = 2*M_PI/100;
   std::uniform_int_distribution<int>     dist_pos(0,N_grid_points-1);
   std::uniform_real_distribution<double> dist_delta_theta(-DELTA_THETA,+DELTA_THETA);
   std::uniform_real_distribution<double> dist_delta_phi(  -DELTA_PHI,  +DELTA_PHI);
   std::uniform_real_distribution<double> dist_p(0,1.0);
   int s_=0;
   for(int s=0;s<N_STEP;++s)
   {
      if((s_%DOWN_STEP)==0){s_=0;kT-=kT_;}
      const int pos = dist_pos(mt_grid);
      const auto stock = random_points.at(pos);
      const std::tuple<theta_<double>,phi_<double> > pnew
         (
            std::get<theta_<double> >(stock).value()+dist_delta_theta(mt_grid),
            std::get<phi_<double>   >(stock).value()+dist_delta_phi(mt_grid)
         );
      const double E_new = E_system(true,pos,pnew);
      const double DE    = E_new-E_current;
      if(DE<=0.0)
      {
         random_points.at(pos)=pnew;
         E_current=E_new;
      }
      else
      {
         if((std::abs(kT)>0.1)&&(std::exp(-DE/kT)<dist_p(mt_grid)))
         {
            random_points.at(pos)=pnew;
            E_current=E_new;
         }
      }
      std::cout<<s<<" "<<E_current<<std::endl;
   }
   return random_points; 
}
