#pragma once
#include <random>
//std::vector<std::tuple<theta_<double>,phi_<double> > > grid_points;

std::mt19937_64 mt_grid(265278465287);

std::vector<std::tuple<theta_<double>,phi_<double> > > generate_random_theta_phi()
{
   std::vector< std::tuple<theta_<double>,phi_<double> > > random_points(N_grid_points);
   std::uniform_real_distribution<double> dist_theta(-0.5*M_PI,+0.5*M_PI);
   std::uniform_real_distribution<double> dist_phi(0.0,2*M_PI);
   for(int i=0;i<N_grid_points;++i)
   {
      auto& [theta,phi] = random_points.at(i);
      theta = theta_<double>(dist_theta(mt_grid));
      phi = phi_<double>(dist_phi(mt_grid));
   }
   const auto S2R = [](std::tuple<theta_<double>,phi_<double> >& p)->Vector3D
   {
      return Vector3D
         (
            std::sin(std::get<theta_<double>>(p).value())*std::cos(std::get<phi_<double> >(p).value()),
            std::sin(std::get<theta_<double>>(p).value())*std::sin(std::get<phi_<double> >(p).value()),
            std::cos(std::get<theta_<double>>(p).value())
         ); 
   };
}

