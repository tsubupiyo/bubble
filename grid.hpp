#pragma once
#include <random>
//std::vector<std::tuple<theta_<double>,phi_<double> > > grid_points;

std::mt19937_64 mt_grid(265278465287);

std::tuple<u_<double>,theta_<double>,phi_<double> > 
RectangleCS2SphericalCS
(const Vector3D& v)
{
   const auto& [x,y,z] = v;
   const auto r = std::hypot(x,y,z);
   const auto sign = [](const auto& a)
   {
      return (a<0)?-1:1;
   };
   return 
   std::tuple<u_<double>,theta_<double>,phi_<double> >
      (u_<double>(r),std::acos(z/r),sign(y)*std::acos(x/std::hypot(x,y)));
}

Vector3D generate_random_point_on_sphere(double R)
{
   std::uniform_real_distribution<double> dist_phi(0,2*M_PI);
   std::uniform_real_distribution<double> dist_z(-R,+R);
   const double f = dist_phi(mt_grid);
   const double z = dist_z(mt_grid);
   const double s = std::sqrt(sqr(R)-sqr(z));
   return 
   {
      s*std::cos(f),
      s*std::sin(f),
      z
   };
}

std::vector<std::tuple<theta_<double>,phi_<double> > > generate_random_theta_phi()
{
   std::vector<std::tuple<theta_<double>,phi_<double> > > res;
   for(int i=0;i<N_grid_points;++i)
   {
      const auto [r, theta, phi] = generate_random_point_on_sphere(1.0);
      res.push_back(std::tuple<theta_<double>,phi_<double> >(theta,phi)); 
   }
   return res;
}
