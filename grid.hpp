#pragma once
#include <random>

std::mt19937_64 mt_grid(265278465287);

#define PATH_TO_GRID_POINTS_SOURCE "GRID_POINTS_SOURCE"
std::vector<std::tuple<theta_<double>,phi_<double> > > generate_random_theta_phi();
std::vector<std::array<size_t, 6> > generate_network(const std::vector<std::tuple<theta_<double>,phi_<double> > >& ps);

inline Vector3D S2R (const std::tuple<theta_<double>,phi_<double> >& p)
{
   return Vector3D
      (
         std::sin(std::get<theta_<double>>(p).value())*std::cos(std::get<phi_<double> >(p).value()),
         std::sin(std::get<theta_<double>>(p).value())*std::sin(std::get<phi_<double> >(p).value()),
         std::cos(std::get<theta_<double>>(p).value())
      ); 
};

std::vector<std::tuple<theta_<double>,phi_<double> > > generate_random_theta_phi()
{
   {//load it
      std::vector< std::tuple<theta_<double>,phi_<double> > > loaded_random_points;
      try
      {
         Getline gl(PATH_TO_GRID_POINTS_SOURCE);
         std::string tmp;
         while(gl.is_open())
         {
            tmp=gl.get();
            std::vector<std::string> vs;
            boost::algorithm::split(vs,tmp,boost::is_any_of(" "));
            loaded_random_points.push_back
               (
                  {
                     theta_<double>(boost::lexical_cast<double> (vs.at(0))),
                     phi_<double>(  boost::lexical_cast<double> (vs.at(1)))
                  }
               );
         }
      }catch(...)
      {
      }
      if(loaded_random_points.size()==N_GRID_POINTS)
      {
         return loaded_random_points;
      }
   }

   std::vector< std::tuple<theta_<double>,phi_<double> > > random_points(N_GRID_POINTS);
   std::uniform_real_distribution<double> dist_theta(0.0,  M_PI);
   std::uniform_real_distribution<double> dist_phi(  0.0,2*M_PI);
   for(size_t i=0;i<N_GRID_POINTS;++i)
   {
      std::get<theta_<double> >(random_points.at(i)).value() = dist_theta(mt_grid);
      std::get<phi_<double>   >(random_points.at(i)).value() = dist_phi(mt_grid);
   }

   const auto LJ = [](const double& distance)->double
   {
      constexpr double EPS = 1.0;
      constexpr double sigma = cexpr_math::sqrt(1.0/N_GRID_POINTS);
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

   constexpr int N_STEP    = 10000*N_GRID_POINTS;
   constexpr int DOWN_STEP = N_STEP/10;
   double kT = 1000.0;
   double kT_=kT/10;
   double E_current = E_system(false);
   constexpr double DELTA_THETA =   M_PI/100;
   constexpr double DELTA_PHI   = 2*M_PI/100;
   std::uniform_int_distribution<int>     dist_pos(0,N_GRID_POINTS-1);
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
      //std::cout<<s<<" "<<E_current<<std::endl;
   }
   std::ofstream ofs(PATH_TO_GRID_POINTS_SOURCE,std::ios::trunc);
   boost::format fmt("%1.15e %1.15e\n");
   for(auto p : random_points)
   {
      ofs<<fmt %(std::get<theta_<double> >(p).value()) %(std::get<phi_<double> >(p).value());
   }
   
   return random_points; 
}

std::vector<std::array<size_t, 6> > generate_network
(
   const std::vector<std::tuple<theta_<double>,phi_<double> > >& ps
)
{
   std::vector<std::array<size_t, 6> > res(ps.size(),std::array<size_t,6>());
   for(size_t i=0,size=ps.size();i<size;++i)
   {
      std::vector<std::tuple<size_t,double> > rank;
      for(size_t nbr=0;nbr<size;++nbr)
      {
         if(i!=nbr)
         {
            rank.push_back(std::tuple<size_t,double>
               (
                  nbr,
                  (S2R(ps.at(i))-S2R(ps.at(nbr))).norm()
               )
            );
         }
      }
      std::sort(rank.begin(),rank.end(),[](const auto& a, const auto& b){return std::get<double>(a)<std::get<double>(b);});
      for(size_t s=0;s<6;++s)
      {
         res.at(i).at(s)=std::get<size_t>(rank.at(s));
      }
   }
   return res;
}
