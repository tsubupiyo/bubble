#ifndef ReadFile_H
#define ReadFile_H

#include <iostream>
#include <cmath>
#include <vector>
#include <set>
#include <thread>
#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/logic/tribool.hpp>
#include <thread>
#include <fstream>

class Getline
{
   public:
      Getline(const std::string& filename)
      {
         ifs.open(filename,std::ios::in);
         back_f=false;
      }
      Getline(){back_f=false;}
      void set(const std::string& filename)
      {
         ifs.open(filename,std::ios::in);
      }
      std::string get()
      {
         if(back_f)
         {
            back_f=false;
            return prev;
         }
         if(!std::getline(ifs,prev))
         {
            ifs.close();
         }
         return prev;
      }
      bool is_open()const
      {
         return ifs.is_open();
      }
      void back()
      {
         if(back_f)
         {
            std::cout<<"dual back, fail"<<std::endl;
            exit(1);
         }
         back_f=true;
      }

   private:
      std::ifstream ifs;
      std::string prev;
      bool back_f;
   public:
      ~Getline()
      {
         if(ifs.is_open())
         {
            ifs.close();
         }
      }
};

//int main()
//{
// Getline test("sample.analysis");
// std::cout<<test.get()<<std::endl;
// test.back();
// test.back();
// std::cout<<test.get()<<std::endl;
// //while(1){
// // const auto hoge=test.get();
// // if(test.is_open())
// // {
// //    std::cout<<hoge<<std::endl;
// // }
// // else
// // {
// //    break;
// // }
// //}
//
// return EXIT_SUCCESS;
//}
#endif
