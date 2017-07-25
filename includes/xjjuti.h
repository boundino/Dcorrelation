#ifndef _XJJUTI_H_
#define _XJJUTI_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <utility>
#include <cstring>
#include <cmath>

namespace xjjuti
{
  template<size_t N, typename T>
  int initarray(T (*array_)[N], T initval_=0);
  template<size_t N>
  int findibin(float (*array_)[N], float element_);
  template<size_t N>
  int findiedge(const float (*array_)[N], float element_);
  template<typename T>
  std::string number_to_string(T param_);
  template<typename T>
  std::string number_remove_zero(T param_);
  void progressbar(int index_, int total_, int morespace_=0);
  void progressbar_summary(int total_);
}

//
template<size_t N, typename T>
int xjjuti::initarray(T (*array_)[N], T initval_/*=0*/)
{
  for(std::size_t i=0;i<(sizeof(*array_)/sizeof((*array_)[0]));i++)
    {
      (*array_)[i]  = initval_;
    }
  return 0;
}

template<size_t N>
int xjjuti::findibin(const float (*array_)[N], float element_)
{
  int ielement = -1;
  for(std::size_t i=0;i<((sizeof(*array_)/sizeof((*array_)[0]))-1) && ielement<0;i++)
    {
      if(element_>=(*array_)[i] && element_<(*array_)[i+1]) ielement = i;
    }
  return ielement;
}

template<size_t N>
int xjjuti::findiedge(const float (*array_)[N], float element_)
{
  int ielement = -1;
  for(std::size_t i=0;i<(sizeof(*array_)/sizeof((*array_)[0])) && ielement<0;i++)
    {
      if(element_==(*array_)[i]) ielement = i;
    }
  return ielement;
}

template<typename T>
std::string xjjuti::number_to_string(T param_)
{
  if(param_<0) return "";
  std::string str = std::to_string(param_);
  std::size_t found = str.find('.');
  if(found==std::string::npos)
    {
      str+="p0";
      return str;
    }
  str.replace(found,1,"p");
  while(*(str.end()-1)=='0' && *(str.end()-2)!='p' && !str.empty()) str.erase(str.end()-1);
  return str;
}

template<typename T>
std::string xjjuti::number_remove_zero(T param_)
{
  if(param_<0) return "";
  std::string str = std::to_string(param_);
  std::size_t found = str.find('.');
  if(found==std::string::npos) return str;
  while(*(str.end()-1)=='0') str.erase(str.end()-1);
  if(*(str.end()-1)=='.') str.erase(str.end()-1);
  return str;
}

void xjjuti::progressbar(int index_, int total_, int morespace_/*=0*/)
{
  std::cout<<std::setiosflags(std::ios::left)<<"  [ \033[1;36m"<<std::setw(10+morespace_)<<index_<<"\033[0m"<<" / "<<std::setw(10+morespace_)<<total_<<" ] "<<"\033[1;36m"<<std::fixed<<std::setprecision(0)<<100.*index_/total_<<"%\033[0m"<<"\r"<<std::flush;
}

void xjjuti::progressbar_summary(int total_)
{
  std::cout<<std::endl<<"  Processed "<<"\033[1;31m"<<total_<<"\033[0m event(s)."<<std::endl;
}


#endif
