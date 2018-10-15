#define _USE_MATH_DEFINES
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <omp.h>
#include <random>
#include <chrono>
#include "RandGen.h"

using namespace std;


RandGen::RandGen()
{ init(); }

void RandGen::init()
{ 
  srand(time(NULL));
  a=1.0f;
  b=10000.0f;
  c=10;
}


void RandGen::SetPar(float x,float y,int z)
{
  a=x;
  b=y;
  c=z;
}


void RandGen::SetPar(float x,float y)
{
  a=x;
  b=y;
}

void RandGen::SetPar(float x)
{ a=x; }

void RandGen::SetPar(int z)
{ c=z; }


float RandGen::Constant(float mean)
{ return mean; }


float RandGen::Uniform(float mean)
{ return  2.0f*mean*RVal(); }


float RandGen::Exponential(float mean)
{ return  -mean*(float)log(RVal()); }


int RandGen::Poisson(float mean) //Special technique required: Box-Muller method...
{
  float sum = 0.0f;
  int i;
  i=-1;
  float z;

  while(sum <=mean)
  {
    z = -(float)log(RVal());
    sum+= z;
    i++;
  }
  return i;
}
    
bool RandGen::PoiST(float mean,float n)  //Poisson Spike Train
{
  float p=mean/n;
  bool T=false;

  if( (float)rand() < (p*(float)RAND_MAX))
    T=true;

  return T;
}



float RandGen::Pareto(float alpha)
{ return 1.0f/(float)pow(RVal(),1.0f/alpha); }


float RandGen::Normal(float mean, float stdev)
{ return mean + stdev*(float)cos(2.0f*M_PI*RVal())*sqrt(-log(RVal())); }


int RandGen::Geometric(float p)
{ return (int)(log(RVal())/log(1.0f-p) - 1.0f); }


float RandGen::Weibull(float scale, float shape)
{ return scale*(float)pow(-(float)log(RVal()),1.0f/shape); }


float RandGen::Erlang(int scale, float shape)  //Special technique required...
{
  float R=1.0f;
  for (int i=1;i<=scale;i++)
    { R = R*RVal();}

  return -shape*(float)log(R);
}

//----------------------------------------------------------------------------

float RandGen::Constant()
{ return a; }


float RandGen::Uniform()
{ return  2.0f*a*RVal(); }


float RandGen::Exponential()
{ return  -a*(float)log(RVal()); }


int RandGen::Poisson() //Special technique required: Box-Muller method...
{
  float sum = 0.0f;
  int i;
  i=-1;
  float z;

  while(sum <=a)
  {
    z = -(float)log(RVal());
    sum+= z;
    i++;
  }
  return i;
}
    
bool RandGen::PoiST()  //Poisson Spike Train
{
  float p=a/b;
  bool T=false;

  if( (float)rand() < (p*(float)RAND_MAX))
    T=true;

  return T;

}



float RandGen::Pareto()
{ return 1.0f/(float)pow(RVal(),1.0f/a); }


float RandGen::Normal()
{ return a + b*(float)cos(2.0f*M_PI*RVal())*sqrt(-log(RVal())); }


int RandGen::Geometric()
{ return (int)(log(RVal())/log(1.0f-a) - 1.0f); }


float RandGen::Weibull()
{ return a*(float)pow(-log(RVal()),1.0f/b); }


float RandGen::Erlang()  //Special technique required...
{
  float R=1.0f;
  for (int i=1;i<=c;i++)
    { R = R*RVal(); }
    
  return -a*(float)log(R);
}

float RandGen::RVal() //Generate a random number between 0 and ~1.
{ return 1e-24f+(float)rand()/(float)RAND_MAX; }

void RandGen::Info()
{
  string s1;
  TexInfo(s1);
  cout<<s1;
}

void RandGen::Info(string& s1)
  { TexInfo(s1); }


void RandGen::TexInfo(string& ostr)
{
  ostringstream ostr_tmp;

  ostr_tmp
      <<"mean is "<<a
      <<" number is "<<b
      <<" scale is "<<c
      <<endl;

  ostr=ostr_tmp.str();
}

