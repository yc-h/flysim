#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <omp.h>
#include <random>
#include <chrono>
#include <ctime>
#include "StiExt.h"


using namespace std;
// global variable: random generator
random_device seed_gen;

StiExt::StiExt()
{
  rdgen.seed(seed_gen());

  normal_distribution<float>::param_type par(MNoiseMean,MNoiseSTD);
  MemNoise.param(par);

  uniform_real_distribution<float>::param_type uni_par(0.0f,1.0f);
  uni_dist.param(uni_par);

  p1=0.0f;
  p2=0.0f;
  p3=0.0f;

}


float StiExt::ISti(float time)
{ return MemNoise(rdgen)+RSin.Irep(time)+RPulse.Irep(time);}

float StiExt::Igauss()
{ return MemNoise(rdgen);}

bool StiExt::R1()
{
  bool T=false;

  if(p1>uni_dist(rdgen))
    T=true;

  return T;
}

bool StiExt::R2()
{
  bool T=false;

  if(p2>uni_dist(rdgen))
    T=true;

  return T;
}

bool StiExt::R3()
{
  bool T=false;

  if(p3>uni_dist(rdgen))
    T=true;

  return T;
}


//----------------------------------------------------------------------
StiExt2::StiExt2()
{
  rdgen.seed(seed_gen());

  normal_distribution<float>::param_type par(MNoiseMean,MNoiseSTD);
  MemNoise.param(par);

  uniform_real_distribution<float>::param_type uni_par(0.0f,1.0f);
  uni_dist.param(uni_par);
}


float StiExt2::ISti(float time)
{ return MemNoise(rdgen)+RSin.Irep(time)+RPulse.Irep(time);}

float StiExt2::Igauss()
{ return MemNoise(rdgen);}

bool StiExt2::Rnd(float p)
{
  bool T=false;

  if(p>uni_dist(rdgen))
    T=true;

  return T;
}

//----------------------------------------------------------------------
// receptor for sinusoidal stimulation

SinGen::SinGen()
{
  amp=0.0*1e-12f;
  freq=0.0f;
  phase=0.0f;
}

void SinGen::SetPar(float a,float f,float p)
{
  amp=a;
  freq=f;
  phase=p;
}


void SinGen::GetPar(float& a,float& f,float& p)
{
  a=amp;
  f=freq;
  p=phase;
}

float SinGen::Irep(float time)
{  return amp*sin(2.0f*(float)M_PI*freq*time+phase); }

void SinGen::TexInfo(string& ostr)
{
  ostringstream ostr_tmp;

  ostr_tmp
      <<"amplitude(A) = "<<amp
      <<" frequency(Hz) = "<<freq
      <<" phase(0~2*PI) = "<<phase
      <<endl;

  ostr=ostr_tmp.str();
}

//----------------------------------------------------------------------
// receptor for pulase stimulation

PulseGen::PulseGen()
{
  amp=0.0*1e-12f;
  freq=0.0f;
  duty=0.0f;
  phase=0.0f;
  TBias=0.0f;
}

void PulseGen::SetPar(float a,float f, float d,float p)
{
  amp=a;
  freq=f;
  duty=d;
  phase=p;

  TBias=(float)(2.0f*M_PI-fmod( phase,(float)(2*M_PI)) ) / (float)(2.0f*M_PI*freq);
}


void PulseGen::GetPar(float& a,float& f,float& d,float& p)
{
  a=amp;
  f=freq;
  d=duty;
  p=phase;
}

float PulseGen::Irep(float time)
{
  float I;

  if(fmod((time + TBias),(1/freq)) <= fmod((time + TBias),(duty/freq)))
    I = amp;
  else
    I = 0.0;

  return I;
}

void PulseGen::TexInfo(string& ostr)
{
  ostringstream ostr_tmp;

  ostr_tmp
      <<"amplitude(A) = "<<amp
      <<" frequency(Hz) = "<<freq
      <<" duty(0~1) = "<<duty
      <<" phase(0~2*PI) = "<<phase
      <<endl;

  ostr=ostr_tmp.str();
}



