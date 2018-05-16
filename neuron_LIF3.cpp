#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <omp.h>
#include <random>
#include <chrono>
#include "neuron.h"

using namespace std;
extern random_device seed_gen;


Neuron_LIF3::Neuron_LIF3()
{
  Irep[0]=0.0f;
  Irep[1]=0.0f;
  Irep[2]=0.0f;
}

Neuron_LIF3::~Neuron_LIF3()
{
  {vector<ConTab2> ().swap(Index);}
  {vector<RV1> ().swap(RepAMPA);}
  {vector<RV1> ().swap(RepGABA);}
  {vector<RV2> ().swap(RepNMDA);}
  {vector<bool> ().swap(AMPA_Sti);}
  {vector<bool> ().swap(GABA_Sti);}
  {vector<bool> ().swap(NMDA_Sti);}
}

void Neuron_LIF3::BgInit()
{
  Isyn=0.0f; //clear Isynapse
  Irep[0]=0.0f;
  Irep[1]=0.0f;
  Irep[2]=0.0f;

  if(!AMPA_Sti.empty())
    AMPA_Sti[AMPA_Sti.size()-1] = R1(); //back ground poisson noise

  if(!GABA_Sti.empty())
    GABA_Sti[GABA_Sti.size()-1] = R2(); //back ground poisson noise

  if(!NMDA_Sti.empty())
    NMDA_Sti[NMDA_Sti.size()-1] = R3(); //back ground poisson noise
}


void Neuron_LIF3::ParInit()
{
  //NeuPar();
  Vth=-50.0f*1e-3f; //-50mV
  El=-70.0f*1e-3f; //-70mV
  gl=2.5f*1e-9f;
  Cm=1.0f*1e-9f;

  Vreset=VRESET;  // mV

  RfPed=21;  // 1.8ms = (21(total) - 3(firing))*0.1m
  SpiDy=0;  // spike is delay 18 DTIME = 1.8ms

  //StiExt();
  p1=0.0f;
  p2=0.0f;
  p3=0.0f;

#ifdef CPP11
  rnd_seed=seed_gen();
  rdgen.seed(rnd_seed);
  normal_distribution<float>::param_type par(MNoiseMean,MNoiseSTD);
  MemNoise.param(par);
#else
  mean=MNoiseMean;
  stddev=MNoiseSTD;
#endif

}

void Neuron_LIF3::VarInit()
{
  //NeuVar();
  Isyn=0.0f;  // current of synapse
  V=VRESTING*1e-3f;  // resting menbrane potential

  Count=0;  // firing tabe index
  Firing=false;  // firing flag, true for firing
  SpiOut=false;  // spike of this neuron, true for spike

  for(unsigned int i=0;i<RepAMPA.size();i++)
  {
    AMPA_Sti[i]=false;
    RepAMPA[i].st=0.0f;
  }

  for(unsigned int i=0;i<RepGABA.size();i++)
  {
    GABA_Sti[i]=false;
    RepGABA[i].st=0.0f;
  }

  for(unsigned int i=0;i<RepNMDA.size();i++)
  {
    NMDA_Sti[i]=false;
    RepNMDA[i].st=0.0f;
    RepNMDA[i].xt=0.0f;
  }
}

void Neuron_LIF3::TexInfo(string& ostr)
{
  ostringstream ostr_tmp;
  string str;

  ostr_tmp
    <<"Vth="<<Vth<<"v"
    <<" El="<<El<<"v"
    <<" g="<<gl<<"S"
    <<" Cm="<<Cm<<"F"
    <<" Vreset="<<Vreset<<"mv"
    <<"\n"
    <<"Axons="<<Index.size()
    <<" AMPAs="<<RepAMPA.size()
    <<" GABAs="<<RepGABA.size()
    <<" NMDAs="<<RepNMDA.size()
    <<"\n";
//----------------------------axion--------------------------------------
  ostr_tmp<<"Axons:"<<"\n";

  for(unsigned int i=0;i<Index.size();i++)
  {
    ostr_tmp
      <<"AxonID="<<i
      <<" NeuID="<<Index[i].NeuID
      <<" AMPA_ID="<<Index[i].AMPA_ID
      <<" GABA_ID="<<Index[i].GABA_ID
      <<" NMDA_ID="<<Index[i].NMDA_ID;

    if((Index[i].Act & AMPA_ACT)!=0)
      ostr_tmp<<" AMPA_Act=1";
    else
      ostr_tmp<<" AMPA_Act=0";

    if((Index[i].Act & GABA_ACT)!=0)
      ostr_tmp<<" GABA_Act=1";
    else
      ostr_tmp<<" GABA_Act=0";

    if((Index[i].Act & NMDA_ACT)!=0)
      ostr_tmp<<" NMDA_Act=1";
    else
      ostr_tmp<<" NMDA_Act=0";

      ostr_tmp<<"\n";
  }

  ostr_tmp<<"\n";

//----------------------------receptors--------------------------------------
  ostr_tmp<<"receptors:"<<"\n";

  ostr_tmp
    <<"AMPA Parameter:"
    <<" tau st="<<ParAMPA.tst
    <<" Reversal potential="<<ParAMPA.Vrev<<"v\n"
    <<"AMPA_ID:g(S)";

  for(unsigned int i=0;i<RepAMPA.size();i++)
  {
    if(i%10==0) ostr_tmp<<"\n";

    ostr_tmp<<i<<":"<<RepAMPA[i].g<<",";
  }
  ostr_tmp<<"\n\n";

  ostr_tmp
    <<"GABA Parameter:"
    <<" tau st="<<ParGABA.tst
    <<" Reversal potential="<<ParGABA.Vrev<<"v\n"
    <<"GABA_ID:g(S)";

  for(unsigned int i=0;i<RepGABA.size();i++)
  {
    if(i%10==0) ostr_tmp<<"\n";

    ostr_tmp<<i<<":"<<RepGABA[i].g<<",";
  }
  ostr_tmp<<"\n\n";


  ostr_tmp
    <<"NMDA Parameter:"
    <<" alpha s="<<ParNMDA.as
    <<" tau st="<<ParNMDA.tst
    <<" tau xt="<<ParNMDA.txt
    <<" mg="<<ParNMDA.mg
    <<" Reversal potential="<<ParNMDA.Vrev<<"v\n"
    <<"NMDA_ID:g(S)";

  for(unsigned int i=0;i<RepNMDA.size();i++)
  {
    if(i%10==0) ostr_tmp<<"\n";
    ostr_tmp<<i<<":"<<RepNMDA[i].g<<",";
  }
  ostr_tmp<<"\n\n";

//----------------------------noise--------------------------------------
  ostr_tmp<<"noise parameters:"<<"\n";

  ostr_tmp
    <<"Membrane injection current:"
    <<"mean="<<MemNoise.mean()
    <<" stddev="<<MemNoise.stddev()
    <<"\n";

  ostr_tmp<<"external stimulation probability of AMPA:"<<p1<<"\n";
  ostr_tmp<<"external stimulation probability of GABA:"<<p2<<"\n";
  ostr_tmp<<"external stimulation probability of NMDA:"<<p3<<"\n";

  RSin.Info(str);
  ostr_tmp<<"FExtSin:"<<str;
  RPulse.Info(str);
  ostr_tmp<<"FExtPulse:"<<str;

  ostr_tmp<<"this neuron connects to "<<Index.size()<<" other neurons, "
      <<"total receptors in this neuron is:"
      <<RepAMPA.size()+RepGABA.size()+RepNMDA.size()
      <<"\n------------------------------------------------------------------\n";

  ostr=ostr_tmp.str();
}

void Neuron_LIF3::TexSInfo(string& ostr)
{
  ostringstream ostr_tmp;
  //string str;

  ostr_tmp
    <<Vth<<" "
    <<El<<" "
    <<gl<<" "
    <<Cm<<" "
    <<Vreset<<" "
    <<RfPed<<" "  // refractory peroid
    <<SpiDy<<" " // spiking delay
    <<Isyn<<" "
    <<V<<" "
    <<Count<<" "
    <<Firing<<" "
    <<SpiOut<<" "
    <<Index.size()<<" "
    <<RepAMPA.size()<<" "
    <<RepGABA.size()<<" "
    <<RepNMDA.size()<<"\n";
//----------------------------axion--------------------------------------
  for(unsigned int i=0;i<Index.size();i++)
  {
    ostr_tmp
      <<Index[i].NeuID<<" "
      <<Index[i].AMPA_ID<<" "
      <<Index[i].GABA_ID<<" "
      <<Index[i].NMDA_ID<<" "
      <<Index[i].Act<<" ";
  }

  ostr_tmp<<"\n";
//----------------------------receptors--------------------------------------
  ostr_tmp
    <<ParAMPA.tst<<" "
    <<ParAMPA.Vrev<<" ";
  for(unsigned int i=0;i<RepAMPA.size();i++)
  {
    ostr_tmp
      <<AMPA_Sti[i]<<" "
      <<RepAMPA[i].g<<" "
      <<RepAMPA[i].st<<" ";
  }
  ostr_tmp<<"\n";

  ostr_tmp
    <<ParGABA.tst<<" "
    <<ParGABA.Vrev<<" ";
  for(unsigned int i=0;i<RepGABA.size();i++)
  {
    ostr_tmp
      <<GABA_Sti[i]<<" "
      <<RepGABA[i].g<<" "
      <<RepGABA[i].st<<" ";
  }
  ostr_tmp<<"\n";

  ostr_tmp
    <<ParNMDA.as<<" "
    <<ParNMDA.tst<<" "
    <<ParNMDA.txt<<" "
    <<ParNMDA.mg<<" "
    <<ParNMDA.Vrev<<" ";
  for(unsigned int i=0;i<RepNMDA.size();i++)
  {
    ostr_tmp
      <<NMDA_Sti[i]<<" "
      <<RepNMDA[i].g<<" "
      <<RepNMDA[i].st<<" "
      <<RepNMDA[i].xt<<" ";
  }
  ostr_tmp<<"\n";

//----------------------------noise--------------------------------------
  ostr_tmp
    <<MemNoise.mean()<<" "
    <<MemNoise.stddev()<<" "
    <<p1<<" "
    <<p2<<" "
    <<p3<<" "
    <<rnd_seed<<"\n";

  ostr=ostr_tmp.str();
}





