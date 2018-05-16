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

Neuron_gp4::Neuron_gp4()
{
  gAMPA=0.0f;
  gGABA=0.0f;
  gNMDA=0.0f;
}

Neuron_gp4::~Neuron_gp4()
{
  {vector<ConTab_gp2> ().swap(Index);}
  {vector<RepVar1_gp> ().swap(RepAMPA);}
  {vector<RepVar1_gp> ().swap(RepGABA);}
  {vector<RepVar1_gp> ().swap(RepNMDA);}
  {vector<NeuVar> ().swap(NeuVm);}
}


void Neuron_gp4::BgInit()
{
  for(unsigned int i=0;i<NeuVm.size();i++)
  { NeuVm[i].Isyn=0.0f;} //clear Isynapse
}

void Neuron_gp4::ParInit()
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

void Neuron_gp4::VarInit()
{

  gAMPA=0.0f;
  gGABA=0.0f;
  gNMDA=0.0f;

  //NeuVar();
  for(unsigned int i=0;i<NeuVm.size();i++)
  {
    NeuVm[i].Isyn=0.0f; //clear Isynapse
    NeuVm[i].V=VRESTING*1e-3f;  // resting menbrane potential
    NeuVm[i].Count=0;  // firing tabe index
    NeuVm[i].Firing=false;  // firing flag, true for firing
    NeuVm[i].SpiOut=false;  // spike of this neuron, true for spike
  }

  for(unsigned int i=0;i<RepAMPA.size();i++)
  { RepAMPA[i].st=0.0f;}

  for(unsigned int i=0;i<RepGABA.size();i++)
  { RepGABA[i].st=0.0f;}

  for(unsigned int i=0;i<RepNMDA.size();i++)
  { RepNMDA[i].st=0.0f; }

}

void Neuron_gp4::TexInfo(string& ostr)
{
  ostringstream ostr_tmp;
  string str;

  ostr_tmp
    <<"population parameters:\n"
    <<"totol neurons="<<NeuVm.size()
    <<" Vth="<<Vth<<"v"
    <<" El="<<El<<"v"
    <<" g="<<gl<<"S"
    <<" Cm="<<Cm<<"F"
    <<" Vreset="<<Vreset<<"mv"
    <<" SelfSti="<<SelfSti
    <<"\n\n";


//----------------------------axion--------------------------------------
  if(Index.size()!=0)
  {
    ostr_tmp<<"post synaptic populaitons:"<<"\n";

    for(unsigned int i=0;i<Index.size();i++)
    {
      ostr_tmp
        <<"Index="<<i
        <<" populaiton ID="<<Index[i].NeuID;

      if(Index[i].AMPA_Act)
      {
        ostr_tmp
          <<" AMPA_MeanG="<<Index[i].AMPA_MeanG
          <<" AMPA_tst="<<Index[i].AMPA_tst;
      }

      if(Index[i].GABA_Act)
      {
        ostr_tmp
          <<" GABA_MeanG="<<Index[i].GABA_MeanG
          <<" GABA_tst="<<Index[i].GABA_tst;
      }

      if(Index[i].NMDA_Act)
      {
        ostr_tmp
          <<" NMDA_MeanG="<<Index[i].NMDA_MeanG
          <<" NMDA_txt="<<Index[i].NMDA_txt
          <<" NMDA_tst="<<Index[i].NMDA_tst
          <<" NMDA_as="<<Index[i].NMDA_as;
      }

      if(Index[i].GAP_Act)
      {
        ostr_tmp
          <<" GAP_MeanG="<<Index[i].GAP_MeanG;
      }

      ostr_tmp<<"\n";
    }
  }
  else
    ostr_tmp<<"post synaptic populaitons:non"<<"\n";

  ostr_tmp<<"\n";

//----------------------------receptors--------------------------------------
  ostr_tmp<<"receptors is each neuron:"<<"\n";

  ostr_tmp
    <<"AMPA Parameter:"
    <<" tau st="<<ParAMPA.tst
    <<" Reversal potential="<<ParAMPA.Vrev<<"v"
    <<" g="<<ParAMPA.g<<"S"
    <<" MeanG="<<ParAMPA.MeanG<<"S"
    <<"\n";

  ostr_tmp
    <<"GABA Parameter:"
    <<" tau st="<<ParGABA.tst
    <<" Reversal potential="<<ParGABA.Vrev<<"v"
    <<" g="<<ParGABA.g<<"S"
    <<" MeanG="<<ParGABA.MeanG<<"S"
    <<"\n";

  ostr_tmp
    <<"NMDA Parameter:"
    <<" alpha s="<<ParNMDA.as
    <<" tau st="<<ParNMDA.tst
    <<" tau xt="<<ParNMDA.txt
    <<" mg="<<ParNMDA.mg
    <<" Reversal potential="<<ParNMDA.Vrev<<"v"
    <<" g="<<ParNMDA.g<<"S"
    <<" MeanG="<<ParNMDA.MeanG<<"S"
    <<"\n\n";

//----------------------------noise--------------------------------------
  ostr_tmp<<"noise parameters:"<<"\n";

#ifdef CPP11
  ostr_tmp
    <<"Membrane injection current:"
    <<"mean="<<MemNoise.mean()
    <<" stddev="<<MemNoise.stddev()
    <<"\n";
#else
  ostr_tmp
    <<"Membrane injection current:"
    <<"mean="<<mean
    <<" stddev="<<stddev
    <<"\n";
#endif

  ostr_tmp<<"external stimulation probability of AMPA:"<<p1<<"\n";
  ostr_tmp<<"external stimulation probability of GABA:"<<p2<<"\n";
  ostr_tmp<<"external stimulation probability of NMDA:"<<p3<<"\n";

  RSin.Info(str);
  ostr_tmp<<"FExtSin:"<<str;
  RPulse.Info(str);
  ostr_tmp<<"FExtPulse:"<<str;

  ostr_tmp<<"this populaiton connects to "<<Index.size()<<" other populaitons, "
      <<"total neurons in this populaiton is:"
      <<NeuVm.size()
      <<"\n------------------------------------------------------------------\n";

  ostr=ostr_tmp.str();
}





