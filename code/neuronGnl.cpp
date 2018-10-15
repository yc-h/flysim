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

NeuronGnl::NeuronGnl()
{ ParInit(DTIME); }

NeuronGnl::~NeuronGnl()
{
  delete [] SynFst1;
  delete [] SynSlo1;
  delete [] SynGlu1;
  delete [] SynCnd1;

  delete [] SFTyp1;
  delete [] SSTyp1;
  delete [] SGTyp1;

  for(unsigned int i =0;i<SynFastSize;i++)
  {
    delete [] (SynFst+i)->PosID;
    delete [] (SynFst+i)->gFast;
  }
  delete [] SynFst;

  for(unsigned int i =0;i<SynSlowSize;i++)
  {
    delete [] (SynSlo+i)->PosID;
    delete [] (SynSlo+i)->gSlow;
  }
  delete [] SynSlo;

  for(unsigned int i =0;i<SynGlutSize;i++)
  {
    delete [] (SynGlu+i)->PosID;
    delete [] (SynGlu+i)->gFast;
    delete [] (SynGlu+i)->gSlow;
  }
  delete [] SynGlu;


  for(unsigned int i =0;i<SynGapSize;i++)
  {
    delete [] (SynCnd+i)->PosID;
    delete [] (SynCnd+i)->gC;
  }
  delete [] SynCnd;

  delete [] Irep;

  delete [] gAMPA;
  delete [] gGABA;
  delete [] gNMDA;
  delete [] gACH;
  delete [] gGCL;

}

void NeuronGnl::BgInit(float dt)
{
  if((ActCtr & AMPA_ACT) !=0)
  {
    if(Rnd(*pb))
      ExtFst->st += 1.0f;
    else
      ExtFst->st -= dt*(ExtFst->st/ParAMPA.tst);

    *gAMPA += ExtFst->g * ExtFst->st;
  }

  if((ActCtr & GABA_ACT) !=0)
  {
    if(Rnd(*(pb+1)))
      (ExtFst+1)->st += 1.0f;
    else
      (ExtFst+1)->st -= dt*((ExtFst+1)->st/ParGABA.tst);

    *gGABA += (ExtFst+1)->g*(ExtFst+1)->st;
  }

  if((ActCtr & NMDA_ACT) !=0)
  {
    if(Rnd(*(pb+2)))
      ExtSlo.xt += 1.0f;
    else
      ExtSlo.xt -= dt*(ExtSlo.xt/ParNMDA.txt);

    ExtSlo.st += dt*( ParNMDA.as*ExtSlo.xt*(1.0f - ExtSlo.st) -  ExtSlo.st/ParNMDA.tst);
    *gNMDA += ExtSlo.g*ExtSlo.st;
  }

  if((ActCtr & ACH_ACT) !=0)
  {
    if(Rnd(*(pb+3)))
      (ExtFst+2)->st += 1.0f;
     else
      (ExtFst+2)->st -= dt*((ExtFst+2)->st/ParACH.tst);

    *gACH += (ExtFst+2)->g*(ExtFst+2)->st;
  }

  if((ActCtr & GCL_ACT) !=0)
  {
    if(Rnd(*(pb+4)))
      (ExtFst+3)->st += 1.0f;
    else
      (ExtFst+3)->st -= dt*((ExtFst+3)->st/ParGCL.tst);

    *gGCL += (ExtFst+3)->g*(ExtFst+3)->st;
  }

  if(IgsEn){ Isyn+= Igauss(); }
}


void NeuronGnl::ParInit(float dt)
{
  ParAMPA.tst=2.0f;
  ParAMPA.Vrev=0.0f;

  ParGABA.tst=2.0f;
  ParGABA.Vrev=-90.0f;
  
  ParNMDA.as=0.6332f; // the alpha st
  ParNMDA.txt=2.0f; // the decay time constant of xt
  ParNMDA.tst=100.0f; //decay time constant for st
  ParNMDA.Vrev=0.0f; // reversal potential
  ParNMDA.mg=1.0f; // the concentration of magnesium

  ParACH.tst=20.0f;
  ParACH.Vrev=0.0f;

  ParGCL.tst=2.0f;
  ParGCL.Vrev=-90.0f;

// ParBase
  Vth=-50.0f; //-50mV
  El=VRESTING; //-70mV
  gl=2.5f;
  Cm=1.0;

  Vreset=VRESET;  // mV
  dVmax = (Vth + 90)/dt;

// ParLIF
  RfPed=(int)1.8f/dt+1;  // 1.8ms = (19(total) - 1(firing))*0.1m
  SpiDy=0;  // spike is delay 18 DTIME = 1.8ms

// ParHH
  Ek=-12.0f + VRESTING;
  Ena=115.0f + VRESTING;

  gk=36.0f;
  gna=120.0f;

// ParStp
  pv=0.6f; // reduce factor of depression factor D
  tD=300.0f; // decay time constant D

  aF=0.0035e-9f; // the alpha of facilitation F
  tF=7000.0f; // the decay time constant of F

  aCap=250.0e-6f; // the alpha of Ca^2+ peak
  tCap=2.0f; // decay time constant  Ca^2+ peak

  aCar=0.1e-6f; // alpha of residual Ca^2+
  tCar=2000.0f; // decay time constant residual Ca^2+

// ParLP
  tLP=50.0f;  // time constant of long term plsiticity factor
  SpiLPDy=0;
  PosISI=1e-5f;
  NegISI=1e-5f;

// synapse variables
  SynFst1=nullptr;
  SynGlu1=nullptr;
  SynSlo1=nullptr;
  SynCnd1=nullptr;

  SynFastSize1=0;
  SynGlutSize1=0;
  SynSlowSize1=0;
  SynGapSize1=0;

  SFTyp1=nullptr;
  SSTyp1=nullptr;
  SGTyp1=nullptr;

  SynFst=nullptr;
  SynGlu=nullptr;
  SynSlo=nullptr;
  SynCnd=nullptr;

  SynFastSize=0;
  SynGlutSize=0;
  SynSlowSize=0;
  SynGapSize=0;

  SFTyp=nullptr;
  SSTyp=nullptr;
  SGTyp=nullptr;
}

void NeuronGnl::VarInit(unsigned int thsize)
{

//Base
  Isyn=0.0f;  // current of synapse
  V=VRESTING;  // menbrane potential

//LIF
  Count=0;  // firing tabe index
  Firing=false;  // firing flag, true for firing
  SpiOut=false;  // spike of this neuron, true for spike

//STP
  D=1.0f;  // depression factor D
  F=1.0f;  // facilitation factor F
  Cap=0.0f*1e-6f; // Ca^2+ peak level
  Car=20.0f*1e-6f;  // residual level of Ca^2+
  //Ca=0.0f;  // Ca^2+ level

//VarHH
  n=0.3176769f;
  m=0.05293248f;
  h=0.5719156f;

//LTP
  LP=0.0f;   // long term plsiticity factor
  FirLP=false;
  SpiLP=false;
  CountLP=0;

//gate variables
  for(unsigned int i=0;i<thsize;i++)
  {
    *(gAMPA+i)=0.0f;
    *(gGABA+i)=0.0f;
    *(gNMDA+i)=0.0f;
    *(gACH+i)=0.0f;
    *(gGCL+i)=0.0f;
  }

//synapses
  if(SynFastSize1!=0)
  {
    for(unsigned int i=0;i<SynFastSize1;i++)
    { (SynFst1+i)->stFast=0.0f; }
  }

  if(SynSlowSize1!=0)
  {
    for(unsigned int i=0;i<SynSlowSize1;i++)
    {
      (SynSlo1+i)->xtSlow=0.0f;
      (SynSlo1+i)->stSlow=0.0f;
    }
  }

  if(SynGlutSize1!=0)
  {
    for(unsigned int i=0;i<SynGlutSize1;i++)
    {
      (SynGlu1+i)->stFast=0.0f;
      (SynGlu1+i)->xtSlow=0.0f;
      (SynGlu1+i)->stSlow=0.0f;
    }
  }


  if(SynFastSize!=0)
  {
    for(unsigned int i=0;i<SynFastSize;i++)
    { (SynFst+i)->stFast=0.0f; }
  }

  if(SynSlowSize!=0)
  {
    for(unsigned int i=0;i<SynSlowSize;i++)
    {
      (SynSlo+i)->xtSlow=0.0f;
      (SynSlo+i)->stSlow=0.0f;
    }
  }

  if(SynGlutSize!=0)
  {
    for(unsigned int i=0;i<SynGlutSize;i++)
    {
      (SynGlu+i)->stFast=0.0f;
      (SynGlu+i)->xtSlow=0.0f;
      (SynGlu+i)->stSlow=0.0f;
    }
  }

//Ext stimulation
  for(unsigned int ix=0;ix<5;ix++)
  { *(Irep+ix)=0.0f; }

  for(unsigned int ix=0;ix<4;ix++)
  { (ExtFst+ix)->st=0.0f; }

  ExtSlo.xt=0.0f;
  ExtSlo.st=0.0f;
}

void NeuronGnl::TexInfo(string& ostr)
{
  ostringstream ostr_tmp;
  string str;

  ostr_tmp
    <<"ActCtr="<<ActCtr
    <<" Vth="<<Vth<<"mv"
    <<" El="<<El<<"mv"
    <<" g="<<gl<<"uS"
    <<" Cm="<<Cm<<"nF"
    <<" Vreset="<<Vreset<<"mv"
    <<" dVmax="<<dVmax<<"mv\n";

  ostr_tmp
    <<"AMPA tst="<<ParAMPA.tst<<"ms"
    <<" Vrev="<<ParAMPA.Vrev<<"mV"
    <<" GABA tst="<<ParGABA.tst<<"ms"
    <<" Vrev="<<ParGABA.Vrev<<"mV\n"

    <<"NMDA tst="<<ParNMDA.tst<<"ms"
    <<" Vrev="<<ParNMDA.Vrev<<"mV"
    <<" txt="<<ParNMDA.txt<<"ms"
    <<" as="<<ParNMDA.as
    <<" mg="<<ParNMDA.mg<<"\n"

    <<"ACH tst="<<ParACH.tst<<"ms"
    <<" Vrev="<<ParACH.Vrev<<"mV"
    <<" GCL tst="<<ParGCL.tst<<"ms"
    <<" Vrev="<<ParGCL.Vrev<<"mV\n";

//LIF
  if((ActCtr&LIF_ACT) !=0)
  {
    ostr_tmp
      <<"LIF:"
      <<" RfPed="<<RfPed
      <<" SpiDy="<<SpiDy<<"\n";
  }

//STP
  if((ActCtr&STP_ACT) !=0)
  {
    ostr_tmp
      <<"STP:"
      <<" pv="<<pv
      <<" tD="<<tD
      <<" aF="<<aF
      <<" tF="<<tF
      <<" aCap="<<aCap
      <<" tCap="<<tCap
      <<" aCar="<<aCar
      <<" tCar="<<tCar<<"\n";
  }

//HH
  if((ActCtr&HH_ACT) !=0)
  {
    ostr_tmp
      <<"HH:"
      <<" Ek="<<Ek
      <<" gk="<<gk
      <<" Ena="<<Ena
      <<" gna="<<gna<<"\n";
  }

//LTP
  if((ActCtr&LTP_ACT) !=0)
  {
    ostr_tmp
      <<"LTP:"
      <<" tLP="<<tLP
      <<" SpiLPDy="<<SpiLPDy
      <<" PosISI="<<PosISI
      <<" NegISI="<<NegISI<<"\n";
  }

//synapses
  if(SynFastSize1!=0)
  {
    ostr_tmp<<"SynFst1=TotalNum;Type,g,PosID;Type,g,PosID;...\n";
    ostr_tmp<<"SynFst1="<<SynFastSize1;
    for(unsigned int i=0;i<SynFastSize1;i++)
    { ostr_tmp<<";"<<(unsigned int)*(SFTyp1+i)<<","<<(SynFst1+i)->gFast<<","<<(SynFst1+i)->PosID; }
    ostr_tmp<<"\n";
  }

  if(SynSlowSize1!=0)
  {
    ostr_tmp<<"SynSlo1=TotalNum;Type,g,PosID;Type,g,PosID;...\n";
    ostr_tmp<<"SynSlo1="<<SynSlowSize1;
    for(unsigned int i=0;i<SynSlowSize1;i++)
    { ostr_tmp<<";"<<(SynSlo1+i)->gSlow<<","<<(SynSlo1+i)->PosID; }
    ostr_tmp<<"\n";
  }

  if(SynGlutSize1!=0)
  {
    ostr_tmp<<"SynGlu1=TotalNum;Type,gFast,gSlow,PosID;Type,gFast,gSlow,PosID;...\n";
    ostr_tmp<<"SynGlu1="<<SynGlutSize1;
    for(unsigned int i=0;i<SynGlutSize1;i++)
    {
      ostr_tmp
        <<";"<<(unsigned int)SGTyp1->A<<","<<(SynGlu1+i)->gFast
        <<","<<(SynGlu1+i)->gSlow
        <<","<<(SynGlu1+i)->PosID;
    }
    ostr_tmp<<"\n";
  }

  if(SynGapSize1!=0)
  {
    ostr_tmp<<"SynCnd1=TotalNum;g,PosID;g,PosID;...\n";
    ostr_tmp<<"SynCnd1="<<SynGapSize1;
    for(unsigned int i=0;i<SynGapSize1;i++)
    { ostr_tmp<<";"<<(SynCnd1+i)->gC<<","<<(SynCnd1+i)->PosID; }
    ostr_tmp<<"\n";
  }

  unsigned int SynSize=0;

  if(SynFastSize!=0)
  {
    if((ActCtr&LTP_ACT) !=0)
      ostr_tmp<<"SynFst=SynNum,TotalNum;Type,PosNum,g,PosID,...,g,PosID;Type,PosNum,g,PosID,...,g,PosID;...\n";
    else
      ostr_tmp<<"SynFst=SynNum,TotalNum;Type,PosNum,g,PosID,...PosID;Type,PosNum,g,PosID,...PosID;...\n";

    for(unsigned int i=0;i<SynFastSize;i++)
    { SynSize += (SynFst+i)->PosIDNum; }

    ostr_tmp<<"SynFst="<<SynFastSize<<","<<SynSize;
    for(unsigned int i=0;i<SynFastSize;i++)
    {
      ostr_tmp<<";"<<(unsigned int)*(SFTyp+i)<<","<<(SynFst+i)->PosIDNum;
      if((ActCtr&LTP_ACT) !=0)
      {
        for(unsigned int j=0;j<(SynFst+i)->PosIDNum;j++)
        { ostr_tmp<<","<<*((SynFst+i)->gFast+j)<<","<<*((SynFst+i)->PosID+j); }
      }
      else
      {
        ostr_tmp<<","<<*((SynFst+i)->gFast);
        for(unsigned int j=0;j<(SynFst+i)->PosIDNum;j++)
        { ostr_tmp<<","<<*((SynFst+i)->PosID+j); }
      }
    }
    ostr_tmp<<"\n";
  }

  SynSize=0;
  if(SynSlowSize!=0)
  {
    if((ActCtr&LTP_ACT) !=0)
      ostr_tmp<<"SynSlo=SynNum,TotalNum;PosNum,g,PosID,...,g,PosID;PosNum,g,PosID,...,g,PosID;...\n";
    else
      ostr_tmp<<"SynSlo=SynNum,TotalNum;PosNum,g,PosID,...PosID;PosNum,g,PosID,...PosID;...\n";

    for(unsigned int i=0;i<SynSlowSize;i++)
    { SynSize += (SynSlo+i)->PosIDNum; }

    ostr_tmp<<"SynSlo="<<SynSlowSize<<","<<SynSize;
    for(unsigned int i=0;i<SynSlowSize;i++)
    {
      ostr_tmp<<";"<<(SynSlo+i)->PosIDNum;
      if((ActCtr&LTP_ACT) !=0)
      {
        for(unsigned int j=0;j<(SynSlo+i)->PosIDNum;j++)
        { ostr_tmp<<","<<*((SynSlo+i)->gSlow+j)<<","<<*((SynSlo+i)->PosID+j); }
      }
      else
      {
        ostr_tmp<<","<<*((SynSlo+i)->gSlow);
        for(unsigned int j=0;j<(SynSlo+i)->PosIDNum;j++)
        { ostr_tmp<<","<<*((SynSlo+i)->PosID+j); }
      }
    }
    ostr_tmp<<"\n";
  }

  SynSize=0;
  if(SynGlutSize!=0)
  {
    if((ActCtr&LTP_ACT) !=0)
      ostr_tmp<<"SynGlu=SynNum,TotalNum;Type,PosNum,gFast,gSlow,PosID,...,gFast,gSlow,PosID;Type,PosNum,gFast,gSlow,PosID,...,gFast,gSlow,PosID;...\n";
    else
      ostr_tmp<<"SynGlu=SynNum,TotalNum;Type,PosNum,gFast,gSlow,PosID,...,PosID;Type,PosNum,gFast,gSlow,PosID,...,PosID;...\n";

    for(unsigned int i=0;i<SynGlutSize;i++)
    { SynSize += (SynGlu+i)->PosIDNum; }

    ostr_tmp<<"SynGlu="<<SynGlutSize<<","<<SynSize;
    for(unsigned int i=0;i<SynGlutSize;i++)
    {
      ostr_tmp<<";"<<(unsigned int)(SGTyp+i)->A<<","<<(SynGlu+i)->PosIDNum;
      if((ActCtr&LTP_ACT) !=0)
      {
        for(unsigned int j=0;j<(SynGlu+i)->PosIDNum;j++)
        {
          ostr_tmp<<","<<*((SynGlu+i)->gFast+j)<<","<<*((SynGlu+i)->gSlow+j)
                  <<","<<*((SynGlu+i)->PosID+j);
        }
      }
      else
      {
        ostr_tmp<<","<<*((SynGlu+i)->gFast)<<","<<*((SynGlu+i)->gSlow);
        for(unsigned int j=0;j<(SynGlu+i)->PosIDNum;j++)
        { ostr_tmp<<","<<*((SynGlu+i)->PosID+j); }
      }
    }
    ostr_tmp<<"\n";
  }

  SynSize=0;
  if(SynGapSize!=0)
  {
    if((ActCtr&LTP_ACT) !=0)
      ostr_tmp<<"SynCnd=SynNum,TotalNum;PosNum,g,PosID,...,g,PosID;PosNum,g,PosID,...,g,PosID;...\n";
    else
      ostr_tmp<<"SynCnd=SynNum,TotalNum;PosNum,g,PosID,...PosID;PosNum,g,PosID,...PosID;...\n";

    for(unsigned int i=0;i<SynGapSize;i++)
    { SynSize += (SynCnd+i)->PosIDNum; }

    ostr_tmp<<"SynCnd="<<SynGapSize<<","<<SynSize;
    for(unsigned int i=0;i<SynGapSize;i++)
    {
      ostr_tmp<<";"<<(SynCnd+i)->PosIDNum;
      if((ActCtr&LTP_ACT) !=0)
      {
        for(unsigned int j=0;j<(SynCnd+i)->PosIDNum;j++)
        { ostr_tmp<<","<<*((SynCnd+i)->gC+j)<<","<<*((SynCnd+i)->PosID+j); }
      }
      else
      {
        ostr_tmp<<*((SynCnd+i)->gC);
        for(unsigned int j=0;j<(SynCnd+i)->PosIDNum;j++)
        { ostr_tmp<<","<<*((SynCnd+i)->PosID+j); }
      }
    }
    ostr_tmp<<"\n";
  }
//----------------------------noise--------------------------------------
  ostr_tmp<<"noise parameters:"<<"\n";

  ostr_tmp
    <<"Injection current:"
    <<"mean="<<Inj.mean()
    <<" stddev="<<Inj.stddev()
    <<"\n";

  ostr_tmp<<"StiExt AMPA:"<<*pb<<"\n";
  ostr_tmp<<"StiExt GABA:"<<*(pb+1)<<"\n";
  ostr_tmp<<"StiExt NMDA:"<<*(pb+2)<<"\n";
  ostr_tmp<<"StiExt Ach:"<<*(pb+3)<<"\n";
  ostr_tmp<<"StiExt GluCl:"<<*(pb+4)<<"\n";

  ostr=ostr_tmp.str();
}

