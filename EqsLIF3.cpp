
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <omp.h>
#include <thread>
#include <mutex>
#include "Eqs.h"

using namespace std;

/*

  float RK4(float *num,unsigned int num_size, void (EqsTyp::*FEqs)(float *,unsigned int)) //Runge Kutta order 4th
  {
    (pObj->*FEqs)( num, num_size);

    for(unsigned int i=0;i<num_size;i++)
    { *(num+i)=*(num+i) + *(k+i)*DTIME/2.0f; }
    (pObj->*FEqs)( num, num_size);

    for(unsigned int i=0;i<num_size;i++)
    { *(num+i)=*(num+i) + *(k+i)*DTIME/2.0f; }
    (pObj->*FEqs)( num, num_size);

    for(unsigned int i=0;i<num_size;i++)
    { *(num+i)=*(num+i) + *(k+i)*DTIME; }
    (pObj->*FEqs)( num, num_size);

    return num + DTIME*( k1 + 2.0f*k2 + 2.0f*k3 + k4 )/6.0f;
  };

void EqsLIF3::FFst(float *num)
{
  for(unsigned int i=0;i<st_size;i++)
  {
    *(k+i)=-*(num+i)/PItr->tst;
  }

}
*/

float EqsLIF3::FFst(float st)
{ return -st/PItr->tst; }

float EqsLIF3::Fxt(float xt)
{ return -xt/PItr2->txt; }

float EqsLIF3::Fst(float st)
{ return ( PItr2->as*VItr2->xt*(1.0f - st)*1e3f - st/PItr2->tst); }

float EqsLIF3::MemPot(float Vm)
{ return ( NItr->Isyn - NItr->gl * (Vm - NItr->El) ) / NItr->Cm ; }

//-----------------------------------------------------------------------
float EqsLIF3::FRK4Irep1Acc()
{
  float GatCnd=0.0f;
  //float y=0.0f;
  //float err=0.0f;

  for(VItr=Rep_begin,SItr=Sti_begin;VItr<Rep_end;VItr++,SItr++)
  {
    if(*SItr)
    {
      SetSti1();
      *SItr=false;
    }
    else
      VItr->st = RK4(VItr->st, &EqsLIF3::FFst);

    GatCnd += (VItr->st * VItr->g);

/*
    //Kahan summation algorithm
    y=(VItr->st * VItr->g)-err;
    t=GatCnd+y;
    err=(t-GatCnd)-y;
    GatCnd=t;
*/
  }

  float Isyn_tmp=GatCnd * ( PItr->Vrev - NItr->V );
  NItr->Isyn += Isyn_tmp;
  return Isyn_tmp;
}

float EqsLIF3::FRK4Irep2Acc()
{
  float GatCnd=0.0f;
  for(VItr2=Rep2_begin,SItr=Sti_begin;VItr2<Rep2_end;VItr2++,SItr++)
  {
    if(*SItr)
    {
      SetSti2();
      *SItr=false;
    }
    else
      VItr2->xt = RK4(VItr2->xt, &EqsLIF3::Fxt);

    VItr2->st = RK4(VItr2->st, &EqsLIF3::Fst);
    GatCnd += (VItr2->st * VItr2->g);
  }

  float Isyn_tmp=GatCnd * ( PItr2->Vrev - NItr->V ) / ( 1.0f + PItr2->mg * (float)exp( -0.062e3f * NItr->V ) / 3.57f );
  NItr->Isyn += Isyn_tmp;
  return Isyn_tmp;
}

float EqsLIF3::FIpvIrep1Acc()
{
  float GatCnd=0.0f;
  for(VItr=Rep_begin,SItr=Sti_begin;VItr<Rep_end;VItr++,SItr++)
  {
    if(*SItr)
    {
      SetSti1();
      *SItr=false;
    }
    else
      VItr->st = IpvEul(VItr->st, &EqsLIF3::FFst);

    GatCnd += (VItr->st * VItr->g);
  }

  float Isyn_tmp=GatCnd * ( PItr->Vrev - NItr->V );
  NItr->Isyn += Isyn_tmp;
  return Isyn_tmp;
}

float EqsLIF3::FIpvIrep2Acc()
{
  float GatCnd=0.0f;
  for(VItr2=Rep2_begin,SItr=Sti_begin;VItr2<Rep2_end;VItr2++,SItr++)
  {
    if(*SItr)
    {
      SetSti2();
      *SItr=false;
    }
    else
      VItr2->xt = IpvEul(VItr2->xt, &EqsLIF3::Fxt);

    VItr2->st = IpvEul(VItr2->st, &EqsLIF3::Fst);

    GatCnd += (VItr2->st * VItr2->g);
  }

  float Isyn_tmp=GatCnd * ( PItr2->Vrev - NItr->V ) / ( 1.0f + PItr2->mg * (float)exp( -0.062e3f * NItr->V ) / 3.57f );
  NItr->Isyn += Isyn_tmp;
  return Isyn_tmp;
}

float EqsLIF3::FEulIrep1Acc()
{
  float GatCnd=0.0f;
  for(VItr=Rep_begin,SItr=Sti_begin;VItr<Rep_end;VItr++,SItr++)
  {
    if(*SItr)
    {
      SetSti1();
      *SItr=false;
    }
    else
      VItr->st = VItr->st - DTIME*(VItr->st/PItr->tst);

    GatCnd += (VItr->st * VItr->g);
  }

  float Isyn_tmp=GatCnd * ( PItr->Vrev - NItr->V );
  NItr->Isyn += Isyn_tmp;
  return Isyn_tmp;
}

float EqsLIF3::FEulIrep2Acc()
{
  float GatCnd=0;
  for(VItr2=Rep2_begin,SItr=Sti_begin;VItr2<Rep2_end;VItr2++,SItr++)
  {
    if(*SItr)
    {
      SetSti2();
      *SItr=false;
    }
    else
      VItr2->xt = VItr2->xt - DTIME*(VItr2->xt / PItr2->txt);

    VItr2->st = VItr2->st + DTIME*(PItr2->as*VItr2->xt*(1.0f - VItr2->st)*1e3f - VItr2->st/PItr2->tst);
    GatCnd += (VItr2->st * VItr2->g);
  }

  float Isyn_tmp=GatCnd * ( PItr2->Vrev - NItr->V ) / ( 1.0f + PItr2->mg * (float)exp( -0.062e3f * NItr->V ) / 3.57f );
  NItr->Isyn += Isyn_tmp;
  return Isyn_tmp;
}

void EqsLIF3::firing()
{
    if(NItr->Count < NItr->RfPed) //Firing
    {

      if(NItr->Count == 0)
        NItr->V=0.0f;
      else
        NItr->V=NItr->Vreset*1e-3f;

      if(NItr->Count==NItr->SpiDy)
        NItr->SpiOut = true;
      else
        NItr->SpiOut = false;

      NItr->Count++;
    }
    else  // silence time is over
    {
      NItr->Count=0;
      NItr->Firing=false;
    }
}

void EqsLIF3::FRK4Vm()
{
  if(NItr->Firing)
  { firing(); }
  else
  {
    NItr->Isyn += NItr->Igauss();
    NItr->V = RK4(NItr->V, &EqsLIF3::MemPot);

    if(NItr->V > NItr->Vth)
    {
      NItr->Firing=true;
      firing();
    }
  }

}

void EqsLIF3::FIpvVm()
{
  if(NItr->Firing)
  { firing(); }
  else
  {
    NItr->Isyn += NItr->Igauss();
    NItr->V = IpvEul(NItr->V, &EqsLIF3::MemPot);

    if(NItr->V > NItr->Vth)
    {
      NItr->Firing=true;
      firing();
    }
  }

}

void EqsLIF3::FEulVm()
{
  if(NItr->Firing)
  { firing(); }
  else
  {
    NItr->Isyn += NItr->Igauss();
    NItr->V = NItr->V + DTIME * ( NItr->Isyn - NItr->gl * (NItr->V - NItr->El) ) / NItr->Cm;

    if(NItr->V > NItr->Vth)
    {
      NItr->Firing=true;
      firing();
    }
  }

}


void EqsLIF3::Spike(unsigned int NID)
{
  NItr=NItr_begin+NID;

  if(NItr->SpiOut)
  {
    for(vector<ConTab2>::iterator CItr=NItr->Index.begin();CItr<NItr->Index.end();CItr++)
    {
      if((CItr->Act & AMPA_ACT)!=0)
        *((NItr_begin + CItr->NeuID)->AMPA_Sti.begin() + CItr->AMPA_ID)=true;

      if((CItr->Act & GABA_ACT)!=0)
        *((NItr_begin + CItr->NeuID)->GABA_Sti.begin() + CItr->GABA_ID)=true;

      if((CItr->Act & NMDA_ACT)!=0)
        *((NItr_begin + CItr->NeuID)->NMDA_Sti.begin() + CItr->NMDA_ID)=true;
    }
  }
}


void EqsLIF3::Activation(unsigned int NID)
{
  NItr=NItr_begin+NID;

  //bcak ground initial
  NItr->BgInit();
  switch(SolType)
  {
    case ACCURATE:
      //synapse current calculation
      //AMBA receptor
      PItr = &NItr->ParAMPA;
      Rep_begin=NItr->RepAMPA.begin();
      Rep_end=NItr->RepAMPA.end();
      Sti_begin=NItr->AMPA_Sti.begin();
      NItr->Irep[0]=FRK4Irep1Acc();

      //GABA receptor
      PItr = &NItr->ParGABA;
      Rep_begin=NItr->RepGABA.begin();
      Rep_end=NItr->RepGABA.end();
      Sti_begin=NItr->GABA_Sti.begin();
      NItr->Irep[1]=FRK4Irep1Acc();

      //NMDA receptor
      PItr2 = &NItr->ParNMDA;
      Rep2_begin=NItr->RepNMDA.begin();
      Rep2_end=NItr->RepNMDA.end();
      Sti_begin=NItr->NMDA_Sti.begin();
      NItr->Irep[2]=FRK4Irep2Acc();

      //membrane potential
      FRK4Vm();
    break;

    case MODERATE:
      //synapse current calculation
      //AMBA receptor
      PItr = &NItr->ParAMPA;
      Rep_begin=NItr->RepAMPA.begin();
      Rep_end=NItr->RepAMPA.end();
      Sti_begin=NItr->AMPA_Sti.begin();
      NItr->Irep[0]=FIpvIrep1Acc();

      //GABA receptor
      PItr = &NItr->ParGABA;
      Rep_begin=NItr->RepGABA.begin();
      Rep_end=NItr->RepGABA.end();
      Sti_begin=NItr->GABA_Sti.begin();
      NItr->Irep[1]=FIpvIrep1Acc();

      //NMDA receptor
      PItr2 = &NItr->ParNMDA;
      Rep2_begin=NItr->RepNMDA.begin();
      Rep2_end=NItr->RepNMDA.end();
      Sti_begin=NItr->NMDA_Sti.begin();
      NItr->Irep[2]=FIpvIrep2Acc();

      //membrane potential
      FIpvVm();
    break;

    case ROUGH:
      //synapse current calculation
      //AMBA receptor
      PItr = &NItr->ParAMPA;
      Rep_begin=NItr->RepAMPA.begin();
      Rep_end=NItr->RepAMPA.end();
      Sti_begin=NItr->AMPA_Sti.begin();
      NItr->Irep[0]=FEulIrep1Acc();

      //GABA receptor
      PItr = &NItr->ParGABA;
      Rep_begin=NItr->RepGABA.begin();
      Rep_end=NItr->RepGABA.end();
      Sti_begin=NItr->GABA_Sti.begin();
      NItr->Irep[1]=FEulIrep1Acc();

      //NMDA receptor
      PItr2 = &NItr->ParNMDA;
      Rep2_begin=NItr->RepNMDA.begin();
      Rep2_end=NItr->RepNMDA.end();
      Sti_begin=NItr->NMDA_Sti.begin();
      NItr->Irep[2]=FEulIrep2Acc();

      //membrane potential
      FEulVm();
    break;

    default:
      //synapse current calculation
      //AMBA receptor
      PItr = &NItr->ParAMPA;
      Rep_begin=NItr->RepAMPA.begin();
      Rep_end=NItr->RepAMPA.end();
      Sti_begin=NItr->AMPA_Sti.begin();
      NItr->Irep[0]=FIpvIrep1Acc();

      //GABA receptor
      PItr = &NItr->ParGABA;
      Rep_begin=NItr->RepGABA.begin();
      Rep_end=NItr->RepGABA.end();
      Sti_begin=NItr->GABA_Sti.begin();
      NItr->Irep[1]=FIpvIrep1Acc();

      //NMDA receptor
      PItr2 = &NItr->ParNMDA;
      Rep2_begin=NItr->RepNMDA.begin();
      Rep2_end=NItr->RepNMDA.end();
      Sti_begin=NItr->NMDA_Sti.begin();
      NItr->Irep[2]=FIpvIrep2Acc();

      //membrane potential
      FIpvVm();
    break;
  }
}


