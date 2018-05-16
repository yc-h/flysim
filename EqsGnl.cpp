
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <omp.h>
#include <thread>
#include <mutex>
#include "Eqs.h"

using namespace std;

float EqsGnl::FFst(float st,float tst)
{ return (-st / tst); }

float EqsGnl::Fst(float st,float tst,float xt,float as)
{ return (as*xt*(1.0f - st) - st/tst);}

float EqsGnl::MemPot(float Vm)
{
  float gAMPA_ACC=0.0f;
  float gGABA_ACC=0.0f;
  float gNMDA_ACC=0.0f;
  float gACH_ACC=0.0f;
  float gGCL_ACC=0.0f;
  for(unsigned int i=0;i<thsize;i++)
  {
    gAMPA_ACC += *(NItr->gAMPA+i);
    gGABA_ACC += *(NItr->gGABA+i);
    gNMDA_ACC += *(NItr->gNMDA+i);
    gACH_ACC += *(NItr->gACH+i);
    gGCL_ACC += *(NItr->gGCL+i);
  }

  *(NItr->Irep+0)=gAMPA_ACC*(NItr->ParAMPA.Vrev - Vm);
  *(NItr->Irep+1)=gGABA_ACC*(NItr->ParGABA.Vrev - Vm);
  *(NItr->Irep+2)=gNMDA_ACC*(NItr->ParNMDA.Vrev - Vm)/( 1.0f + NItr->ParNMDA.mg * (float)exp( -0.062f * Vm ) / 3.57f );
  *(NItr->Irep+3)=gACH_ACC*(NItr->ParACH.Vrev - Vm);
  *(NItr->Irep+4)=gGCL_ACC*(NItr->ParGCL.Vrev - Vm);

  float I =
    + *(NItr->Irep+0)
    + *(NItr->Irep+1)
    + *(NItr->Irep+2)
    + *(NItr->Irep+3)
    + *(NItr->Irep+4)
    + NItr->gl * (NItr->El - Vm)*1e3f
    + NItr->Isyn*1e3f;

  I=1e-3f*I/NItr->Cm;

  if(I>NItr->dVmax)
    I=NItr->dVmax;
  else if(I<-NItr->dVmax)
    I=-NItr->dVmax;

  return I;

}

void EqsGnl::firing()
{
    if(NItr->Count < NItr->RfPed) //Firing
    {

      if((NItr->ActCtr & LIF_ACT) != 0)
      {
        if(NItr->Count == 0)
          NItr->V=0.0f;
        else
          NItr->V=NItr->Vreset;
      }

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


//---------- hogkin huxily model
float EqsGnl::An()
{ return 0.01f*(10.0f-Vx)/(exp((10.0f-Vx)*0.1f)-1.0f); }

float EqsGnl::Bn()
{ return 0.125f*exp(-Vx/80.0f); }

float EqsGnl::Am()
{ return 0.1f*(25.0f-Vx)/(exp((25.0f-Vx)*0.1f)-1.0f); }

float EqsGnl::Bm()
{ return 4.0f*exp(-Vx/18.0f); }

float EqsGnl::Ah()
{ return 0.07f*exp(-Vx/20.0f); }

float EqsGnl::Bh()
{ return 1.0f/(exp((30.0f-Vx)*0.1f)+1.0f); }

float EqsGnl::dx(float x)
{ return Ax*(1.0f-x)-Bx*x; }

float EqsGnl::IHHRK4()
{
  Vx=NItr->V - VRESTING;
  Ax=An();
  Bx=Bn();
  NItr->n=RK4(NItr->n, &EqsGnl::dx);

  Ax=Am();
  Bx=Bm();
  NItr->m=RK4(NItr->m, &EqsGnl::dx);

  Ax=Ah();
  Bx=Bh();
  NItr->h=RK4(NItr->h, &EqsGnl::dx);

  return NItr->gk*pow(NItr->n,4.0f)*(NItr->Ek - NItr->V) + NItr->gna*pow(NItr->m,3.0f)*NItr->h*(NItr->Ena - NItr->V);
}

float EqsGnl::IHHIpv()
{
  Vx=NItr->V - VRESTING;
  Ax=An();
  Bx=Bn();
  NItr->n=IpvEul(NItr->n, &EqsGnl::dx);

  Ax=Am();
  Bx=Bm();
  NItr->m=IpvEul(NItr->m, &EqsGnl::dx);

  Ax=Ah();
  Bx=Bh();
  NItr->h=IpvEul(NItr->h, &EqsGnl::dx);

  return NItr->gk*pow(NItr->n,4.0f)*(NItr->Ek - NItr->V) + NItr->gna*pow(NItr->m,3.0f)*NItr->h*(NItr->Ena - NItr->V);
}

float EqsGnl::IHHEul()
{
  Vx=NItr->V - VRESTING;
  Ax=An();
  Bx=Bn();
  NItr->n += dt*(Ax*(1.0f-NItr->n)-Bx*NItr->n);

  Ax=Am();
  Bx=Bm();
  NItr->m += dt*(Ax*(1.0f-NItr->m)-Bx*NItr->m);

  Ax=Ah();
  Bx=Bh();
  NItr->h += dt*(Ax*(1.0f-NItr->h)-Bx*NItr->h);

  return NItr->gk*pow(NItr->n,4.0f)*(NItr->Ek - NItr->V) + NItr->gna*pow(NItr->m,3.0f)*NItr->h*(NItr->Ena - NItr->V);
}

//-------------long term plasticity
void EqsGnl::BackProp()
{
  if(NItr->FirLP)
  {
    if(NItr->CountLP==NItr->SpiLPDy)
    {
      NItr->SpiLP = true;
      NItr->CountLP=0;
      NItr->FirLP=false;
    }
    else
    {
      NItr->SpiLP = false;
      NItr->CountLP++;
    }
  }
  else
  {  NItr->SpiLP = false; }
}

void EqsGnl::ltpRK4()
{
  if(NItr->SpiLP)
  { NItr->LP += 1.0f; }
  else
  { NItr->LP = RK4Var2(NItr->LP,NItr->tLP,&EqsGnl::FFst); }
}

void EqsGnl::ltpIpv()
{
  if(NItr->SpiLP)
  { NItr->LP +=1.0f; }
  else
  { NItr->LP = IpvEulVar2(NItr->LP,NItr->tLP,&EqsGnl::FFst); }
}

void EqsGnl::ltpEul()
{
  if(NItr->SpiLP)
  { NItr->LP += 1.0f; }
  else
  { NItr->LP -= dt*(NItr->LP/NItr->tLP); }
}

//-------------short term plasticity

float EqsGnl::FCar(float Car)
{ return (CarBase*1e-6f- Car) / NItr->tCar ; }

float EqsGnl::FCap(float Cap)
{ return -Cap/ NItr->tCap ; }

float EqsGnl::FF(float F)
{ return  (NItr->Cap+NItr->Car)*(1.0f-F)* NItr->aF - F/NItr->tF ; }

float EqsGnl::FD(float D)
{ return (1.0f-D)/NItr->tD; }


void EqsGnl::stpRK4()
{
  if(NItr->SpiOut)
  {
    NItr->Car += NItr->aCar;
    NItr->Cap += NItr->aCap;
    NItr->F = RK4(NItr->F, &EqsGnl::FF);
    NItr->D *= ( 1.0f -  NItr->pv * NItr->F);
  }
  else
  {
    NItr->Car = RK4(NItr->Car, &EqsGnl::FCar);
    NItr->Cap = RK4(NItr->Cap, &EqsGnl::FCap);
    NItr->F = RK4(NItr->F, &EqsGnl::FF);
    NItr->D = RK4(NItr->D, &EqsGnl::FD);
  }
}

void EqsGnl::stpIpv()
{
  if(NItr->SpiOut)
  {
    NItr->Car += NItr->aCar;
    NItr->Cap += NItr->aCap;
    NItr->F = IpvEul(NItr->F, &EqsGnl::FF);
    NItr->D *= ( 1.0f -  NItr->pv * NItr->F);
  }
  else
  {
    NItr->Car = IpvEul(NItr->Car, &EqsGnl::FCar);
    NItr->Cap = IpvEul(NItr->Cap, &EqsGnl::FCap);
    NItr->F = IpvEul(NItr->F, &EqsGnl::FF);
    NItr->D = IpvEul(NItr->D, &EqsGnl::FD);
  }
}

void EqsGnl::stpEul()
{
  if(NItr->SpiOut)
  {
    NItr->Car += NItr->aCar;
    NItr->Cap += NItr->aCap;
    NItr->F += dt*( (NItr->Cap+NItr->Car)*(1.0f - NItr->F)*NItr->aF - NItr->F/NItr->tF);
    NItr->D *= ( 1.0f -  NItr->pv * NItr->F);
  }
  else
  {
    NItr->Car += dt*(CarBase*1e-6f - NItr->Car) / NItr->tCar;
    NItr->Cap -= dt*NItr->Cap/ NItr->tCap;
    NItr->F += dt*( (NItr->Cap+NItr->Car)*(1.0f - NItr->F)*NItr->aF - NItr->F/NItr->tF);
    NItr->D += dt*(1.0f - NItr->D)/NItr->tD;
  }
}


void EqsGnl::VmRK4()
{
  if( (NItr->ActCtr & LIF_ACT) != 0 && !NItr->Firing)// LIF model
  { NItr->V = RK4(NItr->V, &EqsGnl::MemPot); }
  else if((NItr->ActCtr & HH_ACT) != 0)
  {
    NItr->Isyn += IHHRK4();
    NItr->V = RK4(NItr->V, &EqsGnl::MemPot);
  }

  if(NItr->Firing)
  { firing(); }
  else
  {
    if(NItr->V > NItr->Vth)
    {
      NItr->Firing=true;
      firing();

      NItr->FirLP=true;
    }
  }

//synapse plasticity
  if((NItr->ActCtr & LTP_ACT) != 0)
  {
    BackProp();
    ltpRK4();
  }

  for(unsigned int i=0;i<thsize;i++)
  {
    *(NItr->gAMPA+i)=0.0f;
    *(NItr->gGABA+i)=0.0f;
    *(NItr->gNMDA+i)=0.0f;
    *(NItr->gACH+i)=0.0f;
    *(NItr->gGCL+i)=0.0f;
  }

  NItr->Isyn=0.0f;
}

void EqsGnl::VmIpv()
{
  if( (NItr->ActCtr & LIF_ACT) != 0 && !NItr->Firing)// LIF model
  { NItr->V = IpvEul(NItr->V, &EqsGnl::MemPot); }
  else if((NItr->ActCtr & HH_ACT) != 0)
  {
    NItr->Isyn += IHHIpv();
    NItr->V = IpvEul(NItr->V, &EqsGnl::MemPot);
  }

  if(NItr->Firing)
  { firing(); }
  else
  {
    if(NItr->V > NItr->Vth)
    {
      NItr->Firing=true;
      firing();

      NItr->FirLP=true;
    }
  }

//synapse plasticity
  if((NItr->ActCtr & LTP_ACT) != 0)
  {
    BackProp();
    ltpIpv();
  }

  for(unsigned int i=0;i<thsize;i++)
  {
    *(NItr->gAMPA+i)=0.0f;
    *(NItr->gGABA+i)=0.0f;
    *(NItr->gNMDA+i)=0.0f;
    *(NItr->gACH+i)=0.0f;
    *(NItr->gGCL+i)=0.0f;
  }

  NItr->Isyn=0.0f;
}

void EqsGnl::VmEul()
{
  if( (NItr->ActCtr & LIF_ACT) != 0 && !NItr->Firing)// LIF model
  { NItr->V += dt * MemPot(NItr->V); }
  else if((NItr->ActCtr & HH_ACT) != 0)
  {
    NItr->Isyn += IHHEul();
    NItr->V += dt * MemPot(NItr->V);
  }

  if(NItr->Firing)
  { firing(); }
  else
  {
    if(NItr->V > NItr->Vth)
    {
      NItr->Firing=true;
      firing();

      NItr->FirLP=true;
    }
  }

//synapse plasticity
  if((NItr->ActCtr & LTP_ACT) != 0)
  {
    BackProp();
    ltpEul();
  }

  for(unsigned int i=0;i<thsize;i++)
  {
    *(NItr->gAMPA+i)=0.0f;
    *(NItr->gGABA+i)=0.0f;
    *(NItr->gNMDA+i)=0.0f;
    *(NItr->gACH+i)=0.0f;
    *(NItr->gGCL+i)=0.0f;
  }

  NItr->Isyn=0.0f;
}


void EqsGnl::Spike(unsigned int thid,unsigned int NID)
{
  float tst=1.0f;
  float NSpi=0.0f;
  float Sti;
  NItr=(NItr_begin+NID);


  if( (NItr->ActCtr & STP_ACT) != 0)
  { Sti=NItr->D*NItr->F; }
  else
  { Sti=1.0f; }

  if(NItr->SpiOut)
    NSpi=1.0f;
  else
    NSpi=0.0f;



  switch(SolType)
  {
    case ACCURATE:
//fast synapse
//SynFastSize1-----------------------------------
      if(NItr->SynFastSize1 != 0)
      {
          for(unsigned int i=0;i<NItr->SynFastSize1;i++)
          {
            switch(*(NItr->SFTyp1+i))
            {
              case KW_AMPA:
                tst=(NItr_begin + (NItr->SynFst1+i)->PosID)->ParAMPA.tst;  // gating variable
              break;

              case KW_GABA:
                tst=(NItr_begin + (NItr->SynFst1+i)->PosID)->ParGABA.tst;  // gating variable
              break;

              case KW_ACH:
                tst=(NItr_begin + (NItr->SynFst1+i)->PosID)->ParACH.tst;  // gating variable
              break;

              case KW_GCL:
                tst=(NItr_begin + (NItr->SynFst1+i)->PosID)->ParGCL.tst;  // gating variable
              break;

              default:;
            }

            (NItr->SynFst1+i)->stFast = NSpi*(Sti + (NItr->SynFst1+i)->stFast) + (1.0f-NSpi)*RK4Var2( (NItr->SynFst1+i)->stFast, tst, &EqsGnl::FFst);
          }
      }

//SynFastSize-----------------------------------
      if(NItr->SynFastSize != 0)
      {
          for(unsigned int i=0;i<NItr->SynFastSize;i++)
          {
            switch(*(NItr->SFTyp+i))
            {
              case KW_AMPA:
                tst=(NItr_begin + *((NItr->SynFst+i)->PosID))->ParAMPA.tst;  // gating variable
              break;

              case KW_GABA:
                tst=(NItr_begin + *((NItr->SynFst+i)->PosID))->ParGABA.tst;  // gating variable
              break;

              case KW_ACH:
                tst=(NItr_begin + *((NItr->SynFst+i)->PosID))->ParACH.tst;  // gating variable
              break;

              case KW_GCL:
                tst=(NItr_begin + *((NItr->SynFst+i)->PosID))->ParGCL.tst;  // gating variable
              break;

              default:;
            }

            (NItr->SynFst+i)->stFast = NSpi*(Sti + (NItr->SynFst+i)->stFast) + (1.0f-NSpi)*RK4Var2( (NItr->SynFst+i)->stFast, tst, &EqsGnl::FFst);
          }
      }

//slow synapse
//SynSlowSize1-----------------------------------
      if(NItr->SynSlowSize1 != 0)
      {
          for(unsigned int i=0;i<NItr->SynSlowSize1;i++)
          {
            (NItr->SynSlo1+i)->xtSlow = NSpi*(Sti + (NItr->SynSlo1+i)->xtSlow) + (1.0f-NSpi)*RK4Var2( (NItr->SynSlo1+i)->xtSlow, (NItr_begin + (NItr->SynSlo1+i)->PosID)->ParNMDA.txt, &EqsGnl::FFst);

            (NItr->SynSlo1+i)->stSlow = RK4Var4( (NItr->SynSlo1+i)->stSlow, (NItr_begin + (NItr->SynSlo1+i)->PosID)->ParNMDA.tst, (NItr->SynSlo1+i)->xtSlow, (NItr_begin + (NItr->SynSlo1+i)->PosID)->ParNMDA.as, &EqsGnl::Fst);
          }
      }
//SynSlowSize-----------------------------------
      if(NItr->SynSlowSize != 0)
      {
          for(unsigned int i=0;i<NItr->SynSlowSize;i++)
          {
            (NItr->SynSlo+i)->xtSlow = NSpi*(Sti + (NItr->SynSlo+i)->xtSlow) + (1.0f-NSpi)*RK4Var2( (NItr->SynSlo+i)->xtSlow, (NItr_begin + *(NItr->SynSlo+i)->PosID)->ParNMDA.txt, &EqsGnl::FFst);

            (NItr->SynSlo+i)->stSlow = RK4Var4( (NItr->SynSlo+i)->stSlow, (NItr_begin + *(NItr->SynSlo+i)->PosID)->ParNMDA.tst, (NItr->SynSlo+i)->xtSlow, (NItr_begin + *(NItr->SynSlo+i)->PosID)->ParNMDA.as, &EqsGnl::Fst);
          }
      }

//glutamate synapse
//SynGlutSize1
      if(NItr->SynGlutSize1 != 0)
      {
          for(unsigned int i=0;i<NItr->SynGlutSize1;i++)
          {
            switch((NItr->SGTyp1+i)->A)
            {
              case KW_AMPA:
                tst=(NItr_begin + (NItr->SynGlu1+i)->PosID)->ParAMPA.tst;  // gating variable
              break;

              case KW_GABA:
                tst=(NItr_begin + (NItr->SynGlu1+i)->PosID)->ParGABA.tst;  // gating variable
              break;

              case KW_ACH:
                tst=(NItr_begin + (NItr->SynGlu1+i)->PosID)->ParACH.tst;  // gating variable
              break;

              case KW_GCL:
                tst=(NItr_begin + (NItr->SynGlu1+i)->PosID)->ParGCL.tst;  // gating variable
              break;

              default:;
            }

            (NItr->SynGlu1+i)->stFast = NSpi*(Sti + (NItr->SynGlu1+i)->stFast) + (1.0f-NSpi)*RK4Var2( (NItr->SynGlu1+i)->stFast, tst, &EqsGnl::FFst);
            (NItr->SynGlu1+i)->xtSlow = NSpi*(Sti + (NItr->SynGlu1+i)->xtSlow) + (1.0f-NSpi)*RK4Var2( (NItr->SynGlu1+i)->xtSlow, (NItr_begin + (NItr->SynGlu1+i)->PosID)->ParNMDA.txt, &EqsGnl::FFst);

            (NItr->SynGlu1+i)->stSlow = RK4Var4( (NItr->SynGlu1+i)->stSlow, (NItr_begin + (NItr->SynGlu1+i)->PosID)->ParNMDA.tst, (NItr->SynGlu1+i)->xtSlow, (NItr_begin + (NItr->SynGlu1+i)->PosID)->ParNMDA.as, &EqsGnl::Fst);
          }
      }

//SynGlutSize
      if(NItr->SynGlutSize != 0)
      {
          for(unsigned int i=0;i<NItr->SynGlutSize;i++)
          {
            switch((NItr->SGTyp+i)->A)
            {
              case KW_AMPA:
                tst=(NItr_begin + *((NItr->SynGlu+i)->PosID))->ParAMPA.tst;  // gating variable
              break;

              case KW_GABA:
                tst=(NItr_begin + *((NItr->SynGlu+i)->PosID))->ParGABA.tst;  // gating variable
              break;

              case KW_ACH:
                tst=(NItr_begin + *((NItr->SynGlu+i)->PosID))->ParACH.tst;  // gating variable
              break;

              case KW_GCL:
                tst=(NItr_begin + *((NItr->SynGlu+i)->PosID))->ParGCL.tst;  // gating variable
              break;

              default:;
            }

            (NItr->SynGlu+i)->stFast = NSpi*(Sti + (NItr->SynGlu+i)->stFast) + (1.0f-NSpi)*RK4Var2( (NItr->SynGlu+i)->stFast, tst, &EqsGnl::FFst);;
            (NItr->SynGlu+i)->xtSlow = NSpi*(Sti + (NItr->SynGlu+i)->xtSlow) + (1.0f-NSpi)*RK4Var2( (NItr->SynGlu+i)->xtSlow, (NItr_begin + *(NItr->SynGlu+i)->PosID)->ParNMDA.txt, &EqsGnl::FFst);

            (NItr->SynGlu+i)->stSlow = RK4Var4( (NItr->SynGlu+i)->stSlow, (NItr_begin + *(NItr->SynGlu+i)->PosID)->ParNMDA.tst, (NItr->SynGlu+i)->xtSlow, (NItr_begin + *(NItr->SynGlu+i)->PosID)->ParNMDA.as, &EqsGnl::Fst);
          }
      }
    break;

    case ROUGH:
//fast synapse
//SynFastSize1-----------------------------------
      if(NItr->SynFastSize1 != 0)
      {
          for(unsigned int i=0;i<NItr->SynFastSize1;i++)
          {
            switch(*(NItr->SFTyp1+i))
            {
              case KW_AMPA:
                tst=(NItr_begin + (NItr->SynFst1+i)->PosID)->ParAMPA.tst;  // gating variable
              break;

              case KW_GABA:
                tst=(NItr_begin + (NItr->SynFst1+i)->PosID)->ParGABA.tst;  // gating variable
              break;

              case KW_ACH:
                tst=(NItr_begin + (NItr->SynFst1+i)->PosID)->ParACH.tst;  // gating variable
              break;

              case KW_GCL:
                tst=(NItr_begin + (NItr->SynFst1+i)->PosID)->ParGCL.tst;  // gating variable
              break;

              default:;
            }

            (NItr->SynFst1+i)->stFast = NSpi*(Sti + (NItr->SynFst1+i)->stFast) + (1.0f-NSpi)*((NItr->SynFst1+i)->stFast - dt*( (NItr->SynFst1+i)->stFast/tst ));

          }
      }

//SynFastSize-----------------------------------
      if(NItr->SynFastSize != 0)
      {
          for(unsigned int i=0;i<NItr->SynFastSize;i++)
          {
            switch(*(NItr->SFTyp+i))
            {
              case KW_AMPA:
                tst=(NItr_begin + *((NItr->SynFst+i)->PosID))->ParAMPA.tst;  // gating variable
              break;

              case KW_GABA:
                tst=(NItr_begin + *((NItr->SynFst+i)->PosID))->ParGABA.tst;  // gating variable
              break;

              case KW_ACH:
                tst=(NItr_begin + *((NItr->SynFst+i)->PosID))->ParACH.tst;  // gating variable
              break;

              case KW_GCL:
                tst=(NItr_begin + *((NItr->SynFst+i)->PosID))->ParGCL.tst;  // gating variable
              break;

              default:;
            }

            (NItr->SynFst+i)->stFast = NSpi*(Sti + (NItr->SynFst+i)->stFast) + (1.0f-NSpi)*((NItr->SynFst+i)->stFast - dt*( (NItr->SynFst+i)->stFast/tst));
          }
      }

//slow synapse
//SynSlowSize1-----------------------------------
      if(NItr->SynSlowSize1 != 0)
      {
          for(unsigned int i=0;i<NItr->SynSlowSize1;i++)
          {
            (NItr->SynSlo1+i)->xtSlow = NSpi*(Sti + (NItr->SynSlo1+i)->xtSlow) + (1.0f-NSpi)*((NItr->SynSlo1+i)->xtSlow - dt*( (NItr->SynSlo1+i)->xtSlow/(NItr_begin + (NItr->SynSlo1+i)->PosID)->ParNMDA.txt ));

            (NItr->SynSlo1+i)->stSlow += dt*( (NItr_begin + (NItr->SynSlo1+i)->PosID)->ParNMDA.as*((NItr->SynSlo1+i)->xtSlow)*(1.0f - (NItr->SynSlo1+i)->stSlow) - (NItr->SynSlo1+i)->stSlow/(NItr_begin + (NItr->SynSlo1+i)->PosID)->ParNMDA.tst);
          }
      }
//SynSlowSize-----------------------------------
      if(NItr->SynSlowSize != 0)
      {
          for(unsigned int i=0;i<NItr->SynSlowSize;i++)
          {
            (NItr->SynSlo+i)->xtSlow = NSpi*(Sti + (NItr->SynSlo+i)->xtSlow) + (1.0f-NSpi)*((NItr->SynSlo+i)->xtSlow - dt*( (NItr->SynSlo+i)->xtSlow/(NItr_begin + *(NItr->SynSlo+i)->PosID)->ParNMDA.txt ));

            (NItr->SynSlo+i)->stSlow += dt*( (NItr_begin + *(NItr->SynSlo+i)->PosID)->ParNMDA.as*(NItr->SynSlo+i)->xtSlow*(1.0f - (NItr->SynSlo+i)->stSlow) - (NItr->SynSlo+i)->stSlow/(NItr_begin + *(NItr->SynSlo+i)->PosID)->ParNMDA.tst);
          }
      }

//glutamate synapse
//SynGlutSize1
      if(NItr->SynGlutSize1 != 0)
      {
          for(unsigned int i=0;i<NItr->SynGlutSize1;i++)
          {
            switch((NItr->SGTyp1+i)->A)
            {
              case KW_AMPA:
                tst=(NItr_begin + (NItr->SynGlu1+i)->PosID)->ParAMPA.tst;  // gating variable
              break;

              case KW_GABA:
                tst=(NItr_begin + (NItr->SynGlu1+i)->PosID)->ParGABA.tst;  // gating variable
              break;

              case KW_ACH:
                tst=(NItr_begin + (NItr->SynGlu1+i)->PosID)->ParACH.tst;  // gating variable
              break;

              case KW_GCL:
                tst=(NItr_begin + (NItr->SynGlu1+i)->PosID)->ParGCL.tst;  // gating variable
              break;

              default:;
            }

            (NItr->SynGlu1+i)->stFast = NSpi*(Sti + (NItr->SynGlu1+i)->stFast) + (1.0f-NSpi)*((NItr->SynGlu1+i)->stFast-dt*( (NItr->SynGlu1+i)->stFast/tst ));
            (NItr->SynGlu1+i)->xtSlow = NSpi*(Sti + (NItr->SynGlu1+i)->xtSlow) + (1.0f-NSpi)*((NItr->SynGlu1+i)->xtSlow-dt*( (NItr->SynGlu1+i)->xtSlow/(NItr_begin + (NItr->SynGlu1+i)->PosID)->ParNMDA.txt ));

            (NItr->SynGlu1+i)->stSlow += dt*(  (NItr_begin + (NItr->SynGlu1+i)->PosID)->ParNMDA.as*((NItr->SynGlu1+i)->xtSlow)*(1.0f - (NItr->SynGlu1+i)->stSlow) - (NItr->SynGlu1+i)->stSlow/(NItr_begin + (NItr->SynGlu1+i)->PosID)->ParNMDA.tst);
          }
      }

//SynGlutSize
      if(NItr->SynGlutSize != 0)
      {
          for(unsigned int i=0;i<NItr->SynGlutSize;i++)
          {
            switch((NItr->SGTyp+i)->A)
            {
              case KW_AMPA:
                tst=(NItr_begin + *((NItr->SynGlu+i)->PosID))->ParAMPA.tst;  // gating variable
              break;

              case KW_GABA:
                tst=(NItr_begin + *((NItr->SynGlu+i)->PosID))->ParGABA.tst;  // gating variable
              break;

              case KW_ACH:
                tst=(NItr_begin + *((NItr->SynGlu+i)->PosID))->ParACH.tst;  // gating variable
              break;

              case KW_GCL:
                tst=(NItr_begin + *((NItr->SynGlu+i)->PosID))->ParGCL.tst;  // gating variable
              break;

              default:;
            }

            (NItr->SynGlu+i)->stFast = NSpi*(Sti + (NItr->SynGlu+i)->stFast) + (1.0f-NSpi)*((NItr->SynGlu+i)->stFast - dt*( (NItr->SynGlu+i)->stFast/tst ));
            (NItr->SynGlu+i)->xtSlow = NSpi*(Sti + (NItr->SynGlu+i)->xtSlow) + (1.0f-NSpi)*((NItr->SynGlu+i)->xtSlow - dt*( (NItr->SynGlu+i)->xtSlow/(NItr_begin + *((NItr->SynGlu+i)->PosID))->ParNMDA.txt  ));

            (NItr->SynGlu+i)->stSlow += dt*(  (NItr_begin + *((NItr->SynGlu+i)->PosID))->ParNMDA.as*((NItr->SynGlu+i)->xtSlow)*(1.0f - (NItr->SynGlu+i)->stSlow) - (NItr->SynGlu+i)->stSlow/(NItr_begin + *((NItr->SynGlu+i)->PosID))->ParNMDA.tst);
          }
      }
    break;

    case MODERATE:
    default:
//fast synapse
//SynFastSize1-----------------------------------
      if(NItr->SynFastSize1 != 0)
      {
          for(unsigned int i=0;i<NItr->SynFastSize1;i++)
          {
            switch(*(NItr->SFTyp1+i))
            {
              case KW_AMPA:
                tst=(NItr_begin + (NItr->SynFst1+i)->PosID)->ParAMPA.tst;  // gating variable
              break;

              case KW_GABA:
                tst=(NItr_begin + (NItr->SynFst1+i)->PosID)->ParGABA.tst;  // gating variable
              break;

              case KW_ACH:
                tst=(NItr_begin + (NItr->SynFst1+i)->PosID)->ParACH.tst;  // gating variable
              break;

              case KW_GCL:
                tst=(NItr_begin + (NItr->SynFst1+i)->PosID)->ParGCL.tst;  // gating variable
              break;

              default:;
            }

            (NItr->SynFst1+i)->stFast = NSpi*(Sti + (NItr->SynFst1+i)->stFast) + (1.0f-NSpi)*IpvEulVar2( (NItr->SynFst1+i)->stFast, tst, &EqsGnl::FFst);
          }
      }

//SynFastSize-----------------------------------
      if(NItr->SynFastSize != 0)
      {
          for(unsigned int i=0;i<NItr->SynFastSize;i++)
          {
            switch(*(NItr->SFTyp+i))
            {
              case KW_AMPA:
                tst=(NItr_begin + *((NItr->SynFst+i)->PosID))->ParAMPA.tst;  // gating variable
              break;

              case KW_GABA:
                tst=(NItr_begin + *((NItr->SynFst+i)->PosID))->ParGABA.tst;  // gating variable
              break;

              case KW_ACH:
                tst=(NItr_begin + *((NItr->SynFst+i)->PosID))->ParACH.tst;  // gating variable
              break;

              case KW_GCL:
                tst=(NItr_begin + *((NItr->SynFst+i)->PosID))->ParGCL.tst;  // gating variable
              break;

              default:;
            }

            (NItr->SynFst+i)->stFast = NSpi*(Sti + (NItr->SynFst+i)->stFast) + (1.0f-NSpi)*IpvEulVar2( (NItr->SynFst+i)->stFast, tst, &EqsGnl::FFst);
          }
      }

//slow synapse
//SynSlowSize1-----------------------------------
      if(NItr->SynSlowSize1 != 0)
      {
          for(unsigned int i=0;i<NItr->SynSlowSize1;i++)
          {
            (NItr->SynSlo1+i)->xtSlow = NSpi*(Sti + (NItr->SynSlo1+i)->xtSlow) + (1.0f-NSpi)*IpvEulVar2( (NItr->SynSlo1+i)->xtSlow, (NItr_begin + (NItr->SynSlo1+i)->PosID)->ParNMDA.txt, &EqsGnl::FFst);

            (NItr->SynSlo1+i)->stSlow = IpvEulVar4( (NItr->SynSlo1+i)->stSlow, (NItr_begin + (NItr->SynSlo1+i)->PosID)->ParNMDA.tst, (NItr->SynSlo1+i)->xtSlow, (NItr_begin + (NItr->SynSlo1+i)->PosID)->ParNMDA.as, &EqsGnl::Fst);
          }
      }
//SynSlowSize-----------------------------------
      if(NItr->SynSlowSize != 0)
      {
          for(unsigned int i=0;i<NItr->SynSlowSize;i++)
          {
            (NItr->SynSlo+i)->xtSlow = NSpi*(Sti + (NItr->SynSlo+i)->xtSlow) + (1.0f-NSpi)*IpvEulVar2( (NItr->SynSlo+i)->xtSlow, (NItr_begin + *(NItr->SynSlo+i)->PosID)->ParNMDA.txt, &EqsGnl::FFst);

            (NItr->SynSlo+i)->stSlow = IpvEulVar4( (NItr->SynSlo+i)->stSlow, (NItr_begin + *(NItr->SynSlo+i)->PosID)->ParNMDA.tst, (NItr->SynSlo+i)->xtSlow, (NItr_begin + *(NItr->SynSlo+i)->PosID)->ParNMDA.as, &EqsGnl::Fst);
          }
      }

//glutamate synapse
//SynGlutSize1
      if(NItr->SynGlutSize1 != 0)
      {
          for(unsigned int i=0;i<NItr->SynGlutSize1;i++)
          {
            switch((NItr->SGTyp1+i)->A)
            {
              case KW_AMPA:
                tst=(NItr_begin + (NItr->SynGlu1+i)->PosID)->ParAMPA.tst;  // gating variable
              break;

              case KW_GABA:
                tst=(NItr_begin + (NItr->SynGlu1+i)->PosID)->ParGABA.tst;  // gating variable
              break;

              case KW_ACH:
                tst=(NItr_begin + (NItr->SynGlu1+i)->PosID)->ParACH.tst;  // gating variable
              break;

              case KW_GCL:
                tst=(NItr_begin + (NItr->SynGlu1+i)->PosID)->ParGCL.tst;  // gating variable
              break;

              default:;
            }

            (NItr->SynGlu1+i)->stFast = NSpi*(Sti + (NItr->SynGlu1+i)->stFast) + (1.0f-NSpi)*IpvEulVar2( (NItr->SynGlu1+i)->stFast, tst, &EqsGnl::FFst);
            (NItr->SynGlu1+i)->xtSlow = NSpi*(Sti + (NItr->SynGlu1+i)->xtSlow) + (1.0f-NSpi)*IpvEulVar2( (NItr->SynGlu1+i)->xtSlow, (NItr_begin + (NItr->SynGlu1+i)->PosID)->ParNMDA.txt, &EqsGnl::FFst);

            (NItr->SynGlu1+i)->stSlow = IpvEulVar4( (NItr->SynGlu1+i)->stSlow, (NItr_begin + (NItr->SynGlu1+i)->PosID)->ParNMDA.tst, (NItr->SynGlu1+i)->xtSlow, (NItr_begin + (NItr->SynGlu1+i)->PosID)->ParNMDA.as, &EqsGnl::Fst);
          }
      }

//SynGlutSize
      if(NItr->SynGlutSize != 0)
      {
          for(unsigned int i=0;i<NItr->SynGlutSize;i++)
          {
            switch((NItr->SGTyp+i)->A)
            {
              case KW_AMPA:
                tst=(NItr_begin + *((NItr->SynGlu+i)->PosID))->ParAMPA.tst;  // gating variable
              break;

              case KW_GABA:
                tst=(NItr_begin + *((NItr->SynGlu+i)->PosID))->ParGABA.tst;  // gating variable
              break;

              case KW_ACH:
                tst=(NItr_begin + *((NItr->SynGlu+i)->PosID))->ParACH.tst;  // gating variable
              break;

              case KW_GCL:
                tst=(NItr_begin + *((NItr->SynGlu+i)->PosID))->ParGCL.tst;  // gating variable
              break;

              default:;
            }

            (NItr->SynGlu+i)->stFast = NSpi*(Sti + (NItr->SynGlu+i)->stFast) + (1.0f-NSpi)*IpvEulVar2( (NItr->SynGlu+i)->stFast, tst, &EqsGnl::FFst);
            (NItr->SynGlu+i)->xtSlow = NSpi*(Sti + (NItr->SynGlu+i)->xtSlow) + (1.0f-NSpi)*IpvEulVar2( (NItr->SynGlu+i)->xtSlow, (NItr_begin + *(NItr->SynGlu+i)->PosID)->ParNMDA.txt, &EqsGnl::FFst);

            (NItr->SynGlu+i)->stSlow = IpvEulVar4( (NItr->SynGlu+i)->stSlow, (NItr_begin + *(NItr->SynGlu+i)->PosID)->ParNMDA.tst, (NItr->SynGlu+i)->xtSlow, (NItr_begin + *(NItr->SynGlu+i)->PosID)->ParNMDA.as, &EqsGnl::Fst);
          }
      }
    break;
  }


//spiking prossing-----------------------------------
//fast synapse
//SynFastSize1-----------------------------------
  float gVar=0.0f;
  float gVar2=0.0f;

  if(NItr->SynFastSize1 != 0)
  {
    for(unsigned int i=0;i<NItr->SynFastSize1;i++)
    {
      gVar = (NItr->SynFst1+i)->gFast * (NItr->SynFst1+i)->stFast;
      switch(*(NItr->SFTyp1+i))
      {
        case KW_AMPA:
        *((NItr_begin + (NItr->SynFst1+i)->PosID)->gAMPA+thid) += gVar;
        break;

        case KW_GABA:
        *((NItr_begin + (NItr->SynFst1+i)->PosID)->gGABA+thid) += gVar;
        break;

        case KW_ACH:
        *((NItr_begin + (NItr->SynFst1+i)->PosID)->gACH+thid) += gVar;
        break;

        case KW_GCL:
        *((NItr_begin + (NItr->SynFst1+i)->PosID)->gGCL+thid) += gVar;
        break;

        default:;
      }
    }
  }

//SynFastSize-----------------------------------
  if(NItr->SynFastSize != 0)
  {
    for(unsigned int i=0;i<NItr->SynFastSize;i++)
    {
      gVar = *((NItr->SynFst+i)->gFast) * (NItr->SynFst+i)->stFast;
      for(unsigned int j=0;j<(NItr->SynFst+i)->PosIDNum;j++) // post link population neurons
      {
        switch(*(NItr->SFTyp+i))
        {
          case KW_AMPA:
          *((NItr_begin + *((NItr->SynFst+i)->PosID+j))->gAMPA+thid) += gVar;
          break;

          case KW_GABA:
          *((NItr_begin + *((NItr->SynFst+i)->PosID+j))->gGABA+thid) += gVar;
          break;

          case KW_ACH:
          *((NItr_begin + *((NItr->SynFst+i)->PosID+j))->gACH+thid) += gVar;
          break;

          case KW_GCL:
          *((NItr_begin + *((NItr->SynFst+i)->PosID+j))->gGCL+thid) += gVar;
          break;

          default:;
        }
      }
    }
  }

//slow synapse
//SynSlowSize1-----------------------------------
  if(NItr->SynSlowSize1 != 0)
  {
    for(unsigned int i=0;i<NItr->SynSlowSize1;i++)
    {
      gVar = (NItr->SynSlo1+i)->gSlow * (NItr->SynSlo1+i)->stSlow;
      *((NItr_begin + (NItr->SynSlo1+i)->PosID)->gNMDA+thid) += gVar;
    }
  }

//SynSlowSize-----------------------------------
  if(NItr->SynSlowSize != 0)
  {
    for(unsigned int i=0;i<NItr->SynSlowSize;i++)
    {
      gVar= *((NItr->SynSlo+i)->gSlow) * (NItr->SynSlo+i)->stSlow;
      for(unsigned int j=0;j<(NItr->SynSlo+i)->PosIDNum;j++) // post link population neurons
      { *((NItr_begin + *((NItr->SynSlo+i)->PosID+j))->gNMDA+thid) += gVar; }
    }
  }

//glutamate synapse
//SynGlutSize1
  if(NItr->SynGlutSize1 != 0)
  {
    for(unsigned int i=0;i<NItr->SynGlutSize1;i++)
    {
      gVar = (NItr->SynGlu1+i)->gFast * (NItr->SynGlu1+i)->stFast;

      if((NItr->SGTyp1+i)->A==KW_AMPA)
      *((NItr_begin + (NItr->SynGlu1+i)->PosID)->gAMPA+thid) += gVar;
      else if((NItr->SGTyp1+i)->A==KW_GABA)
      *((NItr_begin + (NItr->SynGlu1+i)->PosID)->gGABA+thid) += gVar;
      else if((NItr->SGTyp1+i)->A==KW_ACH)
      *((NItr_begin + (NItr->SynGlu1+i)->PosID)->gACH+thid) += gVar;
      else if((NItr->SGTyp1+i)->A==KW_GCL)
      *((NItr_begin + (NItr->SynGlu1+i)->PosID)->gGCL+thid) += gVar;

      gVar = (NItr->SynGlu1+i)->gSlow * (NItr->SynGlu1+i)->stSlow;
      *((NItr_begin + (NItr->SynGlu1+i)->PosID)->gNMDA+thid) += gVar;
    }
  }

//SynGlutSize
  if(NItr->SynGlutSize != 0)
  {
      for(unsigned int i=0;i<NItr->SynGlutSize;i++)
      {
        gVar=*((NItr->SynGlu+i)->gFast) * (NItr->SynGlu+i)->stFast;
        gVar2=*((NItr->SynGlu+i)->gSlow) * (NItr->SynGlu+i)->stSlow;
        for(unsigned int j=0;j<(NItr->SynGlu+i)->PosIDNum;j++) // post link population neurons
        {
          if((NItr->SGTyp+i)->A==KW_AMPA)
          *((NItr_begin + *((NItr->SynGlu+i)->PosID+j))->gAMPA+thid) += gVar;
          else if((NItr->SGTyp+i)->A==KW_GABA)
          *((NItr_begin + *((NItr->SynGlu+i)->PosID+j))->gGABA+thid) += gVar;
          else if((NItr->SGTyp+i)->A==KW_ACH)
          *((NItr_begin + *((NItr->SynGlu+i)->PosID+j))->gACH+thid) += gVar;
          else if((NItr->SGTyp+i)->A==KW_GCL)
          *((NItr_begin + *((NItr->SynGlu+i)->PosID+j))->gGCL+thid) += gVar;

          *((NItr_begin + *((NItr->SynGlu+i)->PosID+j))->gNMDA+thid) += gVar2;
        }
      }
  }

//gap junction synapse
//SynGapSize1
  if(NItr->SynGapSize1 != 0)
  {
      for(unsigned int i=0;i<NItr->SynGapSize1;i++) // pre-neuron synapses
      {
        (NItr_begin + (NItr->SynCnd1+i)->PosID)->Isyn +=  (NItr->SynCnd1+i)->gC * ( NItr->V - (NItr_begin + (NItr->SynCnd1+i)->PosID)->V );
        NItr->Isyn -= (NItr->SynCnd1+i)->gC * ( NItr->V - (NItr_begin + (NItr->SynCnd1+i)->PosID)->V );
      }
  }

//SynGapSize
  if(NItr->SynGapSize != 0)
  {
      for(unsigned int i=0;i<NItr->SynGapSize;i++) // pre-neuron synapses
      {
        gVar = *((NItr->SynCnd+i)->gC);
        for(unsigned int j=0;j<(NItr->SynCnd+i)->PosIDNum;j++) // post link population neurons
        {
          (NItr_begin + *((NItr->SynCnd+i)->PosID+j))->Isyn +=  gVar * ( NItr->V - (NItr_begin + *((NItr->SynCnd+i)->PosID+j))->V );
          NItr->Isyn -= *((NItr->SynCnd+i)->gC+j) * ( NItr->V - (NItr_begin + *((NItr->SynCnd+i)->PosID+j))->V );
        }
      }
  }

//STP update
  if((NItr->ActCtr & STP_ACT) != 0)
  {
    switch(SolType)
    {
      case ACCURATE:
      stpRK4();
      break;

      case ROUGH:
      stpIpv();
      break;

      case MODERATE:
      default:
      stpEul();
      break;
    }
  }

//LTP update
  if((NItr->ActCtr & LTP_ACT) != 0)
  {
//SynFastSize1-----------------------------------
    if(NItr->SynFastSize1 != 0)
    {
      if(NItr->SpiOut) //NegISI
      {
        if(NItr->NegISI>0)  //facilitation
        {
          for(unsigned int i=0;i<NItr->SynFastSize1;i++)
          { (NItr->SynFst1+i)->gFast += (NItr_begin + (NItr->SynFst1+i)->PosID)->LP * NItr->NegISI; } //additive rule
        }
        else //depression
        {
          for(unsigned int i=0;i<NItr->SynFastSize1;i++)
          { (NItr->SynFst1+i)->gFast *= (1.0f + (NItr_begin + (NItr->SynFst1+i)->PosID)->LP * NItr->NegISI); } //multiplicative rule
        }
      }

      if(NItr->PosISI>0)  //PosISI, facilitation
      {
        for(unsigned int i=0;i<NItr->SynFastSize1;i++)
        {
          if((NItr_begin + (NItr->SynFst1+i)->PosID)->SpiOut)
          { (NItr->SynFst1+i)->gFast += NItr->LP * NItr->PosISI; } //additive rule
        }
      }
      else //depression
      {
        for(unsigned int i=0;i<NItr->SynFastSize1;i++)
        {
          if((NItr_begin + (NItr->SynFst1+i)->PosID)->SpiOut)
          { (NItr->SynFst1+i)->gFast *= (1.0f + NItr->LP * NItr->PosISI); } //multiplicative rule
        }
      }
    }
//SynFastSize-----------------------------------
    if(NItr->SynFastSize != 0)
    {
      if(NItr->SpiOut)
      {
        if(NItr->NegISI>0) //facilitation
        {
          for(unsigned int i=0;i<NItr->SynFastSize;i++)
          {
            for(unsigned int j=0;j<(NItr->SynFst+i)->PosIDNum;j++) // post link population neurons
            { *((NItr->SynFst+i)->gFast+j) += (NItr_begin + *((NItr->SynFst+i)->PosID+j))->LP * NItr->NegISI; } //additive rule
          }
        }
        else  //depression
        {
          for(unsigned int i=0;i<NItr->SynFastSize;i++)
          {
            for(unsigned int j=0;j<(NItr->SynFst+i)->PosIDNum;j++) // post link population neurons
            { *((NItr->SynFst+i)->gFast+j) *= (1.0f + (NItr_begin + *((NItr->SynFst+i)->PosID+j))->LP * NItr->NegISI); } //multiplicative rule
          }
        }
      }

      if(NItr->PosISI>0) //facilitation
      {
        for(unsigned int i=0;i<NItr->SynFastSize;i++)
        {
          for(unsigned int j=0;j<(NItr->SynFst+i)->PosIDNum;j++) // post link population neurons
          {
            if((NItr_begin + *((NItr->SynFst+i)->PosID+j))->SpiOut)
            { *((NItr->SynFst+i)->gFast+j) += NItr->LP * NItr->PosISI; }//additive rule
          }
        }
      }
      else //depression
      {
        for(unsigned int i=0;i<NItr->SynFastSize;i++)
        {
          for(unsigned int j=0;j<(NItr->SynFst+i)->PosIDNum;j++) // post link population neurons
          {
            if((NItr_begin + *((NItr->SynFst+i)->PosID+j))->SpiOut)
            { *((NItr->SynFst+i)->gFast+j) *= (1.0f+NItr->LP * NItr->PosISI); }//multiplicative rule
          }
        }
      }
    }

//SynSlowSize1-----------------------------------
    if(NItr->SynSlowSize1 != 0)
    {
      if(NItr->SpiOut) //NegISI
      {
        if(NItr->NegISI>0)  //facilitation
        {
          for(unsigned int i=0;i<NItr->SynSlowSize1;i++)
          { (NItr->SynSlo1+i)->gSlow += (NItr_begin + (NItr->SynSlo1+i)->PosID)->LP * NItr->NegISI; } //additive rule
        }
        else //depression
        {
          for(unsigned int i=0;i<NItr->SynSlowSize1;i++)
          { (NItr->SynSlo1+i)->gSlow *= (1.0f + (NItr_begin + (NItr->SynSlo1+i)->PosID)->LP * NItr->NegISI); } //multiplicative rule
        }
      }

      if(NItr->PosISI>0) //PosISI, facilitation
      {
        for(unsigned int i=0;i<NItr->SynSlowSize1;i++)
        {
          if((NItr_begin + (NItr->SynSlo1+i)->PosID)->SpiOut)
          { (NItr->SynSlo1+i)->gSlow += NItr->LP * NItr->PosISI; } //additive rule
        }
      }
      else //depression
      {
        for(unsigned int i=0;i<NItr->SynSlowSize1;i++)
        {
          if((NItr_begin + (NItr->SynSlo1+i)->PosID)->SpiOut)
          { (NItr->SynSlo1+i)->gSlow *= (1.0f + NItr->LP * NItr->PosISI); } //multiplicative rule
        }
      }
    }

//SynSlowSize-----------------------------------
    if(NItr->SynSlowSize != 0)
    {
      if(NItr->SpiOut)
      {
        if(NItr->NegISI>0) //facilitation
        {
          for(unsigned int i=0;i<NItr->SynSlowSize;i++)
          {
            for(unsigned int j=0;j<(NItr->SynSlo+i)->PosIDNum;j++) // post link population neurons
            { *((NItr->SynSlo+i)->gSlow+j) += (NItr_begin + *((NItr->SynSlo+i)->PosID+j))->LP * NItr->NegISI; } //additive rule
          }
        }
        else //depression
        {
          for(unsigned int i=0;i<NItr->SynSlowSize;i++)
          {
            for(unsigned int j=0;j<(NItr->SynSlo+i)->PosIDNum;j++) // post link population neurons
            { *((NItr->SynSlo+i)->gSlow+j) *= (1.0f + (NItr_begin + *((NItr->SynSlo+i)->PosID+j))->LP * NItr->NegISI); }//multiplicative rule
          }
        }
      }


      if(NItr->PosISI>0) //PosISI, facilitation
      {
        for(unsigned int i=0;i<NItr->SynSlowSize;i++)
        {
          for(unsigned int j=0;j<(NItr->SynSlo+i)->PosIDNum;j++) // post link population neurons
          {
            if((NItr_begin + *((NItr->SynSlo+i)->PosID+j))->SpiOut)
            { *((NItr->SynSlo+i)->gSlow+j) += NItr->LP * NItr->PosISI; }//additive rule
          }
        }
      }
      else  //depression
      {
        for(unsigned int i=0;i<NItr->SynSlowSize;i++)
        {
          for(unsigned int j=0;j<(NItr->SynSlo+i)->PosIDNum;j++) // post link population neurons
          {
            if((NItr_begin + *((NItr->SynSlo+i)->PosID+j))->SpiOut)
            { *((NItr->SynSlo+i)->gSlow+j) *= (1.0f + NItr->LP * NItr->PosISI); }//additive rule
          }
        }
      }
    }

//SynGlutSize1-----------------------------------------
    if(NItr->SynGlutSize1 != 0)
    {
      if(NItr->SpiOut)
      {
        if(NItr->NegISI)
        {
          for(unsigned int i=0;i<NItr->SynGlutSize1;i++)
          {
            (NItr->SynGlu1+i)->gFast += (NItr_begin + (NItr->SynGlu1+i)->PosID)->LP * NItr->NegISI; //additive rule
            (NItr->SynGlu1+i)->gSlow += (NItr_begin + (NItr->SynGlu1+i)->PosID)->LP * NItr->NegISI; //additive rule
          }
        }
        else
        {
          for(unsigned int i=0;i<NItr->SynGlutSize1;i++)
          {
            (NItr->SynGlu1+i)->gFast *= (1.0f + (NItr_begin + (NItr->SynGlu1+i)->PosID)->LP * NItr->NegISI); //multiplicative rule
            (NItr->SynGlu1+i)->gSlow *= (1.0f + (NItr_begin + (NItr->SynGlu1+i)->PosID)->LP * NItr->NegISI); //multiplicative rule
          }
        }
      }

      if(NItr->PosISI)
      {
        for(unsigned int i=0;i<NItr->SynGlutSize1;i++)
        {
          if((NItr_begin + (NItr->SynGlu1+i)->PosID)->SpiOut)
          {
            (NItr->SynGlu1+i)->gFast += NItr->LP * NItr->PosISI; //additive rule
            (NItr->SynGlu1+i)->gSlow += NItr->LP * NItr->PosISI; //additive rule
          }
        }
      }
      else
      {
        for(unsigned int i=0;i<NItr->SynGlutSize1;i++)
        {
          if((NItr_begin + (NItr->SynGlu1+i)->PosID)->SpiOut)
          {
            (NItr->SynGlu1+i)->gFast *= (1.0f + NItr->LP * NItr->PosISI); //additive rule
            (NItr->SynGlu1+i)->gSlow *= (1.0f + NItr->LP * NItr->PosISI); //additive rule
          }
        }
      }
    }

//SynGlutSize-----------------------------------------
    if(NItr->SynGlutSize != 0)
    {
      if(NItr->SpiOut)
      {
        if(NItr->NegISI>0)
        {
          for(unsigned int i=0;i<NItr->SynGlutSize;i++)
          {
            for(unsigned int j=0;j<(NItr->SynGlu+i)->PosIDNum;j++) // post link population neurons
            {
              *((NItr->SynGlu+i)->gFast+j) += (NItr_begin + *((NItr->SynGlu+i)->PosID+j))->LP * NItr->NegISI;//additive rule
              *((NItr->SynGlu+i)->gSlow+j) += (NItr_begin + *((NItr->SynGlu+i)->PosID+j))->LP * NItr->NegISI;//additive rule
            }
          }
        }
        else
        {
          for(unsigned int i=0;i<NItr->SynGlutSize;i++)
          {
            for(unsigned int j=0;j<(NItr->SynGlu+i)->PosIDNum;j++) // post link population neurons
            {
              *((NItr->SynGlu+i)->gFast+j) *= (1.0f + (NItr_begin + *((NItr->SynGlu+i)->PosID+j))->LP * NItr->NegISI);//multiplicative rule
              *((NItr->SynGlu+i)->gSlow+j) *= (1.0f + (NItr_begin + *((NItr->SynGlu+i)->PosID+j))->LP * NItr->NegISI);//multiplicative rule
            }
          }
        }
      }

      if(NItr->PosISI>0)
      {
        for(unsigned int i=0;i<NItr->SynGlutSize;i++)
        {
          for(unsigned int j=0;j<(NItr->SynGlu+i)->PosIDNum;j++) // post link population neurons
          {
            if((NItr_begin + *((NItr->SynGlu+i)->PosID+j))->SpiOut)
            {
              *((NItr->SynGlu+i)->gFast+j) += NItr->LP * NItr->PosISI; //additive rule
              *((NItr->SynGlu+i)->gSlow+j) += NItr->LP * NItr->PosISI; //additive rule
            }
          }
        }
      }
      else
      {
        for(unsigned int i=0;i<NItr->SynGlutSize;i++)
        {
          for(unsigned int j=0;j<(NItr->SynGlu+i)->PosIDNum;j++) // post link population neurons
          {
            if((NItr_begin + *((NItr->SynGlu+i)->PosID+j))->SpiOut)
            {
              *((NItr->SynGlu+i)->gFast+j) *= (1.0f + NItr->LP * NItr->PosISI);//multiplicative rule
              *((NItr->SynGlu+i)->gSlow+j) *= (1.0f + NItr->LP * NItr->PosISI);//multiplicative rule
            }
          }
        }
      }
    }

//SynGapSize1-----------------------------------------
    if(NItr->SynGapSize1 != 0)
    {
      if(NItr->SpiOut)
      {
        if(NItr->NegISI>0)
        {
          for(unsigned int i=0;i<NItr->SynGapSize1;i++)
          { (NItr->SynCnd1+i)->gC += (NItr_begin + (NItr->SynCnd1+i)->PosID)->LP * NItr->NegISI; } //additive rule
        }
        else
        {
          for(unsigned int i=0;i<NItr->SynGapSize1;i++)
          { (NItr->SynCnd1+i)->gC *= (1.0f + (NItr_begin + (NItr->SynCnd1+i)->PosID)->LP * NItr->NegISI); } //multiplicative rule
        }
      }

      if(NItr->PosISI>0)
      {
        for(unsigned int i=0;i<NItr->SynGapSize1;i++) // pre-neuron synapses
        {
          if((NItr_begin + (NItr->SynCnd1+i)->PosID)->SpiOut)
          { (NItr->SynCnd1+i)->gC += NItr->LP * NItr->PosISI; }//additive rule
        }
      }
      else
      {
        for(unsigned int i=0;i<NItr->SynGapSize1;i++) // pre-neuron synapses
        {
          if((NItr_begin + (NItr->SynCnd1+i)->PosID)->SpiOut)
          { (NItr->SynCnd1+i)->gC *= (1.0f + NItr->LP * NItr->PosISI); }//multiplicative rule
        }
      }

    }

//SynGapSize-------------------------------------------
    if(NItr->SynGapSize != 0)
    {
      if(NItr->SpiOut)
      {
        if(NItr->NegISI>0)
        {
          for(unsigned int i=0;i<NItr->SynGapSize;i++)
          {
            for(unsigned int j=0;j<(NItr->SynCnd+i)->PosIDNum;j++) // post link population neurons
            { *((NItr->SynCnd+i)->gC+j) += (NItr_begin + *((NItr->SynCnd+i)->PosID+j))->LP * NItr->NegISI; }//additive rule
          }
        }
        else
        {
          for(unsigned int i=0;i<NItr->SynGapSize;i++)
          {
            for(unsigned int j=0;j<(NItr->SynCnd+i)->PosIDNum;j++) // post link population neurons
            { *((NItr->SynCnd+i)->gC+j) *= (1.0f + (NItr_begin + *((NItr->SynCnd+i)->PosID+j))->LP * NItr->NegISI); }//multiplicative rule
          }
        }
      }

      if(NItr->PosISI>0)
      {
        for(unsigned int i=0;i<NItr->SynGapSize;i++) // pre-neuron synapses
        {
          for(unsigned int j=0;j<(NItr->SynCnd+i)->PosIDNum;j++) // post link population neurons
          {
            if((NItr_begin + *((NItr->SynCnd+i)->PosID+j))->SpiOut)
            { *((NItr->SynCnd+i)->gC+j) += NItr->LP * NItr->PosISI; } //additive rule
          }
        }
      }
      else
      {
        for(unsigned int i=0;i<NItr->SynGapSize;i++) // pre-neuron synapses
        {
          for(unsigned int j=0;j<(NItr->SynCnd+i)->PosIDNum;j++) // post link population neurons
          {
            if((NItr_begin + *((NItr->SynCnd+i)->PosID+j))->SpiOut)
            { *((NItr->SynCnd+i)->gC+j) *= (1.0f + NItr->LP * NItr->PosISI); } //multiplicative rule
          }
        }
      }
    }
  }


}

void EqsGnl::Activation(unsigned int NID)
{
  NItr=(NItr_begin+NID);
  NItr->BgInit(dt);//bcak ground initial

  switch(SolType)
  {
    case ACCURATE:
    VmRK4();
    break;

    case ROUGH:
    VmEul();
    break;

    case MODERATE:
    default:
    VmIpv();
    break;
  }
}

