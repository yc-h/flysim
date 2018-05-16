
#include <string>
#include <cmath>
#include <vector>
#include <omp.h>
#include <thread>
#include <mutex>
#include "Eqs.h"

using namespace std;

float EqsLIF_gp4::FFst(float st)
{ return ( SpiIn - st/PItr->tst); }


float EqsLIF_gp4::Fst(float st)
{ return -st/PItr2->tst; }

float EqsLIF_gp4::FGS1(float st)
{ return -st/tst; }

float EqsLIF_gp4::FGS2(float st)
{ return -st/tst; }

float EqsLIF_gp4::MemPot(float Vm)
{ return ( VmItr->Isyn - NItr->gl * (Vm - NItr->El) ) / NItr->Cm ; }


void EqsLIF_gp4::FRK4Irep1Acc(float GSum,float p)
{
  for(VItr=Rep_begin,VmItr=NItr->NeuVm.begin();VItr<Rep_end;VItr++,VmItr++)
  {
    bool Sti=false;
    if(p>NItr->uni_dist(NItr->rdgen))
      Sti=true;

    if(NItr->SelfSti && VmItr->SpiOut)
    { GSum -= PItr->MeanG;}

    if(Sti)
    { SetSti1(); }
    else
    { VItr->st = RK4(VItr->st, &EqsLIF_gp4::FFst); }

    VmItr->Isyn += (VItr->st*PItr->g + GSum) * ( PItr->Vrev - VmItr->V );
  }
}


void EqsLIF_gp4::FRK4Irep1AccG(float GSum,float p)
{
  for(VItr=Rep_begin,VmItr=NItr->NeuVm.begin();VItr<Rep_end;VItr++,VmItr++)
  {
    bool Sti=false;
    if(p>NItr->uni_dist(NItr->rdgen))
      Sti=true;

    if(NItr->SelfSti && VmItr->SpiOut)
    { GSum -= PItr->MeanG;}

    if(Sti)
    { SetSti1(); }
    else
    { VItr->st = RK4(VItr->st, &EqsLIF_gp4::FFst); }

    VmItr->Isyn += (VItr->st*PItr->g + GSum) * ( PItr->Vrev - VmItr->V );
  }
}



void EqsLIF_gp4::FRK4Irep2Acc(float GSum,float p)
{
  for(VItr=Rep_begin,VmItr=NItr->NeuVm.begin();VItr<Rep_end;VItr++,VmItr++)
  {
    bool Sti=false;
    if(p>NItr->uni_dist(NItr->rdgen))
      Sti=true;

    if(NItr->SelfSti && VmItr->SpiOut)
    { GSum -= PItr2->MeanG;}

    if(Sti)
    { SetSti4(); }
    else
    { VItr->st = RK4(VItr->st, &EqsLIF_gp4::Fst); }

    VmItr->Isyn += (VItr->st*PItr2->g + GSum) * ( PItr2->Vrev - VmItr->V ) / ( 1.0f + 1.0f * PItr2->mg * (float)exp( -0.062e3f * VmItr->V ) / 3.57f );
  }

}


void EqsLIF_gp4::FIpvIrep1Acc(float GSum,float p)
{
  for(VItr=Rep_begin,VmItr=NItr->NeuVm.begin();VItr<Rep_end;VItr++,VmItr++)
  {
    bool Sti=false;
    if(p>NItr->uni_dist(NItr->rdgen))
      Sti=true;

    if(NItr->SelfSti && VmItr->SpiOut)
    { GSum -= PItr->MeanG;}

    if(Sti)
    { SetSti1(); }
    else
    { VItr->st = IpvEul(VItr->st, &EqsLIF_gp4::FFst); }

    VmItr->Isyn += (VItr->st*PItr->g + GSum) * ( PItr->Vrev - VmItr->V );
  }
}


void EqsLIF_gp4::FIpvIrep1AccG(float GSum,float p)
{
  for(VItr=Rep_begin,VmItr=NItr->NeuVm.begin();VItr<Rep_end;VItr++,VmItr++)
  {
    bool Sti=false;
    if(p>NItr->uni_dist(NItr->rdgen))
      Sti=true;

    if(NItr->SelfSti && VmItr->SpiOut)
    { GSum -= PItr->MeanG;}

    if(Sti)
    { SetSti1(); }
    else
    { VItr->st = IpvEul(VItr->st, &EqsLIF_gp4::FFst); }

    VmItr->Isyn += (VItr->st*PItr->g + GSum) * ( PItr->Vrev - VmItr->V );
  }
}



void EqsLIF_gp4::FIpvIrep2Acc(float GSum,float p)
{
  for(VItr=Rep_begin,VmItr=NItr->NeuVm.begin();VItr<Rep_end;VItr++,VmItr++)
  {
    bool Sti=false;
    if(p>NItr->uni_dist(NItr->rdgen))
      Sti=true;

    if(NItr->SelfSti && VmItr->SpiOut)
    { GSum -= PItr2->MeanG;}

    if(Sti)
    { SetSti4(); }
    else
    { VItr->st = IpvEul(VItr->st, &EqsLIF_gp4::Fst); }

    VmItr->Isyn += (VItr->st*PItr2->g + GSum) * ( PItr2->Vrev - VmItr->V ) / ( 1.0f + 1.0f * PItr2->mg * (float)exp( -0.062e3f * VmItr->V ) / 3.57f );
  }

}



void EqsLIF_gp4::FEulIrep1Acc(float GSum,float p)
{
  for(VItr=Rep_begin,VmItr=NItr->NeuVm.begin();VItr<Rep_end;VItr++,VmItr++)
  {
    bool Sti=false;
    if(p>NItr->uni_dist(NItr->rdgen))
      Sti=true;

    if(NItr->SelfSti && VmItr->SpiOut)
    { GSum -= PItr->MeanG;}

    if(Sti)
    { SetSti1(); }
    else
    { VItr->st = VItr->st - DTIME*(VItr->st/PItr->tst); }

    VmItr->Isyn += (VItr->st*PItr->g + GSum) * ( PItr->Vrev - VmItr->V );

  }
}


void EqsLIF_gp4::FEulIrep1AccG(float GSum,float p)
{
  for(VItr=Rep_begin,VmItr=NItr->NeuVm.begin();VItr<Rep_end;VItr++,VmItr++)
  {
    bool Sti=false;
    if(p>NItr->uni_dist(NItr->rdgen))
      Sti=true;

    if(NItr->SelfSti && VmItr->SpiOut)
    { GSum -= PItr->MeanG;}

    if(Sti)
    { SetSti1(); }
    else
    { VItr->st = VItr->st - DTIME*(VItr->st/PItr->tst); }

    VmItr->Isyn += (VItr->st*PItr->g + GSum) * ( PItr->Vrev - VmItr->V );

  }
}



void EqsLIF_gp4::FEulIrep2Acc(float GSum,float p)
{
  for(VItr=Rep_begin,VmItr=NItr->NeuVm.begin();VItr<Rep_end;VItr++,VmItr++)
  {
    bool Sti=false;
    if(p>NItr->uni_dist(NItr->rdgen))
      Sti=true;

    if(NItr->SelfSti && VmItr->SpiOut)
    { GSum -= PItr2->MeanG;}

    if(Sti)
    { SetSti4(); }
    else
    { VItr->st = VItr->st - DTIME*(VItr->st/PItr2->tst); }

    VmItr->Isyn += (VItr->st*PItr2->g + GSum) * ( PItr2->Vrev - VmItr->V ) / ( 1.0f + 1.0f * PItr2->mg * (float)exp( -0.062e3f * VmItr->V ) / 3.57f );

  }
}



void EqsLIF_gp4::firing()
{
    if(VmItr->Count < NItr->RfPed) //Firing
    {

      if(VmItr->Count == 0)
        VmItr->V=0.0f;
      else
        VmItr->V=NItr->Vreset*1e-3f;

      if(VmItr->Count==NItr->SpiDy)
        VmItr->SpiOut=true;
      else
        VmItr->SpiOut=false;

      VmItr->Count++;
    }
    else  // silence time is over
    { 
      VmItr->Count=0;
      VmItr->Firing=false;
    }
}



void EqsLIF_gp4::FRK4Vm()
{
  for(VmItr=NItr->NeuVm.begin();VmItr<NItr->NeuVm.end();VmItr++)
  {
    if(VmItr->Firing)
    { firing(); }
    else
    {
      VmItr->Isyn += NItr->Igauss();
      VmItr->V = RK4(VmItr->V, &EqsLIF_gp4::MemPot);

      if(VmItr->V > NItr->Vth)
      {
        VmItr->Firing=true;
        firing();
      }
    }
  }

  NItr->gAMPA=0.0f;
  NItr->gGABA=0.0f;
  NItr->gNMDA=0.0f;
}


void EqsLIF_gp4::FIpvVm()
{
  for(VmItr=NItr->NeuVm.begin();VmItr<NItr->NeuVm.end();VmItr++)
  {
    if(VmItr->Firing)
    { firing(); }
    else
    {
      VmItr->Isyn += NItr->Igauss();
      VmItr->V = IpvEul(VmItr->V, &EqsLIF_gp4::MemPot);

      if(VmItr->V > NItr->Vth)
      {
        VmItr->Firing=true;
        firing();
      }
    }
  }

  NItr->gAMPA=0.0f;
  NItr->gGABA=0.0f;
  NItr->gNMDA=0.0f;
}



void EqsLIF_gp4::FEulVm()
{
  for(VmItr=NItr->NeuVm.begin();VmItr<NItr->NeuVm.end();VmItr++)
  {
    if(VmItr->Firing)
    { firing(); }
    else
    {
      VmItr->Isyn += NItr->Igauss();
      VmItr->V = VmItr->V + DTIME * ( VmItr->Isyn - NItr->gl * (VmItr->V - NItr->El) ) / NItr->Cm;

      if(VmItr->V > NItr->Vth)
      {
        VmItr->Firing=true;
        firing();
      }
    }
  }

  NItr->gAMPA=0.0f;
  NItr->gGABA=0.0f;
  NItr->gNMDA=0.0f;
}


void EqsLIF_gp4::SpikeRK4()
{
  vector<NeuVar>::iterator VmItr;
  vector<Nlink>::iterator CItr;
  float gAcc;

  for(vector<ConTab_gp2>::iterator IItr=NItr->Index.begin();IItr<NItr->Index.end();IItr++)
  {
    if(IItr->AMPA_Act)
    {
      gAcc=0.0f;
      VmItr=NItr->NeuVm.begin();
      CItr=IItr->Axn.begin();
      tst=IItr->AMPA_tst;
      for(;CItr<IItr->Axn.end();CItr++,VmItr++)
      {
        if(VmItr->SpiOut)
          CItr->AMPA_st = CItr->AMPA_st + 1.0f;
        else
        { CItr->AMPA_st = RK4(CItr->AMPA_st, &EqsLIF_gp4::FGS1); }

        gAcc += CItr->AMPA_st;
      }

      (NItr_begin + IItr->NeuID)->gAMPA += IItr->AMPA_MeanG*gAcc;
    }

    if(IItr->GABA_Act)
    {
      gAcc=0.0f;
      VmItr=NItr->NeuVm.begin();
      CItr=IItr->Axn.begin();
      tst=IItr->GABA_tst;
      for(;CItr<IItr->Axn.end();CItr++,VmItr++)
      {
        if(VmItr->SpiOut)
          CItr->GABA_st = CItr->GABA_st + 1.0f;
        else
        { CItr->GABA_st = RK4(CItr->GABA_st, &EqsLIF_gp4::FGS1); }

        gAcc += CItr->GABA_st;
      }

      (NItr_begin + IItr->NeuID)->gGABA += IItr->GABA_MeanG*gAcc;
    }

    if(IItr->NMDA_Act)
    {
      gAcc=0.0f;
      VmItr=NItr->NeuVm.begin();
      CItr=IItr->Axn.begin();
      tst=IItr->NMDA_tst;
      as=IItr->NMDA_as;
      for(;CItr<IItr->Axn.end();CItr++,VmItr++)
      {
        if(VmItr->SpiOut)
          CItr->NMDA_st = CItr->NMDA_st + IItr->NMDA_as*(1.0f - CItr->NMDA_st);
        else
        { CItr->NMDA_st = RK4(CItr->NMDA_st, &EqsLIF_gp4::FGS2); }

        gAcc += CItr->NMDA_st;
      }

      (NItr_begin + IItr->NeuID)->gNMDA += IItr->NMDA_MeanG*gAcc;
    }

    if(IItr->GAP_Act)
    {
      for(VmItr=NItr->NeuVm.begin();VmItr<NItr->NeuVm.end();VmItr++)
      {
         for(vector<NeuVar>::iterator poVmItr=(NItr_begin + IItr->NeuID)->NeuVm.begin();poVmItr<(NItr_begin + IItr->NeuID)->NeuVm.end();poVmItr++)
         {  poVmItr->Isyn += IItr->GAP_MeanG*(VmItr->V - poVmItr->V);  }
      }
    }

  }
}


void EqsLIF_gp4::SpikeIpv()
{
  vector<NeuVar>::iterator VmItr;
  vector<Nlink>::iterator CItr;
  float gAcc;
  for(vector<ConTab_gp2>::iterator IItr=NItr->Index.begin();IItr<NItr->Index.end();IItr++)
  {
    if(IItr->AMPA_Act)
    {
      gAcc=0.0f;
      VmItr=NItr->NeuVm.begin();
      CItr=IItr->Axn.begin();
      tst=IItr->AMPA_tst;
      for(;CItr<IItr->Axn.end();CItr++,VmItr++)
      {
        if(VmItr->SpiOut)
          CItr->AMPA_st = CItr->AMPA_st + 1.0f;
        else
        { CItr->AMPA_st = IpvEul(CItr->AMPA_st, &EqsLIF_gp4::FGS1); }

        gAcc += CItr->AMPA_st;
      }

      (NItr_begin + IItr->NeuID)->gAMPA += IItr->AMPA_MeanG*gAcc;
    }

    if(IItr->GABA_Act)
    {
      gAcc=0.0f;
      VmItr=NItr->NeuVm.begin();
      CItr=IItr->Axn.begin();
      tst=IItr->GABA_tst;
      for(;CItr<IItr->Axn.end();CItr++,VmItr++)
      {
        if(VmItr->SpiOut)
          CItr->GABA_st = CItr->GABA_st + 1.0f;
        else
        { CItr->GABA_st = IpvEul(CItr->GABA_st, &EqsLIF_gp4::FGS1); }

        gAcc += CItr->GABA_st;
      }

      (NItr_begin + IItr->NeuID)->gGABA += IItr->GABA_MeanG*gAcc;
    }

    if(IItr->NMDA_Act)
    {
      gAcc=0.0f;
      VmItr=NItr->NeuVm.begin();
      CItr=IItr->Axn.begin();
      tst=IItr->NMDA_tst;
      as=IItr->NMDA_as;
      for(;CItr<IItr->Axn.end();CItr++,VmItr++)
      {

        if(VmItr->SpiOut)
          CItr->NMDA_st = CItr->NMDA_st + IItr->NMDA_as*(1.0f - CItr->NMDA_st);
        else
        { CItr->NMDA_st = IpvEul(CItr->NMDA_st, &EqsLIF_gp4::FGS2); }

        gAcc += CItr->NMDA_st;
      }

      (NItr_begin + IItr->NeuID)->gNMDA += IItr->NMDA_MeanG*gAcc;
    }

    if(IItr->GAP_Act)
    {
      for(VmItr=NItr->NeuVm.begin();VmItr<NItr->NeuVm.end();VmItr++)
      {
         for(vector<NeuVar>::iterator poVmItr=(NItr_begin + IItr->NeuID)->NeuVm.begin();poVmItr<(NItr_begin + IItr->NeuID)->NeuVm.end();poVmItr++)
         {  poVmItr->Isyn += IItr->GAP_MeanG*(VmItr->V - poVmItr->V);  }
      }
    }

  }
}


void EqsLIF_gp4::SpikeEul()
{
  vector<NeuVar>::iterator VmItr;
  vector<Nlink>::iterator CItr;
  float gAcc;
  for(vector<ConTab_gp2>::iterator IItr=NItr->Index.begin();IItr<NItr->Index.end();IItr++)
  {
    if(IItr->AMPA_Act)
    {
      gAcc=0.0f;
      VmItr=NItr->NeuVm.begin();
      CItr=IItr->Axn.begin();
      for(;CItr<IItr->Axn.end();CItr++,VmItr++)
      {
        if(VmItr->SpiOut)
          CItr->AMPA_st = CItr->AMPA_st + 1.0f;
        else
          CItr->AMPA_st = CItr->AMPA_st - DTIME*(CItr->AMPA_st/IItr->AMPA_tst);

        gAcc += CItr->AMPA_st;
      }

      (NItr_begin + IItr->NeuID)->gAMPA += IItr->AMPA_MeanG*gAcc;
    }

    if(IItr->GABA_Act)
    {
      gAcc=0.0f;
      VmItr=NItr->NeuVm.begin();
      CItr=IItr->Axn.begin();
      for(;CItr<IItr->Axn.end();CItr++,VmItr++)
      {

        if(VmItr->SpiOut)
          CItr->GABA_st = CItr->GABA_st + 1.0f;
        else
          CItr->GABA_st = CItr->GABA_st - DTIME*(CItr->GABA_st/IItr->GABA_tst);

        gAcc += CItr->GABA_st;
      }

      (NItr_begin + IItr->NeuID)->gGABA += IItr->GABA_MeanG*gAcc;
    }

    if(IItr->NMDA_Act)
    {
      gAcc=0.0f;
      VmItr=NItr->NeuVm.begin();
      CItr=IItr->Axn.begin();
      for(;CItr<IItr->Axn.end();CItr++,VmItr++)
      {
        if(VmItr->SpiOut)
          CItr->NMDA_st = CItr->NMDA_st + IItr->NMDA_as*(1.0f - CItr->NMDA_st);
        else
          CItr->NMDA_st = CItr->NMDA_st - DTIME*(CItr->NMDA_st/IItr->NMDA_tst);

        gAcc += CItr->NMDA_st;
      }

      (NItr_begin + IItr->NeuID)->gNMDA += IItr->NMDA_MeanG*gAcc;
    }

    if(IItr->GAP_Act)
    {
      for(VmItr=NItr->NeuVm.begin();VmItr<NItr->NeuVm.end();VmItr++)
      {
         for(vector<NeuVar>::iterator poVmItr=(NItr_begin + IItr->NeuID)->NeuVm.begin();poVmItr<(NItr_begin + IItr->NeuID)->NeuVm.end();poVmItr++)
         { poVmItr->Isyn += IItr->GAP_MeanG*(VmItr->V - poVmItr->V);  }
      }
    }

  }
}


void EqsLIF_gp4::Spike(unsigned int NID)
{

  NItr=NItr_begin+NID;

  switch(SolType)
  {
    case ACCURATE:
      SpikeRK4();
    break;

    case MODERATE:
      SpikeIpv();
    break;


    case ROUGH:
      SpikeEul();
    break;

    default:
      SpikeIpv();
    break;
  }
}


void EqsLIF_gp4::Activation(unsigned int NID)
{
  NItr=NItr_begin+NID;

  switch(SolType)
  {
    case ACCURATE:
      //synapse current calculation
      //AMBA receptor
      PItr = &NItr->ParAMPA;
      Rep_begin=NItr->RepAMPA.begin();
      Rep_end=NItr->RepAMPA.end();
      FRK4Irep1Acc(NItr->gAMPA,NItr->p1);

      //GABA receptor
      PItr = &NItr->ParGABA;
      Rep_begin=NItr->RepGABA.begin();
      Rep_end=NItr->RepGABA.end();
      FRK4Irep1AccG(NItr->gGABA,NItr->p2);

      //NMDA receptor
      PItr2 = &NItr->ParNMDA;
      Rep_begin=NItr->RepNMDA.begin();
      Rep_end=NItr->RepNMDA.end();
      FRK4Irep2Acc(NItr->gNMDA,NItr->p3);

      //membrane potential
      FRK4Vm();
    break;

    case MODERATE:
      //synapse current calculation
      //AMBA receptor
      PItr = &NItr->ParAMPA;
      Rep_begin=NItr->RepAMPA.begin();
      Rep_end=NItr->RepAMPA.end();
      FIpvIrep1Acc(NItr->gAMPA,NItr->p1);

      //GABA receptor
      PItr = &NItr->ParGABA;
      Rep_begin=NItr->RepGABA.begin();
      Rep_end=NItr->RepGABA.end();
      FIpvIrep1AccG(NItr->gGABA,NItr->p2);

      //NMDA receptor
      PItr2 = &NItr->ParNMDA;
      Rep_begin=NItr->RepNMDA.begin();
      Rep_end=NItr->RepNMDA.end();
      FIpvIrep2Acc(NItr->gNMDA,NItr->p3);

      //membrane potential
      FIpvVm();
    break;


    case ROUGH:
      //synapse current calculation
      //AMBA receptor
      PItr = &NItr->ParAMPA;
      Rep_begin=NItr->RepAMPA.begin();
      Rep_end=NItr->RepAMPA.end();
      FEulIrep1Acc(NItr->gAMPA,NItr->p1);

      //GABA receptor
      PItr = &NItr->ParGABA;
      Rep_begin=NItr->RepGABA.begin();
      Rep_end=NItr->RepGABA.end();
      FEulIrep1AccG(NItr->gGABA,NItr->p2);

      //NMDA receptor
      PItr2 = &NItr->ParNMDA;
      Rep_begin=NItr->RepNMDA.begin();
      Rep_end=NItr->RepNMDA.end();
      FEulIrep2Acc(NItr->gNMDA,NItr->p3);

      //membrane potential
      FEulVm();
    break;

    default:
      //synapse current calculation
      //AMBA receptor
      PItr = &NItr->ParAMPA;
      Rep_begin=NItr->RepAMPA.begin();
      Rep_end=NItr->RepAMPA.end();
      FIpvIrep1Acc(NItr->gAMPA,NItr->p1);

      //GABA receptor
      PItr = &NItr->ParGABA;
      Rep_begin=NItr->RepGABA.begin();
      Rep_end=NItr->RepGABA.end();
      FIpvIrep1AccG(NItr->gGABA,NItr->p2);

      //NMDA receptor
      PItr2 = &NItr->ParNMDA;
      Rep_begin=NItr->RepNMDA.begin();
      Rep_end=NItr->RepNMDA.end();
      FIpvIrep2Acc(NItr->gNMDA,NItr->p3);

      //membrane potential
      FIpvVm();
    break;
  }

  //bcak ground initial
  NItr->BgInit();
}

