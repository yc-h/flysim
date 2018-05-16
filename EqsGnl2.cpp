
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <thread>
#include <atomic>
#include <condition_variable>
#include "Eqs.h"

using namespace std;

//global variables

EqsGnl2::EqsGnl2()
{
  dt=DTIME;
  thsize=1;

  NetSize=0;
  RPsSize=0;
  RVsSize=0;
  NetPBseSize=0;
  NetVBseSize=0;

  NetPLIFSize=0;
  NetVLIFSize=0;

  NetPSTPSize=0;
  NetVSTPSize=0;

  NetPLTPSize=0;
  NetVLTPSize=0;

  SynNMDA1Size=0;
  SynCnd1Size=0;

  AMPA1CntSize=0;
  GABA1CntSize=0;
  NMDA1CntSize=0;
  ACH1CntSize=0;
  GCL1CntSize=0;
  Cnd1CntSize=0;


  PreAMPA1CntSize=0;
  PreGABA1CntSize=0;
  PreACH1CntSize=0;
  PreGCL1CntSize=0;
  PreNMDA1CntSize=0;
  PreCnd1CntSize=0;

  PSynASize=0;
  PSynGSize=0;
  PSynHSize=0;
  PSynLSize=0;
  PSynNSize=0;
  PSynCpSize=0;

  PGatASize=0;
  PGatGSize=0;
  PGatHSize=0;
  PGatLSize=0;
  PGatNSize=0;
  PGapPSize=0;

  PAMPACntSize=0;
  PGABACntSize=0;
  PACHCntSize=0;
  PGCLCntSize=0;
  PNMDACntSize=0;
  PCndPCntSize=0;
  PCndOCntSize=0;

  th=nullptr;
  Nets=nullptr;
  RPs=nullptr;
  RVs=nullptr;
  NetPBse=nullptr;
  NetVBse=nullptr;

  NetPLIF=nullptr;
  NetVLIF=nullptr;

  NetPSTP=nullptr;
  NetVSTP=nullptr;

  NetPLTP=nullptr;
  NetVLTP=nullptr;

  SynFst1=nullptr;
  SynSlo1=nullptr;
  SynCnd1=nullptr;

  NeuMapAMPA1=nullptr;
  NeuMapGABA1=nullptr;
  NeuMapACH1=nullptr;
  NeuMapGCL1=nullptr;
  NeuMapNMDA1=nullptr;
  NeuMapCnd1=nullptr;
  NeuPreMapCnd1=nullptr;

  NeuPreMapAMPA1=nullptr;
  NeuPreMapGABA1=nullptr;
  NeuPreMapACH1=nullptr;
  NeuPreMapGCL1=nullptr;
  NeuPreMapNMDA1=nullptr;
  NeuPreMapCnd1=nullptr;

  PreMapAMPA1=nullptr;
  PreMapGABA1=nullptr;
  PreMapACH1=nullptr;
  PreMapGCL1=nullptr;
  PreMapNMDA1=nullptr;
  PreMapCnd1=nullptr;

  PSynA=nullptr;
  PSynG=nullptr;
  PSynH=nullptr;
  PSynL=nullptr;
  PSynN=nullptr;
  PSynCp=nullptr; //to PreSynapseID
  PSynCo=nullptr; //to PosNeuronID

  PSynFstGatA=nullptr;  //AMPA
  PSynFstGatG=nullptr;  //GABA
  PSynFstGatH=nullptr;  //ACH
  PSynFstGatL=nullptr;  //GCL
  PSynSloGat=nullptr;
  PSynCndP=nullptr; //integrate PreID current to PosID

  NeuMapPAMPA=nullptr;
  NeuMapPGABA=nullptr;
  NeuMapPACH=nullptr;
  NeuMapPGCL=nullptr;
  NeuMapPNMDA=nullptr;

  NeuMapPCndP=nullptr; //integrate PreID current to PosID
  NeuMapPCndO=nullptr; //integrate PosID current to PreID

  MapPAPosID=nullptr;
  MapPGPosID=nullptr;
  MapPHPosID=nullptr;
  MapPLPosID=nullptr;
  MapPNPosID=nullptr;
  MapPCPosID=nullptr;


}

//--------------------membrane equation
void EqsGnl2::VmEul(unsigned int NID)
{
  NeuRepPar2<RP3,RP4> NRP= RPs[NID];
  NeuRepVar2 NRV=RVs[NID];
  ParBase2 NPB=NetPBse[NID];

  if(NetVLIF[NID].NState & FIRE_F)
    firing(NID);
  else
  {
    float Irep[5];
    Irep[0]=NRV.gAMPA*(NRP.ParAMPA.Vrev-NetVBse[NID].V);
    Irep[1]=NRV.gGABA*(NRP.ParGABA.Vrev-NetVBse[NID].V);
    Irep[2]=NRV.gNMDA*(NRP.ParNMDA.Vrev-NetVBse[NID].V)/( 1.0f + NRP.ParNMDA.mg * (float)exp( -0.062f * NetVBse[NID].V ) / 3.57f );
    Irep[3]=NRV.gACH*(NRP.ParACH.Vrev-NetVBse[NID].V);
    Irep[4]=NRV.gGCL*(NRP.ParGCL.Vrev-NetVBse[NID].V);

    float I =
    Irep[0]+Irep[1]+Irep[2]+Irep[3]+Irep[4]+
    NetVBse[NID].Isyn*1e3f+
    NPB.gl*NPB.El*1e3f;

    I=1e-3f*I/NPB.Cm;

    if(I>NPB.dVmax)
    {
      I=NPB.dVmax;
      NetVLIF[NID].NState |= MEMDV_OVF; //active MEMDV_OVF
      NetVLIF[NID].NState &= ~MEMDV_UDF; //deactive MEMDV_UDF
    }
    else if(I<-NPB.dVmax)
    {
      I=-NPB.dVmax;
      NetVLIF[NID].NState &= ~MEMDV_OVF; //deactive MEMDV_OVF
      NetVLIF[NID].NState |= MEMDV_UDF; //active MEMDV_UDF
    }
    else
    {
      NetVLIF[NID].NState &= ~MEMDV_OVF; //deactive MEMDV_OVF
      NetVLIF[NID].NState &= ~MEMDV_UDF; //deactive MEMDV_UDF
    }

    NetVBse[NID].V = NPB.df*NetVBse[NID].V + NPB.ef*I;

    Nets[NID].Irep[0]=Irep[0];
    Nets[NID].Irep[1]=Irep[1];
    Nets[NID].Irep[2]=Irep[2];
    Nets[NID].Irep[3]=Irep[3];
    Nets[NID].Irep[4]=Irep[4];

    if(NetVBse[NID].V > NPB.Vth)
    {
      NetVLIF[NID].NState |= FIRE_F;
      firing(NID);

      if((Nets[NID].ActCtr & LTP_ACT) != 0) NetVLIF[NID].NState |= LP_FIRE_F;
    }
  }

//synapse plasticity
  if((Nets[NID].ActCtr & LTP_ACT) != 0)
  {
    firLP(NID);
    ltpEul(NID);
  }
  NRV.gAMPA=0.0f;
  NRV.gGABA=0.0f;
  NRV.gNMDA=0.0f;
  NRV.gACH=0.0f;
  NRV.gGCL=0.0f;
  NetVBse[NID].Isyn=0.0f;

  RVs[NID]=NRV;
}


//--------------------firing process
void EqsGnl2::firing(unsigned int NID)
{
  if(NetVLIF[NID].Count < NetPLIF[NID].RfPed) //Firing
  {
    if(NetVLIF[NID].Count == 0)
      NetVBse[NID].V=0.0f;
    else
      NetVBse[NID].V=NetPBse[NID].Vreset;

    if(NetVLIF[NID].Count==NetPLIF[NID].SpiDy)
      NetVLIF[NID].NState |= SPIKE_F;
    else
      NetVLIF[NID].NState &= ~SPIKE_F;

    NetVLIF[NID].Count++;
  }
  else  // silence time is over
  {
    NetVLIF[NID].Count=0;
    NetVLIF[NID].NState &= ~FIRE_F;
  }

}

//-------------long term plasticity
void EqsGnl2::firLP(unsigned int NID)
{
  if(NetVLIF[NID].NState & LP_FIRE_F)
  {
    if(NetVLTP[NID].CountLP==NetPLTP[NID].SpiLPDy)
    {
      NetVLIF[NID].NState |= LP_SPIKE_F;
      NetVLIF[NID].NState &= ~LP_FIRE_F;
      NetVLTP[NID].CountLP=0;
    }
    else
    {
      NetVLIF[NID].NState &= ~LP_SPIKE_F;
      NetVLTP[NID].CountLP++;
    }
  }
  else
  {  NetVLIF[NID].NState &= ~LP_SPIKE_F; }
}

void EqsGnl2::ltpEul(unsigned int NID)
{
  if(NetVLIF[NID].NState & LP_SPIKE_F)
    NetVLTP[NID].LP ++;
  else
    NetVLTP[NID].LP = NetPLTP[NID].dfLP*NetVLTP[NID].LP;
}


//-------------background initialization
void EqsGnl2::BgInit(unsigned int thid, unsigned int NID)
{
  unsigned int ActCtr=Nets[NID].ActCtr;

  if((ActCtr & AMPA_ACT)!=0)
  {
    if( BgSt[NID].pb[0]>BgSt[NID].uni_dist(*rdgen) )//AMPA
      BgSt[NID].SynFst1ES[0].st++;
    else
      BgSt[NID].SynFst1ES[0].st *= RPs[NID].ParAMPA.dfst;
    RVs[NID].gAMPA += BgSt[NID].SynFst1ES[0].g * BgSt[NID].SynFst1ES[0].st;
  }

  if((ActCtr & GABA_ACT)!=0)
  {
    if( BgSt[NID].pb[1]>BgSt[NID].uni_dist(*rdgen) )//GABA
      BgSt[NID].SynFst1ES[1].st++;
    else
      BgSt[NID].SynFst1ES[1].st*=RPs[NID].ParGABA.dfst;
    RVs[NID].gGABA += BgSt[NID].SynFst1ES[1].g * BgSt[NID].SynFst1ES[1].st;
  }

  if((ActCtr & ACH_ACT)!=0)
  {
    if( BgSt[NID].pb[2]>BgSt[NID].uni_dist(*rdgen) )//ACH
      BgSt[NID].SynFst1ES[2].st++;
    else
      BgSt[NID].SynFst1ES[2].st*=RPs[NID].ParACH.dfst;
    RVs[NID].gACH += BgSt[NID].SynFst1ES[2].g * BgSt[NID].SynFst1ES[2].st;
  }

  if((ActCtr & GCL_ACT)!=0)
  {
    if( BgSt[NID].pb[3]>BgSt[NID].uni_dist(*rdgen) )//GCL
      BgSt[NID].SynFst1ES[3].st++;
    else
      BgSt[NID].SynFst1ES[3].st*=RPs[NID].ParGCL.dfst;
    RVs[NID].gGCL += BgSt[NID].SynFst1ES[3].g * BgSt[NID].SynFst1ES[3].st;
  }

  if((ActCtr & NMDA_ACT)!=0)
  {
    if( BgSt[NID].pb[4]>BgSt[NID].uni_dist(*rdgen)) //NMDA
      BgSt[NID].SynSlo1ES.xt++;
    else
      BgSt[NID].SynSlo1ES.xt*=RPs[NID].ParNMDA.dfxt;
    BgSt[NID].SynSlo1ES.st = RPs[NID].ParNMDA.dfst*BgSt[NID].SynSlo1ES.st + RPs[NID].ParNMDA.efst*RPs[NID].ParNMDA.as*BgSt[NID].SynSlo1ES.xt*(1.0f - BgSt[NID].SynSlo1ES.st);
    RVs[NID].gNMDA += BgSt[NID].SynSlo1ES.g * BgSt[NID].SynSlo1ES.st;
  }

  if((ActCtr & INJ_ACT)!=0)
  { NetVBse[NID].Isyn += zig[thid].Igauss(BgSt[NID].mean,BgSt[NID].std); }

}

void EqsGnl2::VarInit()
{

//neuron body parameters
  for(unsigned int i=0;i<NetSize;i++)
  {
//Base
    NetVBse[i].Isyn=0.0f;  // current of synapse
    NetVBse[i].V=VRESTING;  // menbrane potential

//LIF
    NetVLIF[i].Count=0;  // firing tabe index
    NetVLIF[i].NState=0;

//STP
    if((Nets[i].ActCtr & STP_ACT) !=0 )
    {
      NetVSTP[i].D=1.0f;  // depression factor D
      NetVSTP[i].F=1.0f;  // facilitation factor F
      NetVSTP[i].Cap=0.0f; // Ca^2+ peak level
      NetVSTP[i].Car=20.0f;  // residual level of Ca^2+
    }
    else if((Nets[i].ActCtr & STD_ACT) !=0 )
      NetVSTP[i].D=1.0f;  // depression factor D

//LTP
    if((Nets[i].ActCtr & LTP_ACT) !=0 )
    {
      NetVLTP[i].LP=0.0f;   // long term plsiticity factor
      NetVLTP[i].CountLP=0;  //back propogation delay counter
    }


//gate variables
    RVs[i].gAMPA=0.0f;
    RVs[i].gGABA=0.0f;
    RVs[i].gNMDA=0.0f;
    RVs[i].gACH=0.0f;
    RVs[i].gGCL=0.0f;

//Ext stimulation
    for(unsigned int ix=0;ix<5;ix++)
    { Nets[i].Irep[ix]=0.0f; }

	  BgSt[i].SynSlo1ES.xt=0.0f;
	  BgSt[i].SynSlo1ES.st=0.0f;

    for(unsigned int j=0;j<4;j++)
    {   BgSt[i].SynFst1ES[j].st=0.0f;   }
  }

//synapses

  for(unsigned int i=0;i<SynGCL1End;i++)
  { SynFst1[i].st=0.0f; }

  for(unsigned int i=0;i<SynNMDA1Size;i++)
  {
    SynSlo1[i].xt=0.0f;
    SynSlo1[i].st=0.0f;
  }

}

//------------------------------post synaptic processing:st*g
void EqsGnl2::ThdSynapse(unsigned int thid)
{
  SynFast2 Fstmp;
  SynSlow2 Slotmp;
  float gTmp;

  for(unsigned int NID=thid;NID<AMPA1CntSize;NID+=thsize)
  {
    gTmp=RVs[SynFst1[NeuMapAMPA1[NID].Beg].PosID].gAMPA;
    for(unsigned int i=NeuMapAMPA1[NID].Beg;i<=NeuMapAMPA1[NID].End;i++)  //AMPA
    {
      Fstmp=SynFst1[i];

      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
        Fstmp.st++;
      else
        Fstmp.st *= RPs[Fstmp.PosID].ParAMPA.dfst;

      gTmp += Fstmp.st * Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapAMPA1[NID].Beg].PosID].gAMPA=gTmp;
  }

  for(unsigned int NID=thid;NID<GABA1CntSize;NID+=thsize)
  {
    gTmp=RVs[SynFst1[NeuMapGABA1[NID].Beg].PosID].gGABA;
    for(unsigned int i=NeuMapGABA1[NID].Beg;i<=NeuMapGABA1[NID].End;i++)  //GABA
    {
      Fstmp=SynFst1[i];

      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
        Fstmp.st++;
      else
        Fstmp.st *= RPs[Fstmp.PosID].ParGABA.dfst;

      gTmp += Fstmp.st * Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapGABA1[NID].Beg].PosID].gGABA=gTmp;
  }

  for(unsigned int NID=thid;NID<ACH1CntSize;NID+=thsize)
  {
    gTmp=RVs[SynFst1[NeuMapACH1[NID].Beg].PosID].gACH;
    for(unsigned int i=NeuMapACH1[NID].Beg;i<=NeuMapACH1[NID].End;i++)  //ACH
    {
      Fstmp=SynFst1[i];

      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
        Fstmp.st++;
      else
        Fstmp.st *= RPs[Fstmp.PosID].ParACH.dfst;

      gTmp += Fstmp.st * Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapACH1[NID].Beg].PosID].gACH=gTmp;
  }

  for(unsigned int NID=thid;NID<GCL1CntSize;NID+=thsize)
  {
    gTmp=RVs[SynFst1[NeuMapGCL1[NID].Beg].PosID].gGCL;
    for(unsigned int i=NeuMapGCL1[NID].Beg;i<=NeuMapGCL1[NID].End;i++)  //GCL
    {
      Fstmp=SynFst1[i];

      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
        Fstmp.st++;
      else
        Fstmp.st *= RPs[Fstmp.PosID].ParGCL.dfst;

      gTmp += Fstmp.st * Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapGCL1[NID].Beg].PosID].gGCL=gTmp;
  }

  for(unsigned int NID=thid;NID<NMDA1CntSize;NID+=thsize)
  {
    gTmp=RVs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].gNMDA;
    for(unsigned int i=NeuMapNMDA1[NID].Beg;i<=NeuMapNMDA1[NID].End;i++)  //NMDA
    {
      Slotmp=SynSlo1[i];
      if(NetVLIF[Slotmp.PreID].NState & SPIKE_F)
        Slotmp.xt++;
      else
        Slotmp.xt *= RPs[Slotmp.PosID].ParNMDA.dfxt;

      Slotmp.st = RPs[Slotmp.PosID].ParNMDA.dfst*Slotmp.st + RPs[Slotmp.PosID].ParNMDA.efst*RPs[Slotmp.PosID].ParNMDA.as*Slotmp.xt*(1.0f-Slotmp.st);

      gTmp += Slotmp.st * Slotmp.g;
      SynSlo1[i].xt=Slotmp.xt;
      SynSlo1[i].st=Slotmp.st;
    }
    RVs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].gNMDA=gTmp;
  }

}


void EqsGnl2::ThdSynapseSTP(unsigned int thid)
{
  SynFast2 Fstmp;
  SynSlow2 Slotmp;
  float gTmp;

  for(unsigned int NID=thid;NID<AMPA1CntSize;NID+=thsize)
  {
    gTmp=RVs[SynFst1[NeuMapAMPA1[NID].Beg].PosID].gAMPA;
    for(unsigned int i=NeuMapAMPA1[NID].Beg;i<=NeuMapAMPA1[NID].End;i++)  //AMPA
    {
      Fstmp=SynFst1[i];

      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
        Fstmp.st += NetVSTP[Fstmp.PreID].D*NetVSTP[Fstmp.PreID].F;
      else
        Fstmp.st *= RPs[Fstmp.PosID].ParAMPA.dfst;

      gTmp += Fstmp.st * Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapAMPA1[NID].Beg].PosID].gAMPA=gTmp;
  }

  for(unsigned int NID=thid;NID<GABA1CntSize;NID+=thsize)
  {
    gTmp=RVs[SynFst1[NeuMapGABA1[NID].Beg].PosID].gGABA;
    for(unsigned int i=NeuMapGABA1[NID].Beg;i<=NeuMapGABA1[NID].End;i++)  //GABA
    {
      Fstmp=SynFst1[i];

      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
        Fstmp.st += NetVSTP[Fstmp.PreID].D*NetVSTP[Fstmp.PreID].F;
      else
        Fstmp.st *= RPs[Fstmp.PosID].ParGABA.dfst;

      gTmp += Fstmp.st * Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapGABA1[NID].Beg].PosID].gGABA=gTmp;
  }

  for(unsigned int NID=thid;NID<ACH1CntSize;NID+=thsize)
  {
    gTmp=RVs[SynFst1[NeuMapACH1[NID].Beg].PosID].gACH;
    for(unsigned int i=NeuMapACH1[NID].Beg;i<=NeuMapACH1[NID].End;i++)  //ACH
    {
      Fstmp=SynFst1[i];

      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
        Fstmp.st += NetVSTP[Fstmp.PreID].D*NetVSTP[Fstmp.PreID].F;
      else
        Fstmp.st *= RPs[Fstmp.PosID].ParACH.dfst;

      gTmp += Fstmp.st * Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapACH1[NID].Beg].PosID].gACH=gTmp;
  }

  for(unsigned int NID=thid;NID<GCL1CntSize;NID+=thsize)
  {
    gTmp=RVs[SynFst1[NeuMapGCL1[NID].Beg].PosID].gGCL;
    for(unsigned int i=NeuMapGCL1[NID].Beg;i<=NeuMapGCL1[NID].End;i++)  //GCL
    {
      Fstmp=SynFst1[i];

      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
        Fstmp.st += NetVSTP[Fstmp.PreID].D*NetVSTP[Fstmp.PreID].F;
      else
        Fstmp.st *= RPs[Fstmp.PosID].ParGCL.dfst;

      gTmp += Fstmp.st * Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapGCL1[NID].Beg].PosID].gGCL=gTmp;
  }

  for(unsigned int NID=thid;NID<NMDA1CntSize;NID+=thsize)
  {
    gTmp=RVs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].gNMDA;
    for(unsigned int i=NeuMapNMDA1[NID].Beg;i<=NeuMapNMDA1[NID].End;i++)  //NMDA
    {
      Slotmp=SynSlo1[i];

      if(NetVLIF[Slotmp.PreID].NState & SPIKE_F)
        Slotmp.xt += NetVSTP[Slotmp.PreID].D*NetVSTP[Slotmp.PreID].F;
      else
        Slotmp.xt *= RPs[Slotmp.PosID].ParNMDA.dfxt;

      Slotmp.st = RPs[Slotmp.PosID].ParNMDA.dfst*Slotmp.st + RPs[Slotmp.PosID].ParNMDA.efst*RPs[Slotmp.PosID].ParNMDA.as*Slotmp.xt*(1.0f-Slotmp.st);

      gTmp += Slotmp.st * Slotmp.g;
      SynSlo1[i].xt=Slotmp.xt;
      SynSlo1[i].st=Slotmp.st;
    }
    RVs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].gNMDA=gTmp;
  }
}

void EqsGnl2::ThdSynapseSTD(unsigned int thid)
{
  SynFast2 Fstmp;
  SynSlow2 Slotmp;
  float gTmp;

  for(unsigned int NID=thid;NID<AMPA1CntSize;NID+=thsize)
  {
    gTmp=RVs[SynFst1[NeuMapAMPA1[NID].Beg].PosID].gAMPA;
    for(unsigned int i=NeuMapAMPA1[NID].Beg;i<=NeuMapAMPA1[NID].End;i++)  //AMPA
    {
      Fstmp=SynFst1[i];

      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
        Fstmp.st += NetVSTP[Fstmp.PreID].D;
      else
        Fstmp.st *= RPs[Fstmp.PosID].ParAMPA.dfst;

      gTmp += Fstmp.st * Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapAMPA1[NID].Beg].PosID].gAMPA=gTmp;
  }

  for(unsigned int NID=thid;NID<GABA1CntSize;NID+=thsize)
  {
    gTmp=RVs[SynFst1[NeuMapGABA1[NID].Beg].PosID].gGABA;
    for(unsigned int i=NeuMapGABA1[NID].Beg;i<=NeuMapGABA1[NID].End;i++)  //GABA
    {
      Fstmp=SynFst1[i];

      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
        Fstmp.st += NetVSTP[Fstmp.PreID].D;
      else
        Fstmp.st *= RPs[Fstmp.PosID].ParGABA.dfst;

      gTmp += Fstmp.st * Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapGABA1[NID].Beg].PosID].gGABA=gTmp;
  }

  for(unsigned int NID=thid;NID<ACH1CntSize;NID+=thsize)
  {
    gTmp=RVs[SynFst1[NeuMapACH1[NID].Beg].PosID].gACH;
    for(unsigned int i=NeuMapACH1[NID].Beg;i<=NeuMapACH1[NID].End;i++)  //ACH
    {
      Fstmp=SynFst1[i];

      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
        Fstmp.st += NetVSTP[Fstmp.PreID].D;
      else
        Fstmp.st *= RPs[Fstmp.PosID].ParACH.dfst;

      gTmp += Fstmp.st * Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapACH1[NID].Beg].PosID].gACH=gTmp;
  }

  for(unsigned int NID=thid;NID<GCL1CntSize;NID+=thsize)
  {
    gTmp=RVs[SynFst1[NeuMapGCL1[NID].Beg].PosID].gGCL;
    for(unsigned int i=NeuMapGCL1[NID].Beg;i<=NeuMapGCL1[NID].End;i++)  //GCL
    {
      Fstmp=SynFst1[i];

      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
        Fstmp.st += NetVSTP[Fstmp.PreID].D;
      else
        Fstmp.st *= RPs[Fstmp.PosID].ParGCL.dfst;

      gTmp += Fstmp.st * Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapGCL1[NID].Beg].PosID].gGCL=gTmp;
  }

  for(unsigned int NID=thid;NID<NMDA1CntSize;NID+=thsize)
  {
    gTmp=RVs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].gNMDA;
    for(unsigned int i=NeuMapNMDA1[NID].Beg;i<=NeuMapNMDA1[NID].End;i++)  //NMDA
    {
      Slotmp=SynSlo1[i];

      if(NetVLIF[Slotmp.PreID].NState & SPIKE_F)
        Slotmp.xt += NetVSTP[Slotmp.PreID].D;
      else
        Slotmp.xt *= RPs[Slotmp.PosID].ParNMDA.dfxt;

      Slotmp.st = RPs[Slotmp.PosID].ParNMDA.dfst*Slotmp.st + RPs[Slotmp.PosID].ParNMDA.efst*RPs[Slotmp.PosID].ParNMDA.as*Slotmp.xt*(1.0f-Slotmp.st);

      gTmp += Slotmp.st * Slotmp.g;
      SynSlo1[i].xt=Slotmp.xt;
      SynSlo1[i].st=Slotmp.st;
    }
    RVs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].gNMDA=gTmp;
  }
}

void EqsGnl2::ThdGapSyn(unsigned int thid)
{
  float gTmp,vTmp;
  for(unsigned int NID=thid;NID<Cnd1CntSize;NID+=thsize)
  {
    gTmp=NetVBse[SynCnd1[NeuMapCnd1[NID].Beg].PosID].Isyn;
    vTmp=NetVBse[SynCnd1[NeuMapCnd1[NID].Beg].PosID].V;
    for(unsigned int i=NeuMapCnd1[NID].Beg;i<=NeuMapCnd1[NID].End;i++)  //gap junction
    { gTmp += SynCnd1[i].g*(NetVBse[SynCnd1[i].PreID].V-vTmp);  }
    NetVBse[SynCnd1[NeuMapCnd1[NID].Beg].PosID].Isyn=gTmp;
  }

  for(unsigned int NID=thid;NID<PreCnd1CntSize;NID+=thsize)
  {
    gTmp=NetVBse[SynCnd1[PreMapCnd1[NeuPreMapCnd1[NID].Beg]].PreID].Isyn;
    vTmp=NetVBse[SynCnd1[PreMapCnd1[NeuPreMapCnd1[NID].Beg]].PreID].V;
    for(unsigned int i=NeuPreMapCnd1[NID].Beg;i<=NeuPreMapCnd1[NID].End;i++)  //gap junction
    { gTmp -= SynCnd1[PreMapCnd1[i]].g*(vTmp-NetVBse[SynCnd1[PreMapCnd1[i]].PosID].V);  }
    NetVBse[SynCnd1[PreMapCnd1[NeuPreMapCnd1[NID].Beg]].PreID].Isyn=gTmp;
  }
}

//--------------------------------------------------
void EqsGnl2::ThdPGatVarSTP(unsigned int thid)
{
  SynFast2 Fstmp;
  SynSlow2 Slotmp;

  for(unsigned int i=thid;i<PGatASize;i+=thsize)  //AMPA
  {
    Fstmp=PSynFstGatA[i];

    if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
      Fstmp.st += NetVSTP[Fstmp.PreID].D*NetVSTP[Fstmp.PreID].F;
    else
      Fstmp.st *= RPs[Fstmp.PosID].ParAMPA.dfst;

    PSynFstGatA[i]=Fstmp;
  }

  for(unsigned int i=thid;i<PGatGSize;i+=thsize)  //GABA
  {
    Fstmp=PSynFstGatG[i];

    if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
      Fstmp.st += NetVSTP[Fstmp.PreID].D*NetVSTP[Fstmp.PreID].F;
    else
      Fstmp.st *= RPs[Fstmp.PosID].ParGABA.dfst;

    PSynFstGatG[i]=Fstmp;
  }

  for(unsigned int i=thid;i<PGatHSize;i+=thsize)  //ACH
  {
    Fstmp=PSynFstGatH[i];

    if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
      Fstmp.st += NetVSTP[Fstmp.PreID].D*NetVSTP[Fstmp.PreID].F;
    else
      Fstmp.st *= RPs[Fstmp.PosID].ParACH.dfst;

    PSynFstGatH[i]=Fstmp;
  }

  for(unsigned int i=thid;i<PGatLSize;i+=thsize)  //GCL
  {
    Fstmp=PSynFstGatL[i];

    if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
      Fstmp.st += NetVSTP[Fstmp.PreID].D*NetVSTP[Fstmp.PreID].F;
    else
      Fstmp.st *= RPs[Fstmp.PosID].ParGCL.dfst;

    PSynFstGatL[i]=Fstmp;
  }

  for(unsigned int i=thid;i<PGatNSize;i+=thsize)  //NMDA
  {
    Slotmp=PSynSloGat[i];

    if(NetVLIF[Slotmp.PreID].NState & SPIKE_F)
      Slotmp.xt += NetVSTP[Slotmp.PreID].D*NetVSTP[Slotmp.PreID].F;
    else
      Slotmp.xt *= RPs[Slotmp.PosID].ParNMDA.dfxt;

    Slotmp.st = RPs[Slotmp.PosID].ParNMDA.dfst*Slotmp.st + RPs[Slotmp.PosID].ParNMDA.efst*RPs[Slotmp.PosID].ParNMDA.as*Slotmp.xt*(1.0f-Slotmp.st);

    PSynSloGat[i]=Slotmp;
  }
}


void EqsGnl2::ThdPGatVarSTD(unsigned int thid)
{
  SynFast2 Fstmp;
  SynSlow2 Slotmp;

  for(unsigned int i=thid;i<PGatASize;i+=thsize)  //AMPA
  {
    Fstmp=PSynFstGatA[i];

    if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
      Fstmp.st += NetVSTP[Fstmp.PreID].D;
    else
      Fstmp.st *= RPs[Fstmp.PosID].ParAMPA.dfst;

    PSynFstGatA[i]=Fstmp;
  }

  for(unsigned int i=thid;i<PGatGSize;i+=thsize)  //GABA
  {
    Fstmp=PSynFstGatG[i];

    if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
      Fstmp.st += NetVSTP[Fstmp.PreID].D;
    else
      Fstmp.st *= RPs[Fstmp.PosID].ParGABA.dfst;

    PSynFstGatG[i]=Fstmp;
  }

  for(unsigned int i=thid;i<PGatHSize;i+=thsize)  //ACH
  {
    Fstmp=PSynFstGatH[i];

    if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
      Fstmp.st += NetVSTP[Fstmp.PreID].D;
    else
      Fstmp.st *= RPs[Fstmp.PosID].ParACH.dfst;

    PSynFstGatH[i]=Fstmp;
  }

  for(unsigned int i=thid;i<PGatLSize;i+=thsize)  //GCL
  {
    Fstmp=PSynFstGatL[i];

    if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
      Fstmp.st += NetVSTP[Fstmp.PreID].D;
    else
      Fstmp.st *= RPs[Fstmp.PosID].ParGCL.dfst;

    PSynFstGatL[i]=Fstmp;
  }

  for(unsigned int i=thid;i<PGatNSize;i+=thsize)  //NMDA
  {
    Slotmp=PSynSloGat[i];

    if(NetVLIF[Slotmp.PreID].NState & SPIKE_F)
      Slotmp.xt += NetVSTP[Slotmp.PreID].D*NetVSTP[Slotmp.PreID].F;
    else
      Slotmp.xt *= RPs[Slotmp.PosID].ParNMDA.dfxt;

    Slotmp.st = RPs[Slotmp.PosID].ParNMDA.dfst*Slotmp.st + RPs[Slotmp.PosID].ParNMDA.efst*RPs[Slotmp.PosID].ParNMDA.as*Slotmp.xt*(1.0f-Slotmp.st);

    PSynSloGat[i]=Slotmp;
  }
}


void EqsGnl2::ThdPGatVar(unsigned int thid)
{
  SynFast2 Fstmp;
  SynSlow2 Slotmp;

  for(unsigned int i=thid;i<PGatASize;i+=thsize)  //AMPA
  {
    Fstmp=PSynFstGatA[i];

    if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
      Fstmp.st++;
    else
      Fstmp.st *= RPs[Fstmp.PosID].ParAMPA.dfst;

    PSynFstGatA[i]=Fstmp;
  }

  for(unsigned int i=thid;i<PGatGSize;i+=thsize)  //GABA
  {
    Fstmp=PSynFstGatG[i];

    if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
      Fstmp.st++;
    else
      Fstmp.st *= RPs[Fstmp.PosID].ParGABA.dfst;

    PSynFstGatG[i]=Fstmp;
  }

  for(unsigned int i=thid;i<PGatHSize;i+=thsize)  //ACH
  {
    Fstmp=PSynFstGatH[i];

    if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
      Fstmp.st++;
    else
      Fstmp.st *= RPs[Fstmp.PosID].ParACH.dfst;

    PSynFstGatH[i]=Fstmp;
  }

  for(unsigned int i=thid;i<PGatLSize;i+=thsize)  //GCL
  {
    Fstmp=PSynFstGatL[i];

    if(NetVLIF[Fstmp.PreID].NState & SPIKE_F)
      Fstmp.st++;
    else
      Fstmp.st *= RPs[Fstmp.PosID].ParGCL.dfst;

    PSynFstGatL[i]=Fstmp;
  }

  for(unsigned int i=thid;i<PGatNSize;i+=thsize)  //NMDA
  {
    Slotmp=PSynSloGat[i];

    if(NetVLIF[Slotmp.PreID].NState & SPIKE_F)
      Slotmp.xt++;
    else
      Slotmp.xt *= RPs[Slotmp.PosID].ParNMDA.dfxt;

    Slotmp.st = RPs[Slotmp.PosID].ParNMDA.dfst*Slotmp.st + RPs[Slotmp.PosID].ParNMDA.efst*RPs[Slotmp.PosID].ParNMDA.as*Slotmp.xt*(1.0f-Slotmp.st);

    PSynSloGat[i]=Slotmp;
  }
}


//------------------------------group connection post synaptic processing:st*g

void EqsGnl2::ThdPSynapse(unsigned int thid)
{
  float gTmp,gTmp2,vTmp;

  for(unsigned int NID=thid;NID<PAMPACntSize;NID+=thsize)
  {
    gTmp=RVs[MapPAPosID[NID]].gAMPA;
    for(unsigned int i=NeuMapPAMPA[NID].Beg;i<=NeuMapPAMPA[NID].End;i++)  //AMPA
    { gTmp += PSynFstGatA[PSynA[i]].st * PSynFstGatA[PSynA[i]].g; }
    RVs[MapPAPosID[NID]].gAMPA=gTmp;
  }

  for(unsigned int NID=thid;NID<PGABACntSize;NID+=thsize)
  {
    gTmp=RVs[MapPGPosID[NID]].gGABA;
    for(unsigned int i=NeuMapPGABA[NID].Beg;i<=NeuMapPGABA[NID].End;i++)  //GABA
    { gTmp += PSynFstGatG[PSynG[i]].st * PSynFstGatG[PSynG[i]].g; }
    RVs[MapPGPosID[NID]].gGABA=gTmp;
  }

  for(unsigned int NID=thid;NID<PACHCntSize;NID+=thsize)
  {
    gTmp=RVs[MapPHPosID[NID]].gACH;
    for(unsigned int i=NeuMapPACH[NID].Beg;i<=NeuMapPACH[NID].End;i++)  //ACH
    { gTmp += PSynFstGatH[PSynH[i]].st * PSynFstGatH[PSynH[i]].g; }
    RVs[MapPHPosID[NID]].gACH=gTmp;
  }

  for(unsigned int NID=thid;NID<PGCLCntSize;NID+=thsize)
  {
    gTmp=RVs[MapPLPosID[NID]].gGCL;
    for(unsigned int i=NeuMapPGCL[NID].Beg;i<=NeuMapPGCL[NID].End;i++)  //GCL
    { gTmp += PSynFstGatL[PSynL[i]].st * PSynFstGatL[PSynL[i]].g; }
    RVs[MapPLPosID[NID]].gGCL=gTmp;
  }

  for(unsigned int NID=thid;NID<PNMDACntSize;NID+=thsize)
  {
    gTmp=RVs[MapPNPosID[NID]].gNMDA;
    for(unsigned int i=NeuMapPNMDA[NID].Beg;i<=NeuMapPNMDA[NID].End;i++)  //NMDA
    { gTmp += PSynSloGat[PSynN[i]].st * PSynSloGat[PSynN[i]].g; }
    RVs[MapPNPosID[NID]].gNMDA=gTmp;
  }


  for(unsigned int NID=thid;NID<PCndPCntSize;NID+=thsize) //integrate PreSID current to PosID
  {
    gTmp=NetVBse[MapPCPosID[NID]].Isyn;
    vTmp=NetVBse[MapPCPosID[NID]].V;
    for(unsigned int i=NeuMapPCndP[NID].Beg;i<=NeuMapPCndP[NID].End;i++)
    { gTmp += PSynCndP[PSynCp[i]].g*(NetVBse[PSynCndP[PSynCp[i]].PreID].V-vTmp); }
    NetVBse[MapPCPosID[NID]].Isyn=gTmp;
  }

  for(unsigned int NID=thid;NID<PCndOCntSize;NID+=thsize) //integrate PosNID current to PreID
  {
    gTmp=NetVBse[PSynCndP[NID].PreID].Isyn;
    gTmp2=PSynCndP[NID].g;
    vTmp=NetVBse[PSynCndP[NID].PreID].V;
    for(unsigned int i=NeuMapPCndO[NID].Beg;i<=NeuMapPCndO[NID].End;i++)
    { gTmp -= gTmp2*(vTmp-NetVBse[PSynCo[i]].V); }
    NetVBse[PSynCndP[NID].PreID].Isyn=gTmp;
  }


}

//-------------Long term plasticity
void EqsGnl2::LTP(unsigned int thid)
{
  for(unsigned int i=thid;i<SynGCL1End;i+=thsize)  //AMPA
  {
    if(NetVLIF[SynFst1[i].PreID].NState & SPIKE_F)
    {
      if(NetPLTP[SynFst1[i].PreID].NegISI>0)  //facilitation
      { SynFst1[i].g +=  NetVLTP[SynFst1[i].PosID].LP * NetPLTP[SynFst1[i].PreID].NegISI; } //additive rule
      else  //depression
      { SynFst1[i].g *= (1.0f + NetVLTP[SynFst1[i].PosID].LP * NetPLTP[SynFst1[i].PreID].NegISI); } //multiplicative rule
    }

    if(NetVLIF[SynFst1[i].PosID].NState & SPIKE_F)
    {
      if(NetPLTP[SynFst1[i].PreID].PosISI>0)  //facilitation
      { SynFst1[i].g +=  NetVLTP[SynFst1[i].PreID].LP * NetPLTP[SynFst1[i].PreID].PosISI; } //additive rule
      else  //depression
      { SynFst1[i].g *= (1.0f + NetVLTP[SynFst1[i].PreID].LP * NetPLTP[SynFst1[i].PreID].PosISI); } //multiplicative rule
    }
  }

  for(unsigned int i=thid;i<SynNMDA1Size;i+=thsize)  //NMDA
  {
    if(NetVLIF[SynSlo1[i].PreID].NState & SPIKE_F)
    {
      if(NetPLTP[SynSlo1[i].PreID].NegISI>0)  //facilitation
      { SynSlo1[i].g +=  NetVLTP[SynSlo1[i].PosID].LP * NetPLTP[SynSlo1[i].PreID].NegISI; } //additive rule
      else  //depression
      { SynSlo1[i].g *= (1.0f + NetVLTP[SynSlo1[i].PosID].LP * NetPLTP[SynSlo1[i].PreID].NegISI); } //multiplicative rule
    }

    if(NetVLIF[SynSlo1[i].PosID].NState & SPIKE_F)
    {
      if(NetPLTP[SynSlo1[i].PreID].PosISI>0)  //facilitation
      { SynSlo1[i].g +=  NetVLTP[SynSlo1[i].PreID].LP * NetPLTP[SynSlo1[i].PreID].PosISI; } //additive rule
      else  //depression
      { SynSlo1[i].g *= (1.0f + NetVLTP[SynSlo1[i].PreID].LP * NetPLTP[SynSlo1[i].PreID].PosISI); } //multiplicative rule
    }
  }
}


//-------------short term plasticity
void EqsGnl2::STP(unsigned int thid)
{
  ParStp2 Pstp;
  VarStp Vstp;
  for(unsigned NID=thid;NID<NetSize;NID+=thsize)
  {
	  Pstp=NetPSTP[NID];
	  Vstp=NetVSTP[NID];

    if(NetVLIF[NID].NState & SPIKE_F)
    {
      Vstp.Car += Pstp.aCar;
      Vstp.Cap += Pstp.aCap;
      Vstp.F += (Vstp.Cap+Vstp.Car)*(1.0f - Vstp.F)*Pstp.aF;
      Vstp.D *= ( 1.0f -  Pstp.pv * Vstp.F);
    }
    else
    {
      Vstp.Car = Pstp.dfCar*Vstp.Car+Pstp.efCar;
      Vstp.Cap = Pstp.dfCap*Vstp.Cap;
      Vstp.F = Pstp.dfF*Vstp.F + Pstp.efF*(Vstp.Cap+Vstp.Car)*(1.0f - Vstp.F)*Pstp.aF;
      Vstp.D = Pstp.dfD*Vstp.D + Pstp.efD;
    }
	  NetVSTP[NID]=Vstp;
  }

}


void EqsGnl2::STD(unsigned int thid)
{
  ParStp2 Pstp;
  for(unsigned NID=thid;NID<NetSize;NID+=thsize)
  {
    Pstp=NetPSTP[NID];

    if(NetVLIF[NID].NState & SPIKE_F)
      NetVSTP[NID].D *= ( 1.0f -  Pstp.pv);
    else
      NetVSTP[NID].D = Pstp.dfD*NetVSTP[NID].D + Pstp.efD;

  }
}
//--------------------------------------------------------------------
void EqsGnl2::EqsOne(unsigned int ActCtr)
{
  Activation(0,ActCtr);
  Spike(0,ActCtr);
}

void EqsGnl2::Activation(unsigned int thid,unsigned int ActCtr)
{
  for(unsigned NID=thid;NID<NetSize;NID+=thsize)
  {
    BgInit(thid,NID); //bcak ground initial
    VmEul(NID);
  }
}

void EqsGnl2::Spike(unsigned int thid,unsigned int ActCtr)
{
  if((ActCtr & STP_ACT)!=0)
    ThdPGatVarSTP(thid);
  else if((ActCtr & STD_ACT)!=0)
    ThdPGatVarSTD(thid);
  else
    ThdPGatVar(thid);

  ThdPSynapse(thid);

  if((ActCtr & STP_ACT)!=0)
    ThdSynapseSTP(thid);
  else if((ActCtr & STD_ACT)!=0)
    ThdSynapseSTD(thid);
  else
    ThdSynapse(thid);

  ThdGapSyn(thid);

  if((ActCtr & STP_ACT)!=0) STP(thid);
  else if((ActCtr & STD_ACT)!=0) STD(thid);

  if((ActCtr & LTP_ACT)!=0) LTP(thid);
}



