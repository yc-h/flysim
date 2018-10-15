
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <thread>
#include <atomic>
#include <condition_variable>
#include "Eqs.h"

using namespace std;

//****************************************************************************80

//global variables

EqsGnl1_1::EqsGnl1_1()
{
  dt=DTIME;
  thsize=1;

  NetSize=0;
  RPsSize=0;
  RVsSize=0;
  RGsSize=0;
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
  //RGs=nullptr;
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
float EqsGnl1_1::dst(float st,float tst) //delta st
{ return (-st / tst); }

float EqsGnl1_1::ndst(float st,float tst,float xt,float as) //NMDA receptor delta st
{ return (as*xt*(1.0f - st) - st/tst);}

float EqsGnl1_1::dD(float D, float tD)
{ return (1.0f-D)/tD; }

float EqsGnl1_1::deltaV(float Vm, unsigned int NID)
{

  NeuRepPar2<RP1,RP2> NRP=RPs[NID];
  NeuRepVar2 NRV=RVs[NID];
  ParBase NPB=NetPBse[NID];
  return
  (
    (
      NRV.gAMPA*(NRP.ParAMPA.Vrev-Vm)
      +NRV.gGABA*(NRP.ParGABA.Vrev-Vm)
      +NRV.gNMDA*(NRP.ParNMDA.Vrev-Vm)/( 1.0f + NRP.ParNMDA.mg * (float)exp( -0.062f * Vm ) / 3.57f )
      +NRV.gACH*(NRP.ParACH.Vrev-Vm)
      +NRV.gGCL*(NRP.ParGCL.Vrev-Vm)
    )*1e-3f
    +NetVBse[NID].Isyn
    +NPB.gl*(NPB.El-Vm)
  )/NPB.Cm;
}

void EqsGnl1_1::MemPot(unsigned int NID)
{
  NeuRepPar2<RP1,RP2> NRP=RPs[NID];
  NeuRepVar2 NRV=RVs[NID];
  ParBase NPB=NetPBse[NID];
  float dv,Vm=NetVBse[NID].V;

  if(NetVLIF[NID].NState & FIRE_F)
    firing(NID);
  else
  {
    switch(SolType)
    {

      case ACCURATE:
      dv = RK4FuI(Vm,NID, &EqsGnl1_1::deltaV);
      break;

      case ROUGH:
      dv = deltaV(Vm,NID);
      break;

      case MODERATE:
      default:
      dv = IpvEulFuI(Vm,NID, &EqsGnl1_1::deltaV);
      break;
    }

    if(dv > NPB.dVmax)
    {
      dv=NPB.dVmax;
      NetVLIF[NID].NState |= MEMDV_OVF; //active MEMDV_OVF
      NetVLIF[NID].NState &= ~MEMDV_UDF; //deactive MEMDV_UDF
    }
    else if(dv < -NPB.dVmax)
    {
      dv=-NPB.dVmax;
      NetVLIF[NID].NState &= ~MEMDV_OVF; //deactive MEMDV_OVF
      NetVLIF[NID].NState |= MEMDV_UDF; //active MEMDV_UDF
    }
    else
    {
      NetVLIF[NID].NState &= ~MEMDV_OVF; //deactive MEMDV_OVF
      NetVLIF[NID].NState &= ~MEMDV_UDF; //deactive MEMDV_UDF
    }
    
    Vm += dt*dv;
    NetVBse[NID].V=Vm;
    if(Vm > NPB.Vth)
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

  Nets[NID].Irep[0]=NRV.gAMPA* (NRP.ParAMPA.Vrev-Vm);
  Nets[NID].Irep[1]=NRV.gGABA* (NRP.ParGABA.Vrev-Vm);
  Nets[NID].Irep[2]=NRV.gNMDA*(NRP.ParNMDA.Vrev-Vm)/( 1.0f + NRP.ParNMDA.mg * (float)exp( -0.062f * Vm ) / 3.57f );
  Nets[NID].Irep[3]=NRV.gACH* (NRP.ParACH.Vrev-Vm);
  Nets[NID].Irep[4]=NRV.gGCL* (NRP.ParGCL.Vrev-Vm);

  NRV.gAMPA=0.0f;
  NRV.gGABA=0.0f;
  NRV.gNMDA=0.0f;
  NRV.gACH=0.0f;
  NRV.gGCL=0.0f;
  NetVBse[NID].Isyn=0.0f;
  RVs[NID]=NRV;
}


//--------------------firing process
void EqsGnl1_1::firing(unsigned int NID)
{
  if(NetVLIF[NID].Count < NetPLIF[NID].RfPed) //Firing
  {
    if(NetVLIF[NID].Count == 0)
      NetVBse[NID].V=0.0f;
    else
      NetVBse[NID].V=NetPBse[NID].Vreset;

    if(NetVLIF[NID].Count==NetPLIF[NID].SpiDy)
    {
      NetVLIF[NID].NState |= SPIKE_F;
    }
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

void EqsGnl1_1::firingHH(unsigned int NID)
{
  if(NetVLIF[NID].Count < NetPLIF[NID].RfPed) //Firing
  {
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
void EqsGnl1_1::firLP(unsigned int NID)
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


//-------------background initialization
void EqsGnl1_1::BgInit(unsigned int thid, unsigned int NID)
{
  unsigned int ActCtr=Nets[NID].ActCtr;
  
  //AMPA
  if((ActCtr & AMPA_ACT)!=0)
  {
    BgSt[NID].SynFst1ES[0].st += ((BgSt[NID].pb[0]>BgSt[NID].uni_dist(*rdgen))?1.0f : -dt*BgSt[NID].SynFst1ES[0].st/RPs[NID].ParAMPA.tst);
    RVs[NID].gAMPA += BgSt[NID].SynFst1ES[0].g*BgSt[NID].SynFst1ES[0].st;
  }

  //GABA
  if((ActCtr & GABA_ACT)!=0)
  {
    BgSt[NID].SynFst1ES[1].st += ((BgSt[NID].pb[1]>BgSt[NID].uni_dist(*rdgen))?1.0f : -dt*BgSt[NID].SynFst1ES[1].st/RPs[NID].ParGABA.tst);
    RVs[NID].gGABA += BgSt[NID].SynFst1ES[1].g*BgSt[NID].SynFst1ES[1].st;
  }

  //ACH
  if((ActCtr & ACH_ACT)!=0)
  {
    BgSt[NID].SynFst1ES[2].st += ((BgSt[NID].pb[2]>BgSt[NID].uni_dist(*rdgen))?1.0f : -dt*BgSt[NID].SynFst1ES[2].st/RPs[NID].ParACH.tst);
    RVs[NID].gACH += BgSt[NID].SynFst1ES[2].g*BgSt[NID].SynFst1ES[2].st;
  }

  //GCL
  if((ActCtr & GCL_ACT)!=0)
  {
    BgSt[NID].SynFst1ES[3].st += ((BgSt[NID].pb[3]>BgSt[NID].uni_dist(*rdgen))?1.0f : -dt*BgSt[NID].SynFst1ES[3].st/RPs[NID].ParGCL.tst);
    RVs[NID].gGCL += BgSt[NID].SynFst1ES[3].g*BgSt[NID].SynFst1ES[3].st;
  }

  //NMDA
  if((ActCtr & NMDA_ACT)!=0)
  {
    BgSt[NID].SynSlo1ES.xt += ((BgSt[NID].pb[4]>BgSt[NID].uni_dist(*rdgen))?1.0f:-dt*(BgSt[NID].SynSlo1ES.xt/RPs[NID].ParNMDA.txt));
    BgSt[NID].SynSlo1ES.st += dt*(RPs[NID].ParNMDA.as*BgSt[NID].SynSlo1ES.xt*(1.0f-BgSt[NID].SynSlo1ES.st)-BgSt[NID].SynSlo1ES.st/RPs[NID].ParNMDA.tst);
    RVs[NID].gNMDA += BgSt[NID].SynSlo1ES.g * BgSt[NID].SynSlo1ES.st;
  }

  if((ActCtr & INJ_ACT)!=0)
  { NetVBse[NID].Isyn += zig[thid].Igauss(BgSt[NID].mean,BgSt[NID].std); }

}

void EqsGnl1_1::VarInit()
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

  for(unsigned int i=0;i<PGatASize;i++)  //AMPA
  { PSynFstGatA[i].st=0.0f; }

  for(unsigned int i=0;i<PGatGSize;i++)  //GABA
  { PSynFstGatG[i].st=0.0f; }

  for(unsigned int i=0;i<PGatHSize;i++)  //ACH
  { PSynFstGatH[i].st=0.0f; }

  for(unsigned int i=0;i<PGatLSize;i++)  //GCL
  { PSynFstGatL[i].st=0.0f; }

  for(unsigned int i=0;i<PGatNSize;i++)  //NMDA
  {
    PSynSloGat[i].xt=0.0f;
    PSynSloGat[i].st=0.0f;
  }
}


void EqsGnl1_1::ThdGapSyn(unsigned int thid)
{
  float gTmp,vTmp;
  for(unsigned int NID=thid;NID<Cnd1CntSize;NID+=thsize)
  {
    gTmp=NetVBse[SynCnd1[NeuMapCnd1[NID].Beg].PosID].Isyn;
    vTmp=NetVBse[SynCnd1[NeuMapCnd1[NID].Beg].PosID].V;
    for(unsigned int i=NeuMapCnd1[NID].Beg;i<=NeuMapCnd1[NID].End;i++)  //integrate PreSID current to PosID
    { gTmp += SynCnd1[i].g*(NetVBse[SynCnd1[i].PreID].V-vTmp);  }
    NetVBse[SynCnd1[NeuMapCnd1[NID].Beg].PosID].Isyn=gTmp;
  }

  for(unsigned int NID=thid;NID<PreCnd1CntSize;NID+=thsize)
  {
    gTmp=NetVBse[SynCnd1[PreMapCnd1[NeuPreMapCnd1[NID].Beg]].PreID].Isyn;
    vTmp=NetVBse[SynCnd1[PreMapCnd1[NeuPreMapCnd1[NID].Beg]].PreID].V;
    for(unsigned int i=NeuPreMapCnd1[NID].Beg;i<=NeuPreMapCnd1[NID].End;i++)  //integrate PosNID current to PreID
    { gTmp -= SynCnd1[PreMapCnd1[i]].g*(vTmp-NetVBse[SynCnd1[PreMapCnd1[i]].PosID].V);  }
    NetVBse[SynCnd1[PreMapCnd1[NeuPreMapCnd1[NID].Beg]].PreID].Isyn=gTmp;
  }
}



//------------------------------group connection post synaptic processing:st*g
void EqsGnl1_1::ThdPSynapse(unsigned int thid)
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


  for(unsigned int NID=thid;NID<PCndPCntSize;NID+=thsize) //gap junction: integrate PreSID current to PosID
  {
    gTmp=NetVBse[MapPCPosID[NID]].Isyn;
    vTmp=NetVBse[MapPCPosID[NID]].V;
    for(unsigned int i=NeuMapPCndP[NID].Beg;i<=NeuMapPCndP[NID].End;i++)
    { gTmp += PSynCndP[PSynCp[i]].g*(NetVBse[PSynCndP[PSynCp[i]].PreID].V-vTmp); }
    NetVBse[MapPCPosID[NID]].Isyn=gTmp;
  }

  for(unsigned int NID=thid;NID<PCndOCntSize;NID+=thsize) //gap junction: integrate PosNID current to PreID
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
void EqsGnl1_1::LTP(unsigned int thid)
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



//--------------------------------------------------------------------
void EqsGnl1_1::EqsOne(unsigned int ActCtr)
{
  Activation(0,ActCtr);
  Spike(0,ActCtr);
}


void EqsGnl1_1::Activation(unsigned int thid,unsigned int ActCtr)
{
  
  if((ActCtr & HH_ACT)!=0) {}
  else if((ActCtr & SOD_ACT)!=0) {}
  else
  {
    for(unsigned NID=thid;NID<NetSize;NID+=thsize)
    {
      BgInit(thid,NID); //bcak ground initial
      MemPot(NID);
    }
  }
}

void EqsGnl1_1::Spike(unsigned int thid,unsigned int ActCtr)
{
  switch(SolType)
  {
    case ACCURATE:
    if((ActCtr & STP_ACT)!=0)
      ThdPGatVarSTPRK4(thid);
    else if((ActCtr & STD_ACT)!=0)
      ThdPGatVarSTDRK4(thid);
    else
      ThdPGatVarRK4(thid);

    ThdPSynapse(thid);

    if((ActCtr & STP_ACT)!=0)
      ThdSynapseSTPRK4(thid);
    else if((ActCtr & STD_ACT)!=0)
      ThdSynapseSTDRK4(thid);
    else
      ThdSynapseRK4(thid);

    ThdGapSyn(thid);

    if((ActCtr & STP_ACT)!=0) STPRK4(thid);
    else if((ActCtr & STD_ACT)!=0) STDRK4(thid);

    if((ActCtr & LTP_ACT)!=0) LTP(thid);
    break;

    case ROUGH:
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
    break;

    case MODERATE:
    default:
    if((ActCtr & STP_ACT)!=0)
      ThdPGatVarSTPIpvEul(thid);
    else if((ActCtr & STD_ACT)!=0)
      ThdPGatVarSTDIpvEul(thid);
    else
      ThdPGatVarIpvEul(thid);

    ThdPSynapse(thid);

    if((ActCtr & STP_ACT)!=0)
      ThdSynapseSTPIpvEul(thid);
    else if((ActCtr & STD_ACT)!=0)
      ThdSynapseSTDIpvEul(thid);
    else
      ThdSynapseIpvEul(thid);

    ThdGapSyn(thid);

    if((ActCtr & STP_ACT)!=0) STPIpvEul(thid);
    else if((ActCtr & STD_ACT)!=0) STDIpvEul(thid);

    if((ActCtr & LTP_ACT)!=0) LTP(thid);

    break;

  }
}



