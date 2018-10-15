
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

EqsGnl2_1::EqsGnl2_1()
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

  th=nullptr;
  Nets=nullptr;
  RPs=nullptr;
  RVs=nullptr;
  RGs=nullptr;
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
}

//--------------------membrane equation
void EqsGnl2_1::VmEul(unsigned int NID)
{
  NeuRepPar2<RP3,RP4> NRP= RPs[NID];
  NeuRepVar2 NRV=RVs[NID];
  NeuGatVar3 NRG=RGs[NID];
  ParBase2 NPB=NetPBse[NID];
  float Vm=NetVBse[NID].V;

  NRG.stA = NRG.stA*NRP.ParAMPA.dfst+NRV.gAMPA;
  NRG.stG = NRG.stG*NRP.ParGABA.dfst+NRV.gGABA;
  NRG.stH = NRG.stH*NRP.ParACH.dfst+NRV.gACH;
  NRG.stL = NRG.stL*NRP.ParGCL.dfst+NRV.gGCL;

  Nets[NID].Irep[0]=NRG.stA*(NRP.ParAMPA.Vrev-Vm);
  Nets[NID].Irep[1]=NRG.stG*(NRP.ParGABA.Vrev-Vm);
  Nets[NID].Irep[2]=NRV.gNMDA*(NRP.ParNMDA.Vrev-Vm)/( 1.0f + NRP.ParNMDA.mg * expApx( -0.062f * Vm ) / 3.57f );
  Nets[NID].Irep[3]=NRG.stH*(NRP.ParACH.Vrev-Vm);
  Nets[NID].Irep[4]=NRG.stL*(NRP.ParGCL.Vrev-Vm);

  if(NetVLIF[NID].NState & FIRE_F)
    firing(NID);
  else
  {
    float dv=(
    (Nets[NID].Irep[0]+Nets[NID].Irep[1]+Nets[NID].Irep[2]+Nets[NID].Irep[3]+Nets[NID].Irep[4])*1e-3f+
    NetVBse[NID].Isyn+
    NPB.gl*NPB.El)/NPB.Cm;

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

    NetVBse[NID].V = NPB.df*Vm + NPB.ef*dv;

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
  RGs[NID]=NRG;
}

void EqsGnl2_1::VmSOD(unsigned int NID)
{
  NeuRepPar2<RP3,RP4> NRP= RPs[NID];
  NeuRepVar2 NRV=RVs[NID];
  NeuGatVar3 NRG=RGs[NID];
  ParBase2 NPB=NetPBse[NID];
  float Vm=NetVBse[NID].V;

  NRG.stA = NRG.stA*NRP.ParAMPA.dfst+NRV.gAMPA;
  NRG.stG = NRG.stG*NRP.ParGABA.dfst+NRV.gGABA;
  NRG.stH = NRG.stH*NRP.ParACH.dfst+NRV.gACH;
  NRG.stL = NRG.stL*NRP.ParGCL.dfst+NRV.gGCL;

  Nets[NID].Irep[0]=NRG.stA*(NRP.ParAMPA.Vrev-Vm);
  Nets[NID].Irep[1]=NRG.stG*(NRP.ParGABA.Vrev-Vm);
  Nets[NID].Irep[2]=NRV.gNMDA*(NRP.ParNMDA.Vrev-Vm)/( 1.0f + NRP.ParNMDA.mg * expApx( -0.062f * Vm ) / 3.57f );
  Nets[NID].Irep[3]=NRG.stH*(NRP.ParACH.Vrev-Vm);
  Nets[NID].Irep[4]=NRG.stL*(NRP.ParGCL.Vrev-Vm);

  float Ax,Bx;
  Ax=NetPHH[NID].AmA*(NetPHH[NID].AmB-Vm)/(exp((NetPHH[NID].AmB-Vm)/NetPHH[NID].AmC)+NetPHH[NID].AmD);  //Am
  Bx=NetPHH[NID].BmA*exp((NetPHH[NID].BmB-Vm)/NetPHH[NID].BmC);  //Bm
  NetVHH[NID].m=HHx3(Ax,Bx,dt,NetVHH[NID].m);

  if(NetVLIF[NID].NState & FIRE_F)
    firing(NID);
  else
  {
    float dv=(
    (Nets[NID].Irep[0]+Nets[NID].Irep[1]+Nets[NID].Irep[2]+Nets[NID].Irep[3]+Nets[NID].Irep[4])*1e-3f+
    NetPHH[NID].gna*pow(NetVHH[NID].m,3.0f)*(NetPHH[NID].Ena-Vm)+
    NetVBse[NID].Isyn+
    NPB.gl*NPB.El)/NPB.Cm;
    
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

    NetVBse[NID].V = NPB.df*Vm + NPB.ef*dv;
    
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
  RGs[NID]=NRG;
}

void EqsGnl2_1::VmHH(unsigned int NID)
{
  NeuRepPar2<RP3,RP4> NRP= RPs[NID];
  NeuRepVar2 NRV=RVs[NID];
  NeuGatVar3 NRG=RGs[NID];
  ParBase2 NPB=NetPBse[NID];
  float Vm=NetVBse[NID].V;

  NRG.stA = NRG.stA*NRP.ParAMPA.dfst+NRV.gAMPA;
  NRG.stG = NRG.stG*NRP.ParGABA.dfst+NRV.gGABA;
  NRG.stH = NRG.stH*NRP.ParACH.dfst+NRV.gACH;
  NRG.stL = NRG.stL*NRP.ParGCL.dfst+NRV.gGCL;

  Nets[NID].Irep[0]=NRG.stA*(NRP.ParAMPA.Vrev-Vm);
  Nets[NID].Irep[1]=NRG.stG*(NRP.ParGABA.Vrev-Vm);
  Nets[NID].Irep[2]=NRV.gNMDA*(NRP.ParNMDA.Vrev-Vm)/( 1.0f + NRP.ParNMDA.mg * expApx( -0.062f * Vm ) / 3.57f );
  Nets[NID].Irep[3]=NRG.stH*(NRP.ParACH.Vrev-Vm);
  Nets[NID].Irep[4]=NRG.stL*(NRP.ParGCL.Vrev-Vm);

  float Ax,Bx,k1;
  Ax=NetPHH[NID].AnA*(NetPHH[NID].AnB-Vm)/(exp((NetPHH[NID].AnB-Vm)/NetPHH[NID].AnC)+NetPHH[NID].AnD);  //An
  Bx=NetPHH[NID].BnA*exp((NetPHH[NID].BnB-Vm)/NetPHH[NID].BnC);  //Bn
  NetVHH[NID].n=HHx3(Ax,Bx,dt,NetVHH[NID].n);

  Ax=NetPHH[NID].AmA*(NetPHH[NID].AmB-Vm)/(exp((NetPHH[NID].AmB-Vm)/NetPHH[NID].AmC)+NetPHH[NID].AmD);  //Am
  Bx=NetPHH[NID].BmA*exp((NetPHH[NID].BmB-Vm)/NetPHH[NID].BmC);  //Bm
  NetVHH[NID].m=HHx3(Ax,Bx,dt,NetVHH[NID].m);

  Ax=NetPHH[NID].AhA*exp((NetPHH[NID].AhB-Vm)/NetPHH[NID].AhC);  //Ah
  Bx=NetPHH[NID].BhA/(exp((NetPHH[NID].BhB-Vm)/NetPHH[NID].BhC)+NetPHH[NID].BhD);  //Bh
  NetVHH[NID].h=HHx3(Ax,Bx,dt,NetVHH[NID].h);

  k1=(2.0f*NPB.Cm*Vm/dt+NetPHH[NID].gk*pow(NetVHH[NID].n,4.0f)*NetPHH[NID].Ek+NetPHH[NID].gna*pow(NetVHH[NID].m,3.0f)*NetVHH[NID].h*NetPHH[NID].Ena+NPB.gl*NPB.El+(Nets[NID].Irep[0]+Nets[NID].Irep[1]+Nets[NID].Irep[2]+Nets[NID].Irep[3]+Nets[NID].Irep[4])*1e-3f+NetVBse[NID].Isyn)/
     (2.0f*NPB.Cm/dt+NetPHH[NID].gk*pow(NetVHH[NID].n,4.0f)+NetPHH[NID].gna*pow(NetVHH[NID].m,3.0f)*NetVHH[NID].h+NPB.gl);
  Vm=2.0f*k1-Vm;
    

  if(NetVLIF[NID].NState & FIRE_F)
    firingHH(NID);
  else
  {

    if(NetVBse[NID].V > NPB.Vth)
    {
      NetVLIF[NID].NState |= FIRE_F;
      firingHH(NID);

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
  RGs[NID]=NRG;
}

//--------------------firing process
void EqsGnl2_1::firing(unsigned int NID)
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

void EqsGnl2_1::firingHH(unsigned int NID)
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
void EqsGnl2_1::firLP(unsigned int NID)
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

void EqsGnl2_1::ltpEul(unsigned int NID)
{ NetVLTP[NID].LP = NetPLTP[NID].dfLP*NetVLTP[NID].LP + ((NetVLIF[NID].NState & LP_SPIKE_F)?1.0f:0.0f); }

//-------------background initialization
void EqsGnl2_1::BgInit(unsigned int thid, unsigned int NID)
{
  unsigned int ActCtr=Nets[NID].ActCtr;

  //AMPA
  if((ActCtr & AMPA_ACT)!=0)
  { if( BgSt[NID].pb[0]>BgSt[NID].uni_dist(*rdgen) ) RVs[NID].gAMPA += BgSt[NID].SynFst1ES[0];}

  //GABA
  if((ActCtr & GABA_ACT)!=0)
  { if( BgSt[NID].pb[1]>BgSt[NID].uni_dist(*rdgen) ) RVs[NID].gGABA += BgSt[NID].SynFst1ES[1];}

  //ACH
  if((ActCtr & ACH_ACT)!=0)
  { if( BgSt[NID].pb[2]>BgSt[NID].uni_dist(*rdgen) ) RVs[NID].gACH += BgSt[NID].SynFst1ES[2];}

  //GCL
  if((ActCtr & GCL_ACT)!=0)
  { if( BgSt[NID].pb[3]>BgSt[NID].uni_dist(*rdgen) ) RVs[NID].gGCL += BgSt[NID].SynFst1ES[3];}

  //NMDA
  if((ActCtr & NMDA_ACT)!=0)
  {
    BgSt[NID].SynSlo1ES.xt = BgSt[NID].SynSlo1ES.xt*RPs[NID].ParNMDA.dfxt + ((BgSt[NID].pb[4]>BgSt[NID].uni_dist(*rdgen))?1.0f:0.0f);
    BgSt[NID].SynSlo1ES.st = RPs[NID].ParNMDA.dfst*BgSt[NID].SynSlo1ES.st + RPs[NID].ParNMDA.efst*RPs[NID].ParNMDA.as*BgSt[NID].SynSlo1ES.xt*(1.0f - BgSt[NID].SynSlo1ES.st);
    RVs[NID].gNMDA += BgSt[NID].SynSlo1ES.g * BgSt[NID].SynSlo1ES.st;
  }


  if((ActCtr & INJ_ACT)!=0)
  { NetVBse[NID].Isyn += zig[thid].Igauss(BgSt[NID].mean,BgSt[NID].std); }

}

void EqsGnl2_1::VarInit()
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

//HH
    if(
        ((Nets[i].ActCtr & HH_ACT) !=0 ) ||
        ((Nets[i].ActCtr & SOD_ACT) !=0 )
      )
    {
      float Ax=0.01f*(10.0f-NetVBse[i].V+VRESTING)/(exp((10.0f-NetVBse[i].V+VRESTING)/10)-1.0f);  //An
      float Bx=0.125f*exp((-NetVBse[i].V+VRESTING)/80.0f);  //Bn
      NetVHH[i].n = Ax/(Ax+Bx);

      Ax=0.1f*(25.0f-NetVBse[i].V+VRESTING)/(exp((25.0f-NetVBse[i].V+VRESTING)/10)-1.0f);  //Am
      Bx=4.0f*exp((-NetVBse[i].V+VRESTING)/18.0f);  //Bm
      NetVHH[i].m = Ax/(Ax+Bx);
    
      Ax=0.07f*exp((-NetVBse[i].V+VRESTING)/20.0f);  //Ah
      Bx=1.0f/(exp((30.0f-NetVBse[i].V+VRESTING)/10)-1.0f);  //Bh
      NetVHH[i].h = Ax/(Ax+Bx);
    }

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
    RGs[i].stA=0.0f;
    RGs[i].stG=0.0f;
    RGs[i].stH=0.0f;
    RGs[i].stL=0.0f;

//Ext stimulation
    for(unsigned int ix=0;ix<5;ix++)
    { Nets[i].Irep[ix]=0.0f; }

      BgSt[i].SynSlo1ES.xt=0.0f;
      BgSt[i].SynSlo1ES.st=0.0f;
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
void EqsGnl2_1::ThdSynapse(unsigned int thid)
{
  SynFast2 Fstmp;
  SynSlow2 Slotmp;
  float gTmp,dfst;

  for(unsigned int NID=thid;NID<AMPA1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    for(unsigned int i=NeuMapAMPA1[NID].Beg;i<=NeuMapAMPA1[NID].End;i++)  //AMPA
    {
      Fstmp=SynFst1[i];
      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F) gTmp += Fstmp.g;
    }
    RVs[SynFst1[NeuMapAMPA1[NID].Beg].PosID].gAMPA+=gTmp;
  }

  for(unsigned int NID=thid;NID<GABA1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    for(unsigned int i=NeuMapGABA1[NID].Beg;i<=NeuMapGABA1[NID].End;i++)  //GABA
    {
      Fstmp=SynFst1[i];
      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F) gTmp += Fstmp.g;
    }
    RVs[SynFst1[NeuMapGABA1[NID].Beg].PosID].gGABA+=gTmp;
  }

  for(unsigned int NID=thid;NID<ACH1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    for(unsigned int i=NeuMapACH1[NID].Beg;i<=NeuMapACH1[NID].End;i++)  //ACH
    {
      Fstmp=SynFst1[i];
      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F) gTmp += Fstmp.g;
    }
    RVs[SynFst1[NeuMapACH1[NID].Beg].PosID].gACH+=gTmp;
  }

  for(unsigned int NID=thid;NID<GCL1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    for(unsigned int i=NeuMapGCL1[NID].Beg;i<=NeuMapGCL1[NID].End;i++)  //GCL
    {
      Fstmp=SynFst1[i];
      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F) gTmp += Fstmp.g;
    }
    RVs[SynFst1[NeuMapGCL1[NID].Beg].PosID].gGCL+=gTmp;
  }

  float dfxt,efst,as;
  for(unsigned int NID=thid;NID<NMDA1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    dfxt=RPs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].ParNMDA.dfxt;
    dfst=RPs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].ParNMDA.dfst;
    efst=RPs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].ParNMDA.efst;
    as  =RPs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].ParNMDA.as;

    for(unsigned int i=NeuMapNMDA1[NID].Beg;i<=NeuMapNMDA1[NID].End;i++)  //NMDA
    {
      Slotmp=SynSlo1[i];
      Slotmp.xt = Slotmp.xt*dfxt + ((NetVLIF[Slotmp.PreID].NState & SPIKE_F)? 1.0f : 0.0f);
      Slotmp.st = dfst*Slotmp.st + efst*as*Slotmp.xt*(1.0f-Slotmp.st);

      gTmp += Slotmp.st * Slotmp.g;
      SynSlo1[i].xt=Slotmp.xt;
      SynSlo1[i].st=Slotmp.st;
    }
    RVs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].gNMDA+=gTmp;
  }

}


void EqsGnl2_1::ThdSynapseSTP(unsigned int thid)
{
  SynFast2 Fstmp;
  SynSlow2 Slotmp;
  float gTmp,dfst;


  for(unsigned int NID=thid;NID<AMPA1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    for(unsigned int i=NeuMapAMPA1[NID].Beg;i<=NeuMapAMPA1[NID].End;i++)  //AMPA
    {
      Fstmp=SynFst1[i];
      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F) gTmp += Fstmp.g*NetVSTP[Fstmp.PreID].D*NetVSTP[Fstmp.PreID].F;
    }
    RVs[SynFst1[NeuMapAMPA1[NID].Beg].PosID].gAMPA+=gTmp;
  }

  for(unsigned int NID=thid;NID<GABA1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    for(unsigned int i=NeuMapGABA1[NID].Beg;i<=NeuMapGABA1[NID].End;i++)  //GABA
    {
      Fstmp=SynFst1[i];
      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F) gTmp += Fstmp.g*NetVSTP[Fstmp.PreID].D*NetVSTP[Fstmp.PreID].F;
    }
    RVs[SynFst1[NeuMapGABA1[NID].Beg].PosID].gGABA+=gTmp;
  }

  for(unsigned int NID=thid;NID<ACH1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    for(unsigned int i=NeuMapACH1[NID].Beg;i<=NeuMapACH1[NID].End;i++)  //ACH
    {
      Fstmp=SynFst1[i];
      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F) gTmp += Fstmp.g*NetVSTP[Fstmp.PreID].D*NetVSTP[Fstmp.PreID].F;
    }
    RVs[SynFst1[NeuMapACH1[NID].Beg].PosID].gACH+=gTmp;
  }

  for(unsigned int NID=thid;NID<GCL1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    for(unsigned int i=NeuMapGCL1[NID].Beg;i<=NeuMapGCL1[NID].End;i++)  //GCL
    {
      Fstmp=SynFst1[i];
      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F) gTmp += Fstmp.g*NetVSTP[Fstmp.PreID].D*NetVSTP[Fstmp.PreID].F;
    }
    RVs[SynFst1[NeuMapGCL1[NID].Beg].PosID].gGCL+=gTmp;
  }

  float dfxt,efst,as;
  for(unsigned int NID=thid;NID<NMDA1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    dfxt=RPs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].ParNMDA.dfxt;
    dfst=RPs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].ParNMDA.dfst;
    efst=RPs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].ParNMDA.efst;
    as  =RPs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].ParNMDA.as;
    
    for(unsigned int i=NeuMapNMDA1[NID].Beg;i<=NeuMapNMDA1[NID].End;i++)  //NMDA
    {
      Slotmp=SynSlo1[i];
      Slotmp.xt = Slotmp.xt*dfxt + ((NetVLIF[Slotmp.PreID].NState & SPIKE_F)? NetVSTP[Slotmp.PreID].D*NetVSTP[Slotmp.PreID].F : 0.0f);
      Slotmp.st = dfst*Slotmp.st + efst*as*Slotmp.xt*(1.0f-Slotmp.st);

      gTmp += Slotmp.st * Slotmp.g;
      SynSlo1[i].xt=Slotmp.xt;
      SynSlo1[i].st=Slotmp.st;
    }
    RVs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].gNMDA+=gTmp;
  }

}

void EqsGnl2_1::ThdSynapseSTD(unsigned int thid)
{
  SynFast2 Fstmp;
  SynSlow2 Slotmp;
  float gTmp,dfst;


  for(unsigned int NID=thid;NID<AMPA1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    for(unsigned int i=NeuMapAMPA1[NID].Beg;i<=NeuMapAMPA1[NID].End;i++)  //AMPA
    {
      Fstmp=SynFst1[i];
      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F) gTmp += Fstmp.g*NetVSTP[Fstmp.PreID].D;
    }
    RVs[SynFst1[NeuMapAMPA1[NID].Beg].PosID].gAMPA+=gTmp;
  }

  for(unsigned int NID=thid;NID<GABA1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    for(unsigned int i=NeuMapGABA1[NID].Beg;i<=NeuMapGABA1[NID].End;i++)  //GABA
    {
      Fstmp=SynFst1[i];
      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F) gTmp += Fstmp.g*NetVSTP[Fstmp.PreID].D;
    }
    RVs[SynFst1[NeuMapGABA1[NID].Beg].PosID].gGABA+=gTmp;
  }

  for(unsigned int NID=thid;NID<ACH1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    for(unsigned int i=NeuMapACH1[NID].Beg;i<=NeuMapACH1[NID].End;i++)  //ACH
    {
      Fstmp=SynFst1[i];
      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F) gTmp += Fstmp.g*NetVSTP[Fstmp.PreID].D;
    }
    RVs[SynFst1[NeuMapACH1[NID].Beg].PosID].gACH+=gTmp;
  }

  for(unsigned int NID=thid;NID<GCL1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    for(unsigned int i=NeuMapGCL1[NID].Beg;i<=NeuMapGCL1[NID].End;i++)  //GCL
    {
      Fstmp=SynFst1[i];
      if(NetVLIF[Fstmp.PreID].NState & SPIKE_F) gTmp += Fstmp.g*NetVSTP[Fstmp.PreID].D;
    }
    RVs[SynFst1[NeuMapGCL1[NID].Beg].PosID].gGCL+=gTmp;
  }

  float dfxt,efst,as;
  for(unsigned int NID=thid;NID<NMDA1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    dfxt=RPs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].ParNMDA.dfxt;
    dfst=RPs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].ParNMDA.dfst;
    efst=RPs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].ParNMDA.efst;
    as  =RPs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].ParNMDA.as;
    
    for(unsigned int i=NeuMapNMDA1[NID].Beg;i<=NeuMapNMDA1[NID].End;i++)  //NMDA
    {
      Slotmp=SynSlo1[i];
      Slotmp.xt = Slotmp.xt*dfxt + ((NetVLIF[Slotmp.PreID].NState & SPIKE_F)? NetVSTP[Slotmp.PreID].D : 0.0f);
      Slotmp.st = dfst*Slotmp.st + efst*as*Slotmp.xt*(1.0f-Slotmp.st);

      gTmp += Slotmp.st * Slotmp.g;
      SynSlo1[i].xt=Slotmp.xt;
      SynSlo1[i].st=Slotmp.st;
    }
    RVs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].gNMDA+=gTmp;
  }

}

void EqsGnl2_1::ThdGapSyn(unsigned int thid)
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


//-------------Long term plasticity
void EqsGnl2_1::LTP(unsigned int thid)
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
void EqsGnl2_1::STP(unsigned int thid)
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

/*
    Vstp.Car = Pstp.dfCar*Vstp.Car+Pstp.efCar +(NetVLIF[NID].NState & SPIKE_F)?Pstp.aCar:0.0f;
    Vstp.Cap = Pstp.dfCap*Vstp.Cap + (NetVLIF[NID].NState & SPIKE_F)?Pstp.aCap:0.0f;
*/
	  NetVSTP[NID]=Vstp;
  }

}

void EqsGnl2_1::STD(unsigned int thid)
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

void EqsGnl2_1::EqsOne(unsigned int ActCtr)
{
  Activation(0,ActCtr);
  Spike(0,ActCtr);
}


void EqsGnl2_1::Activation(unsigned int thid,unsigned int ActCtr)
{
  
  if((ActCtr & HH_ACT)!=0)
  {
    for(unsigned NID=thid;NID<NetSize;NID+=thsize)
    {
      BgInit(thid,NID); //bcak ground initial
      VmHH(NID);
    }
  }
  else if((ActCtr & SOD_ACT)!=0)
  {
    for(unsigned NID=thid;NID<NetSize;NID+=thsize)
    {
      BgInit(thid,NID); //bcak ground initial
      VmSOD(NID);
    }
  }
  else
  {
    for(unsigned NID=thid;NID<NetSize;NID+=thsize)
    {
      BgInit(thid,NID); //bcak ground initial
      VmEul(NID);
    }
  }
}

void EqsGnl2_1::Spike(unsigned int thid,unsigned int ActCtr)
{
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



