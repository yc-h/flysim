
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
void EqsGnl1_1::ltpEul(unsigned int NID)
{ NetVLTP[NID].LP -= dt*(NetVLTP[NID].LP/NetPLTP[NID].tLP) + ((NetVLIF[NID].NState & LP_SPIKE_F)?1.0f:0.0f); }



//------------------------------post synaptic processing:st*g
void EqsGnl1_1::ThdSynapse(unsigned int thid)
{
  SynFast2 Fstmp;
  SynSlow2 Slotmp;
  float gTmp,tst;
  
  for(unsigned int NID=thid;NID<AMPA1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    tst=RPs[SynFst1[NeuMapAMPA1[NID].Beg].PosID].ParAMPA.tst;
    for(unsigned int i=NeuMapAMPA1[NID].Beg;i<=NeuMapAMPA1[NID].End;i++)  //AMPA
    {
      Fstmp=SynFst1[i];
      Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F)? 1.0f:-dt*Fstmp.st/tst);
      gTmp += Fstmp.st*Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapAMPA1[NID].Beg].PosID].gAMPA+=gTmp;
  }

  for(unsigned int NID=thid;NID<GABA1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    tst=RPs[SynFst1[NeuMapGABA1[NID].Beg].PosID].ParGABA.tst;
    for(unsigned int i=NeuMapGABA1[NID].Beg;i<=NeuMapGABA1[NID].End;i++)  //GABA
    {
      Fstmp=SynFst1[i];
      Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F)? 1.0f:-dt*Fstmp.st/tst);
      gTmp += Fstmp.st*Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapGABA1[NID].Beg].PosID].gGABA+=gTmp;
  }

  for(unsigned int NID=thid;NID<ACH1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    tst=RPs[SynFst1[NeuMapACH1[NID].Beg].PosID].ParACH.tst;
    for(unsigned int i=NeuMapACH1[NID].Beg;i<=NeuMapACH1[NID].End;i++)  //ACH
    {
      Fstmp=SynFst1[i];
      Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F)? 1.0f:-dt*Fstmp.st/tst);
      gTmp += Fstmp.st*Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapACH1[NID].Beg].PosID].gACH+=gTmp;
  }

  for(unsigned int NID=thid;NID<GCL1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    tst=RPs[SynFst1[NeuMapGCL1[NID].Beg].PosID].ParGCL.tst;
    for(unsigned int i=NeuMapGCL1[NID].Beg;i<=NeuMapGCL1[NID].End;i++)  //GCL
    {
      Fstmp=SynFst1[i];
      Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F)? 1.0f:-dt*Fstmp.st/tst);
      gTmp += Fstmp.st*Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapGCL1[NID].Beg].PosID].gGCL+=gTmp;
  }

  float txt,as;
  for(unsigned int NID=thid;NID<NMDA1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    txt =RPs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].ParNMDA.txt;
    tst =RPs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].ParNMDA.tst;
    as  =RPs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].ParNMDA.as;

    for(unsigned int i=NeuMapNMDA1[NID].Beg;i<=NeuMapNMDA1[NID].End;i++)  //NMDA
    {
      Slotmp=SynSlo1[i];
      Slotmp.xt += ((NetVLIF[Slotmp.PreID].NState & SPIKE_F)? 1.0f:-dt*(Slotmp.xt/txt));
      Slotmp.st += dt*(as*Slotmp.xt*(1.0f - Slotmp.st) - Slotmp.st/tst);
      gTmp += Slotmp.st * Slotmp.g;
      SynSlo1[i].xt=Slotmp.xt;
      SynSlo1[i].st=Slotmp.st;
    }

    RVs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].gNMDA+=gTmp;
  }

}


void EqsGnl1_1::ThdSynapseSTP(unsigned int thid)
{
  SynFast2 Fstmp;
  SynSlow2 Slotmp;
  float gTmp,tst;


  for(unsigned int NID=thid;NID<AMPA1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    tst=RPs[SynFst1[NeuMapAMPA1[NID].Beg].PosID].ParAMPA.tst;
    for(unsigned int i=NeuMapAMPA1[NID].Beg;i<=NeuMapAMPA1[NID].End;i++)  //AMPA
    {
      Fstmp=SynFst1[i];
      Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F)? NetVSTP[Fstmp.PreID].F*NetVSTP[Fstmp.PreID].D:-dt*Fstmp.st/tst);
      gTmp += Fstmp.st*Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapAMPA1[NID].Beg].PosID].gAMPA+=gTmp;
  }

  for(unsigned int NID=thid;NID<GABA1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    tst=RPs[SynFst1[NeuMapGABA1[NID].Beg].PosID].ParGABA.tst;
    for(unsigned int i=NeuMapGABA1[NID].Beg;i<=NeuMapGABA1[NID].End;i++)  //GABA
    {
      Fstmp=SynFst1[i];
      Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F)? NetVSTP[Fstmp.PreID].F*NetVSTP[Fstmp.PreID].D:-dt*Fstmp.st/tst);
      gTmp += Fstmp.st*Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapGABA1[NID].Beg].PosID].gGABA+=gTmp;
  }

  for(unsigned int NID=thid;NID<ACH1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    tst=RPs[SynFst1[NeuMapACH1[NID].Beg].PosID].ParACH.tst;
    for(unsigned int i=NeuMapACH1[NID].Beg;i<=NeuMapACH1[NID].End;i++)  //ACH
    {
      Fstmp=SynFst1[i];
      Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F)? NetVSTP[Fstmp.PreID].F*NetVSTP[Fstmp.PreID].D:-dt*Fstmp.st/tst);
      gTmp += Fstmp.st*Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapACH1[NID].Beg].PosID].gACH+=gTmp;
  }

  for(unsigned int NID=thid;NID<GCL1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    tst=RPs[SynFst1[NeuMapGCL1[NID].Beg].PosID].ParGCL.tst;
    for(unsigned int i=NeuMapGCL1[NID].Beg;i<=NeuMapGCL1[NID].End;i++)  //GCL
    {
      Fstmp=SynFst1[i];
      Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F)? NetVSTP[Fstmp.PreID].F*NetVSTP[Fstmp.PreID].D:-dt*Fstmp.st/tst);
      gTmp += Fstmp.st*Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapGCL1[NID].Beg].PosID].gGCL+=gTmp;
  }

  float txt,as;
  for(unsigned int NID=thid;NID<NMDA1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    txt =RPs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].ParNMDA.txt;
    tst =RPs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].ParNMDA.tst;
    as  =RPs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].ParNMDA.as;

    for(unsigned int i=NeuMapNMDA1[NID].Beg;i<=NeuMapNMDA1[NID].End;i++)  //NMDA
    {
      Slotmp=SynSlo1[i];
      Slotmp.xt +=((NetVLIF[Slotmp.PreID].NState & SPIKE_F)? NetVSTP[Slotmp.PreID].F*NetVSTP[Slotmp.PreID].D: -dt*(Slotmp.xt/txt));
      Slotmp.st += dt*(as*Slotmp.xt*(1.0f - Slotmp.st) - Slotmp.st/tst);
      gTmp += Slotmp.st * Slotmp.g;
      SynSlo1[i].xt=Slotmp.xt;
      SynSlo1[i].st=Slotmp.st;
    }

    RVs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].gNMDA+=gTmp;
  }

}

void EqsGnl1_1::ThdSynapseSTD(unsigned int thid)
{
  SynFast2 Fstmp;
  SynSlow2 Slotmp;
  float gTmp,tst;


  for(unsigned int NID=thid;NID<AMPA1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    tst=RPs[SynFst1[NeuMapAMPA1[NID].Beg].PosID].ParAMPA.tst;
    for(unsigned int i=NeuMapAMPA1[NID].Beg;i<=NeuMapAMPA1[NID].End;i++)  //AMPA
    {
      Fstmp=SynFst1[i];
      Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F)? NetVSTP[Fstmp.PreID].D:-dt*Fstmp.st/tst);
      gTmp += Fstmp.st*Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapAMPA1[NID].Beg].PosID].gAMPA+=gTmp;
  }

  for(unsigned int NID=thid;NID<GABA1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    tst=RPs[SynFst1[NeuMapGABA1[NID].Beg].PosID].ParGABA.tst;
    for(unsigned int i=NeuMapGABA1[NID].Beg;i<=NeuMapGABA1[NID].End;i++)  //GABA
    {
      Fstmp=SynFst1[i];
      Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F)? NetVSTP[Fstmp.PreID].D:-dt*Fstmp.st/tst);
      gTmp += Fstmp.st*Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapGABA1[NID].Beg].PosID].gGABA+=gTmp;
  }

  for(unsigned int NID=thid;NID<ACH1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    tst=RPs[SynFst1[NeuMapACH1[NID].Beg].PosID].ParACH.tst;
    for(unsigned int i=NeuMapACH1[NID].Beg;i<=NeuMapACH1[NID].End;i++)  //ACH
    {
      Fstmp=SynFst1[i];
      Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F)? NetVSTP[Fstmp.PreID].D:-dt*Fstmp.st/tst);
      gTmp += Fstmp.st*Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapACH1[NID].Beg].PosID].gACH+=gTmp;
  }

  for(unsigned int NID=thid;NID<GCL1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    tst=RPs[SynFst1[NeuMapGCL1[NID].Beg].PosID].ParGCL.tst;
    for(unsigned int i=NeuMapGCL1[NID].Beg;i<=NeuMapGCL1[NID].End;i++)  //GCL
    {
      Fstmp=SynFst1[i];
      Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F)? NetVSTP[Fstmp.PreID].D:-dt*Fstmp.st/tst);
      gTmp += Fstmp.st*Fstmp.g;
      SynFst1[i].st=Fstmp.st;
    }
    RVs[SynFst1[NeuMapGCL1[NID].Beg].PosID].gGCL+=gTmp;
  }

  float txt,as;
  for(unsigned int NID=thid;NID<NMDA1CntSize;NID+=thsize)
  {
    gTmp=0.0f;
    txt =RPs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].ParNMDA.txt;
    tst =RPs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].ParNMDA.tst;
    as  =RPs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].ParNMDA.as;

    for(unsigned int i=NeuMapNMDA1[NID].Beg;i<=NeuMapNMDA1[NID].End;i++)  //NMDA
    {
      Slotmp=SynSlo1[i];
      Slotmp.xt +=((NetVLIF[Slotmp.PreID].NState & SPIKE_F)? NetVSTP[Slotmp.PreID].D: -dt*(Slotmp.xt/txt));
      Slotmp.st += dt*(as*Slotmp.xt*(1.0f - Slotmp.st) - Slotmp.st/tst);
      gTmp += Slotmp.st * Slotmp.g;
      SynSlo1[i].xt=Slotmp.xt;
      SynSlo1[i].st=Slotmp.st;
    }

    RVs[SynSlo1[NeuMapNMDA1[NID].Beg].PosID].gNMDA+=gTmp;
  }

}

//--------------------------------------------------
void EqsGnl1_1::ThdPGatVarSTP(unsigned int thid)
{
  SynFast2 Fstmp;
  SynSlow2 Slotmp;

  for(unsigned int i=thid;i<PGatASize;i+=thsize)  //AMPA
  {
    Fstmp=PSynFstGatA[i];
    Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F) ? NetVSTP[Fstmp.PreID].D*NetVSTP[Fstmp.PreID].F : -dt*Fstmp.st/RPs[Fstmp.PosID].ParAMPA.tst);
    PSynFstGatA[i].st=Fstmp.st;
  }

  for(unsigned int i=thid;i<PGatGSize;i+=thsize)  //GABA
  {
    Fstmp=PSynFstGatG[i];
    Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F) ? NetVSTP[Fstmp.PreID].D*NetVSTP[Fstmp.PreID].F : -dt*Fstmp.st/RPs[Fstmp.PosID].ParGABA.tst);
    PSynFstGatG[i].st=Fstmp.st;
  }

  for(unsigned int i=thid;i<PGatHSize;i+=thsize)  //ACH
  {
    Fstmp=PSynFstGatH[i];
    Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F) ? NetVSTP[Fstmp.PreID].D*NetVSTP[Fstmp.PreID].F : -dt*Fstmp.st/RPs[Fstmp.PosID].ParACH.tst);
    PSynFstGatH[i].st=Fstmp.st;
  }

  for(unsigned int i=thid;i<PGatLSize;i+=thsize)  //GCL
  {
    Fstmp=PSynFstGatL[i];
    Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F) ? NetVSTP[Fstmp.PreID].D*NetVSTP[Fstmp.PreID].F : -dt*Fstmp.st/RPs[Fstmp.PosID].ParGCL.tst);
    PSynFstGatL[i].st=Fstmp.st;
  }

  for(unsigned int i=thid;i<PGatNSize;i+=thsize)  //NMDA
  {
    Slotmp=PSynSloGat[i];
    Slotmp.xt += ((NetVLIF[Slotmp.PreID].NState & SPIKE_F)? NetVSTP[Slotmp.PreID].D*NetVSTP[Slotmp.PreID].F: -dt*(Slotmp.xt/RPs[Slotmp.PosID].ParNMDA.txt));
    Slotmp.st += dt*(RPs[Slotmp.PosID].ParNMDA.as*Slotmp.xt*(1.0f - Slotmp.st) - Slotmp.st/RPs[Slotmp.PosID].ParNMDA.tst);
    PSynSloGat[i].xt=Slotmp.xt;
    PSynSloGat[i].st=Slotmp.st;
  }
}


void EqsGnl1_1::ThdPGatVarSTD(unsigned int thid)
{
  SynFast2 Fstmp;
  SynSlow2 Slotmp;

  for(unsigned int i=thid;i<PGatASize;i+=thsize)  //AMPA
  {
    Fstmp=PSynFstGatA[i];
    Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F) ? NetVSTP[Fstmp.PreID].D : -dt*Fstmp.st/RPs[Fstmp.PosID].ParAMPA.tst);
    PSynFstGatA[i].st=Fstmp.st;
  }

  for(unsigned int i=thid;i<PGatGSize;i+=thsize)  //GABA
  {
    Fstmp=PSynFstGatG[i];
    Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F) ? NetVSTP[Fstmp.PreID].D : -dt*Fstmp.st/RPs[Fstmp.PosID].ParGABA.tst);
    PSynFstGatG[i].st=Fstmp.st;
  }

  for(unsigned int i=thid;i<PGatHSize;i+=thsize)  //ACH
  {
    Fstmp=PSynFstGatH[i];
    Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F) ? NetVSTP[Fstmp.PreID].D : -dt*Fstmp.st/RPs[Fstmp.PosID].ParACH.tst);
    PSynFstGatH[i].st=Fstmp.st;
  }

  for(unsigned int i=thid;i<PGatLSize;i+=thsize)  //GCL
  {
    Fstmp=PSynFstGatL[i];
    Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F) ? NetVSTP[Fstmp.PreID].D : -dt*Fstmp.st/RPs[Fstmp.PosID].ParGCL.tst);
    PSynFstGatL[i].st=Fstmp.st;
  }

  for(unsigned int i=thid;i<PGatNSize;i+=thsize)  //NMDA
  {
    Slotmp=PSynSloGat[i];
    Slotmp.xt += ((NetVLIF[Slotmp.PreID].NState & SPIKE_F)? NetVSTP[Slotmp.PreID].D: -dt*(Slotmp.xt/RPs[Slotmp.PosID].ParNMDA.txt));
    Slotmp.st += dt*(RPs[Slotmp.PosID].ParNMDA.as*Slotmp.xt*(1.0f - Slotmp.st) - Slotmp.st/RPs[Slotmp.PosID].ParNMDA.tst);
    PSynSloGat[i].xt=Slotmp.xt;
    PSynSloGat[i].st=Slotmp.st;
  }
}


void EqsGnl1_1::ThdPGatVar(unsigned int thid)
{
  SynFast2 Fstmp;
  SynSlow2 Slotmp;

  for(unsigned int i=thid;i<PGatASize;i+=thsize)  //AMPA
  {
    Fstmp=PSynFstGatA[i];
    Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F) ? 1.0f : -dt*Fstmp.st/RPs[Fstmp.PosID].ParAMPA.tst);
    PSynFstGatA[i].st=Fstmp.st;
  }

  for(unsigned int i=thid;i<PGatGSize;i+=thsize)  //GABA
  {
    Fstmp=PSynFstGatG[i];
    Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F) ? 1.0f : -dt*Fstmp.st/RPs[Fstmp.PosID].ParGABA.tst);
    PSynFstGatG[i].st=Fstmp.st;
  }

  for(unsigned int i=thid;i<PGatHSize;i+=thsize)  //ACH
  {
    Fstmp=PSynFstGatH[i];
    Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F) ? 1.0f : -dt*Fstmp.st/RPs[Fstmp.PosID].ParACH.tst);
    PSynFstGatH[i].st=Fstmp.st;
  }

  for(unsigned int i=thid;i<PGatLSize;i+=thsize)  //GCL
  {
    Fstmp=PSynFstGatL[i];
    Fstmp.st += ((NetVLIF[Fstmp.PreID].NState & SPIKE_F) ? 1.0f : -dt*Fstmp.st/RPs[Fstmp.PosID].ParGCL.tst);
    PSynFstGatL[i].st=Fstmp.st;
  }

  for(unsigned int i=thid;i<PGatNSize;i+=thsize)  //NMDA
  {
    Slotmp=PSynSloGat[i];
    Slotmp.xt += ((NetVLIF[Slotmp.PreID].NState & SPIKE_F)? 1.0f:-dt*(Slotmp.xt/RPs[Slotmp.PosID].ParNMDA.txt));
    Slotmp.st += dt*(RPs[Slotmp.PosID].ParNMDA.as*Slotmp.xt*(1.0f - Slotmp.st) - Slotmp.st/RPs[Slotmp.PosID].ParNMDA.tst);
    PSynSloGat[i].xt=Slotmp.xt;
    PSynSloGat[i].st=Slotmp.st;
  }
}

//-------------short term plasticity
void EqsGnl1_1::STP(unsigned int thid)
{
  ParStp Pstp;
  VarStp Vstp;
  for(unsigned NID=thid;NID<NetSize;NID+=thsize)
  {
      Pstp=NetPSTP[NID];
      Vstp=NetVSTP[NID];

    if(NetVLIF[NID].NState & SPIKE_F)
    {
      Vstp.Car += Pstp.aCar;
      Vstp.Cap += Pstp.aCap;
      Vstp.F += dt*((Vstp.Cap+Vstp.Car)*(1.0f - Vstp.F)*Pstp.aF-Vstp.F/Pstp.tF);
      Vstp.D *= ( 1.0f -  Pstp.pv * Vstp.F);
    }
    else
    {
      Vstp.Car += dt*(CarBase*1e-6f - Vstp.Car)/Pstp.tCar;
      Vstp.Cap -= dt*Vstp.Cap/Pstp.tCap;
      Vstp.F += dt*((Vstp.Cap+Vstp.Car)*(1.0f - Vstp.F)*Pstp.aF-Vstp.F/Pstp.tF);
      Vstp.D += dt*(1.0f-Vstp.D)/Pstp.tD;
    }

    NetVSTP[NID]=Vstp;
  }

}


void EqsGnl1_1::STD(unsigned int thid)
{
  ParStp Pstp;
  for(unsigned NID=thid;NID<NetSize;NID+=thsize)
  {
      Pstp=NetPSTP[NID];

    if(NetVLIF[NID].NState & SPIKE_F)
      NetVSTP[NID].D *= ( 1.0f -  Pstp.pv);
    else
      NetVSTP[NID].D += dt*(1.0f-NetVSTP[NID].D)/Pstp.tD;
  }

}
