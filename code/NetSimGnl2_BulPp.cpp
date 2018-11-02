#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <stdio.h>
#include <fstream>
#include <cstdlib>
#include <random>
#include <chrono>
#include <unistd.h>
#include <algorithm>
#include <string.h>
#include <condition_variable>
#include <atomic>
#include <thread>
#include <mutex>
#include <cmath>
#include <limits.h>
#include "NetSim.h"



using namespace std;
//--------------setup network population

void NetSimGnl2::SetPp(int thNum,int sop)
{
  dt=dtime;
  thsize=thNum;
  th=new thread[thsize];


//-----struct data initial
//-----------------Confs transcripts to Populations
  vector<IdxSynapse2> IdxSyn;

//-----------------create default neuron network

  NIdx.insert(NIdx.end(),ConfSta.size(),0);
  unsigned int q=0;
  for(unsigned int i=0;i<ConfSta.size();i++)
  {
    NIdx[i]=q;
    q += ConfSta[i].NeuNum;
  }
  TolNeu=q;
  NetSize=q;
  RPsSize=q;
  RVsSize=q;

  NetPBseSize=q;
  NetVBseSize=q;
  NetPLIFSize=q;
  NetVLIFSize=q;
  NetPSTPSize=q;
  NetVSTPSize=q;
  NetPLTPSize=q;
  NetVLTPSize=q;

  Nets = new NeuronGnl2[NetSize];
  RPs = new NeuRepPar2<RP3,RP4>[RPsSize];
  RVs = new NeuRepVar2[RVsSize];

  NetPBse = new ParBase2[NetSize];
  NetVBse = new VarBase[NetSize];

  NetPLIF = new ParLIF[NetSize];
  NetVLIF = new VarLIF2[NetSize];

  NetPSTP = new ParStp2[NetSize];
  NetVSTP = new VarStp[NetSize];

  NetPLTP = new ParLP2[NetSize];
  NetVLTP = new VarLP2[NetSize];

//back ground noise setup
  zig=new ziggurat[thsize];
  BgSt = new BgSti1<RV3,RV4>[NetSize];
  rseeds=new unsigned int[RpSize];
  rdgen = new mt19937;

  if((ActCtr & UDFRDSEED_ACT)!= 0)
    *rseeds=rseed; // rseed from user define command option
  else if(MacRndSeds.size()!=0)
    *rseeds=MacRndSeds[0];
  else
    *rseeds=seed_gen();

  rdgen->seed(*rseeds);
  zig[0].jcong=*rseeds;

  uniform_int_distribution<unsigned int> dist(0,RAND_MAX);

  if(thsize>1)
  {
    for(unsigned int i=1;i<thsize;i++)
    { zig[i].jcong=dist(*rdgen); }
  }

  if(RpSize>1)
  {
    for(unsigned int i=1,s;i<(unsigned int)RpSize;i++)
    {
      if(i<MacRndSeds.size())
        rseeds[i]=MacRndSeds[i];
      else
      {
        bool found=true;
        while(found)
        {
          found=false;
          s=dist(*rdgen);
          for(unsigned int j=0;j<i;j++)
          {
            if( s == rseeds[j] )
            {
              found=true;
              break;
            }
          }
        }

        rseeds[i]=s;
      }
    }
  }


  for(unsigned int i=0;i<NetSize;i++)
  {
    BgSt[i].SynSlo1ES.xt=0.0f;
    BgSt[i].SynSlo1ES.st=0.0f;
    BgSt[i].SynSlo1ES.g=0.0f;

    for(unsigned int j=0;j<4;j++)
    {
      BgSt[i].SynFst1ES[j].st=0.0f;
      BgSt[i].SynFst1ES[j].g=0.0f;
    }

      BgSt[i].mean=MNoiseMean;
      BgSt[i].std=MNoiseSTD;

      if(BgSt[i].mean==0.0f && BgSt[i].std==0.0f)
        Nets[i].ActCtr &= ~INJ_ACT;
      else
        Nets[i].ActCtr |= INJ_ACT;

      uniform_real_distribution<float>::param_type uni_par(0.0f,1.0f);
      BgSt[i].uni_dist.param(uni_par);

    for(unsigned int j=0;j<5;j++)
    { BgSt[i].pb[j]=0.0f; }
  }


//-------------------------
  for(unsigned int i=0;i<ConfSta.size();i++)
  {
    for(unsigned int j=NIdx[i];j<(NIdx[i]+ConfSta[i].NeuNum);j++) // set up each neuron's parameter in every population
    {
      Nets[j].ActCtr=ActCtr;

      //leaky integrate and fire
      NetPBse[j].Vth    = Confs[i].Vth;
      NetPBse[j].El     = Confs[i].El;
      NetPBse[j].gl     = Confs[i].Cm/Confs[i].Taum;
      NetPBse[j].Cm     = Confs[i].Cm;
      NetPBse[j].Vreset = Confs[i].Vreset;
      NetPBse[j].df     = exp(-dt/Confs[i].Taum);
      NetPBse[j].ef     = Confs[i].Taum*(1.0f-NetPBse[j].df);

      NetPLIF[j].RfPed  = Confs[i].RfPed;
      NetPLIF[j].SpiDy  = Confs[i].SpiDy;

      BgSt[i].std=MNoiseSTD*(1.0f+(float)log(1.0f-exp(-0.1f/Confs[i].Taum))/(1.0f-NetPBse[j].df)); //normalize to dt=0.1ms

      //ltp
      NetPLTP[j].dfLP    =  exp(-dt/ConfLtp[i].tLP);
      NetPLTP[j].SpiLPDy = ConfLtp[i].SpiLPDy; //back propogation spike delay
      NetPLTP[j].PosISI  = ConfLtp[i].PosISI;
      NetPLTP[j].NegISI  = ConfLtp[i].NegISI;

      //stp
      NetPSTP[j].pv   = ConfStp[i].pv;
      NetPSTP[j].dfD  = exp(-dt/ConfStp[i].tD);
      NetPSTP[j].efD  = 1.0f-NetPSTP[j].dfD;

      NetPSTP[j].aF   = ConfStp[i].aF;
      NetPSTP[j].dfF  = exp(-dt/ConfStp[i].tF);
      NetPSTP[j].efF  = ConfStp[i].tF*(1.0f-NetPSTP[j].dfF);


      NetPSTP[j].aCap = ConfStp[i].aCap;
      NetPSTP[j].dfCap= exp(-dt/ConfStp[i].tCap);

      NetPSTP[j].aCar = ConfStp[i].aCar;
      NetPSTP[j].dfCar= exp(-dt/ConfStp[i].tCar);
      NetPSTP[j].efCar= (1.0f-NetPSTP[j].dfCar)*CarBase*1e-6f;

      Nets[j].Irep = new float[5];

      for(unsigned int ix=0;ix<5;ix++)
          Nets[j].Irep[ix]=0.0f;

//initial RPs and Nets[j].pb Nets[j].ExtFst.g Nets[j].ExtSlo.g
      RPs[j].ParAMPA.Vrev = 0.0f;
      RPs[j].ParAMPA.dfst = 0.1f;

      RPs[j].ParGABA.Vrev = -90.0f;
      RPs[j].ParGABA.dfst = 0.1f;
      NetPBse[j].dVmax = (NetPBse[j].Vth - RPs[j].ParGABA.Vrev)/NetPBse[j].ef; // assume GABA reversal potential is -90mV

      RPs[j].ParNMDA.as=0.6332f; // the alpha st
      RPs[j].ParNMDA.dfxt = 0.1f;
      RPs[j].ParNMDA.dfst = 0.1f;
      RPs[j].ParNMDA.efst = 0.01f;
      RPs[j].ParNMDA.Vrev = 0.0f;
      RPs[j].ParNMDA.mg=1.0f; // the concentration of magnesium

      RPs[j].ParACH.Vrev = 0.0f;
      RPs[j].ParACH.dfst = 0.1f;

      RPs[j].ParGCL.Vrev = -90.0f;
      RPs[j].ParGCL.dfst = 0.1f;

      for(unsigned int k=0;k<ConfSta[i].PpRepNum;k++)  // parameter in receptors and variable in FExt receptors
      {

        if(PreSynR->NeuRepName==KW_AMPA)
        {
          RPs[j].ParAMPA.Vrev = PreSynR->Vrev;
          RPs[j].ParAMPA.dfst = exp(-dt/PreSynR->Tau);
          BgSt[j].pb[0]=PreSynR->FExt*dt*1e-3f;
          BgSt[j].SynFst1ES[0].g=PreSynR->MnExtEff*PreSynR->MnExtCon;
        }
        else if(PreSynR->NeuRepName==KW_GABA)
        {
          NetPBse[j].dVmax = ((Confs+i)->Vth - PreSynR->Vrev)/NetPBse[j].ef;

          RPs[j].ParGABA.Vrev = PreSynR->Vrev;
          RPs[j].ParGABA.dfst = exp(-dt/PreSynR->Tau);
          BgSt[j].pb[1]=PreSynR->FExt*dt*1e-3f;
          BgSt[j].SynFst1ES[1].g=PreSynR->MnExtEff*PreSynR->MnExtCon;
        }
        else if(PreSynR->NeuRepName==KW_ACH)
        {
          RPs[j].ParACH.Vrev = PreSynR->Vrev;
          RPs[j].ParACH.dfst = exp(-dt/PreSynR->Tau);
          BgSt[j].pb[2]=PreSynR->FExt*dt*1e-3f;
          BgSt[j].SynFst1ES[2].g=PreSynR->MnExtEff*PreSynR->MnExtCon;
        }
        else if(PreSynR->NeuRepName==KW_GCL)
        {
          RPs[j].ParGCL.Vrev = PreSynR->Vrev;
          RPs[j].ParGCL.dfst = exp(-dt/PreSynR->Tau);
          BgSt[j].pb[3]=PreSynR->FExt*dt*1e-3f;
          BgSt[j].SynFst1ES[3].g=PreSynR->MnExtEff*PreSynR->MnExtCon;
        }
        else if(PreSynR->NeuRepName==KW_NMDA)
        {
          RPs[j].ParNMDA.as=0.6332f; // the alpha st
          RPs[j].ParNMDA.dfxt = exp(-dt/PreSynR->TauXT);
          RPs[j].ParNMDA.dfst = exp(-dt/PreSynR->Tau);
          RPs[j].ParNMDA.efst = PreSynR->Tau*(1.0f-RPs[j].ParNMDA.dfst);
          RPs[j].ParNMDA.Vrev = PreSynR->Vrev;
          RPs[j].ParNMDA.mg=1.0f; // the concentration of magnesium

          BgSt[j].pb[4]=PreSynR->FExt*dt*1e-3f;
          BgSt[j].SynSlo1ES.g=PreSynR->MnExtEff*PreSynR->MnExtCon;
        }

      }

    }

  }


  delete [] ConfStp;
  delete [] ConfLtp;
  delete [] ConfHH;


//-----------------------------------------------------------------------------
//population specific

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

  for(unsigned int i=0;i<ConfSta.size();i++) // pre synatic population
  {
    for(unsigned int k=NIdx[i];k<(NIdx[i]+ConfSta[i].NeuNum);k++) // build Fast synapse of j neuron post links
    {
      for(unsigned int j=ConfSta[i].PpLinkBeg,PosIDNum=0;j<ConfSta[i].PpLinkBeg+ConfSta[i].PpLinkNum;j++) // pre synaptic population synapse
      {
        unsigned short Act=GetAct(POSYNREP);
        PosIDNum=(unsigned int)ceil((double)ConfSta[POSYNPP].NeuNum * POSYNCON);

        if((i==POSYNPP) && (ConfSta[POSYNPP].NeuNum==PosIDNum) && !ConfSta[i].SelConEn) // no self connction
        { PosIDNum--; }

        if(PosIDNum<2)  {}
        else if(Act==AMPA_ACT)
        {
          PSynASize+=PosIDNum;
          PGatASize++;
        }
        else if(Act==GABA_ACT)
        {
          PSynGSize+=PosIDNum;
          PGatGSize++;
        }
        else if(Act==ACH_ACT)
        {
          PSynHSize+=PosIDNum;
          PGatHSize++;
        }
        else if(Act==GCL_ACT)
        {
          PSynLSize+=PosIDNum;
          PGatLSize++;
        }
        else if(Act==NMDA_ACT)
        {
          PSynNSize+=PosIDNum;
          PGatNSize++;
        }
        else if(Act==GAP_ACT)
        {
          PSynCpSize+=PosIDNum;
          PGapPSize++;
        }
      }
    }
  }
  SynIdx *TPSynA = new SynIdx[PSynASize];
  SynIdx *TPSynG = new SynIdx[PSynGSize];
  SynIdx *TPSynH = new SynIdx[PSynHSize];
  SynIdx *TPSynL = new SynIdx[PSynLSize];
  SynIdx *TPSynN = new SynIdx[PSynNSize];
  SynIdx *TPSynCp = new SynIdx[PSynCpSize];

  PSynFstGatA=new SynFast2[PGatASize];
  PSynFstGatG=new SynFast2[PGatGSize];
  PSynFstGatH=new SynFast2[PGatHSize];
  PSynFstGatL=new SynFast2[PGatLSize];
  PSynSloGat=new SynSlow2[PGatNSize];
  PSynCndP=new PSynGap[PGapPSize];
  NeuMapPCndO=new NeuMap[PGapPSize];

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

  for(unsigned int i=0;i<ConfSta.size();i++) // pre synatic population
  {
    for(unsigned int k=NIdx[i];k<(NIdx[i]+ConfSta[i].NeuNum);k++) // build Fast synapse of j neuron post links
    {
      for(unsigned int j=ConfSta[i].PpLinkBeg,PosIDNum=0;j<ConfSta[i].PpLinkBeg+ConfSta[i].PpLinkNum;j++) // pre synaptic population synapse
      {
        unsigned short Act=GetAct(POSYNREP);
        PosIDNum=(unsigned int)ceil((double)ConfSta[POSYNPP].NeuNum * POSYNCON);

        if((i==POSYNPP) && (ConfSta[POSYNPP].NeuNum==PosIDNum) && !ConfSta[i].SelConEn) // no self connction
        { PosIDNum--; }

        unsigned int *PosID=PosBuild(k,PosIDNum,NIdx[POSYNPP],ConfSta[POSYNPP].NeuNum,ConfSta[i].SelConEn);
        float g=POSYNEFF * POSYNWGT;

        if(PosIDNum<2)  {}
        else if(Act==AMPA_ACT)
        {
          for(unsigned int kk=0;kk<PosIDNum;kk++)  //copy PosID
          {
            TPSynA[PSynASize+kk].PreID=PGatASize;
            TPSynA[PSynASize+kk].PosID=PosID[kk];
          }
          PSynFstGatA[PGatASize].PreID=k;
          PSynFstGatA[PGatASize].PosID=PosID[0];
          PSynFstGatA[PGatASize].g=g;
          PSynFstGatA[PGatASize].st=0.0f;
          PSynASize+=PosIDNum;
          PGatASize++;
        }
        else if(Act==GABA_ACT)
        {
          for(unsigned int kk=0;kk<PosIDNum;kk++)  //copy PosID
          {
            TPSynG[PSynGSize+kk].PreID=PGatGSize;
            TPSynG[PSynGSize+kk].PosID=PosID[kk];
          }
          PSynFstGatG[PGatGSize].PreID=k;
          PSynFstGatG[PGatGSize].PosID=PosID[0];
          PSynFstGatG[PGatGSize].g=g;
          PSynFstGatG[PGatGSize].st=0.0f;
          PSynGSize+=PosIDNum;
          PGatGSize++;
        }
        else if(Act==ACH_ACT)
        {
          for(unsigned int kk=0;kk<PosIDNum;kk++)  //copy PosID
          {
            TPSynH[PSynHSize+kk].PreID=PGatHSize;
            TPSynH[PSynHSize+kk].PosID=PosID[kk];
          }
          PSynFstGatH[PGatHSize].PreID=k;
          PSynFstGatH[PGatHSize].PosID=PosID[0];
          PSynFstGatH[PGatHSize].g=g;
          PSynFstGatH[PGatHSize].st=0.0f;
          PSynHSize+=PosIDNum;
          PGatHSize++;
        }
        else if(Act==GCL_ACT)
        {
          for(unsigned int kk=0;kk<PosIDNum;kk++)  //copy PosID
          {
            TPSynL[PSynLSize+kk].PreID=PGatLSize;
            TPSynL[PSynLSize+kk].PosID=PosID[kk];
          }
          PSynFstGatL[PGatLSize].PreID=k;
          PSynFstGatL[PGatLSize].PosID=PosID[0];
          PSynFstGatL[PGatLSize].g=g;
          PSynFstGatL[PGatLSize].st=0.0f;
          PSynLSize+=PosIDNum;
          PGatLSize++;
        }
        else if(Act==NMDA_ACT)
        {
          for(unsigned int kk=0;kk<PosIDNum;kk++)  //copy PosID
          {
            TPSynN[PSynNSize+kk].PreID=PGatNSize;
            TPSynN[PSynNSize+kk].PosID=PosID[kk];
          }
          PSynSloGat[PGatNSize].PreID=k;
          PSynSloGat[PGatNSize].PosID=PosID[0];
          PSynSloGat[PGatNSize].g=g;
          PSynSloGat[PGatNSize].xt=0.0f;
          PSynSloGat[PGatNSize].st=0.0f;
          PSynNSize+=PosIDNum;
          PGatNSize++;
        }
        else if(Act==GAP_ACT)
        {
          for(unsigned int kk=0;kk<PosIDNum;kk++)  //copy PosID
          {
            TPSynCp[PSynCpSize+kk].PreID=PGapPSize;
            TPSynCp[PSynCpSize+kk].PosID=PosID[kk];
          }
          NeuMapPCndO[PGapPSize].Beg=PSynCpSize;
          NeuMapPCndO[PGapPSize].End=PSynCpSize+PosIDNum-1;

          PSynCndP[PGapPSize].PreID=k;
          PSynCndP[PGapPSize].g=g;

          PSynCpSize+=PosIDNum;
          PGapPSize++;
        }


        delete [] PosID;
      }

    }
  }
//gap junction: integrate PosID
  PSynCo=new unsigned int[PSynCpSize];
  for(unsigned int i=0;i<PSynCpSize;i++)
  { PSynCo[i]=TPSynCp[i].PosID;}

  sort(TPSynA,TPSynA+PSynASize,[](SynIdx a,SynIdx b){return a.PosID < b.PosID;});
  sort(TPSynG,TPSynG+PSynGSize,[](SynIdx a,SynIdx b){return a.PosID < b.PosID;});
  sort(TPSynH,TPSynH+PSynHSize,[](SynIdx a,SynIdx b){return a.PosID < b.PosID;});
  sort(TPSynL,TPSynL+PSynLSize,[](SynIdx a,SynIdx b){return a.PosID < b.PosID;});
  sort(TPSynN,TPSynN+PSynNSize,[](SynIdx a,SynIdx b){return a.PosID < b.PosID;});
  sort(TPSynCp,TPSynCp+PSynCpSize,[](SynIdx a,SynIdx b){return a.PosID < b.PosID;});

  PSynA = new unsigned int[PSynASize];
  PSynG = new unsigned int[PSynGSize];
  PSynH = new unsigned int[PSynHSize];
  PSynL = new unsigned int[PSynLSize];
  PSynN = new unsigned int[PSynNSize];
  PSynCp = new unsigned int[PSynCpSize];

  for(unsigned int i=0;i<PSynASize;i++)
  {  PSynA[i]=TPSynA[i].PreID;  }

  for(unsigned int i=0;i<PSynGSize;i++)
  {  PSynG[i]=TPSynG[i].PreID;  }

  for(unsigned int i=0;i<PSynHSize;i++)
  {  PSynH[i]=TPSynH[i].PreID;  }

  for(unsigned int i=0;i<PSynLSize;i++)
  {  PSynL[i]=TPSynL[i].PreID;  }

  for(unsigned int i=0;i<PSynNSize;i++)
  {  PSynN[i]=TPSynN[i].PreID;  }

  for(unsigned int i=0;i<PSynCpSize;i++)
  {  PSynCp[i]=TPSynCp[i].PreID;  }


  NeuMapPAMPA=BulBiasDNeuMap<SynIdx>(TPSynA,0,PSynASize,PAMPACntSize);
  NeuMapPGABA=BulBiasDNeuMap<SynIdx>(TPSynG,0,PSynGSize,PGABACntSize);
  NeuMapPACH =BulBiasDNeuMap<SynIdx>(TPSynH,0,PSynHSize, PACHCntSize);
  NeuMapPGCL =BulBiasDNeuMap<SynIdx>(TPSynL,0,PSynLSize, PGCLCntSize);
  NeuMapPNMDA=BulBiasDNeuMap<SynIdx>(TPSynN,0,PSynNSize,PNMDACntSize);
  NeuMapPCndP=BulBiasDNeuMap<SynIdx>(TPSynCp,0,PSynCpSize,PCndPCntSize);

  MapPAPosID = new unsigned int[PAMPACntSize];
  MapPGPosID = new unsigned int[PGABACntSize];
  MapPHPosID = new unsigned int[PACHCntSize];
  MapPLPosID = new unsigned int[PGCLCntSize];
  MapPNPosID = new unsigned int[PNMDACntSize];
  MapPCPosID = new unsigned int[PCndPCntSize];

  for(unsigned int NID=0;NID<PAMPACntSize;NID++)
  { MapPAPosID[NID]=TPSynA[NeuMapPAMPA[NID].Beg].PosID; }

  for(unsigned int NID=0;NID<PGABACntSize;NID++)
  { MapPGPosID[NID]=TPSynG[NeuMapPGABA[NID].Beg].PosID; }

  for(unsigned int NID=0;NID<PACHCntSize;NID++)
  { MapPHPosID[NID]=TPSynH[NeuMapPACH[NID].Beg].PosID; }

  for(unsigned int NID=0;NID<PGCLCntSize;NID++)
  { MapPLPosID[NID]=TPSynL[NeuMapPACH[NID].Beg].PosID; }

  for(unsigned int NID=0;NID<PNMDACntSize;NID++)
  { MapPNPosID[NID]=TPSynN[NeuMapPNMDA[NID].Beg].PosID; }

  for(unsigned int NID=0;NID<PCndPCntSize;NID++)
  { MapPCPosID[NID]=TPSynCp[NeuMapPCndP[NID].Beg].PosID; }

  {
    delete [] TPSynA;
    delete [] TPSynG;
    delete [] TPSynH;
    delete [] TPSynL;
    delete [] TPSynN;
    delete [] TPSynCp;
  };

  TolSynp=PSynASize+PSynGSize+PSynHSize+PSynLSize+PSynNSize+PSynCpSize;


//------------------------------------------------------------------
  unsigned int SynAMPA1Size=0;
  unsigned int SynGABA1Size=0;
  unsigned int SynACH1Size=0;
  unsigned int SynGCL1Size=0;


  for(unsigned int i=0;i<ConfSta.size();i++) // pre synatic population
  {
    cout<<"calculate synapses size:"<<100.0f*(float)(i+1)/ConfSta.size()<<"%        \r"<<flush;

// calculate synapses size
    for(unsigned int k=NIdx[i];k<(NIdx[i]+ConfSta[i].NeuNum);k++) // build Fast synapse of j neuron post links
    {
      for(unsigned int j=ConfSta[i].PpLinkBeg,PosIDNum=0;j<ConfSta[i].PpLinkBeg+ConfSta[i].PpLinkNum;j++) // pre synaptic population synapse
      {

        PosIDNum=(unsigned int)ceil((double)ConfSta[POSYNPP].NeuNum * POSYNCON);
        if(PosIDNum == 1)
        {
          if(GetAct(POSYNREP) == AMPA_ACT) //Fast synapse
            SynAMPA1Size++;

          if(GetAct(POSYNREP) == GABA_ACT)
            SynGABA1Size++;

          if(GetAct(POSYNREP) == ACH_ACT)
            SynACH1Size++;

          if(GetAct(POSYNREP) == GCL_ACT)
            SynGCL1Size++;

          if(GetAct(POSYNREP) == NMDA_ACT) //Slow synapse
            SynNMDA1Size++;

          if((GetAct(POSYNREP) == GAP_ACT) || (GetAct(POSYNREP) == CUBA_ACT))  //Gap junction or current base synapse
            SynCnd1Size++;
        }
      }
    }

  }
  cout<<"\n";

// build up synapses
  SynAMPA1End=0+SynAMPA1Size;
  SynGABA1End=SynAMPA1End+SynGABA1Size;
  SynACH1End=SynGABA1End+SynACH1Size;
  SynGCL1End=SynACH1End+SynGCL1Size;

  if(SynGCL1End!=0)
  {  SynFst1 = new SynFast2[SynGCL1End]; }

  if(SynNMDA1Size!=0)
  {  SynSlo1 = new SynSlow2[SynNMDA1Size]; }

  if(SynCnd1Size!=0)
  {  SynCnd1 = new SynGap2[SynCnd1Size]; }


  unsigned int kAMPA1=0,kGABA1=0,kACH1=0,kGCL1=0,kNMDA1=0,kGAP1=0;
  for(unsigned int i=0;i<ConfSta.size();i++) // pre synatic population
  {
    cout<<"setup network:"<<100.0f*(float)(i+1)/ConfSta.size()<<"%        \r"<<flush;

//setup network
    for(unsigned int k=NIdx[i];k<(NIdx[i]+ConfSta[i].NeuNum);k++) // build Fast synapse of j neuron post links
    {
      for(unsigned int j=ConfSta[i].PpLinkBeg,PosIDNum=0;j<ConfSta[i].PpLinkBeg+ConfSta[i].PpLinkNum;j++) // pre synaptic population synapse
      {
        PosIDNum=(unsigned int)ceil((double)ConfSta[POSYNPP].NeuNum * POSYNCON);
        if(PosIDNum == 1)
        {
          if(GetAct(POSYNREP) == AMPA_ACT) //AMPA synapse
          {
            SynFst1[kAMPA1].PreID=i;
            SynFst1[kAMPA1].PosID=PosBuild1(k,PosIDNum,NIdx[POSYNPP],ConfSta[POSYNPP].NeuNum,ConfSta[i].SelConEn);
            SynFst1[kAMPA1].g=POSYNEFF*POSYNWGT;

            kAMPA1++;
          }
          else if(GetAct(POSYNREP) == GABA_ACT) //GABA synapse
          {
            SynFst1[kGABA1+SynAMPA1End].PreID=i;
            SynFst1[kGABA1+SynAMPA1End].PosID=PosBuild1(k,PosIDNum,NIdx[POSYNPP],ConfSta[POSYNPP].NeuNum,ConfSta[i].SelConEn);
            SynFst1[kGABA1+SynAMPA1End].g=POSYNEFF*POSYNWGT;

            kGABA1++;
          }
          else if(GetAct(POSYNREP) == ACH_ACT) //ACH synapse
          {
            SynFst1[kACH1+SynGABA1End].PreID=i;
            SynFst1[kACH1+SynGABA1End].PosID=PosBuild1(k,PosIDNum,NIdx[POSYNPP],ConfSta[POSYNPP].NeuNum,ConfSta[i].SelConEn);
            SynFst1[kACH1+SynGABA1End].g=POSYNEFF*POSYNWGT;

            kACH1++;
          }
          else if(GetAct(POSYNREP) == GCL_ACT) //GCL synapse
          {
            SynFst1[kGCL1+SynACH1End].PreID=i;
            SynFst1[kGCL1+SynACH1End].PosID=PosBuild1(k,PosIDNum,NIdx[POSYNPP],ConfSta[POSYNPP].NeuNum,ConfSta[i].SelConEn);
            SynFst1[kGCL1+SynACH1End].g=POSYNEFF*POSYNWGT;

            kGCL1++;
          }
          else if(GetAct(POSYNREP) == NMDA_ACT) //NMDA synapse
          {
            SynSlo1[kNMDA1].PreID=i;
            SynSlo1[kNMDA1].PosID=PosBuild1(k,PosIDNum,NIdx[POSYNPP],ConfSta[POSYNPP].NeuNum,ConfSta[i].SelConEn);
            SynSlo1[kNMDA1].g=POSYNEFF*POSYNWGT;

            kNMDA1++;
          }
          else if(GetAct(POSYNREP) == GAP_ACT) //GAP synapse
          {
            SynCnd1[kGAP1].PreID=i;
            SynCnd1[kGAP1].PosID=PosBuild1(k,PosIDNum,NIdx[POSYNPP],ConfSta[POSYNPP].NeuNum,ConfSta[i].SelConEn);
            SynCnd1[kGAP1].g=POSYNEFF*POSYNWGT;

            kGAP1++;
          }
        }
      }
    }

  }

  delete [] PostPp;

  PostSetPp();
  SetOut();
}  //end of SetPp()


void NetSimGnl2::PostSetPp(void)
{
  sort(SynFst1,             SynFst1+SynAMPA1End,[](SynFast2 a,SynFast2 b){return a.PosID < b.PosID;});
  sort(SynFst1+SynAMPA1End, SynFst1+SynGABA1End,[](SynFast2 a,SynFast2 b){return a.PosID < b.PosID;});
  sort(SynFst1+SynGABA1End, SynFst1+SynACH1End, [](SynFast2 a,SynFast2 b){return a.PosID < b.PosID;});
  sort(SynFst1+SynACH1End,  SynFst1+SynGCL1End, [](SynFast2 a,SynFast2 b){return a.PosID < b.PosID;});

  sort(SynSlo1, SynSlo1+SynNMDA1Size,[](SynSlow2 a,SynSlow2 b){return a.PosID < b.PosID;});
  sort(SynCnd1, SynCnd1+SynCnd1Size,[](SynGap2 a,SynGap2 b){return a.PosID < b.PosID;});

  NeuMapAMPA1=BulBiasDNeuMap<SynFast2>(SynFst1,0,          SynAMPA1End,AMPA1CntSize);
  NeuMapGABA1=BulBiasDNeuMap<SynFast2>(SynFst1,SynAMPA1End,SynGABA1End,GABA1CntSize);
  NeuMapACH1 =BulBiasDNeuMap<SynFast2>(SynFst1,SynGABA1End,SynACH1End, ACH1CntSize);
  NeuMapGCL1 =BulBiasDNeuMap<SynFast2>(SynFst1,SynACH1End, SynGCL1End, GCL1CntSize);

  NeuMapNMDA1=BulBiasDNeuMap<SynSlow2>(SynSlo1,0,SynNMDA1Size,NMDA1CntSize);
  NeuMapCnd1=BulBiasDNeuMap<SynGap2>(SynCnd1,0,SynCnd1Size,Cnd1CntSize);

  PreMapCnd1=new unsigned int[SynCnd1Size];
  NeuPreMapCnd1=BulDPreNeuMap<SynGap2>(SynCnd1,0,SynCnd1Size,PreCnd1CntSize,TolNeu,PreMapCnd1);

  TolSynp+=SynGCL1End+SynNMDA1Size+SynCnd1Size;

}


