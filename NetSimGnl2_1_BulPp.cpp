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

void NetSimGnl2_1::SetPp(int thNum,int sop)
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
  RGsSize=q;

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
  RGs = new NeuGatVar3[RGsSize];

  for(unsigned int i=0;i<NetSize;i++)
  { Nets[i].ActCtr=ActCtr; }


  NetPBse = new ParBase2[NetSize];
  NetVBse = new VarBase[NetSize];

  NetPLIF = new ParLIF[NetSize];
  NetVLIF = new VarLIF2[NetSize];

  if((ActCtr & STP_ACT) != 0 || (ActCtr & STD_ACT) != 0)
  {
    NetPSTP = new ParStp2[NetSize];
    NetVSTP = new VarStp[NetSize];
  }

  if((ActCtr & LTP_ACT) != 0)
  {
    NetPLTP = new ParLP2[NetSize];
    NetVLTP = new VarLP2[NetSize];
  }

  if(
      ((ActCtr & HH_ACT) != 0 )||
      ((ActCtr & SOD_ACT) != 0 )
    )
  {
    NetPHH = new ParHH2[NetSize];
    NetVHH = new VarHH[NetSize];
  }



//back ground noise setup
  zig=new ziggurat[thsize];
  BgSt = new BgSti1<float,RV4>[NetSize];
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
      BgSt[i].SynFst1ES[j]=0.0f;
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



//each neuron variable and parameter setup
  for(unsigned int i=0;i<ConfSta.size();i++)
  {
    for(unsigned int j=NIdx[i];j<(NIdx[i]+ConfSta[i].NeuNum);j++) // set up each neuron's parameter in every population
    {
      //Nets[j].ActCtr=ActCtr;

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

      //HH or SODCH
      if(
          ((Nets[j].ActCtr & HH_ACT) != 0) ||
          ((Nets[j].ActCtr & SOD_ACT) != 0)
        )
      {
        NetPHH[j].Ek=ConfHH[i].Ek;
        NetPHH[j].gk=ConfHH[i].gk;
        NetPHH[j].Ena=ConfHH[i].Ena;
        NetPHH[j].gna=ConfHH[i].gna;

        NetPHH[j].AnA=ConfHH[i].AnA;
        NetPHH[j].AnB=ConfHH[i].AnB;
        NetPHH[j].AnC=ConfHH[i].AnC;
        NetPHH[j].AnD=ConfHH[i].AnD;

        NetPHH[j].BnA=ConfHH[i].BnA;
        NetPHH[j].BnB=ConfHH[i].BnB;
        NetPHH[j].BnC=ConfHH[i].BnC;
        NetPHH[j].BnD=ConfHH[i].BnD;

        NetPHH[j].AmA=ConfHH[i].AmA;
        NetPHH[j].AmB=ConfHH[i].AmB;
        NetPHH[j].AmC=ConfHH[i].AmC;
        NetPHH[j].AmD=ConfHH[i].AmD;

        NetPHH[j].BmA=ConfHH[i].BmA;
        NetPHH[j].BmB=ConfHH[i].BmB;
        NetPHH[j].BmC=ConfHH[i].BmC;
        NetPHH[j].BmD=ConfHH[i].BmD;

        NetPHH[j].AhA=ConfHH[i].AhA;
        NetPHH[j].AhB=ConfHH[i].AhB;
        NetPHH[j].AhC=ConfHH[i].AhC;
        NetPHH[j].AhD=ConfHH[i].AhD;

        NetPHH[j].BhA=ConfHH[i].BhA;
        NetPHH[j].BhB=ConfHH[i].BhB;
        NetPHH[j].BhC=ConfHH[i].BhC;
        NetPHH[j].BhD=ConfHH[i].BhD;
      }


      //ltp
      if((Nets[j].ActCtr & LTP_ACT) != 0)
      {
        NetPLTP[j].dfLP    =  exp(-dt/ConfLtp[i].tLP);
        NetPLTP[j].SpiLPDy = ConfLtp[i].SpiLPDy; //back propogation spike delay
        NetPLTP[j].PosISI  = ConfLtp[i].PosISI;
        NetPLTP[j].NegISI  = ConfLtp[i].NegISI;
      }

      //stp,std
      if((Nets[j].ActCtr & STP_ACT) != 0 || (Nets[j].ActCtr & STD_ACT) != 0)
      {
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
      }

      Nets[j].Irep = new float[5];

      for(unsigned int ix=0;ix<5;ix++)
      { Nets[j].Irep[ix]=0.0f; }

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
          if(PreSynR->FExt==0.0f)
            Nets[j].ActCtr &= ~AMPA_ACT;
          else
          {
            Nets[j].ActCtr |= AMPA_ACT;
            BgSt[j].pb[0]=PreSynR->FExt*dt*1e-3f;
          }

          BgSt[j].SynFst1ES[0]=PreSynR->MnExtEff*PreSynR->MnExtCon;
        }
        else if(PreSynR->NeuRepName==KW_GABA)
        {
          NetPBse[j].dVmax = ((Confs+i)->Vth - PreSynR->Vrev)/NetPBse[j].ef;

          RPs[j].ParGABA.Vrev = PreSynR->Vrev;
          RPs[j].ParGABA.dfst = exp(-dt/PreSynR->Tau);
          if(PreSynR->FExt==0.0f)
            Nets[j].ActCtr &= ~GABA_ACT;
          else
          {
            Nets[j].ActCtr |= GABA_ACT;
            BgSt[j].pb[1]=PreSynR->FExt*dt*1e-3f;
          }

          BgSt[j].SynFst1ES[1]=PreSynR->MnExtEff*PreSynR->MnExtCon;
        }
        else if(PreSynR->NeuRepName==KW_ACH)
        {
          RPs[j].ParACH.Vrev = PreSynR->Vrev;
          RPs[j].ParACH.dfst = exp(-dt/PreSynR->Tau);
          if(PreSynR->FExt==0.0f)
            Nets[j].ActCtr &= ~ACH_ACT;
          else
          {
            Nets[j].ActCtr |= ACH_ACT;
            BgSt[j].pb[2]=PreSynR->FExt*dt*1e-3f;
          }

          BgSt[j].SynFst1ES[2]=PreSynR->MnExtEff*PreSynR->MnExtCon;
        }
        else if(PreSynR->NeuRepName==KW_GCL)
        {
          RPs[j].ParGCL.Vrev = PreSynR->Vrev;
          RPs[j].ParGCL.dfst = exp(-dt/PreSynR->Tau);
          if(PreSynR->FExt==0.0f)
            Nets[j].ActCtr &= ~GCL_ACT;
          else
          {
            Nets[j].ActCtr |= GCL_ACT;
            BgSt[j].pb[3]=PreSynR->FExt*dt*1e-3f;
          }

          BgSt[j].SynFst1ES[3]=PreSynR->MnExtEff*PreSynR->MnExtCon;
        }
        else if(PreSynR->NeuRepName==KW_NMDA)
        {
          RPs[j].ParNMDA.as=0.6332f; // the alpha st
          RPs[j].ParNMDA.dfxt = exp(-dt/PreSynR->TauXT);
          RPs[j].ParNMDA.dfst = exp(-dt/PreSynR->Tau);
          RPs[j].ParNMDA.efst = PreSynR->Tau*(1.0f-RPs[j].ParNMDA.dfst);
          RPs[j].ParNMDA.Vrev = PreSynR->Vrev;
          RPs[j].ParNMDA.mg=1.0f; // the concentration of magnesium

          if(PreSynR->FExt==0.0f)
            Nets[j].ActCtr &= ~NMDA_ACT;
          else
          {
            Nets[j].ActCtr |= NMDA_ACT;
            BgSt[j].pb[4]=PreSynR->FExt*dt*1e-3f;
          }
          BgSt[j].SynSlo1ES.g=PreSynR->MnExtEff*PreSynR->MnExtCon;
        }

      }

    }

  }

  delete [] ConfStp;
  delete [] ConfLtp;
  delete [] ConfHH;



//build and flatten population into neuorn level------------------------------------
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
        if((i==POSYNPP) && (ConfSta[POSYNPP].NeuNum==PosIDNum) && !ConfSta[i].SelConEn) // no self connction
        { PosIDNum--; }

        if(PosIDNum == 0){}
        else if(PosIDNum == 1)
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
        else
        {
          if(GetAct(POSYNREP) == AMPA_ACT) //Fast synapse
            SynAMPA1Size+=PosIDNum;

          if(GetAct(POSYNREP) == GABA_ACT)
            SynGABA1Size+=PosIDNum;

          if(GetAct(POSYNREP) == ACH_ACT)
            SynACH1Size+=PosIDNum;

          if(GetAct(POSYNREP) == GCL_ACT)
            SynGCL1Size+=PosIDNum;

          if(GetAct(POSYNREP) == NMDA_ACT) //Slow synapse
            SynNMDA1Size+=PosIDNum;

          if((GetAct(POSYNREP) == GAP_ACT) || (GetAct(POSYNREP) == CUBA_ACT))  //Gap junction or current base synapse
            SynCnd1Size+=PosIDNum;
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

//cout<<"1:"<<SynAMPA1Size<<","<<SynGABA1Size<<","<<SynACH1Size<<","<<SynGCL1Size<<","<<SynNMDA1Size<<endl;


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
        if((i==POSYNPP) && (ConfSta[POSYNPP].NeuNum==PosIDNum) && !ConfSta[i].SelConEn) // no self connction
        { PosIDNum--; }
//cout<<"1.1"<<k<<","<<i<<endl;
        if(PosIDNum == 0){}
        else if(PosIDNum == 1)
        {
          if(GetAct(POSYNREP) == AMPA_ACT) //AMPA synapse
          {
            SynFst1[kAMPA1].PreID=k;
            SynFst1[kAMPA1].PosID=PosBuild1(k,PosIDNum,NIdx[POSYNPP],ConfSta[POSYNPP].NeuNum,ConfSta[i].SelConEn);
            SynFst1[kAMPA1].g=POSYNEFF*POSYNWGT;

            kAMPA1++;
          }
          else if(GetAct(POSYNREP) == GABA_ACT) //GABA synapse
          {
            SynFst1[kGABA1+SynAMPA1End].PreID=k;
            SynFst1[kGABA1+SynAMPA1End].PosID=PosBuild1(k,PosIDNum,NIdx[POSYNPP],ConfSta[POSYNPP].NeuNum,ConfSta[i].SelConEn);
            SynFst1[kGABA1+SynAMPA1End].g=POSYNEFF*POSYNWGT;

            kGABA1++;
          }
          else if(GetAct(POSYNREP) == ACH_ACT) //ACH synapse
          {
            SynFst1[kACH1+SynGABA1End].PreID=k;
            SynFst1[kACH1+SynGABA1End].PosID=PosBuild1(k,PosIDNum,NIdx[POSYNPP],ConfSta[POSYNPP].NeuNum,ConfSta[i].SelConEn);
            SynFst1[kACH1+SynGABA1End].g=POSYNEFF*POSYNWGT;

            kACH1++;
          }
          else if(GetAct(POSYNREP) == GCL_ACT) //GCL synapse
          {
            SynFst1[kGCL1+SynACH1End].PreID=k;
            SynFst1[kGCL1+SynACH1End].PosID=PosBuild1(k,PosIDNum,NIdx[POSYNPP],ConfSta[POSYNPP].NeuNum,ConfSta[i].SelConEn);
            SynFst1[kGCL1+SynACH1End].g=POSYNEFF*POSYNWGT;

            kGCL1++;
          }
          else if(GetAct(POSYNREP) == NMDA_ACT) //NMDA synapse
          {
            SynSlo1[kNMDA1].PreID=k;
            SynSlo1[kNMDA1].PosID=PosBuild1(k,PosIDNum,NIdx[POSYNPP],ConfSta[POSYNPP].NeuNum,ConfSta[i].SelConEn);
            SynSlo1[kNMDA1].g=POSYNEFF*POSYNWGT;

            kNMDA1++;
          }
          else if(GetAct(POSYNREP) == GAP_ACT) //GAP synapse
          {
            SynCnd1[kGAP1].PreID=k;
            SynCnd1[kGAP1].PosID=PosBuild1(k,PosIDNum,NIdx[POSYNPP],ConfSta[POSYNPP].NeuNum,ConfSta[i].SelConEn);
            SynCnd1[kGAP1].g=POSYNEFF*POSYNWGT;

            kGAP1++;
          }
        }
        else //PosIDNum > 1 in ecah population
        {
          if(GetAct(POSYNREP) == AMPA_ACT) //AMPA synapse
          {
            unsigned int *p=PosBuild(k,PosIDNum,NIdx[POSYNPP],ConfSta[POSYNPP].NeuNum,ConfSta[i].SelConEn);
            for(unsigned kk=0;kk<PosIDNum;kk++)
            {
              SynFst1[kAMPA1+kk].PreID=k;
              SynFst1[kAMPA1+kk].PosID=p[kk];
              SynFst1[kAMPA1+kk].g=POSYNEFF*POSYNWGT;
            }
            kAMPA1+=PosIDNum;
            delete [] p;
          }
          else if(GetAct(POSYNREP) == GABA_ACT) //GABA synapse
          {
            unsigned int *p=PosBuild(k,PosIDNum,NIdx[POSYNPP],ConfSta[POSYNPP].NeuNum,ConfSta[i].SelConEn);
            for(unsigned kk=0;kk<PosIDNum;kk++)
            {
              SynFst1[kGABA1+SynAMPA1End+kk].PreID=k;
              SynFst1[kGABA1+SynAMPA1End+kk].PosID=p[kk];
              SynFst1[kGABA1+SynAMPA1End+kk].g=POSYNEFF*POSYNWGT;
            }
            kGABA1+=PosIDNum;
            delete [] p;
          }
          else if(GetAct(POSYNREP) == ACH_ACT) //ACH synapse
          {
            //cout<<"2:"<<k<<","<<PosIDNum<<","<<NIdx[POSYNPP]<<","<<ConfSta[POSYNPP].NeuNum<<endl;

            unsigned int *p=PosBuild(k,PosIDNum,NIdx[POSYNPP],ConfSta[POSYNPP].NeuNum,ConfSta[i].SelConEn);
            for(unsigned kk=0;kk<PosIDNum;kk++)
            {
              SynFst1[kACH1+SynGABA1End+kk].PreID=k;
              SynFst1[kACH1+SynGABA1End+kk].PosID=p[kk];
              SynFst1[kACH1+SynGABA1End+kk].g=POSYNEFF*POSYNWGT;
            }
            kACH1+=PosIDNum;
            delete [] p;
          }
          else if(GetAct(POSYNREP) == GCL_ACT) //GCL synapse
          {
            unsigned int *p=PosBuild(k,PosIDNum,NIdx[POSYNPP],ConfSta[POSYNPP].NeuNum,ConfSta[i].SelConEn);
            for(unsigned kk=0;kk<PosIDNum;kk++)
            {
              SynFst1[kGCL1+SynACH1End+kk].PreID=k;
              SynFst1[kGCL1+SynACH1End+kk].PosID=p[kk];
              SynFst1[kGCL1+SynACH1End+kk].g=POSYNEFF*POSYNWGT;
            }
            kGCL1+=PosIDNum;
            delete [] p;
          }
          else if(GetAct(POSYNREP) == NMDA_ACT) //NMDA synapse
          {
            //cout<<"3:"<<k<<","<<PosIDNum<<","<<NIdx[POSYNPP]<<","<<ConfSta[POSYNPP].NeuNum<<endl;
            unsigned int *p=PosBuild(k,PosIDNum,NIdx[POSYNPP],ConfSta[POSYNPP].NeuNum,ConfSta[i].SelConEn);
            for(unsigned kk=0;kk<PosIDNum;kk++)
            {
              SynSlo1[kNMDA1+kk].PreID=k;
              SynSlo1[kNMDA1+kk].PosID=p[kk];
              SynSlo1[kNMDA1+kk].g=POSYNEFF*POSYNWGT;
            }
            kNMDA1+=PosIDNum;
            delete [] p;
          }
          else if(GetAct(POSYNREP) == GAP_ACT) //GAP synapse
          {
            unsigned int *p=PosBuild(k,PosIDNum,NIdx[POSYNPP],ConfSta[POSYNPP].NeuNum,ConfSta[i].SelConEn);
            for(unsigned kk=0;kk<PosIDNum;kk++)
            {
              SynFst1[kGAP1+kk].PreID=k;
              SynFst1[kGAP1+kk].PosID=p[kk];
              SynFst1[kGAP1+kk].g=POSYNEFF*POSYNWGT;
            }
            kGAP1+=PosIDNum;
            delete [] p;
          }
        }
      }
    }
  }

  delete [] PostPp;

  PostSetPp();
  SetOut();
}  //end of SetPp()


void NetSimGnl2_1::PostSetPp(void)
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


