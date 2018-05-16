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
#include "NetSim.h"



using namespace std;

void NetSimGnl2_1::SInfo(string& ostr) //short information
{
  ostr_tmp
      <<"\ntime step:"<<dt
      <<", total neurons="<<TolNeu
      <<", total synapses="<<TolSynp<<"\n";

  ostr=ostr_tmp.str();
  ostr_tmp.str("");
}


void NetSimGnl2_1::TexInfo(string& ostr)
{
  string str;
//-------show total synapse inforation
  SInfo(str);
  ostr_tmp<<str<<"\n";
  ostr_tmp<<"population information:-------------------------------------------\n";
  for(unsigned int i=0;i<TolNeu;i++)
  {
    (Nets+i)->Info(str);
    ostr_tmp<<"neuron:"<<i<<"\n"<<str;
  }

  ostr_tmp<<"time event information:-------------------------------------------\n";
  for(unsigned int i=0;i<Evts.size();i++)
  {
    Evts[i].Info(str);
    ostr_tmp<<str;
  }

  ostr_tmp<<"output files information:-----------------------------------------\n";
  for(unsigned int i=0;i<OFCtr.size();i++)
  {
    OFCtr[i].Info(str);
    ostr_tmp<<str;
  }

  ostr=ostr_tmp.str();
  ostr_tmp.str("");
}

void NetSimGnl2_1::logDataWrite(const char *logName)
{
  string s1;
  ofstream logdata;
  logdata.open(logName);

//-------show total information with large network
  SInfo(s1);
  logdata<<s1<<"\n";

  if((ActCtr & UDFRDSEED_ACT)!= 0)
    logdata<<"User define random seed enable!\n";
  else
    logdata<<"internal random seed enable!\n";

  logdata<<"repeats="<<RpSize<<",each repeat random seed=";

  for(unsigned int i=0;i<(unsigned int)RpSize;i++)
  { logdata<<" "<<rseeds[i]; }

  logdata<<"\n\npopulation information:-------------------------------------------\n";
  for(unsigned int i=0;i<TolNeu;i++)
  {
    logdata<<"-------------------------------------------\n"
           <<"neuron:"<<i<<"\n";
    logdata
    <<"ActCtr="<<Nets[i].ActCtr
    <<" Vth="<<NetPBse[i].Vth<<"mv"
    <<" El="<<NetPBse[i].El<<"mv"
    <<" g="<<NetPBse[i].gl<<"uS"
    <<" Cm="<<NetPBse[i].Cm<<"nF"
    <<" Vreset="<<NetPBse[i].Vreset<<"mv"
    <<" df="<<NetPBse[i].df
    <<" ef="<<NetPBse[i].ef
    <<" dVmax="<<NetPBse[i].dVmax<<"mv\n";

    logdata
    <<"AMPA df="<<RPs[i].ParAMPA.dfst
    <<" Vrev="<<RPs[i].ParAMPA.Vrev<<"mV"
    <<" GABA df="<<RPs[i].ParGABA.dfst
    <<" Vrev="<<RPs[i].ParGABA.Vrev<<"mV\n"

    <<"NMDA dfxt="<<RPs[i].ParNMDA.dfxt
    <<" dfst="<<RPs[i].ParNMDA.dfst
    <<" efst="<<RPs[i].ParNMDA.efst
    <<" Vrev="<<RPs[i].ParNMDA.Vrev<<"mV"
    <<" as="<<RPs[i].ParNMDA.as
    <<" mg="<<RPs[i].ParNMDA.mg<<"\n"

    <<"ACH df="<<RPs[i].ParACH.dfst
    <<" Vrev="<<RPs[i].ParACH.Vrev<<"mV"
    <<" GCL df="<<RPs[i].ParGCL.dfst
    <<" Vrev="<<RPs[i].ParGCL.Vrev<<"mV\n";

//LIF
    if((Nets[i].ActCtr&LIF_ACT) !=0)
    {
      logdata
      <<"LIF:"
      <<" RfPed="<<NetPLIF[i].RfPed*dt<<"ms"
      <<" SpiDy="<<NetPLIF[i].SpiDy*dt<<"ms\n";
    }

//HH
    if((Nets[i].ActCtr&HH_ACT) !=0)
    {
      logdata
      <<"HH:"
      <<" Ek="<<NetPHH[i].Ek<<" gk="<<NetPHH[i].gk<<" Ena="<<NetPHH[i].Ena<<" gna="<<NetPHH[i].gna<<"\n"
      <<" AnA="<<NetPHH[i].AnA<<" AnB="<<NetPHH[i].AnB<<" AnC="<<NetPHH[i].AnC<<" AnD="<<NetPHH[i].AnD<<"\n"
      <<" BnA="<<NetPHH[i].BnA<<" BnB="<<NetPHH[i].BnB<<" BnC="<<NetPHH[i].BnC<<" BnD="<<NetPHH[i].BnD<<"\n"
      <<" AmA="<<NetPHH[i].AmA<<" AmB="<<NetPHH[i].AmB<<" AmC="<<NetPHH[i].AmC<<" AmD="<<NetPHH[i].AmD<<"\n"
      <<" BmA="<<NetPHH[i].BmA<<" BmB="<<NetPHH[i].BmB<<" BmC="<<NetPHH[i].BmC<<" BmD="<<NetPHH[i].BmD<<"\n"
      <<" AhA="<<NetPHH[i].AhA<<" AhB="<<NetPHH[i].AhB<<" AhC="<<NetPHH[i].AhC<<" AhD="<<NetPHH[i].AhD<<"\n"
      <<" BhA="<<NetPHH[i].BhA<<" BhB="<<NetPHH[i].BhB<<" BhC="<<NetPHH[i].BhC<<" BhD="<<NetPHH[i].BhD<<"\n";
    }

//sodium channel
    if((Nets[i].ActCtr&SOD_ACT) !=0)
    {
      logdata
      <<"Sodium channel:"
      <<" VthSodium="<<NetPHH[i].BmB<<" Ena="<<NetPHH[i].Ena<<" gna="<<NetPHH[i].gna<<"\n";
    }

//STP
    if((Nets[i].ActCtr&STP_ACT) !=0)
    {
      logdata
      <<"STP:"
      <<" pv="<<NetPSTP[i].pv
      <<" dfD="<<NetPSTP[i].dfD
      <<" efD="<<NetPSTP[i].efD
      <<" aF="<<NetPSTP[i].aF
      <<" dfF="<<NetPSTP[i].dfF
      <<" efF="<<NetPSTP[i].efF
      <<" aCap="<<NetPSTP[i].aCap
      <<" dfCap="<<NetPSTP[i].dfCap
      <<" aCar="<<NetPSTP[i].aCar
      <<" dfCar="<<NetPSTP[i].dfCar
      <<" efCar="<<NetPSTP[i].efCar<<"\n";
    }
    else if((Nets[i].ActCtr&STD_ACT) !=0)
    {
      logdata
      <<"STD:"
      <<" pv="<<NetPSTP[i].pv
      <<" dfD="<<NetPSTP[i].dfD
      <<" efD="<<NetPSTP[i].efD<<"\n";
    }

//LTP
    if((Nets[i].ActCtr&LTP_ACT) !=0)
    {
      logdata
      <<"LTP:"
      <<" df="<<NetPLTP[i].dfLP
      <<" SpiLPDy="<<NetPLTP[i].SpiLPDy
      <<" PosISI="<<NetPLTP[i].PosISI
      <<" NegISI="<<NetPLTP[i].NegISI<<"\n";
    }

//background noise

    logdata<<"noise parameters:"<<"\n";

    logdata
      <<"Injection current:"
      <<"mean="<<BgSt[i].mean
      <<" stddev="<<BgSt[i].std
      <<"\n";
    logdata<<"Ext AMPA:prob="<<BgSt[i].pb[0]<<",g="<<BgSt[i].SynFst1ES[0]<<"\n";
    logdata<<"Ext GABA:prob="<<BgSt[i].pb[1]<<",g="<<BgSt[i].SynFst1ES[1]<<"\n";
    logdata<<"Ext Ach:prob="<<BgSt[i].pb[2]<<",g="<<BgSt[i].SynFst1ES[2]<<"\n";
    logdata<<"Ext GluCl:prob="<<BgSt[i].pb[3]<<",g="<<BgSt[i].SynFst1ES[3]<<"\n";
    logdata<<"Ext NMDA:prob="<<BgSt[i].pb[4]<<",g="<<BgSt[i].SynSlo1ES.g<<"\n";

    logdata<<flush;
  }

  logdata<<"\nsynapse connection information----------------------------------\n";


//axoaxonic plsticity
  if(AMPA1CntSize!=0)
  {
    logdata
    <<"\nSynAMPA1:TotalNum,"<<SynAMPA1End<<"\n"
    <<"PosID:Num:PreID,g;PreID,g;...\n";
    for(unsigned int NID=0;NID<AMPA1CntSize;NID++)
    {
      logdata<<SynFst1[NeuMapAMPA1[NID].Beg].PosID<<":"<<NeuMapAMPA1[NID].End-NeuMapAMPA1[NID].Beg+1<<":";
      for(unsigned int i=NeuMapAMPA1[NID].Beg,j=1;i<=NeuMapAMPA1[NID].End;i++,j++)  //AMPA
      {
        logdata<<SynFst1[i].PreID<<","<<SynFst1[i].g<<";";
        if(j%10==0) logdata<<"\n";
      }
      logdata<<endl;
    }
  }

  if(GABA1CntSize!=0)
  {
    logdata
    <<"\nSynGABA1:TotalNum,"<<SynGABA1End-SynAMPA1End<<"\n"
    <<"PosID:Num:PreID,g;PreID,g;...\n";
    for(unsigned int NID=0;NID<GABA1CntSize;NID++)
    {
      logdata<<SynFst1[NeuMapGABA1[NID].Beg].PosID<<":"<<NeuMapGABA1[NID].End-NeuMapGABA1[NID].Beg+1<<":";
      for(unsigned int i=NeuMapGABA1[NID].Beg,j=1;i<=NeuMapGABA1[NID].End;i++,j++)  //GABA
      {
        logdata<<SynFst1[i].PreID<<","<<SynFst1[i].g<<";";
        if(j%10==0) logdata<<"\n";
      }
      logdata<<endl;
    }
  }

  if(ACH1CntSize!=0)
  {
    logdata
    <<"\nSynACH1:TotalNum,"<<SynACH1End-SynGABA1End<<"\n"
    <<"PosID:Num:PreID,g;PreID,g;...\n";
    for(unsigned int NID=0;NID<ACH1CntSize;NID++)
    {
      logdata<<SynFst1[NeuMapACH1[NID].Beg].PosID<<":"<<NeuMapACH1[NID].End-NeuMapACH1[NID].Beg+1<<":";
      for(unsigned int i=NeuMapACH1[NID].Beg,j=1;i<=NeuMapACH1[NID].End;i++,j++)  //ACH
      {
        logdata<<SynFst1[i].PreID<<","<<SynFst1[i].g<<";";
        if(j%10==0) logdata<<"\n";
      }
      logdata<<endl;
    }
  }

  if(GCL1CntSize!=0)
  {
    logdata
    <<"\nSynGCL1:TotalNum,"<<SynGCL1End-SynACH1End<<"\n"
    <<"PosID:Num:PreID,g;PreID,g;...\n";
    for(unsigned int NID=0;NID<GCL1CntSize;NID++)
    {
      logdata<<SynFst1[NeuMapGCL1[NID].Beg].PosID<<":"<<NeuMapGCL1[NID].End-NeuMapGCL1[NID].Beg+1<<":";
      for(unsigned int i=NeuMapGCL1[NID].Beg,j=1;i<=NeuMapGCL1[NID].End;i++,j++)  //GCL
      {
        logdata<<SynFst1[i].PreID<<","<<SynFst1[i].g<<";";
        if(j%10==0) logdata<<"\n";
      }
      logdata<<endl;
    }
  }

  if(NMDA1CntSize!=0)
  {
    logdata
    <<"\nSynNMDA1:TotalNum,"<<SynNMDA1Size<<"\n"
    <<"PosID:Num:PreID,g;PreID,g;...\n";
    for(unsigned int NID=0;NID<NMDA1CntSize;NID++)
    {
      logdata<<SynSlo1[NeuMapNMDA1[NID].Beg].PosID<<":"<<NeuMapNMDA1[NID].End-NeuMapNMDA1[NID].Beg+1<<":";
      for(unsigned int i=NeuMapNMDA1[NID].Beg,j=1;i<=NeuMapNMDA1[NID].End;i++,j++)  //NMDA
      {
        logdata<<SynSlo1[i].PreID<<","<<SynSlo1[i].g<<";";
        if(j%10==0) logdata<<"\n";
      }
      logdata<<endl;
    }
  }

  if(Cnd1CntSize!=0)
  {
    logdata
    <<"\nSynCnd1:TotalNum;"<<SynCnd1Size<<"\n"
    <<"PosID:Num:PreID,g;PreID,g;...\n";
    for(unsigned int NID=0;NID<Cnd1CntSize;NID++)
    {
      logdata<<SynCnd1[NeuMapCnd1[NID].Beg].PosID<<":"<<NeuMapCnd1[NID].End-NeuMapCnd1[NID].Beg<<":";
      for(unsigned int i=NeuMapCnd1[NID].Beg,j=1;i<=NeuMapCnd1[NID].End;i++,j++)  //GAP
      {
        logdata<<SynCnd1[i].PreID<<","<<SynCnd1[i].g<<";";
        if(j%10==0) logdata<<"\n";
      }
      logdata<<endl;
    }
  }


  logdata<<"\ntime event information:-------------------------------------------\n";
  for(unsigned int i=0;i<Evts.size();i++)
  {
    Evts[i].Info(s1);
    logdata<<s1<<flush;
  }

//--output control header discription
  logdata<<"\noutput files information:-----------------------------------------\n";
  for(vector<OutCtr>::size_type i=0;i<OFCtr.size();i++)
  {
    OFCtr[i].Info(s1);
    logdata<<s1;
    switch(OFCtr[i].TypeID)
    {
      case IsynSepSumID:
      logdata<<"1:Time";

      for(unsigned int j=0,kk=2;j<OFCtr[i].ParList.size();j++)
      {
        for(unsigned int k=NIdx[OFCtr[i].ParList[j]];k<(NIdx[OFCtr[i].ParList[j]]+ConfSta[OFCtr[i].ParList[j]].NeuNum);k++)
        {
          logdata<<" NID="<<k;
          logdata<<" "<<kk++<<":AMPA ";
          logdata<<" "<<kk++<<":GABA ";
          logdata<<" "<<kk++<<":NMDA ";
          logdata<<" "<<kk++<<":ACH ";
          logdata<<" "<<kk++<<":GCL ";
        }
      }
      logdata<<"\n";
      break;

      case MemPotID:
      logdata<<"1:Time";

      for(unsigned int j=0,kk=2;j<OFCtr[i].ParList.size();j++)
      {
        for(unsigned int k=NIdx[OFCtr[i].ParList[j]];k<(NIdx[OFCtr[i].ParList[j]]+ConfSta[OFCtr[i].ParList[j]].NeuNum);k++)
        { logdata<<" "<<kk++<<":NID="<<k;}
      }
      logdata<<"\n";
      break;

      case FRateID:
      logdata<<"FiringRateWinodw="<<OFCtr[i].RWin<<"ms PrintStep="<<OFCtr[i].FPnt<<"ms\n";
      logdata<<"1:Time";

        if((OFCtr[i].Status & SUMPARSALL_OF) != 0)
        {
          logdata<<"2:total firing rate";
        }
        else
        {
          for(unsigned int j=0,kk=2;j<OFCtr[i].ParList.size();j++)
          { logdata<<" "<<kk++<<":NID="<<j;}
        }
      logdata<<"\n";
      break;

      case GatVarID:
      if(SynAMPA1End!=0)
      {
        PreMapAMPA1=new unsigned int[SynAMPA1End];
        NeuPreMapAMPA1=BulDPreNeuMap<SynFast2>(SynFst1,0,SynAMPA1End,PreAMPA1CntSize,TolNeu,PreMapAMPA1);
      }

      if((SynGABA1End-SynAMPA1End)!=0)
      {
        PreMapGABA1=new unsigned int[SynGABA1End-SynAMPA1End];
        NeuPreMapGABA1=BulDPreNeuMap<SynFast2>(SynFst1,SynAMPA1End,SynGABA1End,PreGABA1CntSize,TolNeu,PreMapGABA1);
      }

      if((SynACH1End-SynGABA1End)!=0)
      {
        PreMapACH1=new unsigned int[SynACH1End-SynGABA1End];
        NeuPreMapACH1=BulDPreNeuMap<SynFast2>(SynFst1,SynGABA1End,SynACH1End,PreACH1CntSize,TolNeu,PreMapACH1);
      }

      if((SynGCL1End-SynACH1End)!=0)
      {
        PreMapGCL1=new unsigned int[SynGCL1End-SynACH1End];
        NeuPreMapGCL1=BulDPreNeuMap<SynFast2>(SynFst1,SynACH1End,SynGCL1End,PreGCL1CntSize,TolNeu,PreMapGCL1);
      }

      if(SynNMDA1Size!=0)
      {
        PreMapNMDA1=new unsigned int[SynNMDA1Size];
        NeuPreMapNMDA1=BulDPreNeuMap<SynSlow2>(SynSlo1,0,SynNMDA1Size,PreNMDA1CntSize,TolNeu,PreMapNMDA1);
      }

      logdata<<"1:Time ";
      for(unsigned int j=0,kk=2,m=0;j<OFCtr[i].ParList.size();j++)
      {
        m=OFCtr[i].ParList[j];
        for(unsigned int k=NIdx[m];k<(NIdx[m]+ConfSta[m].NeuNum);k++)
        {

          for(unsigned int NID=0;NID<PreAMPA1CntSize;NID++)
          {
            if( SynFst1[PreMapAMPA1[NeuPreMapAMPA1[NID].Beg]].PreID==k)
            {
              for(unsigned int n=NeuPreMapAMPA1[NID].Beg;n<=NeuPreMapAMPA1[NID].End;n++)
              { logdata<<kk++<<":NID="<<k<<"_PosID="<<SynFst1[PreMapAMPA1[n]].PosID<<"_AMAP_st "; }
            }
          }

          for(unsigned int NID=0;NID<PreGABA1CntSize;NID++)
          {
            if( SynFst1[PreMapGABA1[NeuPreMapGABA1[NID].Beg]].PreID==k)
            {
              for(unsigned int n=NeuPreMapGABA1[NID].Beg;n<=NeuPreMapGABA1[NID].End;n++)
              { logdata<<kk++<<":NID="<<k<<"_PosID="<<SynFst1[PreMapGABA1[n]].PosID<<"_GABA_st "; }
            }
          }

          for(unsigned int NID=0;NID<PreACH1CntSize;NID++)
          {
            if( SynFst1[PreMapACH1[NeuPreMapACH1[NID].Beg]].PreID==k)
            {
              for(unsigned int n=NeuPreMapACH1[NID].Beg;n<=NeuPreMapACH1[NID].End;n++)
              { logdata<<kk++<<":NID="<<k<<"_PosID="<<SynFst1[PreMapACH1[n]].PosID<<"_ACH_st "; }
            }
          }

          for(unsigned int NID=0;NID<PreGCL1CntSize;NID++)
          {
            if( SynFst1[PreMapGCL1[NeuPreMapGCL1[NID].Beg]].PreID==k)
            {
              for(unsigned int n=NeuPreMapGCL1[NID].Beg;n<=NeuPreMapGCL1[NID].End;n++)
              { logdata<<kk++<<":NID="<<k<<"_PosID="<<SynFst1[PreMapGCL1[n]].PosID<<"_GCL_st "; }
            }
          }

          for(unsigned int NID=0;NID<PreNMDA1CntSize;NID++)
          {
            if( SynSlo1[PreMapNMDA1[NeuPreMapNMDA1[NID].Beg]].PreID==k)
            {
              for(unsigned int n=NeuPreMapNMDA1[NID].Beg;n<=NeuPreMapNMDA1[NID].End;n++)
              {
                logdata<<kk++<<":NID="<<k<<"_PosID="<<SynSlo1[PreMapNMDA1[n]].PosID<<"_NMDA_xt ";
                logdata<<kk++<<":NID="<<k<<"_PosID="<<SynSlo1[PreMapNMDA1[n]].PosID<<"_NMDA_st ";
              }
            }
          }

        }
      }

      break;

      case SynWgtID:
      if(SynAMPA1End!=0)
      {
        PreMapAMPA1=new unsigned int[SynAMPA1End];
        NeuPreMapAMPA1=BulDPreNeuMap<SynFast2>(SynFst1,0,SynAMPA1End,PreAMPA1CntSize,TolNeu,PreMapAMPA1);
      }

      if((SynGABA1End-SynAMPA1End)!=0)
      {
        PreMapGABA1=new unsigned int[SynGABA1End-SynAMPA1End];
        NeuPreMapGABA1=BulDPreNeuMap<SynFast2>(SynFst1,SynAMPA1End,SynGABA1End,PreGABA1CntSize,TolNeu,PreMapGABA1);
      }

      if((SynACH1End-SynGABA1End)!=0)
      {
        PreMapACH1=new unsigned int[SynACH1End-SynGABA1End];
        NeuPreMapACH1=BulDPreNeuMap<SynFast2>(SynFst1,SynGABA1End,SynACH1End,PreACH1CntSize,TolNeu,PreMapACH1);
      }

      if((SynGCL1End-SynACH1End)!=0)
      {
        PreMapGCL1=new unsigned int[SynGCL1End-SynACH1End];
        NeuPreMapGCL1=BulDPreNeuMap<SynFast2>(SynFst1,SynACH1End,SynGCL1End,PreGCL1CntSize,TolNeu,PreMapGCL1);
      }

      if(SynNMDA1Size!=0)
      {
        PreMapNMDA1=new unsigned int[SynNMDA1Size];
        NeuPreMapNMDA1=BulDPreNeuMap<SynSlow2>(SynSlo1,0,SynNMDA1Size,PreNMDA1CntSize,TolNeu,PreMapNMDA1);
      }

      logdata<<"1:Time ";
      for(unsigned int j=0,kk=2,m=0;j<OFCtr[i].ParList.size();j++)
      {
        m=OFCtr[i].ParList[j];
        for(unsigned int k=NIdx[m];k<(NIdx[m]+ConfSta[m].NeuNum);k++)
        {

          for(unsigned int NID=0;NID<PreAMPA1CntSize;NID++)
          {
            if( SynFst1[PreMapAMPA1[NeuPreMapAMPA1[NID].Beg]].PreID==k)
            {
              for(unsigned int n=NeuPreMapAMPA1[NID].Beg;n<=NeuPreMapAMPA1[NID].End;n++)
              { logdata<<kk++<<":PreNID="<<k<<"_PosID="<<SynFst1[PreMapAMPA1[n]].PosID<<"_AMAP "; }
            }
          }

          for(unsigned int NID=0;NID<PreGABA1CntSize;NID++)
          {
            if( SynFst1[PreMapGABA1[NeuPreMapGABA1[NID].Beg]].PreID==k)
            {
              for(unsigned int n=NeuPreMapGABA1[NID].Beg;n<=NeuPreMapGABA1[NID].End;n++)
              { logdata<<kk++<<":PreNID="<<k<<"_PosID="<<SynFst1[PreMapGABA1[n]].PosID<<"_GABA "; }
            }
          }

          for(unsigned int NID=0;NID<PreACH1CntSize;NID++)
          {
            if( SynFst1[PreMapACH1[NeuPreMapACH1[NID].Beg]].PreID==k)
            {
              for(unsigned int n=NeuPreMapACH1[NID].Beg;n<=NeuPreMapACH1[NID].End;n++)
              { logdata<<kk++<<":PreNID="<<k<<"_PosID="<<SynFst1[PreMapACH1[n]].PosID<<"_ACH "; }
            }
          }

          for(unsigned int NID=0;NID<PreGCL1CntSize;NID++)
          {
            if( SynFst1[PreMapGCL1[NeuPreMapGCL1[NID].Beg]].PreID==k)
            {
              for(unsigned int n=NeuPreMapGCL1[NID].Beg;n<=NeuPreMapGCL1[NID].End;n++)
              { logdata<<kk++<<":PreNID="<<k<<"_PosID="<<SynFst1[PreMapGCL1[n]].PosID<<"_GCL "; }
            }
          }

          for(unsigned int NID=0;NID<PreNMDA1CntSize;NID++)
          {
            if( SynSlo1[PreMapNMDA1[NeuPreMapNMDA1[NID].Beg]].PreID==k)
            {
              for(unsigned int n=NeuPreMapNMDA1[NID].Beg;n<=NeuPreMapNMDA1[NID].End;n++)
              { logdata<<kk++<<":PreNID="<<k<<"_PosID="<<SynSlo1[PreMapNMDA1[n]].PosID<<"_NMDA "; }
            }
          }

        }
      }
      break;

      default:;
    }

    logdata<<endl;
  }

  {string ().swap(s1);}
  logdata.close();

}

void NetSimGnl2_1::FileOut()
{
  string s1;

  for(vector<OutCtr>::size_type i=0;i<OFCtr.size();i++)
  {
    unsigned int kk=0;
    unsigned int n=0;
    switch(OFCtr[i].TypeID)
    {
      case IsynSepSumID:
      ostr_tmp<<EvTime*1e-3f;

      for(unsigned int j=0,m=0;j<OFCtr[i].ParList.size();j++)
      {
        m=OFCtr[i].ParList[j];
        for(unsigned int k=NIdx[m];k<(NIdx[m]+ConfSta[m].NeuNum);k++)
        {
          for(unsigned int m=0;m<5;m++)
          { ostr_tmp<<" "<<Nets[k].Irep[m]*1e-12f; }
        }
      }

      OFiles[i]<<ostr_tmp.str()<<endl;
      ostr_tmp.str("");
      break;

      case MemPotID:
      ostr_tmp<<EvTime*1e-3f;

      for(unsigned int j=0,m=0;j<OFCtr[i].ParList.size();j++)
      {
        m=OFCtr[i].ParList[j];
        for(unsigned int k=NIdx[m];k<(NIdx[m]+ConfSta[m].NeuNum);k++)
        { ostr_tmp<<" "<<NetVBse[k].V*1e-3f; }
      }

      OFiles[i]<<ostr_tmp.str()<<endl;
      ostr_tmp.str("");
      break;

      case dVOvfID:
      for(unsigned int j=0,m=0;j<OFCtr[i].ParList.size();j++)
      {
        m=OFCtr[i].ParList[j];
        for(unsigned int k=NIdx[m];k<(NIdx[m]+ConfSta[m].NeuNum);k++)
        {
          if(NetVLIF[k].NState & MEMDV_OVF)
          {
            if((OFCtr[i].Status & ALIGN_OF)!=0 )
              ostr_tmp<<EvTime*1e-3f<<" "<<n+1<<"\n";
            else
              ostr_tmp<<EvTime*1e-3f<<" "<<k+1<<"\n";

            kk++;
          }
          else if(NetVLIF[k].NState & MEMDV_UDF)
          {

           if((OFCtr[i].Status & ALIGN_OF)!=0 )
             ostr_tmp<<EvTime*1e-3f<<" -"<<n+1<<"\n";
           else
             ostr_tmp<<EvTime*1e-3f<<" -"<<k+1<<"\n";

           kk++;
          }
          n++;
        }
      }
      if(kk!=0)
      {
        OFiles[i]<<ostr_tmp.str()<<flush;
        ostr_tmp.str("");
      }
      break;

      case SpikeID:
      for(unsigned int j=0,m=0;j<OFCtr[i].ParList.size();j++)
      {
        m=OFCtr[i].ParList[j];
        for(unsigned int k=NIdx[m];k<(NIdx[m]+ConfSta[m].NeuNum);k++)
        {
          if((NetVLIF[k].NState & SPIKE_F)!=0)
          {
            if((OFCtr[i].Status & ALIGN_OF)!=0 )
              ostr_tmp<<EvTime*1e-3f<<" "<<n<<"\n";
            else
              ostr_tmp<<EvTime*1e-3f<<" "<<k<<"\n";

            kk++;
          }
          n++;
        }
      }

      if(kk != 0)
      {
        OFiles[i]<<ostr_tmp.str()<<flush;
        ostr_tmp.str("");
      }
      break;

      case FRateID:
      for(unsigned int j=0,m=0;j<OFCtr[i].ParList.size();j++)
      {
        m=OFCtr[i].ParList[j];
        for(unsigned int k=NIdx[m];k<(NIdx[m]+ConfSta[m].NeuNum);k++)
        {
          if(NetVLIF[k].NState & SPIKE_F)
          { OFCtr[i].SpiRWin[OFCtr[i].ParList.size()*OFCtr[i].RWin+j]++; }
        }
      }

      if(Step % OneMsStep==0) //  check per 1 mS
      {
        unsigned long long Step1ms=(unsigned long long)(Step/OneMsStep);
        unsigned int m=Step1ms % OFCtr[i].RWin;

        if(Step1ms%OFCtr[i].FPnt==0) // default is 100 mS , print out to file
        {
          ostr_tmp<<EvTime*1e-3f;

          if((OFCtr[i].Status & SUMPARSALL_OF) != 0)
          {
            unsigned int mm=0,nn=0;
            for(unsigned int j=0;j<OFCtr[i].ParList.size();j++)
            {
              mm+=OFCtr[i].SpiRWin[OFCtr[i].ParList.size()*(OFCtr[i].RWin+1)+j];
              nn+=ConfSta[j].NeuNum;
            }
            ostr_tmp<<" "<<1e3f*mm/OFCtr[i].RWin/nn;

          }
          else
          {
            for(unsigned int j=0;j<OFCtr[i].ParList.size();j++)
            { ostr_tmp<<" "<<(1e3f * OFCtr[i].SpiRWin[OFCtr[i].ParList.size()*(OFCtr[i].RWin+1)+j]/OFCtr[i].RWin/ConfSta[j].NeuNum ); }
          }

          OFiles[i]<<ostr_tmp.str()<<endl;
          ostr_tmp.str("");
        }

        for(unsigned int j=0;j<OFCtr[i].ParList.size();j++) //update average rate in RWin queue
        {
          OFCtr[i].SpiRWin[OFCtr[i].ParList.size()*(OFCtr[i].RWin+1)+j] -= OFCtr[i].SpiRWin[j*(OFCtr[i].RWin)+m]; //clear old data
          OFCtr[i].SpiRWin[OFCtr[i].ParList.size()*(OFCtr[i].RWin+1)+j] += OFCtr[i].SpiRWin[OFCtr[i].ParList.size()*OFCtr[i].RWin+j]; //accumulate new data
          OFCtr[i].SpiRWin[j*(OFCtr[i].RWin)+m]=OFCtr[i].SpiRWin[OFCtr[i].ParList.size()*OFCtr[i].RWin+j];//copy new to overwrite old data
          OFCtr[i].SpiRWin[OFCtr[i].ParList.size()*OFCtr[i].RWin+j]=0; //clear buffer
        }
      }
      break;

      case GatVarID:
      ostr_tmp<<EvTime*1e-3f<<" ";
      for(unsigned int j=0,m=0;j<OFCtr[i].ParList.size();j++)
      {
        m=OFCtr[i].ParList[j];
        for(unsigned int k=NIdx[m];k<(NIdx[m]+ConfSta[m].NeuNum);k++)
        {

          for(unsigned int NID=0;NID<PreAMPA1CntSize;NID++)
          {
            if( SynFst1[PreMapAMPA1[NeuPreMapAMPA1[NID].Beg]].PreID==k)
            {
              for(unsigned int n=NeuPreMapAMPA1[NID].Beg;n<=NeuPreMapAMPA1[NID].End;n++)
              { ostr_tmp<<SynFst1[PreMapAMPA1[n]].st<<" "; }
            }
          }

          for(unsigned int NID=0;NID<PreGABA1CntSize;NID++)
          {
            if( SynFst1[PreMapGABA1[NeuPreMapGABA1[NID].Beg]].PreID==k)
            {
              for(unsigned int n=NeuPreMapGABA1[NID].Beg;n<=NeuPreMapGABA1[NID].End;n++)
              { ostr_tmp<<SynFst1[PreMapGABA1[n]].st<<" "; }
            }
          }

          for(unsigned int NID=0;NID<PreACH1CntSize;NID++)
          {
            if( SynFst1[PreMapACH1[NeuPreMapACH1[NID].Beg]].PreID==k)
            {
              for(unsigned int n=NeuPreMapACH1[NID].Beg;n<=NeuPreMapACH1[NID].End;n++)
              { ostr_tmp<<SynFst1[PreMapACH1[n]].st<<" "; }
            }
          }

          for(unsigned int NID=0;NID<PreGCL1CntSize;NID++)
          {
            if( SynFst1[PreMapGCL1[NeuPreMapGCL1[NID].Beg]].PreID==k)
            {
              for(unsigned int n=NeuPreMapGCL1[NID].Beg;n<=NeuPreMapGCL1[NID].End;n++)
              { ostr_tmp<<SynFst1[PreMapGCL1[n]].st<<" "; }
            }
          }

          for(unsigned int NID=0;NID<PreNMDA1CntSize;NID++)
          {
            if( SynSlo1[PreMapNMDA1[NeuPreMapNMDA1[NID].Beg]].PreID==k)
            {
              for(unsigned int n=NeuPreMapNMDA1[NID].Beg;n<=NeuPreMapNMDA1[NID].End;n++)
              { ostr_tmp<<SynSlo1[PreMapNMDA1[n]].xt<<" "<<SynSlo1[PreMapNMDA1[n]].st<<" "; }
            }
          }

        }
      }
      OFiles[i]<<ostr_tmp.str()<<endl;
      ostr_tmp.str("");
      break;

      case SynWgtID:
      if(Step % (OFCtr[i].FPnt*OneMsStep)==0) // default is 100 mS , print out to file
      {

      ostr_tmp<<EvTime*1e-3f<<" ";
      for(unsigned int j=0,m=0;j<OFCtr[i].ParList.size();j++)
      {
        m=OFCtr[i].ParList[j];
        for(unsigned int k=NIdx[m];k<(NIdx[m]+ConfSta[m].NeuNum);k++)
        {
          for(unsigned int NID=0;NID<PreAMPA1CntSize;NID++)
          {
            if( SynFst1[PreMapAMPA1[NeuPreMapAMPA1[NID].Beg]].PreID==k)
            {
              for(unsigned int n=NeuPreMapAMPA1[NID].Beg;n<=NeuPreMapAMPA1[NID].End;n++)
              { ostr_tmp<<SynFst1[PreMapAMPA1[n]].g<<" "; }
            }
          }

          for(unsigned int NID=0;NID<PreGABA1CntSize;NID++)
          {
            if( SynFst1[PreMapGABA1[NeuPreMapGABA1[NID].Beg]].PreID==k)
            {
              for(unsigned int n=NeuPreMapGABA1[NID].Beg;n<=NeuPreMapGABA1[NID].End;n++)
              { ostr_tmp<<SynFst1[PreMapGABA1[n]].g<<" "; }
            }
          }

          for(unsigned int NID=0;NID<PreACH1CntSize;NID++)
          {
            if( SynFst1[PreMapACH1[NeuPreMapACH1[NID].Beg]].PreID==k)
            {
              for(unsigned int n=NeuPreMapACH1[NID].Beg;n<=NeuPreMapACH1[NID].End;n++)
              { ostr_tmp<<SynFst1[PreMapACH1[n]].g<<" "; }
            }
          }

          for(unsigned int NID=0;NID<PreGCL1CntSize;NID++)
          {
            if( SynFst1[PreMapGCL1[NeuPreMapGCL1[NID].Beg]].PreID==k)
            {
              for(unsigned int n=NeuPreMapGCL1[NID].Beg;n<=NeuPreMapGCL1[NID].End;n++)
              { ostr_tmp<<SynFst1[PreMapGCL1[n]].g<<" "; }
            }
          }

          for(unsigned int NID=0;NID<PreNMDA1CntSize;NID++)
          {
            if( SynSlo1[PreMapNMDA1[NeuPreMapNMDA1[NID].Beg]].PreID==k)
            {
              for(unsigned int n=NeuPreMapNMDA1[NID].Beg;n<=NeuPreMapNMDA1[NID].End;n++)
              { ostr_tmp<<SynSlo1[PreMapNMDA1[n]].xt<<" "<<SynSlo1[PreMapNMDA1[n]].g<<" "; }
            }
          }

        }
      }
      OFiles[i]<<ostr_tmp.str()<<endl;
      ostr_tmp.str("");

      }
      break;

      default:;
    }

  }

  {string ().swap(s1);}
}


