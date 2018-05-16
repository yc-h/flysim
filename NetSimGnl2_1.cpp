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
#include "NetSim.h"



using namespace std;

//global variables
//------------------NetSimGnl
NetSimGnl2_1::NetSimGnl2_1()
{
  TolNeu=0;
  TolSynp=0;
  TolRep=0;
  EvTime=0.0f;
  Step=0;
  Repeat=1;
  RpSize=1;
  dt=DTIME;

  ftokens=0;
  Acnt=0;
  fcnt=0;
  Scnt=0;
}


NetSimGnl2_1::~NetSimGnl2_1()
{
  delete [] Nets;
  {vector<TmEvt> ().swap(NEvts);};
  {vector<OutCtr> ().swap(NOFCtr);};
}


void NetSimGnl2_1::DoEvt(void)
{
  for(unsigned int i=iEvt;i<Evts.size();i++)
  {
    if((Evts[i].TimEv*1e3f>=EvTime) && (Evts[i].TimEv*1e3f<(EvTime+dt)))
    {
      iEvt=i;

      if(Evts[i].Type==EVTEXTFREQ)
      {
        for(unsigned int m=NIdx[Evts[i].Pp];m<(NIdx[Evts[i].Pp]+ConfSta[Evts[i].Pp].NeuNum);m++)
        {
          switch(Evts[i].SubTyp)
          {
            case KW_AMPA:
              if(Evts[i].Var1==0.0f)
                Nets[m].ActCtr &= ~AMPA_ACT;
              else
              {
                Nets[m].ActCtr |= AMPA_ACT;
                BgSt[m].pb[0]=(float)(Evts[i].Var1*dt*1e-3f);
              }
            break;

            case KW_GABA:
              if(Evts[i].Var1==0.0f)
                Nets[m].ActCtr &= ~GABA_ACT;
              else
              {
                Nets[m].ActCtr |= GABA_ACT;
                BgSt[m].pb[1]=(float)(Evts[i].Var1*dt*1e-3f);
              }
            break;

            case KW_ACH:
              if(Evts[i].Var1==0.0f)
                Nets[m].ActCtr &= ~ACH_ACT;
              else
              {
                Nets[m].ActCtr |= ACH_ACT;
                BgSt[m].pb[2]=(float)(Evts[i].Var1*dt*1e-3f);
              }
            break;

            case KW_GCL:
              if(Evts[i].Var1==0.0f)
                Nets[m].ActCtr &= ~GCL_ACT;
              else
              {
                Nets[m].ActCtr |= GCL_ACT;
                BgSt[m].pb[3]=(float)(Evts[i].Var1*dt*1e-3f);
              }
            break;

            case KW_NMDA:
              if(Evts[i].Var1==0.0f)
                Nets[m].ActCtr &= ~NMDA_ACT;
              else
              {
                Nets[m].ActCtr |= NMDA_ACT;
                BgSt[m].pb[4]=(float)(Evts[i].Var1*dt*1e-3f);
              }
            break;

            default:;
          }
        }
      }
      else if(Evts[i].Type==EVTCUNTINJ)
      {
        for(unsigned int m=NIdx[Evts[i].Pp];m<(NIdx[Evts[i].Pp]+ConfSta[Evts[i].Pp].NeuNum);m++)
        {
          BgSt[m].mean=Evts[i].Var1;
          //BgSt[m].std=Evts[i].Var2;
          BgSt[m].std=Evts[i].Var2*(1.0f+(float)log((1.0f-exp(-0.1f/Confs[Evts[i].Pp].Taum))/(1.0f-NetPBse[m].df)));//normalize to dt=0.1ms

          if(Evts[i].Var1==0.0f && Evts[i].Var2==0.0f)
            Nets[m].ActCtr &= ~INJ_ACT;
          else
            Nets[m].ActCtr |= INJ_ACT;
        }
      }

    }
    else if(Evts[i].TimEv*1e3f>=(EvTime+dt))
    { break; }


  }
}

void NetSimGnl2_1::RunInit(int thNum,int sop)
{
// print to screen short population information
  string s1;
  SInfo(s1);
  cout<<s1;

  cout<<"thread number:"<<thNum<<"\n";
// population log file
  cout<<"write log information to network.log , please wait!\n";

  logDataWrite("network.log");
  cout<<"start processing----------------------------------------------------\n";

  {string ().swap(s1);}
}

void NetSimGnl2_1::SimThd(unsigned int thid)
{
  for(unsigned int i=0; (i < Endtime) || NoLmt ; i++) //calculation
  {
    if(thid==0)
    {
      Step=i;
      EvTime=i*dt;
      DoEvt();
      if((i % ReportStep)==0)
      {
        t2=Clock::now();
        if(NoLmt)
        {
          cout
          <<"Expertment time="<<EvTime*1e-3f<<"s "
          <<"real time="<<(chrono::duration_cast<chrono::nanoseconds>(t2-t1).count())*1e-9<<"s "
          <<":infinity End time             \r"
          <<flush;
        }
        else
        {
          cout
          <<"Expertment time="<<EvTime*1e-3f<<"s "
          <<"real time="<<(chrono::duration_cast<chrono::nanoseconds>(t2-t1).count())*1e-9<<"s "
          <<100.0f*(float)(i+1)/Endtime<<"%             \r"
          <<flush;
        }
      }
      Scnt=0;
    }


    {  //sync
      unique_lock<mutex> lck(mtx);
      if(++fcnt==thsize)
      {
        ftokens=1;
        cv.notify_all();
      }
      cv.wait(lck,[](){return ftokens==1?true:false;});
    }

    Activation(thid,ActCtr);

    { //sync
      unique_lock<mutex> lck(mtx);
      if(--fcnt==0)
      {
        ftokens=2;
        cv.notify_all();
      }
      cv.wait(lck,[](){return ftokens==2?true:false;});
    }

    Spike(thid,ActCtr);

    { //sync
      unique_lock<mutex> lck(mtx);
      if(++Scnt==thsize)
      {
        ftokens=5;
        cv.notify_all();
      }
      cv.wait(lck,[](){return ftokens==5?true:false;});
    }

    if(thid==0) FileOut();
  }
}

void NetSimGnl2_1::SimOne(void)
{
  for(unsigned int i=0; (i < Endtime) || NoLmt ; i++)  //calculation
  {
    Step=i;
    EvTime=i*dt;
    DoEvt();
    if((i % ReportStep)==0)
    {
      t2=Clock::now();
      if(NoLmt)
      {
        cout
        <<"Expertment time="<<EvTime*1e-3f<<"s "
        <<"real time="<<(chrono::duration_cast<chrono::nanoseconds>(t2-t1).count())*1e-9<<"s "
        <<":infinity End time             \r"
        <<flush;
      }
      else
      {
        cout
        <<"Expertment time="<<EvTime*1e-3f<<"s "
        <<"real time="<<(chrono::duration_cast<chrono::nanoseconds>(t2-t1).count())*1e-9<<"s "
        <<100.0f*(float)(i+1)/Endtime<<"%             \r"
        <<flush;
      }
    }

    Activation(0,ActCtr);
    Spike(0,ActCtr);

    FileOut();
  }
}

void NetSimGnl2_1::Run()
{
  string s1;
  NoLmt=Evts[Evts.size()-1].NoLmt;
  Endtime = (unsigned long long)floor(Evts[Evts.size()-1].TimEv*1e3f/dt);
  ReportStep=100/dt;
  OneMsStep=(unsigned int)1/dt;

  cout<<setprecision(7)<<showpoint<<flush;

//initial variable in every neuron
  VarInit();

//initial every neuron background noise
  for(unsigned int i=0;i<ConfSta.size();i++)
  {
    for(unsigned int j=NIdx[i];j<(NIdx[i]+ConfSta[i].NeuNum);j++) // set up each neuron's parameter in every population
    {
      for(unsigned int k=0;k<ConfSta[i].PpRepNum;k++)  // parameter in receptors and variable in FExt receptors
      {
        if(PreSynR->NeuRepName==KW_AMPA)
        { BgSt[j].pb[0]=PreSynR->FExt*dt*1e-3f; }
        else if(PreSynR->NeuRepName==KW_GABA)
        { BgSt[j].pb[1]=PreSynR->FExt*dt*1e-3f; }
        else if(PreSynR->NeuRepName==KW_ACH)
        { BgSt[j].pb[2]=PreSynR->FExt*dt*1e-3f; }
        else if(PreSynR->NeuRepName==KW_GCL)
        { BgSt[j].pb[3]=PreSynR->FExt*dt*1e-3f; }
        else if(PreSynR->NeuRepName==KW_NMDA)
        { BgSt[j].pb[4]=PreSynR->FExt*dt*1e-3f; }
      }
    }
  }
// clear output file firing rate buffer
  if(Repeat>1)
  {
    for(vector<OutCtr>::size_type i=0;i<OFCtr.size();i++)
    {
      if(OFCtr[i].SpiRWin.size()!=0)
      {
        for(unsigned int j=0;j<OFCtr[i].SpiRWin.size();j++)
        {OFCtr[i].SpiRWin[j]=0;}
      }
    }

    iEvt=0;
  }


//build up noise
  {
    delete rdgen;
    delete [] zig;
  }

  zig=new ziggurat[thsize];
  rdgen = new mt19937;

  rdgen->seed(rseeds[Repeat-1]);
  zig[0].jcong=rseeds[Repeat-1];

  uniform_int_distribution<unsigned int> dist(0,RAND_MAX);
  if(thsize>1)
  {
    for(unsigned int i=1;i<thsize;i++)
    { zig[i].jcong=dist(*rdgen); }
  }

  cout<<"Iteration:"<<Repeat<<"\n";
  cout<<noshowpoint<<setprecision(7);

  FileOpen();
  if(thsize==1)
  {
    t1=Clock::now();
    SimOne();
  }
  else
  {
    t1=Clock::now();
    for(unsigned int i=0;i<thsize; i++)
    { th[i] = thread(&NetSimGnl2_1::SimThd,this,i); }

    for(unsigned int i=0;i<thsize;i++)
    {th[i].join();}
  }

  t2=Clock::now();
  cout
  <<"Expertment time="<<EvTime*1e-3f<<"s "
  <<"real time="<<(chrono::duration_cast<chrono::nanoseconds>(t2-t1).count())*1e-9<<"s 100%                 "<<endl;

  FileClose();

}


