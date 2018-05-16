#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <stdio.h>
#include <fstream>
#include <cstdlib>
#include <omp.h>
#include <random>
#include <chrono>
#include "NetSim.h"


using namespace std;

//------------------NetSim_gp4
NetSim_gp4::NetSim_gp4()
{
  TolNeu=0;
  TolSynp=0;
  TolRep=0;
  EvTime=0.0f;
  Step=0;
  Repeat=1;
}

NetSim_gp4::~NetSim_gp4()
{ {vector<Neuron_gp4> ().swap(Nets);} }

void NetSim_gp4::SInfo(string& ostr) //short information
{
  ostringstream ostr_tmp;
  float dt=DTIME;
  string SType;
//-------show total synapse inforation
  switch(SolType)
  {
    case ACCURATE:
    SType="accurate";
    break;

    case MODERATE:
    SType="moderate";
    break;

    case ROUGH:
    SType="rough";
    break;

    default:
    SType="moderate";
    break;
  }

  ostr_tmp
      <<"\nneuron model:LIF_GP4"
      <<", solver:"<<SType
      <<", total neurons="<<TolNeu
      <<", total synapses="<<TolSynp
      <<", total receptors="<<TolRep
      <<", time step="<<dt
      <<"\n";

  ostr=ostr_tmp.str();
}

void NetSim_gp4::TexInfo(string& ostr)
{
  ostringstream ostr_tmp;
  string str;
//-------show total synapse inforation
  SInfo(str);
  ostr_tmp<<str<<"\n";
  ostr_tmp<<"population information:-------------------------------------------\n";
  for(unsigned int i=0;i<Nets.size();i++)
  {
    Nets[i].Info(str);
    ostr_tmp<<"population:"<<i<<"\n"<<str;
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
}

void NetSim_gp4::SetPp(int thNum,int sop)
{
  unsigned int a,c;
  bool IdxNew;
  vector<unsigned int> IdxNum;
  int state=IdxIns;

  Neuron_gp4 iniNeu;
  NeuVar iniVm;
  ConTab_gp2 iniIdx;
  Nlink iniNlink;
  RepVar1_gp iniVar1;

//-----------------create default neuron network


  NIdx.insert(NIdx.end(),ConfSta.size(),0);
  Nets.insert(Nets.end(),ConfSta.size(),iniNeu); // build total neurons

  unsigned int m=0;
  for(unsigned int i=0;i<ConfSta.size();i++)
  {
    Nets[i].RepAMPA.insert(Nets[i].RepAMPA.end(),ConfSta[i].NeuNum,iniVar1);
    Nets[i].RepGABA.insert(Nets[i].RepGABA.end(),ConfSta[i].NeuNum,iniVar1);
    Nets[i].RepNMDA.insert(Nets[i].RepNMDA.end(),ConfSta[i].NeuNum,iniVar1);
    Nets[i].NeuVm.insert(Nets[i].NeuVm.end(),ConfSta[i].NeuNum,iniVm);

    NIdx[i]=m;
    m += ConfSta[i].NeuNum;
  }
  TolNeu=m;


  for(unsigned int i=0;i<ConfSta.size();i++)  //initialize parameters in each population
  {
    //leaky integrate and fire
    Nets[i].Vth    = (Confs+i)->Vth*1e-3f;
    Nets[i].El     = (Confs+i)->El*1e-3f;
    Nets[i].gl     = (Confs+i)->Cm/(Confs+i)->Taum *1e-6f;
    Nets[i].Cm     = (Confs+i)->Cm*1e-9f;
    Nets[i].Vreset = (Confs+i)->Vreset;
    Nets[i].RfPed  = (Confs+i)->RfPed;
    Nets[i].SpiDy  = (Confs+i)->SpiDy;

    for(unsigned int k=0;k<ConfSta[i].PpRepNum;k++)  // parameters of receptors in the population
    {
      if(PreSynR->NeuRepName==KW_AMPA)
      {
        Nets[i].ParAMPA.g=PreSynR->MnExtEff*1e-9f*PreSynR->MnExtCon;  // FExt Efficacy
        Nets[i].ParAMPA.tst  = PreSynR->Tau*1e-3f;
        Nets[i].ParAMPA.Vrev = PreSynR->Vrev*1e-3f;
        Nets[i].p1=PreSynR->FExt*DTIME;
      }
      else if(PreSynR->NeuRepName==KW_GABA)
      {
        Nets[i].ParGABA.g=PreSynR->MnExtEff*1e-9f*PreSynR->MnExtCon;  // FExt Efficacy
        Nets[i].ParGABA.tst  = PreSynR->Tau*1e-3f;
        Nets[i].ParGABA.Vrev = PreSynR->Vrev*1e-3f;
        Nets[i].p2=PreSynR->FExt*DTIME;
      }
      else if(PreSynR->NeuRepName==KW_NMDA)
      {
        Nets[i].ParNMDA.g=PreSynR->MnExtEff*1e-9f*PreSynR->MnExtCon;  // FExt Efficacy
        Nets[i].ParNMDA.tst  = PreSynR->Tau*1e-3f;
        Nets[i].ParNMDA.Vrev = PreSynR->Vrev*1e-3f;
        Nets[i].p3=PreSynR->FExt*DTIME;
      }
    }

    Nets[i].rnd_seed=seed_gen();
    Nets[i].rdgen.seed(Nets[i].rnd_seed);
  }//initialize end

//-----------------------------------------------------------------------------
//-----------------Confs transcripts to Nets

  for(unsigned int i=0;i<ConfSta.size();i++) // figure out pre population
  {
    cout<<"setup network:"<<100.0f*(float)(i+1)/ConfSta.size()<<"%        \r"<<flush;
    a=0;
    if(ConfSta[i].PpLinkNum!=0)
    {
      IdxNum.insert(IdxNum.end(),ConfSta[i].PpLinkNum,0);
      state=IdxIns;
      a=0;
      for(unsigned int j=0;j<ConfSta[i].PpLinkNum;j++) // pre synaptic population synapse
      {
        IdxNew=true;
        for(c=0;c<j;c++)
        {
          if( PoSynPp == ((Confs+i)->PpLink+c)->TarPpName)
          {
            IdxNew=false;
            break;
          }

        }

        switch(state)
        {
          case IdxIns:
          if(!IdxNew)
          {
            a=j-1;
            state=IdxFiPre;
          }
          break;

          case IdxFiPre:
          if(!IdxNew)
          { state=IdxFiPre; }
          else
          { state=IdxRe;}
          break;

          default:;
        }

        if(j==0)
        {
          IdxNum[j]=1;
          Nets[i].Index.push_back(iniIdx);
          Nets[i].Index[Nets[i].Index.size()-1].Axn.insert(Nets[i].Index[Nets[i].Index.size()-1].Axn.begin(),ConfSta[PoSynPp].NeuNum,iniNlink);
        }
        else if(state==IdxIns)
        {
          IdxNum[j]=IdxNum[j-1]+1;
          Nets[i].Index.push_back(iniIdx);
          Nets[i].Index[Nets[i].Index.size()-1].Axn.insert(Nets[i].Index[Nets[i].Index.size()-1].Axn.begin(),ConfSta[PoSynPp].NeuNum,iniNlink);
        }
        else if(state==IdxFiPre)
        { IdxNum[j]=IdxNum[c]; }
        else if(state==IdxRe)
        {
          IdxNum[j]=IdxNum[a]+1;

          Nets[i].Index.push_back(iniIdx);
          Nets[i].Index[Nets[i].Index.size()-1].Axn.insert(Nets[i].Index[Nets[i].Index.size()-1].Axn.begin(),ConfSta[PoSynPp].NeuNum,iniNlink);

          state=IdxIns;
        }

      }

      a=0;
      for(unsigned int k=0;k<IdxNum.size();k++)
      { if(a<IdxNum[k]) a=IdxNum[k]; }

    }

    TolSynp += a;
    TolRep += ConfSta[i].NeuNum*3;

    for(unsigned int j=0;j<ConfSta[i].PpLinkNum;j++)// pre synaptic population synapse
    {

      if(i == PoSynPp) // self-stimulation
      { //setup MeanG for subtraction

        Nets[i].SelfSti=true;

        if(PoSyn->TarRep == KW_AMPA)
        { Nets[i].ParAMPA.MeanG = PoSyn->MeanEff*1e-9f*PoSyn->weight; }
        else if(PoSyn->TarRep == KW_GABA)
        { Nets[i].ParGABA.MeanG = PoSyn->MeanEff*1e-9f*PoSyn->weight; }
        else if(PoSyn->TarRep == KW_NMDA)
        { Nets[i].ParNMDA.MeanG = PoSyn->MeanEff*1e-9f*PoSyn->weight; }
      }

      Nets[i].Index[IdxNum[j]-1].NeuID=PoSynPp;
      for(unsigned int k=0;k<ConfSta[PoSynPp].PpRepNum;k++)  // parameters of receptors in the population
      {
        if((ConfSta[PoSynPp].PpRep+k)->NeuRepName==KW_AMPA)
        { Nets[i].Index[IdxNum[j]-1].AMPA_tst=(ConfSta[PoSynPp].PpRep+k)->Tau*1e-3f; }
        else if((ConfSta[PoSynPp].PpRep+k)->NeuRepName==KW_GABA)
        { Nets[i].Index[IdxNum[j]-1].GABA_tst=(ConfSta[PoSynPp].PpRep+k)->Tau*1e-3f; }
        else if((ConfSta[PoSynPp].PpRep+k)->NeuRepName==KW_NMDA)
        { Nets[i].Index[IdxNum[j]-1].NMDA_tst=(ConfSta[PoSynPp].PpRep+k)->Tau*1e-3f; }
      }

      if(PoSyn->TarRep == KW_AMPA)  // initial and insert post neuron's receptors "connection" properties
      {
        Nets[i].Index[IdxNum[j]-1].AMPA_Act=true;
        Nets[i].Index[IdxNum[j]-1].AMPA_MeanG=PoSyn->MeanEff*1e-9f * PoSyn->weight;
      }
      else if(PoSyn->TarRep == KW_GABA)
      {
        Nets[i].Index[IdxNum[j]-1].GABA_Act=true;
        Nets[i].Index[IdxNum[j]-1].GABA_MeanG=PoSyn->MeanEff*1e-9f * PoSyn->weight;
      }
      else if(PoSyn->TarRep == KW_NMDA)
      {
        Nets[i].Index[IdxNum[j]-1].NMDA_Act=true;
        Nets[i].Index[IdxNum[j]-1].NMDA_MeanG=PoSyn->MeanEff*1e-9f * PoSyn->weight;
      }
      else if(PoSyn->TarRep == KW_GAP)
      {
        Nets[i].Index[IdxNum[j]-1].GAP_Act=true;
        Nets[i].Index[IdxNum[j]-1].GAP_MeanG=PoSyn->MeanEff*1e-9f * PoSyn->weight;
      }

    } // end of pre synaptic population synapse

    {vector<unsigned int> ().swap(IdxNum);}
    delete [] (Confs+i)->PpLink;

  } // end of pre synatic population

  delete [] Confs;
  cout<<"\n";

  SetOut();

}  //end of SetPp()


void NetSim_gp4::DoEvt()
{
  for(unsigned int i=iEvt;i<Evts.size();i++)
  {
    if((Evts[i].TimEv>=EvTime) && (Evts[i].TimEv<(EvTime+DTIME)))
    {
      iEvt=i;

      if(Evts[i].Type==EVTEXTFREQ)
      {
        if(Evts[i].SubTyp==KW_AMPA)
        { Nets[Evts[i].Pp].p1=(float)(Evts[i].Var1*DTIME);}
        else if(Evts[i].SubTyp==KW_GABA)
        { Nets[Evts[i].Pp].p2=(float)(Evts[i].Var1*DTIME);}
        else if(Evts[i].SubTyp==KW_NMDA)
        { Nets[Evts[i].Pp].p3=(float)(Evts[i].Var1*DTIME);}
        else if(Evts[i].SubTyp==KW_SINE)
        { Nets[Evts[i].Pp].RSin.SetPar(FExtSinAmp,Evts[i].Var1,0.0f); }
      }
      else if(Evts[i].Type==EVTCUNTINJ)
      {
        normal_distribution<float>::param_type par(Evts[i].Var1*1e-9f,Evts[i].Var2*1e-9f);
        Nets[Evts[i].Pp].MemNoise.param(par);
      }

    } //  at event time end
  }
}


void NetSim_gp4::RunInit(int thid,int sop)
{
  SolType=sop;
  iniEqsThd(Nets,(unsigned int)thid);

// print to screen short population information
  string s1;
  SInfo(s1);
  cout<<s1;

  cout<<"thread nember:"<<thid<<"\n";
// population log file
  cout<<"write log information to network.log , please wait!\n";

  ofstream logdata;
  logdata.open("network.log");

//-------show total synapse inforation
  SInfo(s1);
  logdata<<s1<<"\n";
  logdata<<"population information:-------------------------------------------\n";
  for(unsigned int i=0;i<Nets.size();i++)
  {
    Nets[i].Info(s1);
    logdata<<"population:"<<i<<"\n"<<s1<<flush;
  }

  logdata<<"time event information:-------------------------------------------\n";
  for(unsigned int i=0;i<Evts.size();i++)
  {
    Evts[i].Info(s1);
    logdata<<s1<<flush;
  }

  logdata<<"output files information:-----------------------------------------\n";
  for(unsigned int i=0;i<OFCtr.size();i++)
  {
    OFCtr[i].Info(s1);
    logdata<<s1<<flush;
  }

  logdata.close();
  cout<<"start processing----------------------------------------------------\n";
}

void NetSim_gp4::Run()
{
  clock_t t1,t2;
  int Endtime = (int)floor(Evts[Evts.size()-1].TimEv/DTIME);

//iniital every neuron
  for(unsigned int i=0;i<Nets.size();i++)
  { Nets[i].VarInit();  }

//set every neuron to default external stimulation
  for(unsigned int i=0;i<ConfSta.size();i++)
  {
      for(unsigned int k=0;k<ConfSta[i].PpRepNum;k++)  // parameter in receptors and variable in FExt receptors
      {
        if(PreSynR->NeuRepName==KW_AMPA)
        { Nets[i].p1=PreSynR->FExt*DTIME; }
        else if(PreSynR->NeuRepName==KW_GABA)
        { Nets[i].p2=PreSynR->FExt*DTIME; }
        else if(PreSynR->NeuRepName==KW_NMDA)
        { Nets[i].p3=PreSynR->FExt*DTIME; }
      }
#ifdef CPP11
      normal_distribution<float>::param_type par(MNoiseMean,MNoiseSTD);
      Nets[i].MemNoise.param(par);
#else
      Nets[i].mean=MNoiseMean;
      Nets[i].stddev=MNoiseSTD;
#endif

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

  cout<<"Iteration:"<<Repeat<<"\n";

  FileOpen();

//write gating variable print out information
  for(vector<OutCtr>::size_type i=0;i<OFCtr.size();i++)
  {
    if(GatVarID == OFCtr[i].TypeID)
    {
      *(OFiles+i)<<"output format:\n"
               <<"time:1";

      for(unsigned int j=0,jj=2;j<OFCtr[i].ParList.size();j++,jj++)
      {
        for(unsigned int k=0;k<Nets[OFCtr[i].ParList[j]].NeuVm.size();k++)
        {

          //*(OFiles+i)<<" Population="<<OFCtr[i].PName[j];
          *(OFiles+i)<<" Population="<<OFCtr[i].ParList[j];
          *(OFiles+i)<<" NeuronID="<<k;
          *(OFiles+i)<<" External_AMPA_st:"<<jj++;
          *(OFiles+i)<<" External_GABA_st:"<<jj++;
          *(OFiles+i)<<" External_NMDA_st:"<<jj++;

          for(unsigned int m=0;m<Nets.size();m++)
          {
            for(unsigned int n=0;n<Nets[m].Index.size();n++)
            {
              if(Nets[m].Index[n].NeuID==OFCtr[i].ParList[j])
              {
                for(unsigned int p=0;p<Nets[m].Index[n].Axn.size();p++)
                {
                  if(Nets[m].Index[n].AMPA_Act)
                    *(OFiles+i)<<" AMAP_st:"<<jj++;

                  if(Nets[m].Index[n].GABA_Act)
                    *(OFiles+i)<<" GABA_st:"<<jj++;

                  if(Nets[m].Index[n].NMDA_Act)
                    *(OFiles+i)<<" NMDA_st:"<<jj++;
                }
              }
            }
          }

          *(OFiles+i)<<",";
        }
      }

      *(OFiles+i)<<"\n";
    }

  }
//--------------------------------------------------------------------
  cout<<setprecision(4)<<showpoint;
  t1=clock();

  for(Step=0; Step < Endtime; Step++)  //calculation
  {
    EvTime=Step*DTIME;
    EqsTime=Step*DTIME;
    DoEvt();


    if((Step%ReportClock)==0)
    {
      t2=clock();
      cout
        <<"Expertment time="<<setw(4)<<EvTime<<"s "
        <<"real time="<<setw(4)<<(float)(t2-t1)/CLOCKS_PER_SEC<<"s "
        <<setw(4)<<100.0f*(float)(Step+1)/Endtime<<"%             \r"
        <<flush;
    }

    RunThd();
    //Equations[0].SpikeIpv();

    FileOut();
  }

  t2=clock();
  cout
    <<"Expertment time="<<setw(4)<<EvTime<<"s "
    <<"real time="<<setw(4)<<(float)(t2-t1)/CLOCKS_PER_SEC<<"s "
    <<setprecision(8)
    <<"total steps:"<<Step<<"              \n"
    <<"Expertment End!\n\n";

  FileClose();

}

void NetSim_gp4::FileOut()
{
  for(vector<OutCtr>::size_type i=0;i<OFCtr.size();i++)
  {
    switch(OFCtr[i].TypeID)
    {
      case MemPotID:
      *(OFiles+i)<<setprecision(7);
      *(OFiles+i)<<EvTime;
      *(OFiles+i)<<setprecision(6);

      for(unsigned int j=0;j<OFCtr[i].ParList.size();j++)
      {
        for(unsigned int k=0;k<Nets[OFCtr[i].ParList[j]].NeuVm.size();k++)
        { *(OFiles+i)<<" "<<Nets[OFCtr[i].ParList[j]].NeuVm[k].V; }
      }

      *(OFiles+i)<<endl;
      break;

      case SpikeID:
      *(OFiles+i)<<setprecision(7)<<noshowpoint;

      for(unsigned int j=0;j<OFCtr[i].ParList.size();j++)
      {
        for(unsigned int k=0;k<Nets[OFCtr[i].ParList[j]].NeuVm.size();k++)
          if(Nets[OFCtr[i].ParList[j]].NeuVm[k].SpiOut)
          { *(OFiles+i)<<EvTime<<" "<<NIdx[OFCtr[i].ParList[j]]+k<<"\n"; }
      }

      *(OFiles+i)<<flush;
      break;

      case FRateID:
      for(unsigned int j=0;j<OFCtr[i].ParList.size();j++)
      {
        for(unsigned int k=0;k<Nets[j].NeuVm.size();k++)
        {
          if(Nets[j].NeuVm[k].SpiOut)
          { OFCtr[i].SpiRWin[OFCtr[i].ParList.size()*OFCtr[i].RWin+j]++; }
        }
      }

      if(Step%OneMilliSecond==0) //  check per 1 mS
      {
        int Step1ms=(int)(Step/OneMilliSecond);
        int m=Step1ms%OFCtr[i].RWin;

        if(Step1ms%OFCtr[i].FPnt==0) // default is 100 mS , print out to file
        {
          *(OFiles+i)<<setprecision(7);
          *(OFiles+i)<<EvTime;
          *(OFiles+i)<<setprecision(6);

          for(unsigned int j=0;j<OFCtr[i].ParList.size();j++)
          { *(OFiles+i)<<" "<<(1e3f * OFCtr[i].SpiRWin[OFCtr[i].ParList.size()*(OFCtr[i].RWin+1)+j]/OFCtr[i].RWin/ConfSta[j].NeuNum ); }

          *(OFiles+i)<<endl;
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
      *(OFiles+i)<<setprecision(7);
      *(OFiles+i)<<EvTime;
      *(OFiles+i)<<setprecision(6);

      for(unsigned int j=0;j<OFCtr[i].ParList.size();j++)
      {
        for(unsigned int k=0;k<Nets[OFCtr[i].ParList[j]].NeuVm.size();k++)
        {
          *(OFiles+i)<<" "<<Nets[OFCtr[i].ParList[j]].RepAMPA[k].st;
          *(OFiles+i)<<" "<<Nets[OFCtr[i].ParList[j]].RepGABA[k].st;
          *(OFiles+i)<<" "<<Nets[OFCtr[i].ParList[j]].RepNMDA[k].st;

          for(unsigned int m=0;m<Nets.size();m++)
          {
            for(unsigned int n=0;n<Nets[m].Index.size();n++)
            {
              if(Nets[m].Index[n].NeuID==OFCtr[i].ParList[j])
              {
                for(unsigned int p=0;p<Nets[m].Index[n].Axn.size();p++)
                {
                  if(Nets[m].Index[n].AMPA_Act)
                    *(OFiles+i)<<" "<<Nets[m].Index[n].Axn[p].AMPA_st;

                  if(Nets[m].Index[n].GABA_Act)
                    *(OFiles+i)<<" "<<Nets[m].Index[n].Axn[p].GABA_st;

                  if(Nets[m].Index[n].NMDA_Act)
                    *(OFiles+i)<<" "<<Nets[m].Index[n].Axn[p].NMDA_st;
                }
              }
            }
          }

        }
      }

      *(OFiles+i)<<"\n";
      break;

      default:;
    }

  }

}

