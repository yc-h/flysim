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
#include <unistd.h>
#include "NetSim.h"

using namespace std;

//------------------NetSim_LIF3
NetSim_LIF3::NetSim_LIF3()
{
  TolNeu=0;
  TolSynp=0;
  TolRep=0;
  EvTime=0.0f;
  Step=0;
  Repeat=1;
}

NetSim_LIF3::~NetSim_LIF3()
{ {vector<Neuron_LIF3> ().swap(Nets);} }

void NetSim_LIF3::SInfo(string& ostr) //short information
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
      <<"\nneuron model:LIF3"
      <<", solver:"<<SType
      <<", total neurons="<<TolNeu
      <<", total synapses="<<TolSynp
      <<", total receptors="<<TolRep
      <<", time step="<<dt
      <<"\n";

  ostr=ostr_tmp.str();
}

void NetSim_LIF3::TexInfo(string& ostr)
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
}

//--------------setup network population

void NetSim_LIF3::SetPp(int thNum,int sop)
{
  unsigned int a,c;
  bool IdxNew;
  vector<unsigned int> IdxNum;
  int state=IdxIns;

  Neuron_LIF3 iniNeu;
  ConTab2 iniIdx;
  RV1 iniVar1;
  RV2 iniVar2;
//-----struct data initial
  iniNeu.ParAMPA.tst=2.0e-3f;
  iniNeu.ParAMPA.Vrev=0.0f;

  iniNeu.ParGABA.tst=2.0e-3f;
  iniNeu.ParGABA.Vrev=0.0f;

  iniNeu.ParNMDA.as=0.6332f; // the alpha st
  iniNeu.ParNMDA.txt=2.0e-3f; // the decay time constant of xt
  iniNeu.ParNMDA.tst=100.0e-3f; //decay time constant for st
  iniNeu.ParNMDA.Vrev=0.0e-3f; // reversal potential
  iniNeu.ParNMDA.mg=1.0f; // the concentration of magnesium

  iniIdx.NeuID=0;
  iniIdx.AMPA_ID=0;
  iniIdx.GABA_ID=0;
  iniIdx.NMDA_ID=0;
  iniIdx.Act=0;

  iniVar1.g=0.11e-9f; // synaptic conductance
  iniVar1.st=0.0f;

  iniVar2.g=0.11e-9f; // synaptic conductance
  iniVar2.st=0.0f;
  iniVar2.xt=0.0f;
//-----------------create default neuron network

  NIdx.insert(NIdx.end(),ConfSta.size(),0);
  unsigned int q=0;
  for(unsigned int i=0;i<ConfSta.size();i++)
  {
    NIdx[i]=q;
    q += ConfSta[i].NeuNum;
  }
  TolNeu=q;

  Nets.insert(Nets.end(),TolNeu,iniNeu); // build total neurons

  for(unsigned int i=0;i<ConfSta.size();i++)
  {
    for(unsigned int j=NIdx[i];j<(NIdx[i]+ConfSta[i].NeuNum);j++) // set up each neuron's parameter in every population
    {
      //leaky integrate and fire
      Nets[j].Vth    = (Confs+i)->Vth*1e-3f;
      Nets[j].El     = (Confs+i)->El*1e-3f;
      Nets[j].gl     = (Confs+i)->Cm/(Confs+i)->Taum *1e-6f;
      Nets[j].Cm     = (Confs+i)->Cm*1e-9f;
      Nets[j].Vreset = (Confs+i)->Vreset;
      Nets[j].RfPed  = (Confs+i)->RfPed;
      Nets[j].SpiDy  = (Confs+i)->SpiDy;

      for(unsigned int k=0;k<ConfSta[i].PpRepNum;k++)  // parameter in receptors and variable in FExt receptors
      {

        if(PreSynR->NeuRepName==KW_AMPA)
        {
          Nets[j].ParAMPA.tst  = PreSynR->Tau*1e-3f;
          Nets[j].ParAMPA.Vrev = PreSynR->Vrev*1e-3f;
          Nets[j].p1=PreSynR->FExt*DTIME;
        }
        else if(PreSynR->NeuRepName==KW_GABA)
        {
          Nets[j].ParGABA.tst  = PreSynR->Tau*1e-3f;
          Nets[j].ParGABA.Vrev = PreSynR->Vrev*1e-3f;
          Nets[j].p2=PreSynR->FExt*DTIME;
        }
        else if(PreSynR->NeuRepName==KW_NMDA)
        {
          Nets[j].ParNMDA.tst  = PreSynR->Tau*1e-3f;
          Nets[j].ParNMDA.Vrev = PreSynR->Vrev*1e-3f;
          Nets[j].p3=PreSynR->FExt*DTIME;
        }

      }

      Nets[j].rnd_seed=seed_gen();
      Nets[j].rdgen.seed(Nets[j].rnd_seed);
    }

  }

//-----------------------------------------------------------------------------
//-----------------Confs transcripts to Nets
  for(unsigned int i=0;i<ConfSta.size();i++) // pre synatic population
  {
    cout<<"setup network:"<<100.0f*(float)(i+1)/ConfSta.size()<<"%        \r"<<flush;

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
        { IdxNum[j]=ConfSta[PoSynPp].NeuNum; }
        else if(state==IdxIns)
        { IdxNum[j]=IdxNum[j-1]+ConfSta[PoSynPp].NeuNum; }
        else if(state==IdxFiPre)
        { IdxNum[j]=IdxNum[c]; }
        else if(state==IdxRe)
        {
          IdxNum[j]=IdxNum[a]+ConfSta[PoSynPp].NeuNum;
          state=IdxIns;
        }
      }

      a=0;
      for(unsigned int k=0;k<IdxNum.size();k++)
      { if(a<IdxNum[k]) a=IdxNum[k]; }

      for(unsigned int k=0;k<ConfSta[i].NeuNum;k++)
      { Nets[NIdx[i]+k].Index.insert(Nets[NIdx[i]+k].Index.end(),a,iniIdx); }

      TolSynp += ConfSta[i].NeuNum*a;
    }



    for(unsigned int j=0;j<ConfSta[i].PpLinkNum;j++) // pre synaptic population synapse
    {
      // initial and insert post neuron's receptors
      for(unsigned int n=0;n<ConfSta[PoSynPp].NeuNum;n++) // post synaptic population
      {
        // define parameters from synapse
        iniVar1.g = PoSyn->MeanEff*1e-9f * PoSyn->weight;
        iniVar2.g = PoSyn->MeanEff*1e-9f * PoSyn->weight;

        if(PoSyn->TarRep == KW_AMPA)  // initial and insert post neuron's receptors
        {
          Nets[NIdx[PoSynPp]+n].RepAMPA.insert(Nets[NIdx[PoSynPp]+n].RepAMPA.end(),ConfSta[i].NeuNum,iniVar1);
          Nets[NIdx[PoSynPp]+n].AMPA_Sti.insert(Nets[NIdx[PoSynPp]+n].AMPA_Sti.end(),ConfSta[i].NeuNum,false);
        }
        else if(PoSyn->TarRep == KW_GABA)
        {
          Nets[NIdx[PoSynPp]+n].RepGABA.insert(Nets[NIdx[PoSynPp]+n].RepGABA.end(),ConfSta[i].NeuNum,iniVar1);
          Nets[NIdx[PoSynPp]+n].GABA_Sti.insert(Nets[NIdx[PoSynPp]+n].GABA_Sti.end(),ConfSta[i].NeuNum,false);
        }
        else if(PoSyn->TarRep == KW_NMDA)
        {
          Nets[NIdx[PoSynPp]+n].RepNMDA.insert(Nets[NIdx[PoSynPp]+n].RepNMDA.end(),ConfSta[i].NeuNum,iniVar2);
          Nets[NIdx[PoSynPp]+n].NMDA_Sti.insert(Nets[NIdx[PoSynPp]+n].NMDA_Sti.end(),ConfSta[i].NeuNum,false);
        }

        TolRep += ConfSta[i].NeuNum;
      }

      // initail pre neurons connection table
      // must be built after receptors built for using "Nets[NIdx[k]+n].RepAMPA.size()+n"
      // build connection table of each neuron
      for(unsigned int m=0;m<ConfSta[i].NeuNum;m++) // pre synaptic population
      {
        for(unsigned int n=0;n<ConfSta[PoSynPp].NeuNum;n++) // post synaptic population
        {
          Nets[NIdx[i]+m].Index[IdxBgn+n].NeuID=NIdx[PoSynPp]+n;

          if(PoSyn->TarRep == KW_AMPA)  // initial and insert post neuron's receptors "connection" properties
          {
            Nets[NIdx[i]+m].Index[IdxBgn+n].AMPA_ID=(unsigned short)(Nets[NIdx[PoSynPp]+n].RepAMPA.size()-ConfSta[i].NeuNum+m);

            if(NIdx[i]+m != NIdx[PoSynPp]+n) // prevent neuron self-link
              Nets[NIdx[i]+m].Index[IdxBgn+n].Act = Nets[NIdx[i]+m].Index[IdxBgn+n].Act | AMPA_ACT;
          }
          else if(PoSyn->TarRep == KW_GABA)
          {
            Nets[NIdx[i]+m].Index[IdxBgn+n].GABA_ID=(unsigned short)(Nets[NIdx[PoSynPp]+n].RepGABA.size()-ConfSta[i].NeuNum+m);

            if(NIdx[i]+m != NIdx[PoSynPp]+n) // prevent neuron self-link
              Nets[NIdx[i]+m].Index[IdxBgn+n].Act = Nets[NIdx[i]+m].Index[IdxBgn+n].Act | GABA_ACT;
          }
          else if(PoSyn->TarRep == KW_NMDA)
          {
            Nets[NIdx[i]+m].Index[IdxBgn+n].NMDA_ID=(unsigned short)(Nets[NIdx[PoSynPp]+n].RepNMDA.size()-ConfSta[i].NeuNum+m);

            if(NIdx[i]+m != NIdx[PoSynPp]+n) // prevent neuron self-link
              Nets[NIdx[i]+m].Index[IdxBgn+n].Act = Nets[NIdx[i]+m].Index[IdxBgn+n].Act | NMDA_ACT;
          }
        }
      }
    } // end of pre synaptic population synapse

    {vector<unsigned int> ().swap(IdxNum);}
    delete [] (Confs+i)->PpLink;

  } // end of pre synatic population
  delete [] Confs;
  cout<<"\n";


//--------------------build FExt receptors
  for(unsigned int i=0;i<ConfSta.size();i++)
  {
    for(unsigned int j=NIdx[i];j<(NIdx[i]+ConfSta[i].NeuNum);j++) // neurons in each population
    {
      for(unsigned int k=0;k<ConfSta[i].PpRepNum;k++)  // variable in FExt receptors
      {
        iniVar1.g=PreSynR->MnExtEff*1e-9f*PreSynR->MnExtCon;
        iniVar2.g=PreSynR->MnExtEff*1e-9f*PreSynR->MnExtCon;

        if(PreSynR->NeuRepName==KW_AMPA)
        {
          Nets[j].RepAMPA.push_back(iniVar1);
          Nets[j].AMPA_Sti.push_back(false);
        }
        else if(PreSynR->NeuRepName==KW_GABA)
        {
          Nets[j].RepGABA.push_back(iniVar1);
          Nets[j].GABA_Sti.push_back(false);
        }
        else if(PreSynR->NeuRepName==KW_NMDA)
        {
          Nets[j].RepNMDA.push_back(iniVar2);
          Nets[j].NMDA_Sti.push_back(false);
        }
      }
    }
  }

  SetOut();
}  //end of SetPp()

void NetSim_LIF3::DoEvt()
{
  for(unsigned int i=iEvt;i<Evts.size();i++)
  {
    if((Evts[i].TimEv>=EvTime) && (Evts[i].TimEv<(EvTime+DTIME)))
    {
      iEvt=i;

      if(Evts[i].Type=="ChangeExtFreq")
      {
        for(unsigned int m=NIdx[Evts[i].Pp];m<(NIdx[Evts[i].Pp]+ConfSta[Evts[i].Pp].NeuNum);m++)
        {
          if(Evts[i].SubTyp==KW_AMPA)
          { Nets[m].p1=(float)(Evts[i].Var1*DTIME); }
          else if(Evts[i].SubTyp==KW_GABA)
          { Nets[m].p2=(float)(Evts[i].Var1*DTIME); }
          else if(Evts[i].SubTyp==KW_NMDA)
          { Nets[m].p3=(float)(Evts[i].Var1*DTIME); }
          else if(Evts[i].SubTyp==KW_SINE)
          { Nets[m].RSin.SetPar(FExtSinAmp,Evts[i].Var1,0.0f); }
        }
      }
      else if(Evts[i].Type=="ChangeMembraneNoise")
      {
        for(unsigned int m=NIdx[Evts[i].Pp];m<(NIdx[Evts[i].Pp]+ConfSta[Evts[i].Pp].NeuNum);m++)
        {
          normal_distribution<float>::param_type par(Evts[i].Var1*1e-9f,Evts[i].Var2*1e-9f);
          Nets[m].MemNoise.param(par);
        }
      }
      else if(Evts[i].Type=="Resume")
      {
        ifstream fin;
        unsigned int IdxSize;
        unsigned int AMPSize;
        unsigned int GABSize;
        unsigned int NMDSize;
        bool Stitmp;
        string s1;
        fin.open(Evts[i].Label.c_str());

        if(NIdx[Evts[i].Pp]!=0) // seek position of first neuron needed to resume
        {
          unsigned int k;
          bool finded=false;
          while(!finded)
          {
             fin>>s1;
             if(s1=="NID")
             {
                fin>>k;
                if(k==(NIdx[Evts[i].Pp]-1))
                {
                    while(!finded)
                    {
                      fin>>s1;
                      if(s1=="EndNID")
                      {  finded=true; }
                    }
                }
             }
          }
        }

        if(Evts[i].SubTyp==RSM_SYNCOND)
        {
          for(unsigned int m=NIdx[Evts[i].Pp];m<(NIdx[Evts[i].Pp]+ConfSta[Evts[i].Pp].NeuNum);m++)
          {
            fin>>s1>>s1
                >>s1>>s1>>s1>>s1>>s1>>s1>>s1>>s1>>s1>>s1>>s1>>s1
                >>IdxSize
                >>AMPSize
                >>GABSize
                >>NMDSize;
//----------------------------axion--------------------------------------
            for(unsigned int i=0;i<IdxSize;i++)
            {  fin>>s1>>s1>>s1>>s1>>s1; }
//----------------------------receptors--------------------------------------
            fin>>s1>>s1;
            for(unsigned int i=0;i<AMPSize;i++)
            {
              fin>>s1
                 >>Nets[m].RepAMPA[i].g
                 >>s1;
            }

            fin>>s1>>s1;
            for(unsigned int i=0;i<GABSize;i++)
            {
              fin>>s1
                 >>Nets[m].RepGABA[i].g
                 >>s1;
            }

            fin>>s1>>s1>>s1>>s1>>s1;
            for(unsigned int i=0;i<NMDSize;i++)
            {
              fin>>s1
                 >>Nets[m].RepNMDA[i].g
                 >>s1
                 >>s1;
            }
//----------------------------noise--------------------------------------
            fin>>s1>>s1>>s1>>s1>>s1>>s1
                >>s1;
          }
        }
        else if(Evts[i].SubTyp==RSM_GATVAR)
        {
          for(unsigned int m=NIdx[Evts[i].Pp];m<(NIdx[Evts[i].Pp]+ConfSta[Evts[i].Pp].NeuNum);m++)
          {
            fin>>s1>>s1
                >>s1>>s1>>s1>>s1>>s1>>s1>>s1>>s1>>s1>>s1>>s1>>s1
                >>IdxSize
                >>AMPSize
                >>GABSize
                >>NMDSize;
//----------------------------axion--------------------------------------
            for(unsigned int i=0;i<IdxSize;i++)
            {  fin>>s1>>s1>>s1>>s1>>s1;  }
//----------------------------receptors--------------------------------------
            fin>>s1>>s1;
            for(unsigned int i=0;i<AMPSize;i++)
            {
              fin>>s1
                 >>s1
                 >>Nets[m].RepAMPA[i].st;
            }

            fin>>s1>>s1;
            for(unsigned int i=0;i<GABSize;i++)
            {
              fin>>s1
                 >>s1
                 >>Nets[m].RepGABA[i].st;
            }

            fin>>s1>>s1>>s1>>s1>>s1;
            for(unsigned int i=0;i<NMDSize;i++)
            {
              fin>>s1
                 >>s1
                 >>Nets[m].RepNMDA[i].st
                 >>Nets[m].RepNMDA[i].xt;
            }
//----------------------------noise--------------------------------------
            fin>>s1>>s1>>s1>>s1>>s1>>s1
                >>s1;
          }
        }
        else if(Evts[i].SubTyp==RSM_MEMPOT)
        {
          for(unsigned int m=NIdx[Evts[i].Pp];m<(NIdx[Evts[i].Pp]+ConfSta[Evts[i].Pp].NeuNum);m++)
          {
            fin>>s1>>s1
                >>s1>>s1>>s1>>s1>>s1>>s1>>s1>>s1
                >>Nets[m].V
                >>s1>>s1>>s1
                >>IdxSize
                >>AMPSize
                >>GABSize
                >>NMDSize;
//----------------------------axion--------------------------------------
            for(unsigned int i=0;i<IdxSize;i++)
            {  fin>>s1>>s1>>s1>>s1>>s1; }
//----------------------------receptors--------------------------------------
            fin>>s1>>s1;
            for(unsigned int i=0;i<AMPSize;i++)
            {  fin>>s1>>s1>>s1; }

            fin>>s1>>s1;
            for(unsigned int i=0;i<GABSize;i++)
            {  fin>>s1>>s1>>s1; }

            fin>>s1>>s1>>s1>>s1>>s1;
            for(unsigned int i=0;i<NMDSize;i++)
            {  fin>>s1>>s1>>s1>>s1; }
//----------------------------noise--------------------------------------
            fin>>s1>>s1>>s1>>s1>>s1>>s1
                >>s1;
          }
        }
        else if(Evts[i].SubTyp==RSM_RNDSED)
        {
          for(unsigned int m=NIdx[Evts[i].Pp];m<(NIdx[Evts[i].Pp]+ConfSta[Evts[i].Pp].NeuNum);m++)
          {
            fin>>s1>>s1
                >>s1>>s1>>s1>>s1>>s1>>s1>>s1>>s1>>s1>>s1>>s1>>s1
                >>IdxSize
                >>AMPSize
                >>GABSize
                >>NMDSize;
//----------------------------axion--------------------------------------
            for(unsigned int i=0;i<IdxSize;i++)
            {  fin>>s1>>s1>>s1>>s1>>s1; }
//----------------------------receptors--------------------------------------
            fin>>s1>>s1;
            for(unsigned int i=0;i<AMPSize;i++)
            {  fin>>s1>>s1>>s1; }

            fin>>s1>>s1;
            for(unsigned int i=0;i<GABSize;i++)
            {  fin>>s1>>s1>>s1; }

            fin>>s1>>s1>>s1>>s1>>s1;
            for(unsigned int i=0;i<NMDSize;i++)
            {  fin>>s1>>s1>>s1>>s1; }
//----------------------------noise--------------------------------------
            fin>>s1>>s1>>s1>>s1>>s1
                >>Nets[m].rnd_seed
                >>s1;

            Nets[m].rdgen.seed(Nets[m].rnd_seed);
          }
        }
        else if(Evts[i].SubTyp==RSM_ALL)
        {
          for(unsigned int m=NIdx[Evts[i].Pp];m<(NIdx[Evts[i].Pp]+ConfSta[Evts[i].Pp].NeuNum);m++)
          {
            fin>>s1>>s1
                >>Nets[m].Vth
                >>Nets[m].El
                >>Nets[m].gl
                >>Nets[m].Cm
                >>Nets[m].Vreset
                >>Nets[m].RfPed
                >>Nets[m].SpiDy
                >>Nets[m].Isyn
                >>Nets[m].V
                >>Nets[m].Count
                >>Nets[m].Firing
                >>Nets[m].SpiOut
                >>IdxSize
                >>AMPSize
                >>GABSize
                >>NMDSize;
//----------------------------axion--------------------------------------
            for(unsigned int i=0;i<IdxSize;i++)
            {
              fin>>Nets[m].Index[i].NeuID
                 >>Nets[m].Index[i].AMPA_ID
                 >>Nets[m].Index[i].GABA_ID
                 >>Nets[m].Index[i].NMDA_ID
                 >>Nets[m].Index[i].Act;
            }
//----------------------------receptors--------------------------------------

            fin>>Nets[m].ParAMPA.tst
                >>Nets[m].ParAMPA.Vrev;

            for(unsigned int i=0;i<AMPSize;i++)
            {
              fin>>Stitmp
                 >>Nets[m].RepAMPA[i].g
                 >>Nets[m].RepAMPA[i].st;
              Nets[m].AMPA_Sti[i]=Stitmp;
            }

            fin>>Nets[m].ParGABA.tst
                >>Nets[m].ParGABA.Vrev;
            for(unsigned int i=0;i<GABSize;i++)
            {
              fin>>Stitmp
                 >>Nets[m].RepGABA[i].g
                 >>Nets[m].RepGABA[i].st;
              Nets[m].GABA_Sti[i]=Stitmp;
            }

            fin>>Nets[m].ParNMDA.as
                >>Nets[m].ParNMDA.tst
                >>Nets[m].ParNMDA.txt
                >>Nets[m].ParNMDA.mg
                >>Nets[m].ParNMDA.Vrev;
            for(unsigned int i=0;i<NMDSize;i++)
            {
              fin>>Stitmp
                 >>Nets[m].RepNMDA[i].g
                 >>Nets[m].RepNMDA[i].st
                 >>Nets[m].RepNMDA[i].xt;
              Nets[m].NMDA_Sti[i]=Stitmp;
            }
//----------------------------noise--------------------------------------
            float Me,Std;

            fin>>Me
                >>Std
                >>Nets[m].p1
                >>Nets[m].p2
                >>Nets[m].p3
                >>Nets[m].rnd_seed
                >>s1;

            normal_distribution<float>::param_type par(Me,Std*1e-9f);
            Nets[m].MemNoise.param(par);

            Nets[m].rdgen.seed(Nets[m].rnd_seed);

          }
        }  //end of ALL

        fin.close();
      } // end of resume

    } //  at event time end
  }
}

void NetSim_LIF3::RunInit(int thid,int sop)
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

//-------show total information with large network
  SInfo(s1);
  logdata<<s1<<"\n";
  logdata<<"population information:-------------------------------------------\n";
  for(unsigned int i=0;i<Nets.size();i++)
  {
    Nets[i].Info(s1);
    logdata<<"neuron:"<<i<<"\n"<<s1<<flush;
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

  {string ().swap(s1);}
  logdata.close();

  cout<<"start processing----------------------------------------------------\n";
}

void NetSim_LIF3::Runsteps(int steps)
{
  for(int Step_end=Step+steps; Step < Step_end; Step++)  //calculation
  {
    EvTime=Step*DTIME;
    EqsTime=Step*DTIME;

    RunActSteps(0,Nets.size());
  }
}

void NetSim_LIF3::Run()
{
  clock_t t1,t2;
  int Endtime = (int)floor(Evts[Evts.size()-1].TimEv/DTIME);

//iniital every neuron
  for(unsigned int i=0;i<Nets.size();i++)
  { Nets[i].VarInit();  }

//set every neuron to default external stimulation
  for(unsigned int i=0;i<ConfSta.size();i++)
  {
    for(unsigned int j=NIdx[i];j<(NIdx[i]+ConfSta[i].NeuNum);j++) // set up each neuron's parameter in every population
    {
      for(unsigned int k=0;k<ConfSta[i].PpRepNum;k++)  // parameter in receptors and variable in FExt receptors
      {
        if(PreSynR->NeuRepName==KW_AMPA)
        { Nets[j].p1=PreSynR->FExt*DTIME; }
        else if(PreSynR->NeuRepName==KW_GABA)
        { Nets[j].p2=PreSynR->FExt*DTIME; }
        else if(PreSynR->NeuRepName==KW_NMDA)
        { Nets[j].p3=PreSynR->FExt*DTIME; }
      }
#ifdef CPP11
      normal_distribution<float>::param_type par(MNoiseMean,MNoiseSTD);
      Nets[j].MemNoise.param(par);
#else
      Nets[j].mean=MNoiseMean;
      Nets[j].stddev=MNoiseSTD;
#endif
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

  cout<<"Iteration:"<<Repeat<<"\n";

  FileOpen();

//write gating variable print out information
  for(vector<OutCtr>::size_type i=0;i<OFCtr.size();i++)
  {
    if(GatVarID == OFCtr[i].TypeID)
    {
      ostr_tmp<<"output format:\n"
               <<"time:1";

      for(unsigned int j=0,jj=2;j<OFCtr[i].ParList.size();j++,jj++)
      {
        for(unsigned int k=NIdx[OFCtr[i].ParList[j]];k<(NIdx[OFCtr[i].ParList[j]]+ConfSta[OFCtr[i].ParList[j]].NeuNum);k++)
        {
          ostr_tmp<<" Population="<<OFCtr[i].ParList[j];
          ostr_tmp<<" NeuronID="<<k;

          if(Nets[k].RepAMPA.size()!=0)
          {
            for(unsigned int m=0;m<Nets[k].RepAMPA.size()-1;m++)
            { ostr_tmp<<" AMPA_st:"<<jj++; }

            ostr_tmp<<" External_AMPA_st:"<<jj++;
          }

          if(Nets[k].RepGABA.size()!=0)
          {
            for(unsigned int m=0;m<Nets[k].RepGABA.size()-1;m++)
            { ostr_tmp<<" GABA_st:"<<jj++; }

            ostr_tmp<<" External_GABA_st:"<<jj++;
          }

          if(Nets[k].RepNMDA.size()!=0)
          {
            for(unsigned int m=0;m<Nets[k].RepNMDA.size()-1;m++)
            {
              ostr_tmp<<" NMDA_xt:"<<jj++;
              ostr_tmp<<" NMDA_st:"<<jj++;
            }

            ostr_tmp<<" External_NMDA_xt:"<<jj++;
            ostr_tmp<<" External_NMDA_st:"<<jj++;
          }

          ostr_tmp<<",";
        }

      }

      ostr_tmp<<"\n";
    }
    else if(IsynSepSumID==OFCtr[i].TypeID)
    {
      ostr_tmp<<"output format:\n"
               <<"time:1";

      for(unsigned int j=0,jj=2;j<OFCtr[i].ParList.size();j++,jj++)
      {
        for(unsigned int k=NIdx[OFCtr[i].ParList[j]];k<(NIdx[OFCtr[i].ParList[j]]+ConfSta[OFCtr[i].ParList[j]].NeuNum);k++)
        {
          ostr_tmp<<" Population="<<OFCtr[i].ParList[j];
          ostr_tmp<<" NeuronID="<<k;

          ostr_tmp<<" I_AMPA:"<<jj++;
          ostr_tmp<<" I_GABA:"<<jj++;
          ostr_tmp<<" I_NMDA:"<<jj++;
        }
      }
    }
    if(SaveID == OFCtr[i].TypeID)
    {
        //ostr_tmp<<conf;
    }

    *(OFiles+i)<<ostr_tmp.str();
    ostr_tmp.str("");
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

void NetSim_LIF3::FileOut()
{
  string s1;
  for(vector<OutCtr>::size_type i=0;i<OFCtr.size();i++)
  {
    switch(OFCtr[i].TypeID)
    {
      case IsynSepSumID:
      ostr_tmp<<setprecision(7);
      ostr_tmp<<EvTime;
      ostr_tmp<<setprecision(6);

      for(unsigned int j=0;j<OFCtr[i].ParList.size();j++)
      {
        for(unsigned int k=NIdx[OFCtr[i].ParList[j]];k<(NIdx[OFCtr[i].ParList[j]]+ConfSta[OFCtr[i].ParList[j]].NeuNum);k++)
        {
          for(unsigned int m=0;m<3;m++)
          { ostr_tmp<<" "<<Nets[k].Irep[m]; }
        }
      }

      *(OFiles+i)<<ostr_tmp.str()<<endl;
      ostr_tmp.str("");
      break;

      case MemPotID:
      ostr_tmp<<setprecision(7);
      ostr_tmp<<EvTime;
      ostr_tmp<<setprecision(6);

      for(unsigned int j=0;j<OFCtr[i].ParList.size();j++)
      {
        for(unsigned int k=NIdx[OFCtr[i].ParList[j]];k<(NIdx[OFCtr[i].ParList[j]]+ConfSta[OFCtr[i].ParList[j]].NeuNum);k++)
        { ostr_tmp<<" "<<Nets[k].V; }
      }

      *(OFiles+i)<<ostr_tmp.str()<<endl;
      ostr_tmp.str("");
      break;

      case SpikeID:
      ostr_tmp<<setprecision(7)<<noshowpoint;

      for(unsigned int j=0;j<OFCtr[i].ParList.size();j++)
      {
        for(unsigned int k=NIdx[OFCtr[i].ParList[j]];k<(NIdx[OFCtr[i].ParList[j]]+ConfSta[OFCtr[i].ParList[j]].NeuNum);k++)
        {
          if(Nets[k].SpiOut)
          { ostr_tmp<<EvTime<<" "<<k<<"\n"; }
        }
      }

      s1=ostr_tmp.str();
      if(s1.size()!=0)
      {
        *(OFiles+i)<<ostr_tmp.str()<<flush;
        ostr_tmp.str("");
      }
      break;

      case FRateID:
      for(unsigned int j=0;j<OFCtr[i].ParList.size();j++)
      {
        for(unsigned int k=NIdx[OFCtr[i].ParList[j]];k<(NIdx[OFCtr[i].ParList[j]]+ConfSta[OFCtr[i].ParList[j]].NeuNum);k++)
        {
          if(Nets[k].SpiOut)
          { OFCtr[i].SpiRWin[OFCtr[i].ParList.size()*OFCtr[i].RWin+j]++; }
        }
      }

      if(Step%OneMilliSecond==0) //  check per 1 mS
      {
        int Step1ms=(int)(Step/OneMilliSecond);
        int m=Step1ms%OFCtr[i].RWin;

        if(Step1ms%OFCtr[i].FPnt==0) // default is 100 mS , print out to file
        {
          ostr_tmp<<setprecision(7);
          ostr_tmp<<EvTime;
          ostr_tmp<<setprecision(6);

          for(unsigned int j=0;j<OFCtr[i].ParList.size();j++)
          { ostr_tmp<<" "<<(1e3f * OFCtr[i].SpiRWin[OFCtr[i].ParList.size()*(OFCtr[i].RWin+1)+j]/OFCtr[i].RWin/ConfSta[j].NeuNum ); }

          *(OFiles+i)<<ostr_tmp.str()<<endl;
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
      ostr_tmp<<setprecision(7);
      ostr_tmp<<EvTime;
      ostr_tmp<<setprecision(6);

      for(unsigned int j=0;j<OFCtr[i].ParList.size();j++)
      {
        for(unsigned int k=NIdx[OFCtr[i].ParList[j]];k<(NIdx[OFCtr[i].ParList[j]]+ConfSta[OFCtr[i].ParList[j]].NeuNum);k++)
        {
          for(unsigned int m=0;m<Nets[k].RepAMPA.size();m++)
          { ostr_tmp<<" "<<Nets[k].RepAMPA[m].st; }

          for(unsigned int m=0;m<Nets[k].RepGABA.size();m++)
          { ostr_tmp<<" "<<Nets[k].RepGABA[m].st; }

          for(unsigned int m=0;m<Nets[k].RepNMDA.size();m++)
          {
            ostr_tmp<<" "<<Nets[k].RepNMDA[m].xt;
            ostr_tmp<<" "<<Nets[k].RepNMDA[m].st;
          }
        }
      }

      *(OFiles+i)<<ostr_tmp.str()<<"\n";
      ostr_tmp.str("");
      break;


      case SaveID:
      if((OFCtr[i].TimEv>=EvTime)&&(OFCtr[i].TimEv<(EvTime+DTIME)))
      {
        for(unsigned int j=0;j<OFCtr[i].ParList.size();j++)
        {
          for(unsigned int k=NIdx[OFCtr[i].ParList[j]];k<(NIdx[OFCtr[i].ParList[j]]+ConfSta[OFCtr[i].ParList[j]].NeuNum);k++)
          {
            Nets[k].TexSInfo(s1);
            ostr_tmp<<"NID "<<k<<"\n"<<s1<<"EndNID"<<endl;

            *(OFiles+i)<<ostr_tmp.str();
            ostr_tmp.str("");
          }
        }
      }
      break;


      default:;
    }

  }

  {string ().swap(s1);}

}


