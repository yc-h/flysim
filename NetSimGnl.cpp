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
#include "NetSim.h"


using namespace std;
extern random_device seed_gen;

default_random_engine udf_reng;
uniform_int_distribution<unsigned int> udf_gen(0,RAND_MAX);//RAND_MAX


//------------------NetSimGnl
NetSimGnl::NetSimGnl()
{
  TolNeu=0;
  TolSynp=0;
  TolRep=0;
  EvTime=0.0f;
  Step=0;
  Repeat=1;
  dt=DTIME;

#ifdef RABBITMQ
  MQini();
  CMQini();

  thrbMQ=new thread;
  *thrbMQ=thread(&NetSimGnl::MQListen,this);
#endif

}

NetSimGnl::~NetSimGnl()
{
  delete [] Nets;
  {vector<TmEvt> ().swap(NEvts);};
  {vector<OutCtr> ().swap(NOFCtr);};
#ifdef RABBITMQ
  MQListop=true;
  MQclose(MQconn);
  MQclose(CMQconn);
  //thrbMQ->join();
  delete thrbMQ;
#endif
}

void NetSimGnl::SInfo(string& ostr) //short information
{
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
      <<"\nsolver:"<<SType
      <<", time step:"<<dt
      <<", total neurons="<<TolNeu
      <<", total synapses="<<TolSynp
      <<", total receptors="<<TolRep
      <<"\n";

  ostr=ostr_tmp.str();
  ostr_tmp.str("");
}


void NetSimGnl::TexInfo(string& ostr)
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

//--------------setup network population

void NetSimGnl::SetPp(int thNum,int sop)
{
  bool IdxNew;

  if((ActCtr & UDFRDSEED_ACT)!= 0)
    rdgen.seed(rseed);
  else
    rdgen.seed(seed_gen());


//-----struct data initial
//-----------------Confs transcripts to Populations
  vector<IdxSynapse> IdxSyn;
  IdxSynapse iniIdxSyn;
  iniIdxSyn.PosPp=0;
  iniIdxSyn.PosIDNum=0;
  iniIdxSyn.gFast=0.0f;
  iniIdxSyn.gSlow=0.0f;
  iniIdxSyn.Act=0;

  SynFast iniSynFast;
  iniSynFast.PosIDNum=0;
  iniSynFast.PosID=nullptr;
  iniSynFast.gFast=nullptr;
  iniSynFast.stFast=0.0f;

  SynSlow iniSynSlow;
  iniSynSlow.PosIDNum=0;
  iniSynSlow.PosID=nullptr;
  iniSynSlow.gSlow=nullptr;
  iniSynSlow.xtSlow=0.0f;
  iniSynSlow.stSlow=0.0f;

  SynGlut iniSynGlut;
  iniSynGlut.PosIDNum=0;
  iniSynGlut.PosID=nullptr;
  iniSynGlut.gFast=nullptr;
  iniSynGlut.stFast=0.0f;
  iniSynGlut.gSlow=nullptr;
  iniSynGlut.xtSlow=0.0f;
  iniSynGlut.stSlow=0.0f;

  SynGap iniSynCond;
  iniSynCond.PosIDNum=0;
  iniSynCond.PosID=nullptr;
  iniSynCond.gC=nullptr;

//-----------------create default neuron network
  dt=dtime;

  NIdx.insert(NIdx.end(),ConfSta.size(),0);
  unsigned int q=0;
  for(unsigned int i=0;i<ConfSta.size();i++)
  {
    NIdx[i]=q;
    q += ConfSta[i].NeuNum;
  }
  TolNeu=q;

  Nets = new NeuronGnl[TolNeu];

  for(unsigned int i=0;i<ConfSta.size();i++)
  {
    for(unsigned int j=NIdx[i];j<(NIdx[i]+ConfSta[i].NeuNum);j++) // set up each neuron's parameter in every population
    {
      (Nets+j)->ActCtr=ActCtr;
      (Nets+j)->ParInit(dt);

      //leaky integrate and fire
      (Nets+j)->Vth    = (Confs+i)->Vth;
      (Nets+j)->El     = (Confs+i)->El;

      if((Confs+i)->Taum > MIN_TAU)
        (Nets+j)->gl     = (Confs+i)->Cm/(Confs+i)->Taum;
      else
        (Nets+j)->gl     = (Confs+i)->Cm/MIN_TAU;

      (Nets+j)->Cm     = (Confs+i)->Cm;
      (Nets+j)->Vreset = (Confs+i)->Vreset;

      (Nets+j)->RfPed  = (Confs+i)->RfPed;
      (Nets+j)->SpiDy  = (Confs+i)->SpiDy;

      //ltp
      if(ConfLtp[i].tLP>MIN_TAU)
        (Nets+j)->tLP = ConfLtp[i].tLP;  // time constant of long term plsiticity factor
      else
        (Nets+j)->tLP = MIN_TAU;  // time constant of long term plsiticity factor

      (Nets+j)->SpiLPDy= ConfLtp[i].SpiLPDy; //back propogation spike delay
      (Nets+j)->PosISI  = ConfLtp[i].PosISI;
      (Nets+j)->NegISI  = ConfLtp[i].NegISI;

      //stp
      (Nets+j)->pv   = ConfStp[i].pv;

      if(ConfStp[i].tD>MIN_TAU)
        (Nets+j)->tD = ConfStp[i].tD;
      else
        (Nets+j)->tD = MIN_TAU;

      (Nets+j)->aF   = ConfStp[i].aF*1e-9f;

      if(ConfStp[i].tF>MIN_TAU)
        (Nets+j)->tF   = ConfStp[i].tF;
      else
        (Nets+j)->tF   = MIN_TAU;

      (Nets+j)->aCap = ConfStp[i].aCap*1e-6f;

      if(ConfStp[i].tCap>MIN_TAU)
        (Nets+j)->tCap = ConfStp[i].tCap;
      else
        (Nets+j)->tCap = MIN_TAU;

      (Nets+j)->aCar = ConfStp[i].aCar*1e-6f;

      if(ConfStp[i].tCar>MIN_TAU)
        (Nets+j)->tCar = ConfStp[i].tCar;
      else
        (Nets+j)->tCar = MIN_TAU;

      (Nets+j)->Irep = new float[5];
      (Nets+j)->pb = new float[5];
      (Nets+j)->ExtFst = new RV1[4];

      for(unsigned int ix=0;ix<5;ix++)
      {
        *(Nets[j].Irep+ix)=0.0f;
        *(Nets[j].pb+ix)=0.0f;
      }

      for(unsigned int ix=0;ix<4;ix++)
      {
        (Nets[j].ExtFst+ix)->st=0.0f;
        (Nets[j].ExtFst+ix)->g=0.0f;
      }

      Nets[j].ExtSlo.xt=0.0f;
      Nets[j].ExtSlo.st=0.0f;
      Nets[j].ExtSlo.g=0.0f;


      for(unsigned int k=0;k<ConfSta[i].PpRepNum;k++)  // parameter in receptors and variable in FExt receptors
      {

        if(PreSynR->NeuRepName==KW_AMPA)
        {
          if(PreSynR->Tau> MIN_TAU)
            Nets[j].ParAMPA.tst = PreSynR->Tau;
          else
            Nets[j].ParAMPA.tst = MIN_TAU;

          Nets[j].ParAMPA.Vrev = PreSynR->Vrev;
          *(Nets[j].pb)=PreSynR->FExt*dt*1e-3f;
          Nets[j].ExtFst->g=PreSynR->MnExtEff*PreSynR->MnExtCon;
        }
        else if(PreSynR->NeuRepName==KW_GABA)
        {
          Nets[j].dVmax = ((Confs+i)->Vth - PreSynR->Vrev)/dt;

          if(PreSynR->Tau> MIN_TAU)
            Nets[j].ParGABA.tst  = PreSynR->Tau;
          else
            Nets[j].ParGABA.tst = MIN_TAU;

          Nets[j].ParGABA.Vrev = PreSynR->Vrev;
          *(Nets[j].pb+1)=PreSynR->FExt*dt*1e-3f;
          (Nets[j].ExtFst+1)->g=PreSynR->MnExtEff*PreSynR->MnExtCon;
        }
        else if(PreSynR->NeuRepName==KW_NMDA)
        {
          if(PreSynR->TauXT> MIN_TAU)
            Nets[j].ParNMDA.txt  = PreSynR->TauXT;
          else
            Nets[j].ParNMDA.txt = MIN_TAU;

          if(PreSynR->Tau> MIN_TAU)
            Nets[j].ParNMDA.tst  = PreSynR->Tau;
          else
            Nets[j].ParNMDA.tst = MIN_TAU;

          Nets[j].ParNMDA.Vrev = PreSynR->Vrev;
          *(Nets[j].pb+2)=PreSynR->FExt*dt*1e-3f;
          Nets[j].ExtSlo.g=PreSynR->MnExtEff*PreSynR->MnExtCon;
        }
        else if(PreSynR->NeuRepName==KW_ACH)
        {
          if(PreSynR->Tau> MIN_TAU)
            Nets[j].ParACH.tst  = PreSynR->Tau;
          else
            Nets[j].ParACH.tst = MIN_TAU;

          Nets[j].ParACH.Vrev = PreSynR->Vrev;
          *(Nets[j].pb+3)=PreSynR->FExt*dt*1e-3f;
          (Nets[j].ExtFst+2)->g=PreSynR->MnExtEff*PreSynR->MnExtCon;
        }
        else if(PreSynR->NeuRepName==KW_GCL)
        {
          if(PreSynR->Tau> MIN_TAU)
            Nets[j].ParGCL.tst  = PreSynR->Tau;
          else
            Nets[j].ParGCL.tst = MIN_TAU;

          Nets[j].ParGCL.Vrev = PreSynR->Vrev;
          *(Nets[j].pb+4)=PreSynR->FExt*dt*1e-3f;
          (Nets[j].ExtFst+3)->g=PreSynR->MnExtEff*PreSynR->MnExtCon;
        }

      }

      (Nets+j)->rdgen=&rdgen;
    }

  }

  delete [] ConfStp;
  delete [] ConfLtp;
  delete [] ConfHH;

//-----------------------------------------------------------------------------
//-----------------Confs transcripts to Nets

  for(unsigned int i=0;i<ConfSta.size();i++) // pre synatic population
  {
    cout<<"setup network:"<<100.0f*(float)(i+1)/ConfSta.size()<<"%        \r"<<flush;

//synapse type setup
    for(unsigned int j=0;j<ConfSta[i].PpLinkNum;j++) // pre synaptic population synapse
    {
      IdxNew=true;
      unsigned int k;
      unsigned short SynAct = GetAct(PoSyn->TarRep);
      for(k=0;k<IdxSyn.size();k++)
      {
        if(IdxSyn[k].PosPp==PoSynPp)
        {
          IdxNew=false;
          break;
        }
      }

      if(IdxNew)
      {
        IdxSyn.push_back(iniIdxSyn);
        IdxSyn[IdxSyn.size()-1].PosPp = PoSynPp;
        IdxSyn[IdxSyn.size()-1].PosIDNum = (unsigned int)ceil((double)ConfSta[PoSynPp].NeuNum * PoSyn->Connectivity);
        IdxSyn[IdxSyn.size()-1].Act = SynAct;

        if((i==PoSynPp) && (ConfSta[PoSynPp].NeuNum==IdxSyn[IdxSyn.size()-1].PosIDNum) && !ConfSta[i].SelConEn) // no self connction
        { IdxSyn[IdxSyn.size()-1].PosIDNum--; }

        if(IdxSyn[IdxSyn.size()-1].PosIDNum==0 && !ConfSta[i].SelConEn) // no self connction
          IdxSyn.erase(IdxSyn.end()-1);
        else  // self connction
          TolSynp += IdxSyn[IdxSyn.size()-1].PosIDNum;
      }
      else
        IdxSyn[k].Act |= SynAct;

      if(SynAct==NMDA_ACT)
        IdxSyn[k].gSlow = PoSyn->MeanEff*PoSyn->weight;
      else
        IdxSyn[k].gFast = PoSyn->MeanEff*PoSyn->weight;

    }

// build up synapses

    for(unsigned int j=NIdx[i];j<(NIdx[i]+ConfSta[i].NeuNum);j++) // set up each neuron's parameter in every population
    {
      for(unsigned int k=0;k<IdxSyn.size();k++) //synapse statisitic
      {
         //Fast synapse
        if((IdxSyn[k].Act == AMPA_ACT) || (IdxSyn[k].Act == GABA_ACT) || (IdxSyn[k].Act == ACH_ACT) || (IdxSyn[k].Act == GCL_ACT))
        {
          if(IdxSyn[k].PosIDNum == 1)
          { (Nets+j)->SynFastSize1++; }
          else
          { (Nets+j)->SynFastSize++; }
        }

        if(IdxSyn[k].Act == NMDA_ACT)  //Slow synapse
        {
          if(IdxSyn[k].PosIDNum == 1)
          { (Nets+j)->SynSlowSize1++; }
          else
          { (Nets+j)->SynSlowSize++; }
        }

        if(IdxSyn[k].Act == (AMPA_ACT | NMDA_ACT))  //Glutamate synapse
        {
          if(IdxSyn[k].PosIDNum == 1)
          { (Nets+j)->SynGlutSize1++; }
          else
          { (Nets+j)->SynGlutSize++; }
        }

        if((IdxSyn[k].Act == GAP_ACT) || (IdxSyn[k].Act == CUBA_ACT))  //Gap junction or current base synapse
        {
          if(IdxSyn[k].PosIDNum == 1)
          { (Nets+j)->SynGapSize1++; }
          else
          { (Nets+j)->SynGapSize++; }
        }
      }

      if((Nets+j)->SynFastSize1!=0)
      {
        (Nets+j)->SynFst1 = new SynFast1[(Nets+j)->SynFastSize1]; // build all Fast1 synapse
        (Nets+j)->SFTyp1 = new unsigned char[(Nets+j)->SynFastSize1]; // build all Fast synapse type
      }

      if((Nets+j)->SynFastSize!=0)
      {
        (Nets+j)->SynFst = new SynFast[(Nets+j)->SynFastSize]; // build all Fast synapse
        (Nets+j)->SFTyp = new unsigned char[(Nets+j)->SynFastSize]; // build all Fast synapse type
      }

      if((Nets+j)->SynSlowSize1!=0)
      {
        (Nets+j)->SynSlo1 = new SynSlow1[(Nets+j)->SynSlowSize1]; // build Slow1 synapse
        (Nets+j)->SSTyp1 = new unsigned char[(Nets+j)->SynSlowSize1]; // build all Slow synapse type
      }

      if((Nets+j)->SynSlowSize!=0)
      {
        (Nets+j)->SynSlo = new SynSlow[(Nets+j)->SynSlowSize]; // build Slow synapse
        (Nets+j)->SSTyp = new unsigned char[(Nets+j)->SynSlowSize]; // build all Slow synapse type
      }

      if((Nets+j)->SynGlutSize1!=0)
      {
        (Nets+j)->SynGlu1 = new SynGlut1[(Nets+j)->SynGlutSize1]; // build glutamate1 synapse
        (Nets+j)->SGTyp1 = new SynGlutType[(Nets+j)->SynGlutSize1]; // build all glutamate synapse type
      }

      if((Nets+j)->SynGlutSize!=0)
      {
        (Nets+j)->SynGlu = new SynGlut[(Nets+j)->SynGlutSize]; // build glutamate synapse
        (Nets+j)->SGTyp = new SynGlutType[(Nets+j)->SynGlutSize]; // build all glutamate synapse type
      }

      if((Nets+j)->SynGapSize1!=0)
      { (Nets+j)->SynCnd1 = new SynGap1[(Nets+j)->SynGapSize1]; }// build all gap1 junction synapse

      if((Nets+j)->SynGapSize!=0)
      { (Nets+j)->SynCnd = new SynGap[(Nets+j)->SynGapSize]; }// build all gap junction synapse

//Fast synapse--------------------------
      unsigned int kk1=0,kk=0;
      for(unsigned int k=0;k<IdxSyn.size();k++) // build Fast synapse of j neuron post links
      {
        if(IdxSyn[k].PosIDNum == 1)
        {
          if(IdxSyn[k].Act == AMPA_ACT)
          { *((Nets+j)->SFTyp1+kk1)=KW_AMPA; }
          else if(IdxSyn[k].Act == GABA_ACT)
          { *((Nets+j)->SFTyp1+kk1)=KW_GABA; }
          else if(IdxSyn[k].Act == ACH_ACT)
          { *((Nets+j)->SFTyp1+kk1)=KW_ACH; }
          else if(IdxSyn[k].Act == GCL_ACT)
          { *((Nets+j)->SFTyp1+kk1)=KW_GCL; }

          if((IdxSyn[k].Act == AMPA_ACT) || (IdxSyn[k].Act == GABA_ACT) || (IdxSyn[k].Act == ACH_ACT) || (IdxSyn[k].Act == GCL_ACT))
          {
            TolRep += IdxSyn[k].PosIDNum;
            ((Nets+j)->SynFst1+kk1)->PosID=PosBuild1(j,IdxSyn[k].PosIDNum,NIdx[IdxSyn[k].PosPp],ConfSta[IdxSyn[k].PosPp].NeuNum,ConfSta[i].SelConEn);
            ((Nets+j)->SynFst1+kk1)->gFast=IdxSyn[k].gFast;
            kk1++;
          }
        }
        else
        {
          if(IdxSyn[k].Act == AMPA_ACT)
          { *((Nets+j)->SFTyp+kk)=KW_AMPA; }
          else if(IdxSyn[k].Act == GABA_ACT)
          { *((Nets+j)->SFTyp+kk)=KW_GABA; }
          else if(IdxSyn[k].Act == ACH_ACT)
          { *((Nets+j)->SFTyp+kk)=KW_ACH; }
          else if(IdxSyn[k].Act == GCL_ACT)
          { *((Nets+j)->SFTyp+kk)=KW_GCL; }

          if((IdxSyn[k].Act == AMPA_ACT) || (IdxSyn[k].Act == GABA_ACT) || (IdxSyn[k].Act == ACH_ACT) || (IdxSyn[k].Act == GCL_ACT))
          {
            TolRep += IdxSyn[k].PosIDNum;
            ((Nets+j)->SynFst+kk)->PosIDNum=IdxSyn[k].PosIDNum;
            ((Nets+j)->SynFst+kk)->PosID=PosBuild(j,IdxSyn[k].PosIDNum,NIdx[IdxSyn[k].PosPp],ConfSta[IdxSyn[k].PosPp].NeuNum,ConfSta[i].SelConEn);

            if(((Nets+j)->ActCtr & LTP_ACT) != 0)
            {
              ((Nets+j)->SynFst+kk)->gFast= new float[IdxSyn[k].PosIDNum];
              for(unsigned int jj=0;jj<IdxSyn[k].PosIDNum;jj++)
              { *(((Nets+j)->SynFst+kk)->gFast+jj)=IdxSyn[k].gFast; }
            }
            else
            {
              ((Nets+j)->SynFst+kk)->gFast= new float;
              *(((Nets+j)->SynFst+kk)->gFast)=IdxSyn[k].gFast;
            }

            kk++;
          }
        }
      }

//NMDA synapse--------------------------
      kk1=0;
      kk=0;
      for(unsigned int k=0;k<IdxSyn.size();k++) // build post neuron links in each Slow synapse
      {
        if(IdxSyn[k].PosIDNum == 1)
        {
          if(IdxSyn[k].Act == NMDA_ACT)
          {
            *((Nets+j)->SSTyp1+kk1)=KW_NMDA;

            TolRep += IdxSyn[k].PosIDNum;
            ((Nets+j)->SynSlo1+kk1)->PosID=PosBuild1(j,IdxSyn[k].PosIDNum,NIdx[IdxSyn[k].PosPp],ConfSta[IdxSyn[k].PosPp].NeuNum,ConfSta[i].SelConEn);
            ((Nets+j)->SynSlo1+kk1)->gSlow=IdxSyn[k].gSlow;
            kk1++;
          }
        }
        else
        {
          if(IdxSyn[k].Act == NMDA_ACT)
          {
            *((Nets+j)->SSTyp+kk)=KW_NMDA;

            TolRep += IdxSyn[k].PosIDNum;
            ((Nets+j)->SynSlo+kk)->PosIDNum=IdxSyn[k].PosIDNum;
            ((Nets+j)->SynSlo+kk)->PosID=PosBuild(j,IdxSyn[k].PosIDNum,NIdx[IdxSyn[k].PosPp],ConfSta[IdxSyn[k].PosPp].NeuNum,ConfSta[i].SelConEn);

            if(((Nets+j)->ActCtr & LTP_ACT) != 0)
            {
              ((Nets+j)->SynSlo+kk)->gSlow= new float[IdxSyn[k].PosIDNum];
              for(unsigned int jj=0;jj<IdxSyn[k].PosIDNum;jj++)
              { *(((Nets+j)->SynSlo+kk)->gSlow+jj)=IdxSyn[k].gSlow; }
            }
            else
            {
              ((Nets+j)->SynSlo+kk)->gSlow= new float;
              *(((Nets+j)->SynSlo+kk)->gSlow)=IdxSyn[k].gSlow;
            }

            kk++;
          }
        }
      }

//glutamate synapse------------------------------
      kk1=0;
      kk=0;
      for(unsigned int k=0;k<IdxSyn.size();k++) // build post neuron links in each glutamate synapse
      {
        if(IdxSyn[k].PosIDNum == 1)
        {
          if(IdxSyn[k].Act == (AMPA_ACT | NMDA_ACT))
          {
            ((Nets+j)->SGTyp1+kk1)->A=KW_AMPA;
            ((Nets+j)->SGTyp1+kk1)->N=KW_NMDA;

            TolRep += 2*IdxSyn[k].PosIDNum;
            ((Nets+j)->SynGlu1+kk1)->PosID=PosBuild1(j,IdxSyn[k].PosIDNum,NIdx[IdxSyn[k].PosPp],ConfSta[IdxSyn[k].PosPp].NeuNum,ConfSta[i].SelConEn);
            ((Nets+j)->SynGlu1+kk1)->gFast=IdxSyn[k].gFast;
            ((Nets+j)->SynGlu1+kk1)->gSlow=IdxSyn[k].gSlow;
            kk1++;
          }
        }
        else
        {
          if(IdxSyn[k].Act == (AMPA_ACT | NMDA_ACT))
          {
            ((Nets+j)->SGTyp+kk)->A=KW_AMPA;
            ((Nets+j)->SGTyp+kk)->N=KW_NMDA;

            TolRep += 2*IdxSyn[k].PosIDNum;
            ((Nets+j)->SynGlu+kk)->PosIDNum=IdxSyn[k].PosIDNum;
            ((Nets+j)->SynGlu+kk)->PosID=PosBuild(j,IdxSyn[k].PosIDNum,NIdx[IdxSyn[k].PosPp],ConfSta[IdxSyn[k].PosPp].NeuNum,ConfSta[i].SelConEn);

            if(((Nets+j)->ActCtr & LTP_ACT) != 0)
            {
              ((Nets+j)->SynGlu+kk)->gFast= new float[IdxSyn[k].PosIDNum];
              ((Nets+j)->SynGlu+kk)->gSlow= new float[IdxSyn[k].PosIDNum];

              for(unsigned int jj=0;jj<IdxSyn[k].PosIDNum;jj++)
              {
                *(((Nets+j)->SynGlu+kk)->gFast+jj)=IdxSyn[k].gFast;
                *(((Nets+j)->SynGlu+kk)->gSlow+jj)=IdxSyn[k].gSlow;
              }
            }
            else
            {
              ((Nets+j)->SynGlu+kk)->gFast= new float;
              ((Nets+j)->SynGlu+kk)->gSlow= new float;

              *(((Nets+j)->SynGlu+kk)->gFast)=IdxSyn[k].gFast;
              *(((Nets+j)->SynGlu+kk)->gSlow)=IdxSyn[k].gSlow;
            }

            kk++;
          }
        }
      }

//gap junction synapse---------------------------------
      kk1=0;
      kk=0;
      for(unsigned int k=0;k<IdxSyn.size();k++) // build post neuron links in each gap junction synapse
      {
        if(IdxSyn[k].PosIDNum == 1)
        {
          if((IdxSyn[k].Act == GAP_ACT) || (IdxSyn[k].Act == CUBA_ACT))
          {
            TolRep += IdxSyn[k].PosIDNum;
            ((Nets+j)->SynCnd1+kk1)->PosID=PosBuild1(j,IdxSyn[k].PosIDNum,NIdx[IdxSyn[k].PosPp],ConfSta[IdxSyn[k].PosPp].NeuNum,ConfSta[i].SelConEn);
            ((Nets+j)->SynCnd1+kk1)->gC=IdxSyn[k].gFast*1e-3f;
            kk1++;
          }
        }
        else
        {
          if((IdxSyn[k].Act == GAP_ACT) || (IdxSyn[k].Act == CUBA_ACT))
          {
            TolRep += IdxSyn[k].PosIDNum;
            ((Nets+j)->SynCnd+kk)->PosIDNum=IdxSyn[k].PosIDNum;
            ((Nets+j)->SynCnd+kk)->PosID=PosBuild(j,IdxSyn[k].PosIDNum,NIdx[IdxSyn[k].PosPp],ConfSta[IdxSyn[k].PosPp].NeuNum,ConfSta[i].SelConEn);

            if(((Nets+j)->ActCtr & LTP_ACT) != 0)
            {
              ((Nets+j)->SynCnd+kk)->gC= new float[IdxSyn[k].PosIDNum];
              for(unsigned int jj=0;jj<IdxSyn[k].PosIDNum;jj++)
              { *(((Nets+j)->SynCnd+kk)->gC+jj)=IdxSyn[k].gFast*1e-3f; }
            }
            else
            {
              ((Nets+j)->SynCnd+kk)->gC= new float;
              *(((Nets+j)->SynCnd+kk)->gC)=IdxSyn[k].gFast*1e-3f;
            }

            kk++;
          }
        }
      }
    }
    {vector<IdxSynapse> ().swap(IdxSyn);}
    delete [] (Confs+i)->PpLink;
  }

  SetOut();
}  //end of SetPp()


void NetSimGnl::DoNowEvt()
{
  for(unsigned int i=0;i<NEvts.size();i++)
  {
    if(NEvts[i].Type==EVTEXTFREQ)
    {
      for(unsigned int m=NIdx[NEvts[i].Pp];m<(NIdx[NEvts[i].Pp]+ConfSta[NEvts[i].Pp].NeuNum);m++)
      {
        switch(NEvts[i].SubTyp)
        {
          case KW_AMPA:
            *((Nets+m)->pb)=(float)(NEvts[i].Var1*dt*1e-3f);
          break;

          case KW_GABA:
            *((Nets+m)->pb+1)=(float)(NEvts[i].Var1*dt*1e-3f);
          break;

          case KW_NMDA:
            *((Nets+m)->pb+2)=(float)(NEvts[i].Var1*dt*1e-3f);
          break;

          case KW_ACH:
            *((Nets+m)->pb+3)=(float)(NEvts[i].Var1*dt*1e-3f);
          break;

          case KW_GCL:
            *((Nets+m)->pb+4)=(float)(NEvts[i].Var1*dt*1e-3f);
          break;

          default:;
        }
      }
    }
    else if(NEvts[i].Type==EVTCUNTINJ)
    {
      for(unsigned int m=NIdx[NEvts[i].Pp];m<(NIdx[NEvts[i].Pp]+ConfSta[NEvts[i].Pp].NeuNum);m++)
      {
        normal_distribution<float>::param_type par(NEvts[i].Var1,NEvts[i].Var2);
        (Nets+m)->Inj.param(par);
        (Nets+m)->SetIgsEn();
      }
    }
  }

  {vector<TmEvt> ().swap(NEvts);}
  {vector<OutCtr> ().swap(NOFCtr);}
}


#ifdef RABBITMQ
void NetSimGnl::NowEvtPrs(char *ch) // Now Event parse
{
  char *pch;
  string s1;

  if(strcmp(ch,"")!=0)
    pch = strtok (ch," :;,'=\t\r\n");
  else
    pch=ch;

  if(strcmp(pch,"CMD")==0)
  {
    pch = strtok (NULL," :;,'=\t\r\n");
    if(strcmp(pch,"UsrDefSti1")==0)
    {
      cout<<"\nNow ";
      string NevtsName = "UsrDefSti1_30.pro";
      parse(NevtsName,NEvts,NOFCtr); // for now network.pro

      pch = strtok (NULL," :;,'=\t\r\n");
      float sti=(float)atof(pch);
      cout<<"sti="<<sti<<endl;

      ostr2_tmp<<"{\"type\":\"MG\",\"text\":\"Ack:UsrDefSti1="<<sti<<"\"}";
      s1=ostr2_tmp.str();
      MQSend(s1);
      ostr2_tmp.str("");

      for(unsigned int i=0;i<NEvts.size();i++)
      { NEvts[i].Var1=sti; }
    }
  }
  strcpy(ch,"");
}
#endif

void NetSimGnl::DoEvt()
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
              *((Nets+m)->pb)=(float)(Evts[i].Var1*dt*1e-3f);
            break;

            case KW_GABA:
              *((Nets+m)->pb+1)=(float)(Evts[i].Var1*dt*1e-3f);
            break;

            case KW_NMDA:
              *((Nets+m)->pb+2)=(float)(Evts[i].Var1*dt*1e-3f);
            break;

            case KW_ACH:
              *((Nets+m)->pb+3)=(float)(Evts[i].Var1*dt*1e-3f);
            break;

            case KW_GCL:
              *((Nets+m)->pb+4)=(float)(Evts[i].Var1*dt*1e-3f);
            break;

            default:;
          }
        }
      }
      else if(Evts[i].Type==EVTCUNTINJ)
      {
        for(unsigned int m=NIdx[Evts[i].Pp];m<(NIdx[Evts[i].Pp]+ConfSta[Evts[i].Pp].NeuNum);m++)
        {
          normal_distribution<float>::param_type par(Evts[i].Var1,Evts[i].Var2);
          (Nets+m)->Inj.param(par);
          (Nets+m)->SetIgsEn();
        }
      }

    }
    else if(Evts[i].TimEv*1e3f>=(EvTime+dt))
    { break; }


  }
}

void NetSimGnl::RunInit(int thid,int sop)
{
  float TolLoad=0.0f;
  float AvgLoad=0.0f;

// membrane load calculation-------------------------------
  AvgLoad=(float)TolNeu/thid;
  ActNum = new unsigned int[thid];

  for(unsigned int i=0,k=0;i<TolNeu;i++)
  {
    TolLoad += 1.0f;

    if(TolLoad >= (float)AvgLoad*(k+1) )
    {
      *(ActNum+k)=i;
      k++;
    }
  }

  *(ActNum+thid-1)=TolNeu;

// synapse load calculation-------------------------------
  TolLoad=0.0f;
  for(unsigned int i=0;i<TolNeu;i++)
  {
    TolLoad +=
      (Nets+i)->SynFastSize1*1.2f
      + (Nets+i)->SynGlutSize1*3.2f
      + (Nets+i)->SynSlowSize1*2.0f;

    for(unsigned int j=0;j<(Nets+i)->SynFastSize;j++)
    { TolLoad += ((Nets+i)->SynFst+j)->PosIDNum*1.22f; }

    for(unsigned int j=0;j<(Nets+i)->SynGlutSize;j++)
    { TolLoad += ((Nets+i)->SynGlu+j)->PosIDNum*3.22f; }

    for(unsigned int j=0;j<(Nets+i)->SynSlowSize;j++)
    { TolLoad += ((Nets+i)->SynSlo+j)->PosIDNum*2.02f; }
  }

  AvgLoad=(float)TolLoad/thid;
  SpiNum = new unsigned int[thid];

  TolLoad=0.0f;
  for(unsigned int i=0,k=0;i<TolNeu;i++)
  {
    TolLoad +=
      (Nets+i)->SynFastSize1*1.2f
      + (Nets+i)->SynGlutSize1*3.2f
      + (Nets+i)->SynSlowSize1*2.0f;

    for(unsigned int j=0;j<(Nets+i)->SynFastSize;j++)
    { TolLoad += ((Nets+i)->SynFst+j)->PosIDNum*1.22f; }

    for(unsigned int j=0;j<(Nets+i)->SynGlutSize;j++)
    { TolLoad += ((Nets+i)->SynGlu+j)->PosIDNum*3.22f; }

    for(unsigned int j=0;j<(Nets+i)->SynSlowSize;j++)
    { TolLoad += ((Nets+i)->SynSlo+j)->PosIDNum*2.02f; }

    if(TolLoad >= (float)AvgLoad*(k+1) )
    {
      *(SpiNum+k)=i;
      k++;
    }
  }

  *(SpiNum+thid-1)=TolNeu;

  //dt=dtime;
  SolType=sop;
  iniEqsThd(Nets,thid);

  for(unsigned int i=0;i<TolNeu;i++)
  {
    (Nets+i)->gAMPA=new float[thid];
    (Nets+i)->gGABA=new float[thid];
    (Nets+i)->gNMDA=new float[thid];
    (Nets+i)->gACH=new float[thid];
    (Nets+i)->gGCL=new float[thid];

    for(unsigned int j=0;j<(unsigned int)thid;j++)
    {
      *((Nets+i)->gAMPA+j)=0.0f;
      *((Nets+i)->gGABA+j)=0.0f;
      *((Nets+i)->gNMDA+j)=0.0f;
      *((Nets+i)->gACH+j)=0.0f;
      *((Nets+i)->gGCL+j)=0.0f;
    }
  }
  
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
  if((ActCtr & UDFRDSEED_ACT)!= 0)
    logdata<<"User define random seed enable!\n";
  else
    logdata<<"internal random seed enable!\n";

  logdata<<"\npopulation information:-------------------------------------------\n";
  for(unsigned int i=0;i<TolNeu;i++)
  {
    (Nets+i)->Info(s1);
    logdata<<"-------------------------------------------\n"
           <<"neuron:"<<i<<"\n"<<s1<<flush;
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
      logdata<<"1:Time";

      for(unsigned int j=0,kk=2;j<OFCtr[i].ParList.size();j++)
      {
        for(unsigned int k=NIdx[OFCtr[i].ParList[j]];k<(NIdx[OFCtr[i].ParList[j]]+ConfSta[OFCtr[i].ParList[j]].NeuNum);k++)
        {
          logdata<<" NID="<<k;
          if((Nets+k)->SynFastSize1!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynFastSize1;m++)
            { logdata<<" "<<kk++<<":stFast="<<m; }
          }

          if((Nets+k)->SynSlowSize1!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynSlowSize1;m++)
            {
              logdata<<" "<<kk++<<":xtSlow="<<m;
              logdata<<" "<<kk++<<":stSlow="<<m;
            }
          }

          if((Nets+k)->SynGlutSize1!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynGlutSize1;m++)
            {
              logdata<<" "<<kk++<<":stFast="<<m;
              logdata<<" "<<kk++<<":xtSlow="<<m;
              logdata<<" "<<kk++<<":stSlow="<<m;
            }
          }

          if((Nets+k)->SynFastSize!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynFastSize;m++)
            { logdata<<" "<<kk++<<":stFast="<<m; }
          }

          if((Nets+k)->SynSlowSize!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynSlowSize;m++)
            {
              logdata<<" "<<kk++<<":xtSlow="<<m;
              logdata<<" "<<kk++<<":stSlow="<<m;
            }
          }

          if((Nets+k)->SynGlutSize!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynGlutSize;m++)
            {
              logdata<<" "<<kk++<<":stFast="<<m;
              logdata<<" "<<kk++<<":xtSlow="<<m;
              logdata<<" "<<kk++<<":stSlow="<<m;
            }
          }
        }
      }
      logdata<<"\n";
      break;

      case SynWgtID:
      logdata<<"PrintStep="<<OFCtr[i].FPnt<<"ms\n";
      logdata<<"1:Time";
      for(unsigned int j=0,kk=2;j<OFCtr[i].ParList.size();j++)
      {
        for(unsigned int k=NIdx[OFCtr[i].ParList[j]];k<(NIdx[OFCtr[i].ParList[j]]+ConfSta[OFCtr[i].ParList[j]].NeuNum);k++)
        {
          logdata<<" NID="<<k;
          if((Nets+k)->SynFastSize1!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynFastSize1;m++)
            { logdata<<" "<<kk++<<":gFast="<<m; }
          }

          if((Nets+k)->SynSlowSize1!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynSlowSize1;m++)
            { logdata<<" "<<kk++<<":gSlow="<<m; }
          }

          if((Nets+k)->SynGlutSize1!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynGlutSize1;m++)
            {
              logdata<<" "<<kk++<<":gFast="<<m;
              logdata<<" "<<kk++<<":gSlow="<<m;
            }
          }

          if(((Nets+k)->ActCtr & LTP_ACT) != 0)
          {
          if((Nets+k)->SynFastSize!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynFastSize;m++)
            {
              logdata<<" SynFast:"<<m;
              for(unsigned int n=0;n<((Nets+k)->SynFst+m)->PosIDNum;n++)
              { logdata<<" "<<kk++<<":gFast="<<n;}
            }
          }

          if((Nets+k)->SynSlowSize!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynSlowSize;m++)
            {
              logdata<<" SynSlow:"<<m;
              for(unsigned int n=0;n<((Nets+k)->SynSlo+m)->PosIDNum;n++)
              { logdata<<" "<<kk++<<":gSlow="<<n;}
            }
          }

          if((Nets+k)->SynGlutSize!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynGlutSize;m++)
            {
              logdata<<" SynGlut:"<<m;
              for(unsigned int n=0;n<((Nets+k)->SynGlu+m)->PosIDNum;n++)
                {
                logdata<<" "<<kk++<<":gFast="<<n;
                logdata<<" "<<kk++<<":gSlow="<<n;
                }

            }
          }
          }
          else
          {
          if((Nets+k)->SynFastSize!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynFastSize;m++)
            { logdata<<" "<<kk++<<":gFast="<<m; }
          }

          if((Nets+k)->SynSlowSize!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynSlowSize;m++)
            { logdata<<" "<<kk++<<":gSlow="<<m; }
          }

          if((Nets+k)->SynGlutSize!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynGlutSize;m++)
            {
              logdata<<" "<<kk++<<":gFast="<<m;
              logdata<<" "<<kk++<<":gSlow="<<m;
            }
          }
          }
        }
      }
      break;

      default:;
    }

  logdata<<endl;
  }
//-------------------------------------------------
  {string ().swap(s1);}
  logdata.close();

  cout<<"start processing----------------------------------------------------\n";


}

void NetSimGnl::Run()
{
  //clock_t t1,t2;
  typedef std::chrono::high_resolution_clock Clock;
  auto t1=Clock::now();
  auto t2=Clock::now();

  string s1;
  bool NoLmt=Evts[Evts.size()-1].NoLmt;
  unsigned long long Endtime = (unsigned long long)floor(Evts[Evts.size()-1].TimEv*1e3f/dt);
  unsigned int ReportStep=100/dt;
  OneMsStep=(unsigned int)1/dt;

  cout<<setprecision(7)<<showpoint<<flush;

//iniital variable in every neuron 
  for(unsigned int i=0;i<TolNeu;i++)
  { (Nets+i)->VarInit(thsize); }
//iniital every neuron
  for(unsigned int i=0;i<ConfSta.size();i++)
  {
    for(unsigned int j=NIdx[i];j<(NIdx[i]+ConfSta[i].NeuNum);j++) // set up each neuron's parameter in every population
    {
      for(unsigned int k=0;k<ConfSta[i].PpRepNum;k++)  // parameter in receptors and variable in FExt receptors
      {
        if(PreSynR->NeuRepName==KW_AMPA)
        { *(Nets[j].pb)=PreSynR->FExt*dt*1e-3f; }
        else if(PreSynR->NeuRepName==KW_GABA)
        { *(Nets[j].pb+1)=PreSynR->FExt*dt*1e-3f; }
        else if(PreSynR->NeuRepName==KW_NMDA)
        { *(Nets[j].pb+2)=PreSynR->FExt*dt*1e-3f; }
        else if(PreSynR->NeuRepName==KW_ACH)
        { *(Nets[j].pb+3)=PreSynR->FExt*dt*1e-3f; }
        else if(PreSynR->NeuRepName==KW_GCL)
        { *(Nets[j].pb+4)=PreSynR->FExt*dt*1e-3f; }
      }
    }
  }
//iniital every neuron LTP/LTD
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
  cout<<noshowpoint<<setprecision(7);

  FileOpen();

  t1=Clock::now();
  for(Step=0; (Step < Endtime) || NoLmt ; Step++)  //calculation
  {
    EvTime=Step*dt;

#ifdef RABBITMQ
    mtx.lock();
    NowEvtPrs(chbuff);
    DoNowEvt();
    mtx.unlock();
#endif

    DoEvt();
    if((Step % ReportStep)==0)
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
        <<100.0f*(float)(Step+1)/Endtime<<"%             \r"
        <<flush;
      }
#ifdef RABBITMQ
      if(NoLmt)
      {
        ostr2_tmp
        <<"{\"type\":\"MG\",\"text\":\""<<EvTime*1e-3f
        <<":"<<(chrono::duration_cast<chrono::nanoseconds>(t2-t1).count())*1e-9
        <<":infinity\""<<"}";
      }
      else
      {
        ostr2_tmp<<"{\"type\":\"MG\",\"text\":\""<<EvTime*1e-3f
        <<":"<<(chrono::duration_cast<chrono::nanoseconds>(t2-t1).count())*1e-9
        <<":"<<100.0f*(float)(Step+1)/Endtime<<"\%\"}";
      }
      s1=ostr2_tmp.str();
      MQSend(s1);
      ostr2_tmp.str("");
#endif
    }

    RunThd();
    FileOut();
  }

  t2=Clock::now();

  cout
  <<"Expertment time="<<EvTime*1e-3f<<"s "
  <<"real time="<<(chrono::duration_cast<chrono::nanoseconds>(t2-t1).count())*1e-9<<"s 100%                 "<<endl;

  FileClose();

}

void NetSimGnl::FileOut()
{
  unsigned int kk=0;
  string s1;
  istringstream sin;

#ifdef RABBITMQ
  float EvTime2=0.0f;
#endif

  for(vector<OutCtr>::size_type i=0;i<OFCtr.size();i++)
  {
    switch(OFCtr[i].TypeID)
    {
      case IsynSepSumID:
      ostr_tmp<<EvTime*1e-3f;

      for(unsigned int j=0;j<OFCtr[i].ParList.size();j++)
      {
        for(unsigned int k=NIdx[OFCtr[i].ParList[j]];k<(NIdx[OFCtr[i].ParList[j]]+ConfSta[OFCtr[i].ParList[j]].NeuNum);k++)
        {
          for(unsigned int m=0;m<5;m++)
          { ostr_tmp<<" "<<*((Nets+k)->Irep+m)*1e-12f; }
        }
      }

      *(OFiles+i)<<ostr_tmp.str()<<endl;
      ostr_tmp.str("");
      break;

      case MemPotID:
      ostr_tmp<<EvTime*1e-3f;

      for(unsigned int j=0;j<OFCtr[i].ParList.size();j++)
      {
        for(unsigned int k=NIdx[OFCtr[i].ParList[j]];k<(NIdx[OFCtr[i].ParList[j]]+ConfSta[OFCtr[i].ParList[j]].NeuNum);k++)
        { ostr_tmp<<" "<<(Nets+k)->V*1e-3f; }
      }

      *(OFiles+i)<<ostr_tmp.str()<<endl;
      ostr_tmp.str("");
      break;

      case SpikeID:

      for(unsigned int j=0;j<OFCtr[i].ParList.size();j++)
      {
        for(unsigned int k=NIdx[OFCtr[i].ParList[j]];k<(NIdx[OFCtr[i].ParList[j]]+ConfSta[OFCtr[i].ParList[j]].NeuNum);k++)
        {
          if((Nets+k)->SpiOut)
          {
            ostr_tmp<<EvTime*1e-3f<<" "<<k<<"\n";
            kk++;
          }
        }
      }

      if(kk != 0)
      {

#ifdef RABBITMQ
        unsigned int kkk=0;
        sin.str(ostr_tmp.str());
        if(kk==1)
        {
          sin>>EvTime2;
          sin>>kkk;
          ostr2_tmp<<"{\"type\":\"SP\",\"time\":"<<EvTime2<<",\"neuronID\":"<<kkk<<"}";
        }
        else // kk >=2
        {
          ostr2_tmp<<"[";
          for(unsigned int k=0;k<kk-1;k++)
          {
            sin>>EvTime2;
            sin>>kkk;
            ostr2_tmp<<"{\"type\":\"SP\",\"time\":"<<EvTime2<<",\"neuronID\":"<<kkk<<"}";
            ostr2_tmp<<",";
          }
          sin>>EvTime2;
          sin>>kkk;
          ostr2_tmp<<"{\"type\":\"SP\",\"time\":"<<EvTime2<<",\"neuronID\":"<<kkk<<"}";
          ostr2_tmp<<"]";
        }
        s1=ostr2_tmp.str();
        MQSend(s1);
        ostr2_tmp.str("");
#endif
        *(OFiles+i)<<ostr_tmp.str()<<flush;
        ostr_tmp.str("");
      }
      break;

      case FRateID:
      for(unsigned int j=0;j<OFCtr[i].ParList.size();j++)
      {
        for(unsigned int k=NIdx[OFCtr[i].ParList[j]];k<(NIdx[OFCtr[i].ParList[j]]+ConfSta[OFCtr[i].ParList[j]].NeuNum);k++)
        {
          if((Nets+k)->SpiOut)
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

#ifdef RABBITMQ
            ostr2_tmp<<"{\"type\":\"FR\",\"time\":"<<EvTime*1e-3f<<",\"FiringRate\":"<<1e3f*mm/OFCtr[i].RWin/nn<<"}";
#endif

          }
          else
          {
            for(unsigned int j=0;j<OFCtr[i].ParList.size();j++)
            { ostr_tmp<<" "<<(1e3f * OFCtr[i].SpiRWin[OFCtr[i].ParList.size()*(OFCtr[i].RWin+1)+j]/OFCtr[i].RWin/ConfSta[j].NeuNum ); }
          }

#ifdef RABBITMQ
          s1=ostr2_tmp.str();
          MQSend(s1);
          ostr2_tmp.str("");
#endif

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
      ostr_tmp<<EvTime*1e-3f;

      for(unsigned int j=0;j<OFCtr[i].ParList.size();j++)
      {
        for(unsigned int k=NIdx[OFCtr[i].ParList[j]];k<(NIdx[OFCtr[i].ParList[j]]+ConfSta[OFCtr[i].ParList[j]].NeuNum);k++)
        {
          if((Nets+k)->SynFastSize1!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynFastSize1;m++)
            { ostr_tmp<<" "<<((Nets+k)->SynFst1+m)->stFast; }
          }

          if((Nets+k)->SynSlowSize1!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynSlowSize1;m++)
            {
              ostr_tmp<<" "<<((Nets+k)->SynSlo1+m)->xtSlow;
              ostr_tmp<<" "<<((Nets+k)->SynSlo1+m)->stSlow;
            }
          }

          if((Nets+k)->SynGlutSize1!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynGlutSize1;m++)
            {
              ostr_tmp<<" "<<((Nets+k)->SynGlu1+m)->stFast;
              ostr_tmp<<" "<<((Nets+k)->SynGlu1+m)->xtSlow;
              ostr_tmp<<" "<<((Nets+k)->SynGlu1+m)->stSlow;
            }
          }

          if((Nets+k)->SynFastSize!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynFastSize;m++)
            { ostr_tmp<<" "<<((Nets+k)->SynFst+m)->stFast; }
          }

          if((Nets+k)->SynSlowSize!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynSlowSize;m++)
            {
              ostr_tmp<<" "<<((Nets+k)->SynSlo+m)->xtSlow;
              ostr_tmp<<" "<<((Nets+k)->SynSlo+m)->stSlow;
            }
          }

          if((Nets+k)->SynGlutSize!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynGlutSize;m++)
            {
              ostr_tmp<<" "<<((Nets+k)->SynGlu+m)->stFast;
              ostr_tmp<<" "<<((Nets+k)->SynGlu+m)->xtSlow;
              ostr_tmp<<" "<<((Nets+k)->SynGlu+m)->stSlow;
            }
          }
        }
      }

      *(OFiles+i)<<ostr_tmp.str()<<"\n";
      ostr_tmp.str("");
      break;

      case SynWgtID:
      if(Step % (OFCtr[i].FPnt*OneMsStep)==0) // default is 100 mS , print out to file
      {

      ostr_tmp<<EvTime*1e-3f;
      for(unsigned int j=0;j<OFCtr[i].ParList.size();j++)
      {
        for(unsigned int k=NIdx[OFCtr[i].ParList[j]];k<(NIdx[OFCtr[i].ParList[j]]+ConfSta[OFCtr[i].ParList[j]].NeuNum);k++)
        {
          if((Nets+k)->SynFastSize1!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynFastSize1;m++)
            { ostr_tmp<<" "<<((Nets+k)->SynFst1+m)->gFast; }
          }

          if((Nets+k)->SynSlowSize1!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynSlowSize1;m++)
            { ostr_tmp<<" "<<((Nets+k)->SynSlo1+m)->gSlow; }
          }

          if((Nets+k)->SynGlutSize1!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynGlutSize1;m++)
            {
              ostr_tmp<<" "<<((Nets+k)->SynGlu1+m)->gFast;
              ostr_tmp<<" "<<((Nets+k)->SynGlu1+m)->gSlow;
            }
          }

          if(((Nets+k)->ActCtr & LTP_ACT) != 0)
          {
          if((Nets+k)->SynFastSize!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynFastSize;m++)
            {
              for(unsigned int n=0;n<((Nets+k)->SynFst+m)->PosIDNum;n++)
              { ostr_tmp<<" "<<(((Nets+k)->SynFst+m)->gFast+n);}
            }
          }

          if((Nets+k)->SynSlowSize!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynSlowSize;m++)
            {
              for(unsigned int n=0;n<((Nets+k)->SynSlo+m)->PosIDNum;n++)
              { ostr_tmp<<" "<<(((Nets+k)->SynSlo+m)->gSlow+n);}
            }
          }

          if((Nets+k)->SynGlutSize!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynGlutSize;m++)
            {
              for(unsigned int n=0;n<((Nets+k)->SynGlu+m)->PosIDNum;n++)
                {
              ostr_tmp<<" "<<(((Nets+k)->SynGlu+m)->gFast+n);
              ostr_tmp<<" "<<(((Nets+k)->SynGlu+m)->gSlow+n);
                }

            }
          }
          }
          else
          {
          if((Nets+k)->SynFastSize!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynFastSize;m++)
            { ostr_tmp<<" "<<((Nets+k)->SynFst+m)->gFast; }
          }

          if((Nets+k)->SynSlowSize!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynSlowSize;m++)
            { ostr_tmp<<" "<<((Nets+k)->SynSlo+m)->gSlow; }
          }

          if((Nets+k)->SynGlutSize!=0)
          {
            for(unsigned int m=0;m<(Nets+k)->SynGlutSize;m++)
            {
              ostr_tmp<<" "<<((Nets+k)->SynGlu+m)->gFast;
              ostr_tmp<<" "<<((Nets+k)->SynGlu+m)->gSlow;
            }
          }
          }
        }
      }

      *(OFiles+i)<<ostr_tmp.str()<<"\n";
      ostr_tmp.str("");
      }
      break;

      default:;
    }

  }

  {string ().swap(s1);}

}


