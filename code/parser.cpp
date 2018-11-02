#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include <random>
#include <unistd.h>
#include <limits.h>
#include <algorithm>
#include <ctime>
#include <cfloat>
#include <limits.h>
#include "parser.h"

using namespace std;


Parser::Parser()
{
  CommMem=false;
  CommSpi=false;
  CommFR=false;
  OutMem="Non";
  OutSpi="Non";
  OutFR="Non";
  delimiters = " :;,'=\t\r\n";

  state = STATE_Initial;
  RSeed=time(NULL);
}

Parser::~Parser()
{
  for(unsigned int i=0;i<ConfSta.size();i++)
  { delete [] ConfSta[i].PpRep; }

  {vector<NetConfSta> ().swap(ConfSta);}

  for(unsigned int i=0;i<tokens.size();i++)
  {{string ().swap(tokens[i]);}}

  {vector<string> ().swap(tokens);}

  for(unsigned int i=0;i<GpExt.size();i++)
  {{vector<unsigned int> ().swap(GpExt[i].Pp);}}
  {vector<GroupExtent> ().swap(GpExt);}
}
//-------------------------------------------------network.conf statsitic
void Parser::CSta(const string& CmdFileName,float dt)
{
  string::size_type startpos;
  string::size_type pos;
  string s1,str;

  char ch[65536];
  ifstream fin;

  fin.open(CmdFileName.c_str());
  if(fin)
  {
    fin.get(ch,65536,EOF);
    str = ch;
  }
  else
  { cout<<".conf open file error!"<<endl; }


  startpos = 0;
  pos = str.find_first_of(delimiters, startpos);

  cout<<"Configuration file statsitic"<<endl;
  while( (pos != string::npos) || (startpos != string::npos) )
  {
    {string ().swap(s1);}
    s1 = str.substr(startpos, pos - startpos);
    startpos = str.find_first_not_of(delimiters, pos);
    pos = str.find_first_of(delimiters, startpos);

    switch(state)
    {
      case STATE_Initial:
      if(s1 == "NeuralPopulation")
      {
        NetConfSta NewSta;
        NewSta.NeuNum=0;
        NewSta.PpRepNum=0;
        NewSta.PpLinkNum=0;
        NewSta.PpLinkBeg=0;
        NewSta.SelConEn=false;
        ConfSta.push_back(NewSta);

        state = STATE_NeuPopuName;
      }
      break;

//-----Neural Population

      case STATE_NeuPopuName:
      ConfSta[ConfSta.size()-1].NeuPpName=s1;
      state = STATE_NeuPopu;
      break;

      case STATE_NeuPopu:
      if(s1 == "Receptor")
      {
        ConfSta[ConfSta.size()-1].PpRepNum++;  //Receptor number increase
        state = STATE_Receptor;
      }
      else if(s1 == "TargetPopulation")
      {
        ConfSta[ConfSta.size()-1].PpLinkNum++;  //post synaptic number increase
        state = STATE_NeuCon;
      }
      else if(s1 == "SelfConnection")
      { state = STATE_SelfCon; }
      else if(s1 == "EndNeuralPopulation") // return to STATE_Initial
        state = STATE_Initial;
      break;

//-------------- Receptor end

      case STATE_Receptor:
      if(s1 == "EndReceptor")
        state = STATE_NeuPopu;
      break;
//-------------- Pupulation link
      case STATE_NeuCon:
      if(s1 == "EndTargetPopulation")
        state = STATE_NeuPopu;
      break;

//-------------- Pupulation link
      case STATE_SelfCon:
        if(s1=="true")
          ConfSta[ConfSta.size()-1].SelConEn=true;
        else
          ConfSta[ConfSta.size()-1].SelConEn=false;

        state = STATE_NeuPopu;
      break;


      default:;
    }

    if(fin.eof()){} //end of file, do not thing
    else if(startpos == string::npos)
    {
      fin.get(ch,65536,EOF);

      {string ().swap(str);}
      str = ch;

      startpos = 0;
      pos = str.find_first_of(delimiters, startpos);
    }
    else if(pos == string::npos)
    {
      s1=str.substr(startpos,str.size()-1);
      fin.get(ch,65536,EOF);

      {string ().swap(str);}
      str = s1+ch;

      startpos = 0;
      pos = str.find_first_of(delimiters, startpos);
    }

  }

  {string ().swap(s1);}
  {string ().swap(str);}
  fin.close();  // close file
  fin.clear();  // close file

  for(unsigned int i=0;i<ConfSta.size();i++)
  { ConfSta[i].PpRep = new NeuRep[ConfSta[i].PpRepNum]; }

// initializaiton
  unsigned int Pp=0;
  for(unsigned int i=0;i<ConfSta.size();i++)
  {
    for(unsigned int j=0;j<ConfSta[i].PpRepNum;j++)
    {
      ConfSta[i].PpRep[j].NeuRepName = KW_AMPA;
      ConfSta[i].PpRep[j].Tau = 5.0f;
      ConfSta[i].PpRep[j].TauXT = 2.0f;
      ConfSta[i].PpRep[j].Vrev = 0.0f;
      ConfSta[i].PpRep[j].FExt = 0.0f;
      ConfSta[i].PpRep[j].MnExtEff = 2.1f;
      ConfSta[i].PpRep[j].MnExtCon = 1.0f;
    }

    ConfSta[i].PpLinkBeg = Pp;
    Pp += ConfSta[i].PpLinkNum;
  }
  PostPp = new NeuCon[Pp];


  Confs = new NetConf[ConfSta.size()];

  for(unsigned int i=0;i<ConfSta.size();i++)
  {

// membrance parameter
    Confs[i].Taum=2.0f;
    Confs[i].Vth=-50.0f; //-50mV
    Confs[i].El=VRESTING; //-70mV
    Confs[i].gl=2.5f;
    Confs[i].Cm=1.0f;
    Confs[i].Vreset=VRESET;  // mV
    Confs[i].RfPed=1.8f/dt+dt;  // dt=0.1ms: 1.8ms = 19(total)*dt - 1(firing))*dt
    Confs[i].SpiDy=0.0f;  // spike is delay 18 DTIME = 1.8ms

    for(unsigned int j=ConfSta[i].PpLinkBeg;j<ConfSta[i].PpLinkBeg+ConfSta[i].PpLinkNum;j++)
    {
      
      PostPp[j].TarPpName = 0;//UINT_MAX
      PostPp[j].TarRep = KW_AMPA;//UCHAR_MAX
      PostPp[j].Connectivity =1.0f;
      PostPp[j].MeanEff =0.0f;
      PostPp[j].weight =1.0f;
    }
  }

// short-term plasticity
  ConfStp = new ParStp[ConfSta.size()];
  for(unsigned int i=0;i<ConfSta.size();i++)
  {
    ConfStp[i].pv=0.6f; // reduce factor of depression factor D
    ConfStp[i].tD=300.0f; // decay time constant D

    ConfStp[i].aF=0.0035f; // the alpha of facilitation F
    ConfStp[i].tF=7000.0f; // the decay time constant of F

    ConfStp[i].aCap=250.0f; // the alpha of Ca^2+ peak
    ConfStp[i].tCap=2.0f; // decay time constant  Ca^2+ peak

    ConfStp[i].aCar=0.1f; // alpha of residual Ca^2+
    ConfStp[i].tCar=2000.0f; // decay time constant residual Ca^2+
  }

//long-term plasticity
  ConfLtp = new ParLP[ConfSta.size()];
  for(unsigned int i=0;i<ConfSta.size();i++)
  {
    ConfLtp[i].tLP=50.0f;  // time constant of long term plsiticity factor
    ConfLtp[i].SpiLPDy=0.0f; //back propogation spike delay
    ConfLtp[i].PosISI=1e-5f;
    ConfLtp[i].NegISI=1e-5f;
  }

//HH
  ConfHH = new ParHH2[ConfSta.size()];
  for(unsigned int i=0;i<ConfSta.size();i++)
  {
    ConfHH[i].Ek=0.0f;//-12.0f+VRESTING;
    ConfHH[i].gk=0.0f;//36.0f;
    ConfHH[i].Ena=0.0f;//120.0f+VRESTING;
    ConfHH[i].gna=0.0f;//120.0f;

    ConfHH[i].AnA=0.01f;
    ConfHH[i].AnB=10.0f+VRESTING;
    ConfHH[i].AnC=10.0f;
    ConfHH[i].AnD=-1.0f;

    ConfHH[i].BnA=0.125f;
    ConfHH[i].BnB=0.0f+VRESTING;
    ConfHH[i].BnC=80.0f;
    ConfHH[i].BnD=0.0f;

    ConfHH[i].AmA=0.1f;
    ConfHH[i].AmB=25.0f+VRESTING;
    ConfHH[i].AmC=10.0f;
    ConfHH[i].AmD=-1.0f;

    ConfHH[i].BmA=4.0f;
    ConfHH[i].BmB=0.0f+VRESTING;
    ConfHH[i].BmC=18.0f;
    ConfHH[i].BmD=0.0f;

    ConfHH[i].AhA=0.07f;
    ConfHH[i].AhB=0.0f+VRESTING;
    ConfHH[i].AhC=20.0f;
    ConfHH[i].AhD=0.0f;

    ConfHH[i].BhA=1.0f;
    ConfHH[i].BhB=30.0f+VRESTING;
    ConfHH[i].BhC=10.0f;
    ConfHH[i].BhD=1.0f;
  }
}

//-------------------------------------------------network.conf
void Parser::parse(const string& CmdFileName,float dt)
{
  string::size_type startpos;
  string::size_type pos;
  string s1,str;

  char ch[65536];
  ifstream fin;
  unsigned int KwTmp=KW_AMPA;

  unsigned int Ni=0;
  unsigned int Ri=0;
  unsigned int Li=0;

  float ConLmt=0.0f;

  fin.open(CmdFileName.c_str());
  if(fin)
  {
    fin.get(ch,65536,EOF);
    str = ch;
  }
  else
  { cout<<".conf open file error!"<<endl; }

  startpos = 0;
  pos = str.find_first_of(delimiters, startpos);

  while( (pos != string::npos) || (startpos != string::npos) )
  {
    {string ().swap(s1);}
    s1 = str.substr(startpos, pos - startpos);
    startpos = str.find_first_not_of(delimiters, pos);
    pos = str.find_first_of(delimiters, startpos);

    switch(state)
    {
      case STATE_Initial:
      if(s1 == "NeuralPopulation")
        state = STATE_NeuPopuName;
      break;

//-----Neural Population

      case STATE_NeuPopuName:
      cout<<"Configuration file setup:"<<100.0f*(float)(Ni+1)/ConfSta.size()<<"%        \r"<<flush;
      state = STATE_NeuPopu;
      break;

      case STATE_NeuPopu:
      if(s1 == "N")
        state = STATE_N;
      else if(s1 == "Taum")  //---- leaky integrate and fire model
        state = STATE_Taum;
      else if(s1 == "Threshold")
        state = STATE_Vth;
      else if(s1 == "LeakyConductance")
        state = STATE_gl;
      else if(s1 == "C")
        state = STATE_Cm;
      else if(s1 == "RevPot")  // legacy mode
        state = STATE_Vresting;
      else if(s1 == "RestPot")
        state = STATE_Vresting;
      else if(s1 == "ResetPot")
        state = STATE_Vreset;
      else if(s1 == "RefractoryPeriod")
        state = STATE_RfPed;
      else if(s1 == "SpikeDly")
        state = STATE_SpiDy;

      else if(s1 == "STP_pv") //stp
        state = STATE_STP_pv;
      else if(s1 == "STP_tD")
        state = STATE_STP_tD;
      else if(s1 == "STP_aF")
        state = STATE_STP_aF;
      else if(s1 == "STP_tF")
        state = STATE_STP_tF;
      else if(s1 == "STP_aCap")
        state = STATE_STP_aCap;
      else if(s1 == "STP_tCap")
        state = STATE_STP_tCap;
      else if(s1 == "STP_aCar")
        state = STATE_STP_aCar;
      else if(s1 == "STP_tCar")
        state = STATE_STP_tCar;

      else if(s1 == "LTP_tLP") //ltp
        state = STATE_LTP_tLP;
      else if(s1 == "LTP_SpiLPDy")
        state = STATE_LTP_SpiLPDy;
      else if(s1 == "LTP_PosISI")
        state = STATE_LTP_PosISI;
      else if(s1 == "LTP_NegISI")
        state = STATE_LTP_NegISI;

      else if(s1 == "VthSodium") //sodium channel model
        state = STATE_VthSod;
      else if(s1 == "HH_Ek") //HH model
        state = STATE_HH_Ek;
      else if(s1 == "HH_gk")
        state = STATE_HH_gk;
      else if(s1 == "HH_Ena")
        state = STATE_HH_Ena;
      else if(s1 == "HH_gna")
        state = STATE_HH_gna;

      else if(s1 == "HH_AnA")
        state = STATE_HH_AnA;
      else if(s1 == "HH_AnB")
        state = STATE_HH_AnB;
      else if(s1 == "HH_AnC")
        state = STATE_HH_AnC;
      else if(s1 == "HH_AnD")
        state = STATE_HH_AnD;
      else if(s1 == "HH_BnA")
        state = STATE_HH_BnA;
      else if(s1 == "HH_BnB")
        state = STATE_HH_BnB;
      else if(s1 == "HH_BnC")
        state = STATE_HH_BnC;
      else if(s1 == "HH_BnD")
        state = STATE_HH_BnD;

      else if(s1 == "HH_AmA")
        state = STATE_HH_AmA;
      else if(s1 == "HH_AmB")
        state = STATE_HH_AmB;
      else if(s1 == "HH_AmC")
        state = STATE_HH_AmC;
      else if(s1 == "HH_AmD")
        state = STATE_HH_AmD;
      else if(s1 == "HH_BmA")
        state = STATE_HH_BmA;
      else if(s1 == "HH_BmB")
        state = STATE_HH_BmB;
      else if(s1 == "HH_BmC")
        state = STATE_HH_BmC;
      else if(s1 == "HH_BmD")
        state = STATE_HH_BmD;

      else if(s1 == "HH_AhA")
        state = STATE_HH_AhA;
      else if(s1 == "HH_AhB")
        state = STATE_HH_AhB;
      else if(s1 == "HH_AhC")
        state = STATE_HH_AhC;
      else if(s1 == "HH_AhD")
        state = STATE_HH_AhD;
      else if(s1 == "HH_BhA")
        state = STATE_HH_BhA;
      else if(s1 == "HH_BhB")
        state = STATE_HH_BhB;
      else if(s1 == "HH_BhC")
        state = STATE_HH_BhC;
      else if(s1 == "HH_BhD")
        state = STATE_HH_BhD;

      else if(s1 == "Receptor")
        state = STATE_RepName;
      else if(s1 == "TargetPopulation")
        state = STATE_NeuConName;
      else if(s1 == "EndNeuralPopulation") // return to STATE_Initial
      {
        Ni++;
        Ri=0;
        Li=0;
        state = STATE_Initial;
      }
      break;

//--------------neuron basic parameter

      case STATE_N:
      ConfSta[Ni].NeuNum=stoi(s1);
      state = STATE_NeuPopu;
      break;

//--------- leaky integrate and fire
      case STATE_Taum:
      Confs[Ni].Taum=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_Vth:
      Confs[Ni].Vth=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_gl:
      Confs[Ni].gl=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_Cm:
      Confs[Ni].Cm=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_Vreset:
      Confs[Ni].Vreset=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_Vresting:
      Confs[Ni].El=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_RfPed:
      Confs[Ni].RfPed=floor(stof(s1)/dt)+dt;
      state = STATE_NeuPopu;
      break;

      case STATE_SpiDy:
      Confs[Ni].SpiDy=floor(stof(s1)/dt);
      state = STATE_NeuPopu;
      break;

//------------stp
      case STATE_STP_pv:
      ConfStp[Ni].pv=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_STP_tD:
      ConfStp[Ni].tD=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_STP_aF:
      ConfStp[Ni].aF=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_STP_tF:
      ConfStp[Ni].tF=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_STP_aCap:
      ConfStp[Ni].aCap=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_STP_tCap:
      ConfStp[Ni].tCap=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_STP_aCar:
      ConfStp[Ni].aCar=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_STP_tCar:
      ConfStp[Ni].tCar=stof(s1);
      state = STATE_NeuPopu;
      break;

//------------ltp
      case STATE_LTP_tLP:
      ConfLtp[Ni].tLP=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_LTP_SpiLPDy:
      ConfLtp[Ni].SpiLPDy=floor(stof(s1)/dt);
      state = STATE_NeuPopu;
      break;

      case STATE_LTP_PosISI:
      ConfLtp[Ni].PosISI=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_LTP_NegISI:
      ConfLtp[Ni].NegISI=stof(s1);
      state = STATE_NeuPopu;
      break;



//------------sodium channel
      case STATE_VthSod:
      ConfHH[Ni].AmB=25.0f+stof(s1);
      ConfHH[Ni].BmB=stof(s1);
      state = STATE_NeuPopu;
      break;

//------------HH
      case STATE_HH_Ek:
      ConfHH[Ni].Ek=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_gk:
      ConfHH[Ni].gk=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_Ena:
      ConfHH[Ni].Ena=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_gna:
      ConfHH[Ni].gna=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_AnA:
      ConfHH[Ni].AnA=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_AnB:
      ConfHH[Ni].AnB=10.0f+stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_AnC:
      ConfHH[Ni].AnC=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_AnD:
      ConfHH[Ni].AnD=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_AmA:
      ConfHH[Ni].AmA=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_AmB:
      ConfHH[Ni].AmB=25.0f+stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_AmC:
      ConfHH[Ni].AmC=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_AmD:
      ConfHH[Ni].AmD=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_AhA:
      ConfHH[Ni].AhA=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_AhB:
      ConfHH[Ni].AhB=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_AhC:
      ConfHH[Ni].AhC=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_AhD:
      ConfHH[Ni].AhD=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_BnA:
      ConfHH[Ni].BnA=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_BnB:
      ConfHH[Ni].BnB=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_BnC:
      ConfHH[Ni].BnC=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_BnD:
      ConfHH[Ni].BnD=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_BmA:
      ConfHH[Ni].BmA=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_BmB:
      ConfHH[Ni].BmB=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_BmC:
      ConfHH[Ni].BmC=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_BmD:
      ConfHH[Ni].BmD=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_BhA:
      ConfHH[Ni].BhA=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_BhB:
      ConfHH[Ni].BhB=30.0f+stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_BhC:
      ConfHH[Ni].BhC=stof(s1);
      state = STATE_NeuPopu;
      break;

      case STATE_HH_BhD:
      ConfHH[Ni].BhD=stof(s1);
      state = STATE_NeuPopu;
      break;

//-------------- Receptor name

      case STATE_RepName:
      if(s1=="AMPA")
        KwTmp=KW_AMPA;
      else if(s1=="GABA")
        KwTmp=KW_GABA;
      else if(s1=="NMDA")
        KwTmp=KW_NMDA;
      else if(s1=="GAP")
        KwTmp=KW_GAP;
      else if(s1=="Sine")
        KwTmp=KW_SINE;
      else if(s1=="Ach")
        KwTmp=KW_ACH;
      else if(s1=="GluCl")
        KwTmp=KW_GCL;

      ConfSta[Ni].PpRep[Ri].NeuRepName=KwTmp;
      state = STATE_Receptor;
      break;
//-------------- Receptor parameter

      case STATE_Receptor:
      if(s1=="Tau")
        state = STATE_RepTauST;
      else if(s1=="TauXT")
        state = STATE_RepTauXT;
      else if(s1=="RevPot")
        state = STATE_RepVrv;
      else if(s1=="FreqExt")
        state = STATE_RepFExt;
      else if(s1=="MeanExtEff")
        state = STATE_RepMnExtEff;
      else if(s1=="MeanExtCon")
        state = STATE_RepMnExtCon;
      else if(s1 == "EndReceptor")
      {
        Ri++;
        state = STATE_NeuPopu;
      }

      break;

      case STATE_RepTauST:
      ConfSta[Ni].PpRep[Ri].Tau=stof(s1);
      state = STATE_Receptor;
      break;

      case STATE_RepTauXT:
      ConfSta[Ni].PpRep[Ri].TauXT=stof(s1);
      state = STATE_Receptor;
      break;

      case STATE_RepVrv:
      ConfSta[Ni].PpRep[Ri].Vrev=stof(s1);
      state = STATE_Receptor;
      break;

      case STATE_RepFExt:
      ConfSta[Ni].PpRep[Ri].FExt=stof(s1);
      state = STATE_Receptor;
      break;

      case STATE_RepMnExtEff:
      ConfSta[Ni].PpRep[Ri].MnExtEff=stof(s1);
      state = STATE_Receptor;
      break;

      case STATE_RepMnExtCon:
      ConfSta[Ni].PpRep[Ri].MnExtCon=stof(s1);
      state = STATE_Receptor;
      break;

//-------------- Pupulation link name

      case STATE_NeuConName:
      for(unsigned int i=0;i<ConfSta.size();i++)
      {
        if(ConfSta[i].NeuPpName==s1)
        {
          //Confs[Ni].PpLink[Li].TarPpName=i;
          PostPp[ConfSta[Ni].PpLinkBeg+Li].TarPpName=i;
          break;
        }
      }
      state = STATE_NeuCon;
      break;
//-------------- Pupulation link
      case STATE_NeuCon:
      if(s1=="Connectivity")
        state = STATE_NeuCon_Conct;
      else if(s1=="TargetReceptor")
        state = STATE_NeuCon_TarRep;
      else if(s1=="MeanEff")
        state = STATE_NeuCon_MnEff;
      else if(s1=="weight")
        state = STATE_NeuCon_weight;
      else if(s1 == "EndTargetPopulation")
      {
        Li++;
        state = STATE_NeuPopu;
      }
      break;

      case STATE_NeuCon_Conct:
      ConLmt=stof(s1);
      if(ConLmt > 1.0f){ ConLmt = 1.0f; }

      PostPp[ConfSta[Ni].PpLinkBeg+Li].Connectivity=ConLmt;
      state = STATE_NeuCon;
      break;

      case STATE_NeuCon_TarRep:

      if(s1=="AMPA")
        KwTmp=KW_AMPA;
      else if(s1=="GABA")
        KwTmp=KW_GABA;
      else if(s1=="NMDA")
        KwTmp=KW_NMDA;
      else if(s1=="GAP")
        KwTmp=KW_GAP;
      else if(s1=="Sine")
        KwTmp=KW_SINE;
      else if(s1=="Ach")
        KwTmp=KW_ACH;
      else if(s1=="GluCl")
        KwTmp=KW_GCL;

      PostPp[ConfSta[Ni].PpLinkBeg+Li].TarRep=KwTmp;
      state = STATE_NeuCon;
      break;

      case STATE_NeuCon_MnEff:
      PostPp[ConfSta[Ni].PpLinkBeg+Li].MeanEff=stof(s1);
      state = STATE_NeuCon;
      break;

      case STATE_NeuCon_weight:
      PostPp[ConfSta[Ni].PpLinkBeg+Li].weight=stof(s1);
      state = STATE_NeuCon;
      break;

      default:;
    }

    if(fin.eof()){} //end of file, do not thing
    else if(startpos == string::npos)
    {
      fin.get(ch,65536,EOF);

      {string ().swap(str);}
      str = ch;

      startpos = 0;
      pos = str.find_first_of(delimiters, startpos);
    }
    else if(pos == string::npos)
    {
      s1=str.substr(startpos,str.size()-1);
      fin.get(ch,65536,EOF);

      {string ().swap(str);}
      str = s1+ch;

      startpos = 0;
      pos = str.find_first_of(delimiters, startpos);
    }

  }

  cout<<endl;
  {string ().swap(s1);}
  {string ().swap(str);}

  fin.close();  // close file
  fin.clear();  // close file
}


//------------------- define Macro in network.pro
void Parser::DefMacro(const string& CmdFileName)
{
  GroupExtent NewGpExt;

  string::size_type startpos;
  string::size_type pos;
  string s1;

  char ch[65536];
  ifstream fin;

  fin.open(CmdFileName.c_str());
  if(fin)
  {
    fin.get(ch,65536,EOF);
    str = ch;
  }
  else
  { cout<<".pro open file error!"<<endl; }

  startpos = 0;
  pos = str.find_first_of(delimiters, startpos);

  cout<<"protocal Macros setup"<<endl;
  while( (pos != string::npos) || (startpos != string::npos) )
  {
    {string ().swap(s1);}
    s1 = str.substr(startpos, pos - startpos);
    startpos = str.find_first_not_of(delimiters, pos);
    pos = str.find_first_of(delimiters, startpos);

    switch(state)
    {
      case STATE_Initial:
      if(s1=="DefineMacro")
        state = STATE_DefMac;
      break;

//Define Macro-----------------------
      case STATE_DefMac:
      if(s1=="GroupName")
        state = STATE_MacGpName;
      else if(s1=="GroupMembers")
        state = STATE_MacGpMbs;
      else if(s1=="RandomSeeds")
        state = STATE_MacRndSeds;
      else if(s1=="EndDefineMacro")
        state = STATE_Initial;

      break;

      case STATE_MacGpName:
      NewGpExt.GpName=s1;
      state = STATE_DefMac;
      break;

      case STATE_MacRndSeds:
      if(s1=="EndRandomSeeds")
        state = STATE_DefMac;
      else
      {
        if( (unsigned int)stoi(s1)  <= RAND_MAX)  MacRndSeds.push_back((unsigned int)stoi(s1));
        state = STATE_MacRndSeds;
      }
      break;

      case STATE_MacGpMbs:
      if(s1=="EndGroupMembers")
      {
        GpExt.push_back(NewGpExt);
        NewGpExt.Pp.clear();
        state = STATE_DefMac;
      }
      else
      {
        for(unsigned int i=0;i<ConfSta.size();i++)
        {
          if(ConfSta[i].NeuPpName==s1)
          {
            NewGpExt.Pp.push_back(i);
            break;
          }
        }
        state = STATE_MacGpMbs;
      }
      break;

      default:;
    }

    if(fin.eof()){} //end of file, do not thing
    else if(startpos == string::npos)
    {
      fin.get(ch,65536,EOF);

      {string ().swap(str);}
      str = ch;

      startpos = 0;
      pos = str.find_first_of(delimiters, startpos);
    }
    else if(pos == string::npos)
    {
      s1=str.substr(startpos,str.size()-1);
      fin.get(ch,65536,EOF);

      {string ().swap(str);}
      str = s1+ch;

      startpos = 0;
      pos = str.find_first_of(delimiters, startpos);
    }
  }

  {string ().swap(s1);}
  {string ().swap(str);}
  fin.close();  // close file
  fin.clear();  // close file
}

//---------------------------------------------------network.pro
void Parser::parse(const string& CmdFileName,vector<TmEvt>& Evts,vector<OutCtr>& OutFilesCtr)
{
  TmEvt NewEvt;
  vector<TmEvt> UnSolEvts;

  OutCtr NewOutF;

  string::size_type startpos;
  string::size_type pos;
  string s1;

  char ch[65536];
  ifstream fin;
  unsigned int KwTmp=KW_AMPA;

  fin.open(CmdFileName.c_str());
  if(fin)
  {
    fin.get(ch,65536,EOF);
    str = ch;
  }
  else
  { cout<<".pro open file error!"<<endl; }

  startpos = 0;
  pos = str.find_first_of(delimiters, startpos);

  cout<<"protocal file setup"<<endl;
  while( (pos != string::npos) || (startpos != string::npos) )
  {
    {string ().swap(s1);}
    s1 = str.substr(startpos, pos - startpos);
    startpos = str.find_first_not_of(delimiters, pos);
    pos = str.find_first_of(delimiters, startpos);

    switch(state)
    {
      case STATE_Initial:
      if(s1 == "EventTime")
        state = STATE_TimEv;
      else if(s1=="OutControl")
        state = STATE_OutCtr;
      break;

      case STATE_TimEv:
      if(s1=="NoLimit")
      {
        NewEvt.NoLmt = true;
        NewEvt.TimEv=FLT_MAX;
      }
      else
      {
        NewEvt.NoLmt = false;
        NewEvt.TimEv=stof(s1)*1e-3f;
      }

      state = STATE_TimEvDef;
      break;

      case STATE_TimEvDef:
      if(s1=="Type")
        state = STATE_TimTyp;
      else if(s1=="Label" || s1=="InFile")
        state = STATE_TimLbl;
      else if(s1=="Population")
        state = STATE_TimPp;
      else if(s1=="Receptor" || s1=="SubType")
        state = STATE_TimSubTyp;
      else if(s1=="FreqExt" || s1=="GaussMean")
        state = STATE_TimVar1;
      else if(s1=="GaussSTD")
        state = STATE_TimVar2;
      else if(s1=="EndEvent") // write back to event sequence
      {
        UnSolEvts.push_back(NewEvt);
        NewEvt.init();
        state = STATE_Initial;
      }
      break;

      case STATE_TimTyp:
      if(s1=="ChangeExtFreq")
        NewEvt.Type=EVTEXTFREQ;
      else if((s1=="ChangeMembraneNoise") || (s1=="CurrentInjec"))
        NewEvt.Type=EVTCUNTINJ;

      state = STATE_TimEvDef;
      break;

      case STATE_TimLbl:
      //NewEvt.Label=s1;
      state = STATE_TimEvDef;
      break;

      case STATE_TimPp:
      // single population
      for(unsigned int i=0;i<ConfSta.size();i++)
      {
        if(ConfSta[i].NeuPpName==s1)
        {
          NewEvt.Pp=i;
          break;
        }
      }

      // Group extention
      for(unsigned int i=0;i<GpExt.size();i++)
      {
        if(GpExt[i].GpName==s1)
        {
          NewEvt.Pp=i+ConfSta.size();
          break;
        }
      }
      if(s1=="AllPopulation")
        NewEvt.Pp=UINT_MAX;

      state = STATE_TimEvDef;
      break;

      case STATE_TimSubTyp:
      if(s1=="AMPA")
        KwTmp=KW_AMPA;
      else if(s1=="GABA")
        KwTmp=KW_GABA;
      else if(s1=="NMDA")
        KwTmp=KW_NMDA;
      else if(s1=="Ach")
        KwTmp=KW_ACH;
      else if(s1=="GluCl")
        KwTmp=KW_GCL;
      else if(s1=="GAP")
        KwTmp=KW_GAP;
      else if(s1=="Sine")
        KwTmp=KW_SINE;
      else if(s1=="ResumeSynCdt")
        KwTmp=RSM_SYNCOND;
      else if(s1=="ResumeGatingVariable")
        KwTmp=RSM_GATVAR;
      else if(s1=="ResumeMemPot")
        KwTmp=RSM_MEMPOT;
      else if(s1=="ResumeRndSeed")
        KwTmp=RSM_RNDSED;
      else if(s1=="ResumeAll")
        KwTmp=RSM_ALL;

      NewEvt.SubTyp=KwTmp;
      state = STATE_TimEvDef;
      break;

      case STATE_TimVar1:
      NewEvt.Var1=stof(s1);
      state = STATE_TimEvDef;
      break;

      case STATE_TimVar2:
      NewEvt.Var2=stof(s1);
      state = STATE_TimEvDef;
      break;

//output file control-----------------------
      case STATE_OutCtr:
      if(s1=="FileName")
        state = STATE_OutName;
      else if(s1=="Type")
        state = STATE_OutType;
      else if(s1=="FiringRateWindow")
        state = STATE_OutRWin;
      else if(s1=="PrintStep")
        state = STATE_OutFPnt;
      else if(s1=="SumAllPopulation")
        state = STATE_OutSumAllPp;
      else if(s1=="Align")
        state = STATE_OutAlign;
      else if(s1=="EventTime")
        state = STATE_OutTimEv;
      else if(s1=="population")
        state = STATE_OutPopu;
      else if(s1=="EndOutputFile") // write back to event sequence
      {
        if(NewOutF.Type != "Non")
        { OutFilesCtr.push_back(NewOutF); }

        NewOutF.init();
        NewOutF.ParList.clear();
        state = STATE_OutCtr;
      }
      else if(s1=="EndOutControl")
        state = STATE_Initial;
      break;

      case STATE_OutName:
      NewOutF.FileName=s1;
      state = STATE_OutCtr;
      break;

      case STATE_OutType:
      if((s1=="MemPot")&&(CommMem==false))
      {
        NewOutF.Type=s1;
        NewOutF.TypeID = MemPotID;
      }
      else if((s1=="Spike")&&(CommSpi==false))
      {
        NewOutF.Type=s1;
        NewOutF.TypeID = SpikeID;
      }
      else if((s1=="FiringRate")&&(CommFR==false))
      {
        NewOutF.Type=s1;
        NewOutF.TypeID = FRateID;
      }
      else if(s1=="GatingVariable")
      {
        NewOutF.Type=s1;
        NewOutF.TypeID = GatVarID;
      }
      else if(s1=="IsynSeparatedSum")
      {
        NewOutF.Type=s1;
        NewOutF.TypeID = IsynSepSumID;
      }
      else if(s1=="SaveAllStates")
      {
        NewOutF.Type=s1;
        NewOutF.TypeID = SaveID;
      }
      else if(s1=="SynapticWeight")
      {
        NewOutF.Type=s1;
        NewOutF.TypeID = SynWgtID;
      }
      else if(s1=="DeltaMemPotLimit")
      {
        NewOutF.Type=s1;
        NewOutF.TypeID = dVOvfID;
      }

      state = STATE_OutCtr;
      break;

      case STATE_OutRWin:
      NewOutF.RWin=ceil(stof(s1));
      state = STATE_OutCtr;
      break;

      case STATE_OutFPnt:
      NewOutF.FPnt=ceil(stof(s1));
      state = STATE_OutCtr;
      break;

      case STATE_OutSumAllPp:
      if(s1=="true")
        NewOutF.Status |= SUMPARSALL_OF;
      else
        NewOutF.Status &= ~SUMPARSALL_OF;

      state = STATE_OutCtr;
      break;
        
      case STATE_OutAlign:
      if(s1=="true")
        NewOutF.Status |= ALIGN_OF;
      else
        NewOutF.Status &= ~ALIGN_OF;

      state = STATE_OutCtr;
      break;

      case STATE_OutTimEv:
      NewOutF.TimEv=stof(s1)*1e-3f;
      state = STATE_OutCtr;
      break;


      case STATE_OutPopu:
      if(s1=="AllPopulation")
      {
        for(unsigned int i=0;i<ConfSta.size();i++)
        { NewOutF.ParList.push_back(i); }
      }
      else
      {
      // single population
        for(unsigned int i=0;i<ConfSta.size();i++)
        {
          if( ConfSta[i].NeuPpName==s1)
          {
            NewOutF.ParList.push_back(i);
            break;
          }
        }

      // Group extention
        for(unsigned int i=0;i<GpExt.size();i++)
        {
          if(GpExt[i].GpName==s1)
          {
            NewOutF.ParList.insert(NewOutF.ParList.end(),GpExt[i].Pp.begin(),GpExt[i].Pp.end());
            break;
          }
        }

      }
      state = STATE_OutCtr;
      break;

      default:;
    }

    if(fin.eof()){} //end of file, do not thing
    else if(startpos == string::npos)
    {
      fin.get(ch,65536,EOF);

      {string ().swap(str);}
      str = ch;

      startpos = 0;
      pos = str.find_first_of(delimiters, startpos);
    }
    else if(pos == string::npos)
    {
      s1=str.substr(startpos,str.size()-1);
      fin.get(ch,65536,EOF);

      {string ().swap(str);}
      str = s1+ch;

      startpos = 0;
      pos = str.find_first_of(delimiters, startpos);
    }
  }

//time events post process of Macro extension----------------

  for(unsigned int i=0,EvSz=0;i<UnSolEvts.size();i++)
  {
    if(UnSolEvts[i].Pp<ConfSta.size())
      Evts.push_back(UnSolEvts[i]);
    else if(UnSolEvts[i].Pp==UINT_MAX)
    {
      EvSz=Evts.size();
      Evts.insert(Evts.end(),ConfSta.size(),UnSolEvts[i]);
      for(unsigned int m=0;m<ConfSta.size();m++)
      { Evts[EvSz+m].Pp=m; }
    }
    else  //group extension
    {
      for(unsigned int j=0;j<GpExt.size();j++)
      {
        if(UnSolEvts[i].Pp==(ConfSta.size()+j)) // UnSolEvts[i] indicate to GpExt[j]
        {
          EvSz=Evts.size();
          Evts.insert(Evts.end(),GpExt[j].Pp.size(),UnSolEvts[i]);
          for(unsigned int m=0;m<GpExt[j].Pp.size();m++)
          { Evts[EvSz+m].Pp=GpExt[j].Pp[m]; }
        }
      }
    }
  }

  //sort time events for non-reduntent DoEvt()
  sort(Evts.begin(),Evts.end(),[](TmEvt a,TmEvt b){ return a.TimEv < b.TimEv; });

  for(unsigned int i=0;i<GpExt.size();i++)
  {{vector<unsigned int> ().swap(GpExt[i].Pp);}}
  {vector<GroupExtent> ().swap(GpExt);}
  {vector<TmEvt> ().swap(UnSolEvts);}

  {string ().swap(s1);}
  {string ().swap(str);}
  fin.close();  // close file
  fin.clear();  // close file

}

//-----------------------------------------------------------------------
void Parser::Tokenize()
{
  string::size_type startpos = 0;
  string::size_type pos = str.find_first_of(delimiters, startpos);

  while ( (pos != string::npos) || (startpos != string::npos) )
  {
    tokens.push_back( str.substr(startpos, pos - startpos) );
    startpos = str.find_first_not_of(delimiters, pos);
    pos = str.find_first_of(delimiters, startpos);
  }
}

//-----------------------------------------------------------------------
void Parser::Comm(vector<OutCtr>& OutFilesCtr)
{
  if(CommMem)// membrane potential
  {
    OutCtr NewOutF;
    NewOutF.Type="MemPot";
    NewOutF.TypeID=MemPotID;
    NewOutF.FileName=OutMem;

    for(unsigned int i=0;i<ConfSta.size();i++)
    { NewOutF.ParList.push_back(i); }

    OutFilesCtr.push_back(NewOutF);
  }

  if(CommSpi)// spikes
  {
    OutCtr NewOutF;
    NewOutF.Type="Spike";
    NewOutF.TypeID=SpikeID;
    NewOutF.FileName=OutSpi;

    for(unsigned int i=0;i<ConfSta.size();i++)
    { NewOutF.ParList.push_back(i); }

    OutFilesCtr.push_back(NewOutF);
  }

  if(CommFR) // firing rate
  {
    OutCtr NewOutF;
    NewOutF.Type="FiringRate";
    NewOutF.TypeID=FRateID;
    NewOutF.FileName=OutFR;

    for(unsigned int i=0;i<ConfSta.size();i++)
    { NewOutF.ParList.push_back(i); }

    OutFilesCtr.push_back(NewOutF);
  }
}

//-----------------------------------------------------------------------
unsigned short Parser::GetAct(unsigned char f)
{
  if(f==KW_AMPA)
    return AMPA_ACT;
  else if(f==KW_GABA)
    return GABA_ACT;
  else if(f==KW_NMDA)
    return NMDA_ACT;
  else if(f==KW_GAP)
    return GAP_ACT;
  else if(f==KW_CUBA)
    return CUBA_ACT;
  else if(f==KW_ACH)
    return ACH_ACT;
  else if(f==KW_GCL)
    return GCL_ACT;
  else
    return 0;
}

//-----------------------------------------------------------------------
void Parser::ConfsInfo(const string& CmdFileName)
{
  ofstream ostr_tmp;
  ostr_tmp.open(CmdFileName);

  cout<<"write configuration information to ConfsInfo.log, please wait!\n";
  ostr_tmp<<"configuration information:----------------------------------------\n";
  for(unsigned int i=0;i<ConfSta.size();i++)
  {
    ostr_tmp
    <<"population:"<<i<<"\n"
    <<"NeuPopuName:"<<ConfSta[i].NeuPpName
    <<", N="<<ConfSta[i].NeuNum
    <<" Taum="<<Confs[i].Taum

    <<" Vth="<<Confs[i].Vth
    <<" El="<<Confs[i].El
    <<" g="<<Confs[i].gl
    <<" Cm="<<Confs[i].Cm
    <<" Vreset="<<Confs[i].Vreset
    <<"\n";
//HH
    ostr_tmp
    <<"Ek="<<ConfHH[i].Ek<<" gk="<<ConfHH[i].gk<<" Ena="<<ConfHH[i].Ena<<" gna="<<ConfHH[i].gna<<"\n"
    <<"AnA="<<ConfHH[i].AnA<<" AnB="<<ConfHH[i].AnB<<" AnC="<<ConfHH[i].AnC<<" AnD="<<ConfHH[i].AnD<<"\n"
    <<"BnA="<<ConfHH[i].BnA<<" BnB="<<ConfHH[i].BnB<<" BnC="<<ConfHH[i].BnC<<" BnD="<<ConfHH[i].BnD<<"\n"
    <<"AmA="<<ConfHH[i].AmA<<" AmB="<<ConfHH[i].AmB<<" AmC="<<ConfHH[i].AmC<<" AmD="<<ConfHH[i].AmD<<"\n"
    <<"BmA="<<ConfHH[i].BmA<<" BmB="<<ConfHH[i].BmB<<" BmC="<<ConfHH[i].BmC<<" BmD="<<ConfHH[i].BmD<<"\n"
    <<"AhA="<<ConfHH[i].AhA<<" AhB="<<ConfHH[i].AhB<<" AhC="<<ConfHH[i].AhC<<" AhD="<<ConfHH[i].AhD<<"\n"
    <<"BhA="<<ConfHH[i].BhA<<" BhB="<<ConfHH[i].BhB<<" BhC="<<ConfHH[i].BhC<<" BhD="<<ConfHH[i].BhD<<"\n";

//STP
    ostr_tmp
    <<"pv="<<ConfStp[i].pv
    <<" tD="<<ConfStp[i].tD
    <<" aF="<<ConfStp[i].aF
    <<" tF="<<ConfStp[i].tF
    <<" aCap="<<ConfStp[i].aCap
    <<" tCap="<<ConfStp[i].tCap
    <<" aCar="<<ConfStp[i].aCar
    <<" tCar="<<ConfStp[i].tCar
    <<"\n";

//LTP
    ostr_tmp
    <<"tLP="<<ConfLtp[i].tLP
    <<" SpiLPDy="<<ConfLtp[i].SpiLPDy
    <<" PosISI="<<ConfLtp[i].PosISI
    <<" NegISI="<<ConfLtp[i].NegISI
    <<"\n";


    for(unsigned int k=0;k<ConfSta[i].PpRepNum;k++)
    {

      ostr_tmp<<"NeuRepName:";
      switch(PreSynR->NeuRepName)
      {
        case KW_AMPA:
        ostr_tmp<<"AMPA";
        break;

        case KW_GABA:
        ostr_tmp<<"GABA";
        break;

        case KW_NMDA:
        ostr_tmp<<"NMDA";
        break;

        case KW_ACH:
        ostr_tmp<<"Ach";
        break;

        case KW_GCL:
        ostr_tmp<<"GluCl";
        break;

        default:;
      }

      ostr_tmp<<", Tau="<<PreSynR->Tau;

      if(PreSynR->NeuRepName==KW_NMDA)
        ostr_tmp<<", TauXT="<<PreSynR->TauXT;

      ostr_tmp
        <<" Vrv="<<PreSynR->Vrev
        <<" FExt="<<PreSynR->FExt
        <<" MnExtEff="<<PreSynR->MnExtEff
        <<" MnExtCon="<<PreSynR->MnExtCon
        <<"\n";
    }
    ostr_tmp<<"\n";

    for(unsigned int j=ConfSta[i].PpLinkBeg;j<ConfSta[i].PpLinkBeg+ConfSta[i].PpLinkNum;j++)
    {
      ostr_tmp
        <<"TarPopuName:"<<ConfSta[POSYNPP].NeuPpName
        <<" TarReceptor:";

      switch(POSYNREP)
      {
        case KW_AMPA:
        ostr_tmp<<"AMPA";
        break;

        case KW_GABA:
        ostr_tmp<<"GABA";
        break;

        case KW_NMDA:
        ostr_tmp<<"NMDA";
        break;

        case KW_ACH:
        ostr_tmp<<"Ach";
        break;

        case KW_GCL:
        ostr_tmp<<"GluCl";
        break;

        case KW_GAP:
        ostr_tmp<<"GAP";
        break;

        default:;
      }

      ostr_tmp
        <<", Connectivity="<<POSYNCON
        <<" MeanEff="<<POSYNEFF
        <<" weight="<<POSYNWGT
        <<"\n";
    }

    ostr_tmp<<"------------------------------------------------------------------"<<endl;

  }

  ostr_tmp.close();

}

