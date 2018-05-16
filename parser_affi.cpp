#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <random>
#include <omp.h>
#include "parser.h"


using namespace std;

TmEvt::TmEvt(){ init(); }

void TmEvt::init()
{
  TimEv=0.0f;
  Type=0;
  //Label="Non";
  Pp=0;
  SubTyp=KW_AMPA; //target receptor

  Var1=0.0f;
  Var2=0.0f;
  NoLmt=false;
}

void TmEvt::TexInfo(string& ostr)
{
  ostringstream ostr_tmp;

  ostr_tmp<<"Event Time="<<TimEv;

  if(Type==EVTEXTFREQ)
    ostr_tmp<<", type:ChangeExtFreq";
  else if(Type==EVTCUNTINJ)
    ostr_tmp<<", type:CurrentInjec";

  ostr_tmp
      <<", Target population:"<<Pp
      <<", Receptor:"<<(unsigned int)SubTyp
      <<", Variable1="<<Var1
      <<", Variable2="<<Var2
      <<"\n";
  ostr=ostr_tmp.str();
}

//-----------------------------------------------------------


OutCtr::OutCtr(){ init();}

void OutCtr::init()
{
  Type="Non";
  FileName="Non";
  TypeID=0;
  RWin=50; // 50ms
  FPnt=100; // 100ms
  TimEv=0.0f;
  //SumPars=false;
  Status=0;
}

OutCtr::~OutCtr()
{
  {string ().swap(Type);}
  {string ().swap(FileName);}
  {vector<unsigned int> ().swap(ParList);}
  {vector<unsigned int> ().swap(SpiRWin);}
}

void OutCtr::TexInfo(string& ostr)
{
  ostringstream ostr_tmp;
  ostr_tmp
      <<"Output File Name:"<<FileName
      <<", Type:"<<Type<<"\n"
      <<"Population:";

  for(unsigned int i=0;i<ParList.size();i++)
  {
    ostr_tmp<<ParList[i]<<", ";
    if(i%10==0) ostr_tmp<<"\n";
  }
  ostr_tmp<<"\n";

  ostr=ostr_tmp.str();
}

