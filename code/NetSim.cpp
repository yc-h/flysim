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
NetSim::NetSim()
{
  EvTime=0.0f;
  Repeat=1;
  iEvt=0;
  dtime=0.1f;
  rseed=0;
}

NetSim::~NetSim()
{
  {vector<OutCtr> ().swap(OFCtr);}
  {vector<TmEvt> ().swap(Evts);}
  delete [] OFiles;
}

void NetSim::SetOut()
{
  OFiles = new ofstream[OFCtr.size()];

  for(unsigned int i=0;i<OFCtr.size();i++)
  {
// firing rate---------------------------------------------
    if(OFCtr[i].TypeID==FRateID)
    { OFCtr[i].SpiRWin.insert(OFCtr[i].SpiRWin.begin(),(OFCtr[i].ParList.size()*(OFCtr[i].RWin+2)),0); }
  }
}


void NetSim::FileOpen()
{

  if(Repeat>1)
  {
    string s1;
    for(unsigned int i=0;i<OFCtr.size();i++)
    {
      s1=OFCtr[i].FileName+"_"+to_string((long long int)Repeat);
      (OFiles+i)->open(s1.c_str());
    }
  }
  else
  {
    for(unsigned int i=0;i<OFCtr.size();i++)
    { (OFiles+i)->open(OFCtr[i].FileName.c_str()); }
  }

  for(unsigned int i=0;i<OFCtr.size();i++)
  { *(OFiles+i)<<noshowpoint<<setprecision(7); }

}

void NetSim::FileClose()
{
  for(unsigned int i=0;i<OFCtr.size();i++)
  {
    (OFiles+i)->close();
    (OFiles+i)->clear();
  }
}

