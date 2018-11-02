#include <iostream>
#include <cmath>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <unistd.h>
#include <thread>
#include <mutex>
#include "NetSim.h"

#include <fstream>
#include <string.h>
#include <sys/socket.h> /*用來定址 address family*/
#include <sys/types.h>
#include <sys/stat.h>
#include <netinet/in.h> /*用於in addr_t Address*/
#include <arpa/inet.h>
#include <netdb.h>

using namespace std;

NetSim *GetNewSim(int,float);


int main (int argc, char * const argv[]) {

//input file
  string conf="network.conf";
  string pro="network.pro";

//output file
  string OutMem="MemPot.dat";
  string OutSpi="Spikes.dat";
  string OutFR="FRate.dat";
  bool OutMem_en=false;
  bool OutSpi_en=false;
  bool OutFR_en=false;

  bool STP=false;
  bool STD=false;
  bool LTP=false;
  bool SOD=false;
  bool rseed_en=false;
  unsigned int rseed=0;

//-----------------------------
  NetSim *FlySim;  // network pointer
  string s1;
  int Rp=1;
  int state=STATE_SimName;
  int sop=ROUGH;
  int thsize=1;
  OutCtr NewOutF;
  int nmodel=GNL2;
  float dt=0.1f;

  for(int i=0;i<argc;i++)
  {
    s1=argv[i];
    switch(state)
    {
      case STATE_SimName:
      state=STATE_command;
      break;

      case STATE_command:
      if(s1=="-pro")
        state=STATE_pro;
      else if(s1=="-conf")
        state=STATE_conf;
      else if(s1=="-om")
        state=STATE_OutMem;
      else if(s1=="-os")
        state=STATE_OutSpi;
      else if(s1=="-or")
        state=STATE_OutFR;
      else if(s1=="-rp")
        state=STATE_Repeat;
      else if(s1=="-t")
        state=STATE_thread;
      else if(s1=="-s")
        state=STATE_SolverType;
      else if(s1=="-nmodel")
        state=STATE_NeuModel;
      else if(s1=="-daemon")
        state=STATE_Daemon;
      else if(s1=="-dt")
        state=STATE_DTIME;
      else if(s1=="-udfsed")
        state=STATE_UDFRSEED;
      else if(s1=="-STP")
        STP=true;
      else if(s1=="-STD")
        STD=true;
      else if(s1=="-LTP")
        LTP=true;
      else
      {  // help file
        cout
        <<"neuron and neural network simulator designed by HYC in ccLo Lab @2018\n\n"
        <<"flysim -option parameters\n\n"
        <<"examples:\n"
        <<"-pro network.pro      # read protocal file: default=network.pro\n"
        <<"-conf network.conf    # read configuration file: default=network.conf\n"
        <<"-om my_membrane.dat   # for batch operation: output membrane potential file\n"
        <<"-os my_spike.dat      # for batch operation: output spikes file\n"
        <<"-or my_rate.dat       # for batch operation: output firing rate: default rate window=50ms, print out=100ms\n"
        <<"-rp 1                 # set repeat times: default=1\n"
        <<"-t 4                  # set multithreading: default=1\n"
        <<"-dt 0.1               # time step(default=0.1ms)\n"
        <<"-STP                  # use short term plasticity synapse\n"
        <<"-STD                  # use short term depression synapse and this option is disabled when -STP used\n"
        <<"-LTP                  # use long term plasticity synapse(STDP)\n"
        <<"-nmodel GNL2          # neuron model:\n"
        <<"                      # GNL2(default): classical leaky integrate and fire model\n";

        i=argc;
        state=STATE_help;
      }
      break;

      case STATE_pro:
      pro=s1;
      state=STATE_command;
      break;

      case STATE_conf:
      conf=s1;
      state=STATE_command;
      break;

      case STATE_OutMem:
      OutMem=s1;
      OutMem_en=true;
      state=STATE_command;
      break;

      case STATE_OutSpi:
      OutSpi=s1;
      OutSpi_en=true;
      state=STATE_command;
      break;

      case STATE_OutFR:
      OutFR=s1;
      OutFR_en=true;
      state=STATE_command;
      break;

      case STATE_Repeat:
      Rp=atoi(s1.c_str());
      state=STATE_command;
      break;

      case STATE_thread:
      thsize=atoi(s1.c_str());
      state=STATE_command;
      break;

      case STATE_SolverType:
      if(s1=="accurate")
        sop=ACCURATE;
      else if(s1=="moderate")
        sop=MODERATE;
      else if(s1=="rough")
        sop=ROUGH;
      else
        sop=ROUGH;

      state=STATE_command;
      break;



      case STATE_NeuModel:
      if(s1=="LIF")
        nmodel=LIF;
      else if(s1=="LIF_GP")
        nmodel=LIF_GP4;
      else if(s1=="HH")
        nmodel=HH;
      else if(s1=="sim06")
        nmodel=GNL;
      else if(s1=="sim07")
        nmodel=GNL2;
      else if(s1=="GNL2")
        nmodel=GNL2;
      else
        nmodel=LIF_GP4;

      state=STATE_command;
      break;

      case STATE_DTIME:
      dt=(float)atof(s1.c_str());
      state=STATE_command;
      break;

      case STATE_UDFRSEED:
      rseed=(unsigned int )atoi(s1.c_str());
      rseed_en=true;
      state=STATE_command;
      break;

      default:
      i=argc;
      state=STATE_help;
      break;
    }

  }

//--------------------------------------------------------------------------

  if(state!=STATE_help)
  {
      FlySim=GetNewSim(nmodel,dt); //new a NetSim

      if(rseed_en)
        FlySim->ActCtr |= UDFRDSEED_ACT;
      else
        FlySim->ActCtr &= ~UDFRDSEED_ACT;

      FlySim->RpSize=Rp;
      FlySim->rseed=rseed;

      if(SOD && (nmodel == LIF)) FlySim->ActCtr |= SOD_ACT;
      if(STP) FlySim->ActCtr |= STP_ACT;
      if(STD) FlySim->ActCtr |= STD_ACT;
      if(LTP) FlySim->ActCtr |= LTP_ACT;

      FlySim->CSta(conf,dt);  // read in netwrok configuration
      FlySim->parse(conf,dt);  // build netwrok

      FlySim->CommMem=OutMem_en;
      FlySim->CommSpi=OutSpi_en;
      FlySim->CommFR=OutFR_en;

      FlySim->OutMem=OutMem;
      FlySim->OutSpi=OutSpi;
      FlySim->OutFR=OutFR;

      FlySim->Comm(FlySim->OFCtr); // command line option for network.pro
      FlySim->DefMacro(pro);
      FlySim->parse(pro,FlySim->Evts,FlySim->OFCtr);  // read time events, and output file events
      FlySim->ConfsInfo("ConfsInfo.log");
      FlySim->SetPp(thsize,sop);  // setup population
      FlySim->RunInit(thsize,sop);

      for(FlySim->Repeat=1;FlySim->Repeat<=Rp;FlySim->Repeat++)
      { FlySim->Run(); }
  } //state end

  return 0;
}


NetSim *GetNewSim(int nmodel,float dt)
{
  NetSim *Sim;  // network pointer
  Sim=new NetSimGnl2;
  Sim->ActCtr=AMPA_ACT | GABA_ACT | NMDA_ACT | ACH_ACT | GCL_ACT | GAP_ACT | LIF_ACT;
  Sim->dtime=dt;
  return Sim;
}


