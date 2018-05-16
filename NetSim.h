#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <omp.h>
#include <algorithm>
#include <mutex>
#include <string.h>
#include <atomic>
#include <thread>
#include <chrono>
#include "parser.h"

#ifdef RABBITMQ
#include <amqp_tcp_socket.h>
#include <amqp.h>
#include <amqp_framing.h>
#endif


using namespace std;
//extern random_device seed_gen;
typedef std::chrono::high_resolution_clock Clock;

static mutex mtx;
static condition_variable cv;
static atomic<unsigned int> ftokens;
static atomic<unsigned int> Acnt;
static atomic<unsigned int> fcnt;
static atomic<unsigned int> Scnt;

static auto t1=Clock::now();
static auto t2=Clock::now();


#ifdef RABBITMQ
class rbmq
{
  public:
  rbmq()
  {
    chbuff=new char[128];
    strcpy(chbuff,"");
    MQListop=true;
  };
  ~rbmq(){delete chbuff;};

  void MQini();
  amqp_socket_t *MQsocket = NULL;
  amqp_connection_state_t MQconn;
  char const *MQhost="localhost";
  int MQport;
  amqp_basic_properties_t props;

  void CMQini();
  amqp_socket_t *CMQsocket = NULL;
  amqp_connection_state_t CMQconn;
  char const *CMQhost="localhost";
  int CMQport;

  char const *EXCHANGE_NAME="flysimMQ";
  char const *EXCHANGE_NAME2="flysimCMQ";
  char const *exchangetype="fanout"; // direct, topic, headers and fanout
  char const *routingkey="test"; //for sending
  char const *bindingkey="test"; //for listening

  void MQSend(string &);
  void MQListen();
  void MQclose(amqp_connection_state_t &);
  char *chbuff;
  bool MQListop;

  mutex mtx;
};
#endif


class NetSim : public Parser
{
  public:
  NetSim();
  ~NetSim();

  vector<TmEvt> Evts;
  vector<OutCtr> OFCtr;
  float EvTime;

  ofstream *OFiles;
  void SetOut(void);
  void FileOpen(void);
  void FileClose(void);

  int Repeat;  //repeat times
  int RpSize;  //total repeat times

  unsigned int iEvt;
  unsigned int ActCtr;
  float dtime;
  unsigned int rseed;

  virtual void DoEvt(void)=0;
  virtual void SetPp(int,int)=0;
  virtual void RunInit(int,int)=0;
  virtual void Run(void)=0;
  virtual void FileOut(void)=0;

};


class NetSim_LIF3 : public NetSim, public Information, public EqsThd<EqsLIF3,Neuron_LIF3>
{
  public:
  NetSim_LIF3();
  ~NetSim_LIF3();

  vector<Neuron_LIF3> Nets;

  void DoEvt(void);
  void SetPp(int,int);
  void RunInit(int,int);
  void Run(void);
  void Runsteps(int);
  void FileOut(void);

  //for Confs
  vector<unsigned int> NIdx;
  unsigned int TolNeu;
  unsigned int TolSynp;
  unsigned int TolRep;

  int Step;

  private:
  void SInfo(string&);
  void TexInfo(string&);
  ostringstream ostr_tmp;

};

class NetSim_gp4 : public NetSim, public Information, public EqsThd<EqsLIF_gp4,Neuron_gp4>
{
  public:
  NetSim_gp4();
  ~NetSim_gp4();

  int SolType; // solver type
  vector<Neuron_gp4> Nets;

  void DoEvt(void);
  void SetPp(int,int);
  void RunInit(int,int);
  void Run(void);

  private:
  void SInfo(string&);
  void TexInfo(string&);

  int Step;
  void Activation(void);
  void FileOut(void);

  //for Confs
  vector<unsigned int> NIdx;
  unsigned int TolNeu;
  unsigned int TolSynp;
  unsigned int TolRep;

};


#ifdef RABBITMQ
class NetSimGnl : public NetSim, public EqsGnlThd, public Information, public rbmq
#else
class NetSimGnl : public NetSim, public EqsGnlThd, public Information
#endif
{
  public:
  NetSimGnl();
  ~NetSimGnl();

  vector<TmEvt> NEvts;
  vector<OutCtr> NOFCtr;
  NeuronGnl *Nets;

  void DoEvt(void);
  void SetPp(int,int);
  void RunInit(int,int);
  void Run(void);
  void FileOut(void);

  void DoNowEvt(void);

#ifdef RABBITMQ
  void NowEvtPrs(char *); // Now Event parse
  thread *thrbMQ;
#endif

  //for Confs
  vector<unsigned int> NIdx;
  unsigned int TolNeu;
  unsigned int TolSynp;
  unsigned int TolRep;

  unsigned long long Step;
  unsigned int OneMsStep;
  ostringstream ostr_tmp;
  ostringstream ostr2_tmp;

  unsigned int PosBuild1(unsigned int j,unsigned int PosIDNum,unsigned int PosIDStar,unsigned int PosNum, bool SelfConEn)
  {
    unsigned int PosID;
    if( PosNum == 1)
      PosID=PosIDStar;
    else
    {
      vector<unsigned int> ConAry(PosNum,0);
      bool SelfCon=false;
      unsigned int jStart=j-PosIDStar;

      for(unsigned int jj=0;jj<PosNum;jj++)
      {
        ConAry[jj]=PosIDStar+jj;
        if( (jStart==jj) && !SelfConEn ) { SelfCon=true; }

      }

      if(SelfCon) { ConAry.erase(ConAry.begin()+jStart); } //self connection

      shuffle(ConAry.begin(), ConAry.end(),*(Nets[j].rdgen));
      PosID=ConAry[0];
      {vector<unsigned int> ().swap(ConAry);}
    }

    return PosID;
  }

  unsigned int* PosBuild(unsigned int j,unsigned int PosIDNum,unsigned int PosIDStar,unsigned int PosNum, bool SelfConEn)
  {
    vector<unsigned int> ConAry(PosNum,0);
    bool SelfCon=false;
    unsigned int *PosID = new unsigned int[PosIDNum];
    unsigned int jStart=j-PosIDStar;

    for(unsigned int jj=0;jj<PosNum;jj++)
    {
      ConAry[jj]=PosIDStar+jj;
      if( (jStart==jj) && !SelfConEn ) { SelfCon=true; }
    }

    if(SelfCon) { ConAry.erase(ConAry.begin()+jStart); } //self connection

    shuffle(ConAry.begin(), ConAry.end(),*(Nets[j].rdgen));
    sort(ConAry.begin(),(ConAry.begin()+PosIDNum));

    for(unsigned int jj=0;jj<PosIDNum;jj++)
    { *(PosID+jj)=ConAry[jj]; }

    {vector<unsigned int> ().swap(ConAry);}
    return PosID;
  }

  private:
  void SInfo(string&);
  void TexInfo(string&);

};


class NetSimGnl1_1 : public NetSim, public EqsGnl1_1, public Information
{
  public:
  NetSimGnl1_1();
  ~NetSimGnl1_1();

  vector<TmEvt> NEvts;
  vector<OutCtr> NOFCtr;


  void SetPp(int,int);
  void DoEvt(void);
  void RunInit(int,int);
  void Run(void);
  void SimThd(unsigned int);
  void SimOne(void);
  void FileOut(void);
  //void DoNowEvt(void);
  void PostSetPp(void);
  void logDataWrite(const char *);

  bool NoLmt;
  unsigned long long Endtime;
  unsigned int ReportStep;

  //for Confs
  vector<unsigned int> NIdx;
  unsigned int TolNeu;
  unsigned int TolSynp;
  unsigned int TolRep;

  unsigned long long Step;
  unsigned int OneMsStep;
  ostringstream ostr_tmp;

//post neuron build up
  unsigned int PosBuild1(unsigned int j,unsigned int PosIDNum,unsigned int PosIDStar,unsigned int PosNum, bool SelfConEn)
  {
    unsigned int PosID;
    if( PosNum == 1)
      PosID=PosIDStar;
    else
    {
      vector<unsigned int> ConAry(PosNum,0);
      bool SelfCon=false;
      unsigned int jStart=j-PosIDStar;

      for(unsigned int jj=0;jj<PosNum;jj++)
      {
        ConAry[jj]=PosIDStar+jj;
        if( (jStart==jj) && !SelfConEn ) { SelfCon=true; }

      }

      if(SelfCon) { ConAry.erase(ConAry.begin()+jStart); } //self connection

      shuffle(ConAry.begin(), ConAry.end(),*rdgen);
      PosID=ConAry[0];
      {vector<unsigned int> ().swap(ConAry);}
    }

    return PosID;
  }

  unsigned int* PosBuild(unsigned int j,unsigned int PosIDNum,unsigned int PosIDStar,unsigned int PosNum, bool SelfConEn)
  {
    vector<unsigned int> ConAry(PosNum,0);
    bool SelfCon=false;
    unsigned int *PosID = new unsigned int[PosIDNum];
    unsigned int jStart=j-PosIDStar;

    for(unsigned int jj=0;jj<PosNum;jj++)
    {
      ConAry[jj]=PosIDStar+jj;
      if( (jStart==jj) && !SelfConEn ) { SelfCon=true; }
    }

    if(SelfCon) { ConAry.erase(ConAry.begin()+jStart); } //self connection

    shuffle(ConAry.begin(), ConAry.end(),*rdgen);
    sort(ConAry.begin(),(ConAry.begin()+PosIDNum));

    for(unsigned int jj=0;jj<PosIDNum;jj++)
    { *(PosID+jj)=ConAry[jj]; }

    {vector<unsigned int> ().swap(ConAry);}
    return PosID;
  }

  private:
  void SInfo(string&);
  void TexInfo(string&);

};


class NetSimGnl2 : public NetSim, public EqsGnl2, public Information
{
  public:
  NetSimGnl2();
  ~NetSimGnl2();

  vector<TmEvt> NEvts;
  vector<OutCtr> NOFCtr;


  void SetPp(int,int);
  void DoEvt(void);
  void RunInit(int,int);
  void Run(void);
  void SimThd(unsigned int);
  void SimOne(void);
  void FileOut(void);
  //void DoNowEvt(void);
  void PostSetPp(void);
  void logDataWrite(const char *);

  bool NoLmt;
  unsigned long long Endtime;
  unsigned int ReportStep;

  //for Confs
  vector<unsigned int> NIdx;
  unsigned int TolNeu;
  unsigned int TolSynp;
  unsigned int TolRep;

  unsigned long long Step;
  unsigned int OneMsStep;
  ostringstream ostr_tmp;

//post neuron build up
  unsigned int PosBuild1(unsigned int j,unsigned int PosIDNum,unsigned int PosIDStar,unsigned int PosNum, bool SelfConEn)
  {
    unsigned int PosID;
    if( PosNum == 1)
      PosID=PosIDStar;
    else
    {
      vector<unsigned int> ConAry(PosNum,0);
      bool SelfCon=false;
      unsigned int jStart=j-PosIDStar;

      for(unsigned int jj=0;jj<PosNum;jj++)
      {
        ConAry[jj]=PosIDStar+jj;
        if( (jStart==jj) && !SelfConEn ) { SelfCon=true; }

      }

      if(SelfCon) { ConAry.erase(ConAry.begin()+jStart); } //self connection

      shuffle(ConAry.begin(), ConAry.end(),*rdgen);
      PosID=ConAry[0];
      {vector<unsigned int> ().swap(ConAry);}
    }

    return PosID;
  }

  unsigned int* PosBuild(unsigned int j,unsigned int PosIDNum,unsigned int PosIDStar,unsigned int PosNum, bool SelfConEn)
  {
    vector<unsigned int> ConAry(PosNum,0);
    bool SelfCon=false;
    unsigned int *PosID = new unsigned int[PosIDNum];
    unsigned int jStart=j-PosIDStar;

    for(unsigned int jj=0;jj<PosNum;jj++)
    {
      ConAry[jj]=PosIDStar+jj;
      if( (jStart==jj) && !SelfConEn ) { SelfCon=true; }
    }

    if(SelfCon) { ConAry.erase(ConAry.begin()+jStart); } //self connection

    shuffle(ConAry.begin(), ConAry.end(),*rdgen);
    sort(ConAry.begin(),(ConAry.begin()+PosIDNum));

    for(unsigned int jj=0;jj<PosIDNum;jj++)
    { *(PosID+jj)=ConAry[jj]; }

    {vector<unsigned int> ().swap(ConAry);}
    return PosID;
  }

  private:
  void SInfo(string&);
  void TexInfo(string&);

};



class NetSimGnl2_1 : public NetSim, public EqsGnl2_1, public Information
{
  public:
  NetSimGnl2_1();
  ~NetSimGnl2_1();

  vector<TmEvt> NEvts;
  vector<OutCtr> NOFCtr;


  void SetPp(int,int);
  void DoEvt(void);
  void RunInit(int,int);
  void Run(void);
  void SimThd(unsigned int);
  void SimOne(void);
  void FileOut(void);
  //void DoNowEvt(void);
  void PostSetPp(void);
  void logDataWrite(const char *);

  bool NoLmt;
  unsigned long long Endtime;
  unsigned int ReportStep;

  //for Confs
  vector<unsigned int> NIdx;
  unsigned int TolNeu;
  unsigned int TolSynp;
  unsigned int TolRep;

  unsigned long long Step;
  unsigned int OneMsStep;
  ostringstream ostr_tmp;

//post neuron build up
  unsigned int PosBuild1(unsigned int NID,unsigned int PosLinkNum,unsigned int PosIDStar,unsigned int PosPpNum, bool SelfConEn)
  {
    unsigned int PosID;
    if( PosPpNum == 1)
      PosID=PosIDStar;
    else
    {
      vector<unsigned int> ConAry(PosPpNum,0);
      bool SelfCon=false;

      for(unsigned int jj=PosIDStar;jj<(PosIDStar+PosPpNum);jj++)
      {
        ConAry[jj-PosIDStar]=jj;
        if( (NID==jj) && !SelfConEn ) { SelfCon=true; }

      }

      if(SelfCon) { ConAry.erase(ConAry.begin()+NID-PosIDStar); } //self connection

      shuffle(ConAry.begin(), ConAry.end(),*rdgen);
      PosID=ConAry[0];
      {vector<unsigned int> ().swap(ConAry);}
    }

    return PosID;
  }

  unsigned int* PosBuild(unsigned int NID,unsigned int PosLinkNum,unsigned int PosIDStar,unsigned int PosPpNum, bool SelfConEn)
  {
    vector<unsigned int> ConAry(PosPpNum,0);
    bool SelfCon=false;
    unsigned int *PosID = new unsigned int[PosLinkNum];

    for(unsigned int jj=PosIDStar;jj<(PosIDStar+PosPpNum);jj++)
    {
      ConAry[jj-PosIDStar]=jj;
      if( (NID==jj) && !SelfConEn ) { SelfCon=true; }
    }

    if(SelfCon) { ConAry.erase(ConAry.begin()+NID-PosIDStar); } //self connection

    shuffle(ConAry.begin(), ConAry.end(),*rdgen);
    sort(ConAry.begin(),(ConAry.begin()+PosLinkNum));

    for(unsigned int jj=0;jj<PosLinkNum;jj++)
    { PosID[jj]=ConAry[jj]; }

    {vector<unsigned int> ().swap(ConAry);}
    return PosID;
  }

  private:
  void SInfo(string&);
  void TexInfo(string&);

};


class NetSimLEGWorm : public Information
{
  public:
  NetSimLEGWorm();
  ~NetSimLEGWorm();

  vector<NetSimGnl2 *> Worms;

  unsigned int AddWorm(string&,int);
  void DelWorm(unsigned int);
  void DoEvt(string&,unsigned int);
  void RunWorm(int);

  void DaemonMode(int);
  int CltReq(int);
  void Ack(int);

  int TSteps;
  private:
  void TexInfo(string&);
  void LdEst(int);  //load estimate
  long long CID;
};


