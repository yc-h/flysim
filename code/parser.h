#include <string>
#include <vector>
#include <random>
#include <memory>
#include "Eqs.h"



struct NeuCon  // link to target population feature
{
  unsigned int TarPpName;
  unsigned char TarRep;

  float Connectivity;
  float MeanEff;
  float weight;

};

struct NeuRep  // own population feature
{
  unsigned char NeuRepName;

  float Tau; // tst
  float TauXT;
  float Vrev; // reversal potentail of receptor
  float FExt;
  float MnExtEff;
  float MnExtCon;

};

struct NetConfSta
{
  string NeuPpName;
  unsigned int NeuNum; //number of neurons
  unsigned int PpRepNum;

  unsigned int PpLinkBeg;
  unsigned int PpLinkNum;

  unsigned int TriPpLinkNum;
  bool SelConEn;
  NeuRep *PpRep;
};

//struct NetConf : public NeuPar, public ParStp, public ParLP
struct NetConf : public NeuPar
{

// leaky integrate and fire mode
  float Taum; //membrane time constant = g/c

//calcium plasticity
/*
  float Vrsting_ca;
  float Alpha_ca;
  float Tau_ca;
  float Eff_ca;

  float STDeP;
  float STDeTau;
*/
  //unsigned int PpLinkBeg;
  NeuCon *PpLink;
};

class TmEvt : public Information
{
  public:
  TmEvt();
  ~TmEvt(){};

  void init();

  float TimEv;  // event time
  bool NoLmt;  // no limitation
  unsigned int Type;  // event type
  //string Type;  // event type
  //string Label;  // event label
  unsigned int Pp; //target population name
  unsigned char SubTyp; //sub type

  float Var1;
  float Var2;

  private:
  void TexInfo(string&);
};

struct GroupExtent
{
  string GpName;
  vector<unsigned int> Pp; //target population name
};


class OutCtr : public Information
{
  public:
  OutCtr();
  ~OutCtr();

  void init();
  unsigned int Status; //reorder,SumPars

  string Type;
  int TypeID;
  string FileName;
  vector<unsigned int> ParList;

  unsigned int RWin;
  unsigned int FPnt;
  vector<unsigned int> SpiRWin;

  float TimEv;  // event time
  //bool SumPars;
  private:
  void TexInfo(string&);
};



class Parser
{
  public:
  Parser();
  ~Parser();

  string OutMem;
  string OutSpi;
  string OutFR;

  vector<NetConfSta> ConfSta;
  NetConf *Confs;
  ParHH2 *ConfHH;
  ParStp *ConfStp;
  ParLP *ConfLtp;

  NeuCon *PostPp;

  vector<GroupExtent> GpExt;
  vector<unsigned int> MacRndSeds;

  void CSta(const string&,float);  // for network.conf staticitic
  void parse(const string&,float);  // for network.conf
  void parse(const string&,vector<TmEvt>&,vector<OutCtr>&); // for network.pro
  void DefMacro(const string&); // for network.pro


  bool CommMem; // command line option for membrane potential
  bool CommSpi; // command line option for spike raster
  bool CommFR; // command line option for firing rate
  void Comm(vector<OutCtr>&); // command line option for network.pro

  void ConfsInfo(const string&);

  unsigned short GetAct(unsigned char);

  unsigned int RSeed;

  private:
  int state;
  string str;
  vector<string> tokens;
  string delimiters;

  void Tokenize();

};




