#include <vector>
#include <string>
#include <random>
#include "neuron_affi.h"


using namespace std;




//----------------------------------------------------------------------------------------
class Neuron_LIF3 : public NeuDatBase<RP1,RV1,RP2,RV2,ConTab2>, public Information, public NeuPar, public NeuVar, public StiExt
{
  public:
  Neuron_LIF3();
  virtual ~Neuron_LIF3();

  float Irep[3];
  void BgInit();  // back ground initial
  void ParInit(); // Parameter initial
  void VarInit(); // variable initial
  void TexSInfo(string&); // short information of this neuron for save/resume

  private:
  void TexInfo(string&); // information of this neuron

};

class Neuron_gp4 : public NeuGPBase<RepPar1_gp2,RepVar1_gp,RepPar2_gp2,RepVar1_gp,ConTab_gp2>, public Information, public NeuPar_gp, public StiExt
{
  public:
  Neuron_gp4();
  virtual ~Neuron_gp4();

  void BgInit();  // back ground initial
  void ParInit(); // Parameter initial
  void VarInit(); // variable initial

  private:
  void TexInfo(string&); // information of this neuron

};

class NeuronGnl :
public Information, public StiExt3<RV1,RV2>,
public ParBase, public ParLIF, public ParHH, public ParStp, public ParLP,
public VarBase, public VarLIF, public VarHH, public VarStp, public VarLP,
public NeuRepPar<RP1,RP2>
{
  public:
  NeuronGnl();
  virtual ~NeuronGnl();

  void BgInit(float);  // back ground initial
  void ParInit(float); // Parameter initial
  void VarInit(unsigned int); // variable initial


  unsigned int ActCtr;

  float *Irep;

//synapses
  SynFast1 *SynFst1;
  SynGlut1 *SynGlu1;
  SynSlow1 *SynSlo1;
  SynGap1 *SynCnd1;

  unsigned int SynFastSize1;
  unsigned int SynGlutSize1;
  unsigned int SynSlowSize1;
  unsigned int SynGapSize1;

  unsigned char *SFTyp1;
  unsigned char *SSTyp1;
  SynGlutType *SGTyp1;

  SynFast *SynFst;
  SynGlut *SynGlu;
  SynSlow *SynSlo;
  SynGap  *SynCnd;

  unsigned int SynFastSize;
  unsigned int SynGlutSize;
  unsigned int SynSlowSize;
  unsigned int SynGapSize;

  unsigned char *SFTyp;
  unsigned char *SSTyp;
  SynGlutType *SGTyp;

  private:
  void TexInfo(string&); // information of this neuron
};



//class NeuronGnl2 :
//public Information, public StiExt3<RV3,RV4>
class NeuronGnl2 :
public Information
{
  public:
  NeuronGnl2(){};
  ~NeuronGnl2(){};

  void VarInit(float); // Parameter initial

  unsigned int ActCtr;
  float *Irep;

  private:
  void TexInfo(string&); // information of this neuron
};



//class NeuronGnl2_1:
class NeuronGnl2_1 :
public Information
{
  public:
  NeuronGnl2_1(){};
  ~NeuronGnl2_1(){};

  void VarInit(float); // Parameter initial

  unsigned int ActCtr;
  float *Irep;

  private:
  void TexInfo(string&); // information of this neuron
};
