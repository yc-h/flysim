
#include <vector>
#include "StiExt.h"

using namespace std;

class ConTab
{
  public:
  ConTab();
  ~ConTab(){};

  unsigned int NeuID;  //post neuron ID

  unsigned short AMPA_ID; // AMPA receptor ID on post neuron's dendrites
  unsigned short GABA_ID; // GABA receptor ID on post neuron's dendrites
  unsigned short NMDA_ID; // NMDA receptor ID on post neuron's dendrites

  bool AMPA_Act; //AMPA receptor active flag, true for active
  bool GABA_Act; //GABA receptor active flag, true for active
  bool NMDA_Act; //NMDA receptor active flag, true for active

};

class ConTab_gp
{
  public:
  ConTab_gp();
  ~ConTab_gp(){};

  unsigned int NeuID;  //post neuron ID

  float AMPA_MeanG;  // pre-synapse to post-synapse connection mean efficacy
  float GABA_MeanG;  // pre-synapse to post-synapse connection mean efficacy
  float NMDA_MeanG;  // pre-synapse to post-synapse connection mean efficacy

  bool AMPA_Act; //AMPA receptor active flag, true for active
  bool GABA_Act; //GABA receptor active flag, true for active
  bool NMDA_Act; //NMDA receptor active flag, true for active

};


class Nlink
{
  public:
  Nlink();
  ~Nlink(){};

  float AMPA_st,AMPA_tst;  // pre-synapse to post-synapse connection mean efficacy
  float GABA_st,GABA_tst;  // pre-synapse to post-synapse connection mean efficacy
  float NMDA_xt,NMDA_st,NMDA_txt,NMDA_tst,NMDA_as;  // pre-synapse to post-synapse connection mean efficacy

};


class ConTab_gp2 : public ConTab_gp
{
  public:
  ConTab_gp2();
  ~ConTab_gp2(){};

  vector<Nlink> Axn;

  float AMPA_tst;  // pre-synapse to post-synapse connection mean efficacy
  float GABA_tst;  // pre-synapse to post-synapse connection mean efficacy
  float NMDA_txt,NMDA_tst,NMDA_as;  // pre-synapse to post-synapse connection mean efficacy

  float GAP_MeanG;  // pre-synapse to post-synapse connection mean efficacy
  bool GAP_Act; //GAP junction active flag, true for active

};

struct ConTab2
{
  unsigned int NeuID;  //post-Synaptic neuron ID
  unsigned short AMPA_ID; // AMPA receptor ID on post neuron's dendrites
  unsigned short GABA_ID; // GABA receptor ID on post neuron's dendrites
  unsigned short NMDA_ID; // NMDA receptor ID on post neuron's dendrites
  unsigned short Act; // receptors active flag
};

struct ConTab3
{
  unsigned int NeuID;  //post-Synaptic neuron ID
  unsigned short AMPA_ID; // AMPA receptor ID on post neuron's dendrites
  unsigned short ACH_ID; // Ach receptor ID on post neuron's dendrites
  unsigned short GABA_ID; // GABA receptor ID on post neuron's dendrites
  unsigned short GCL_ID; // Glu Cl receptor ID on post neuron's dendrites
  unsigned short NMDA_ID; // NMDA receptor ID on post neuron's dendrites
  unsigned short GAP_ID; // NMDA receptor ID on post neuron's dendrites
  unsigned short Act; // receptors active flag
};

//-------------------------------neuron parameters

struct ParBase
{
// membrane parameter
  float Vth;  //threthod voltage
  float El;  //reversal potential
  float gl; // menbrane conductance
  float Cm; // menbrane capacitance

  float Vreset; //reset potential
  float dVmax;
};

struct ParBase2 : public ParBase
{
  float df;  // decay factor:exp(-dt/taum)
  float ef;  // error factor:taum*(1-exp(-dt/taum))
};

struct ParLIF
{
  int RfPed;  // refractory peroid
  int SpiDy; // spiking delay
};

struct ParStp
{

  float pv; // reduce factor of depression factor D
  float tD; // decay time constant D

  float aF; // the alpha of facilitation F
  float tF; // the decay time constant F

  float aCap; // the alpha of Ca^2+ peak
  float tCap; // decay time constant Ca^2+ peak

  float aCar; // alpha of Ca^2+ residual
  float tCar; // decay time constant Ca^2+ residual
};

struct ParStp2
{

  float pv; // reduce factor of depression factor D
  float dfD; // decay factor D
  float efD; // error factor D

  float aF; // the alpha of facilitation F
  float dfF; // decay factor F
  float efF; // error factor F

  float aCap; // the alpha of Ca^2+ peak
  float dfCap; // decay factor  Ca^2+ peak

  float aCar; // alpha of Ca^2+ residual
  float dfCar; // decay factor Ca^2+ residual
  float efCar; // error factor Ca^2+ residual
};

struct ParHH
{
  float Ek; //reversal potential potassium-ion
  float gk; // menbrane conductance
  float Ena; //reversal potential sodium-ion
  float gna; // menbrane conductance
};

struct ParHH2
{
  float Ek; //reversal potential potassium-ion
  float gk; // menbrane conductance
  float Ena; //reversal potential sodium-ion
  float gna; // menbrane conductance

  float AnA;
  float AnB;
  float AnC;
  float AnD;

  float BnA;
  float BnB;
  float BnC;
  float BnD;

  float AmA;
  float AmB;
  float AmC;
  float AmD;

  float BmA;
  float BmB;
  float BmC;
  float BmD;

  float AhA;
  float AhB;
  float AhC;
  float AhD;

  float BhA;
  float BhB;
  float BhC;
  float BhD;

};

struct ParLP
{
// long term plsiticity
  float tLP; // time constant of long term plsiticity factor
  int SpiLPDy; //long term plsiticity propogation spike delay

  float PosISI; // +additive or -multiplicative rule of positive interspike intervals
  float NegISI; // +additive or -multiplicative rule of negative interspike intervals
};

struct ParLP2
{
// long term plsiticity
  float dfLP;
  int SpiLPDy; //long term plsiticity propogation spike delay

  float PosISI; // +additive or -multiplicative rule of positive interspike intervals
  float NegISI; // +additive or -multiplicative rule of negative interspike intervals
};

struct ParAxA
{
// AxoAxonic dynsmic
  float tAxA; // time constant of AxoAxonic plsiticity factor
  int SpiAxADy; //AxoAxonic dynsmic spike delay
};

class NeuPar : public ParBase, public ParLIF
{
  public:
  NeuPar();
  ~NeuPar(){};
};

class NeuPar_stp : public NeuPar, public ParStp
{
  public:
  NeuPar_stp();
  ~NeuPar_stp(){};

};

class NeuPar_HH : public NeuPar, public ParHH
{
  public:
  NeuPar_HH();
  ~NeuPar_HH(){};

};

class NeuPar_HST : public NeuPar_stp, public ParHH
{
  public:
  NeuPar_HST();
  ~NeuPar_HST(){};
};

class NeuPar_bp : public NeuPar
{
  public:
  NeuPar_bp();
  ~NeuPar_bp(){};

  float tx;  // taum factor
  float ax;
  float vx;

};

class NeuPar_uni : public NeuPar_bp
{
  public:
  NeuPar_uni();
  ~NeuPar_uni(){};

  float ta;  // taum factor
  float aa;
  float va;

};

class NeuPar_gp : public NeuPar
{
  public:
  NeuPar_gp();
  ~NeuPar_gp(){};

  bool SelfSti; //self-stimilation flag

};

struct NeuPar_LP : public ParBase, public ParLIF, public ParLP
{};


//------------neuron varables
struct VarBase
{
  float Isyn;  // current of synapse
  float V;  // menbrane potential
};

struct VarLIF
{
  int Count;  // firing table index
  bool Firing;  // firing flag, true for firing
  bool SpiOut;  // spike of this neuron, true for spike
};

struct VarLIF2
{
  int Count;  // firing table index
  unsigned int NState;
};

struct VarStp
{
// short-term plasticity
  float D;  // depression factor
  float F;  // facilitataion factor
  float Cap; // Ca^2+ peak
  float Car;  // Ca^2+ residual
};

struct VarHH
{
//Hodgkin-Huxley model
  float n;
  float m;
  float h;
};

struct VarLP
{
  float LP; // long term plsiticity factor
  bool FirLP;  // firing flag
  bool SpiLP; // delayed spike
  int CountLP;  // delay counter
};

struct VarLP2 //FirLP & SpiLP are combined into NSatae
{
  float LP; // long term plsiticity factor
  int CountLP;  // delay counter
};

struct VarAxA
{
  float AxA; // AxoAxonic plsiticity factor
  bool FirAxA;  //AxoAxonic dynsmic spike flag
  bool SpiAxA; //AxoAxonic dynsmic delayed spike
  int CountAxA;  //AxoAxonic dynsmic delay counter
};


class NeuVar : public VarBase, public VarLIF
{
  public:
  NeuVar();
  ~NeuVar(){};
};

class NeuVar_stp : public NeuVar , public VarStp
{
  public:
  NeuVar_stp();
  ~NeuVar_stp(){};
};

class NeuVar_HH : public NeuVar , public VarHH
{
  public:
  NeuVar_HH();
  ~NeuVar_HH(){};
};

class NeuVar_HST : public NeuVar , public VarHH , public VarStp
{
  public:
  NeuVar_HST();
  ~NeuVar_HST(){};
};


class NeuVar_bp : public NeuVar
{
  public:
  NeuVar_bp();
  ~NeuVar_bp(){};

  float Vd;
};

class NeuVar_uni : public NeuVar_bp
{
  public:
  NeuVar_uni();
  ~NeuVar_uni(){};

  float Va;
};


struct NeuVar_LP : public VarBase, public VarLIF, public VarLP
{};


//--------------------------------------------------------------------
struct NVar : public VarBase, public VarLIF
{};

struct NVar_ltp : public NVar, public VarLP
{};

struct NPar : public ParBase, public ParLIF
{};

struct NPar_ltp : public NPar, public ParLP
{};

//------------------------------------------------------------------------


struct SynCon
{
  unsigned int PosIDNum;
  unsigned int *PosID;
};
struct VarFast
{
  float *gFast; // synaptic conductance
  float stFast;  // gating variable
};

struct VarSlow
{
  float *gSlow; // synaptic conductance
  float xtSlow; // gating variable
  float stSlow; // gating variable
};


struct SynFast : public SynCon, public VarFast
{};

struct SynSlow : public SynCon, public VarSlow
{};

struct SynGlut : public SynCon, public VarFast , public VarSlow
{};

struct SynGlutType
{ unsigned char A,N; };

struct SynGap : public SynCon
{ float *gC; }; // Gap junciton synaptic conductance or current basd model resistent

struct SynFast1
{
  unsigned int PosID;
  float gFast; // synaptic conductance
  float stFast;  // gating variable
};

struct SynSlow1
{
  unsigned int PosID;
  float gSlow; // synaptic conductance
  float xtSlow; // gating variable
  float stSlow; // gating variable
};

struct SynGlut1
{
  unsigned int PosID;
  float gFast; // synaptic conductance
  float stFast;  // gating variable
  float gSlow; // synaptic conductance
  float xtSlow; // gating variable
  float stSlow; // gating variable
};

struct SynGap1
{
  unsigned int PosID;
  float gC; // synaptic conductance
};

struct IdxSynapse
{
  unsigned int PosPp;
  unsigned int PosIDNum;
  float gFast;
  float gSlow;
  unsigned short Act;
};


//------- multithread Post ID remap indics

struct NeuMap
{
  unsigned int Beg;
  unsigned int End;
};

struct NeuMap2
{
  unsigned int Beg;
  unsigned int End;
  unsigned int PosID;
};
//------------------------------------------------------------------------

struct SynIdx
{
  unsigned int PreID;
  unsigned int PosID;
};


struct SynFast2 : public SynIdx
{
  float g; // synaptic conductance
  float st;  // gating variable
};

struct SynSlow2 : public SynIdx
{
  float g; // synaptic conductance
  float xt; // gating variable
  float st; // gating variable
};

struct SynGap2 : public SynIdx
{ float g; };// synaptic conductance

//----population model
//-------AxoAnonal synapse

struct PSynIdx
{  unsigned int PreID;};

struct PSynFast2 : public PSynIdx
{ float st; }; // gating variable


struct PSynSlow2 : public PSynIdx
{
  float xt; // gating variable
  float st; // gating variable
};

struct PSynGap
{
  unsigned int PreID;
  float g; // gating variable
};


struct IdxSynapse2
{
  unsigned int PosPp;
  unsigned int PosIDNum;
  float g;
  unsigned short Act;
};


//---------------------------------------------------
NeuDatTmp
struct NeuDatBase
{
  vector<ConT> Index; // index of post neurons

// receptor properties
  RepP ParAMPA;
  RepP ParGABA;
  RepP2 ParNMDA;

  vector<RepV> RepAMPA; // AMPA receptors in this neuron
  vector<RepV> RepGABA; // GABA receptors in this neuron
  vector<RepV2> RepNMDA; // NMDA receptors in this neuron

  vector<bool> AMPA_Sti; // AMPA stimulation index of post neurons
  vector<bool> GABA_Sti; // GABA stimulation index of post neurons
  vector<bool> NMDA_Sti; // NMDA stimulation index of post neurons
};


NeuDatTmp
struct NeuDatBase2 : public NeuDatGen
{
  RP1 ParACH;
  RP1 ParGCL;

  vector<RV1> RepACH; // ACH receptors in this neuron
  vector<RV1> RepGCL; // GCL receptors in this neuron
  vector<float> RepGAP; // GAP receptors in this neuron

  vector<bool> ACH_Sti; // ACH stimulation index of post neurons
  vector<bool> GCL_Sti; // GCL stimulation index of post neurons
};

NeuDatTmp
struct NeuGPBase : public NeuDatGen
{
  vector<NeuVar> NeuVm; // each neuron's membrane potential

  float gAMPA;
  float gGABA;
  float gNMDA;

};


NeuParTmp
struct NeuRepPar
{
// receptor properties
  RepP ParAMPA;
  RepP ParGABA;
  RepP2 ParNMDA;
  RepP ParACH;
  RepP ParGCL;

//population gating variable sumation
  float *gAMPA;
  float *gGABA;
  float *gNMDA;
  float *gACH;
  float *gGCL;
};

NeuParTmp
struct NeuRepPar2
{
// receptor properties
  RepP ParAMPA;
  RepP ParGABA;
  RepP2 ParNMDA;
  RepP ParACH;
  RepP ParGCL;
};

struct NeuRepVar2
{
//population gating variable sumation
  float gAMPA;
  float gGABA;
  float gNMDA;
  float gACH;
  float gGCL;
};

struct NeuRepVar3
{
//population conductance sumation
  float gA;
  float gG;
  float gH;
  float gL;
};
  
struct NeuGatVar3
{
//population gating variable sumation
  float stA;
  float stG;
  float stH;
  float stL;
};





