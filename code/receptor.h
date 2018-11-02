
#include <string>

using namespace std;

class RepPar1
{
  public:
  RepPar1();
  ~RepPar1(){};

  // parameters
  float tst; //tau decay time constant st
  float Vrev; // reversal potential

};

class RepPar2 : public RepPar1
{
  public:
  RepPar2();
  ~RepPar2(){};

  float as; // the alpha of st
  float txt; // the decay time constant of xt
  float mg; // the concentration of magnesium
};

class RepVar1
{
  public:
  RepVar1();
  ~RepVar1(){};

  float g; // synaptic conductance
  bool Sti; // stimulate input
  float st;

};

class RepVar2 : public RepVar1
{
  public:
  RepVar2();
  ~RepVar2(){};

  float xt;

};

//---------stp

class RepVar1_stp : public RepVar1
{
  public:
  RepVar1_stp();
  ~RepVar1_stp(){};

  float ax; // the alpha of xt

};

class RepVar2_stp : public RepVar2
{
  public:
  RepVar2_stp();
  ~RepVar2_stp(){};

  float ax; // the alpha of xt

};

//---------group model------------------------------------

class RepPar1_gp : public RepPar1
{
  public:
  RepPar1_gp();
  ~RepPar1_gp(){};

  float g; // Maxium Externel stimulaiton synaptic conductance: MeanEff*MeanCon
  float GTotal; // totel conductance
  float MeanG; //synaptic conductance: MeanEff, for self-stimilation subtraction
};


class RepPar2_gp : public RepPar2
{
  public:
  RepPar2_gp();
  ~RepPar2_gp(){};

  float g; // Maxium Externel stimulaiton synaptic conductance: MeanEff*MeanCon
  float GTotal; // totel conductance
  float MeanG; //synaptic conductance: MeanEff, for self-stimilation subtraction
};

class RepPar1_gp2 : public RepPar1
{
  public:
  RepPar1_gp2();
  ~RepPar1_gp2(){};

  float g; // Maxium Externel stimulaiton synaptic conductance: MeanEff*MeanCon
  float MeanG; //synaptic conductance: MeanEff, for self-stimilation subtraction
};


class RepPar2_gp2 : public RepPar2
{
  public:
  RepPar2_gp2();
  ~RepPar2_gp2(){};

  float g; // Maxium Externel stimulaiton synaptic conductance: MeanEff*MeanCon
  float MeanG; //synaptic conductance: MeanEff, for self-stimilation subtraction
};


class RepVar1_gp
{
  public:
  RepVar1_gp();
  ~RepVar1_gp(){};

  float st;

};

class RepVar2_gp : public RepVar1_gp
{
  public:
  RepVar2_gp();
  ~RepVar2_gp(){};

  float xt;

};

//-----------compact LIF

struct RP1
{
  float tst; //tau decay time constant st
  float Vrev; // reversal potential
};

struct RP2 : public RP1
{
  float as; // the alpha of st
  float txt; // the decay time constant of xt
  float mg; // the concentration of magnesium
};

struct RV1
{
  float g; // synaptic conductance
  float st;
};

struct RV2 :  public RV1
{ float xt; };



//---------------------
struct RP3
{
  float Vrev; // reversal potential
  float dfst; // decay factor: st
};

struct RP4
{
  float Vrev; // reversal potential
  float dfst; // decay factor: st
  float efst; // error factor: st

  float as; // the alpha of st
  float dfxt; // decay factor: xt
  float mg; // the concentration of magnesium
};

struct RV3
{
  float g; // synaptic conductance
  float st;
};

struct RV4 :  public RV3
{ float xt; };

