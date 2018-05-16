
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <omp.h>
#include <thread>
#include <mutex>
#include <atomic>
#include <pthread.h>
#include <random>
#include <condition_variable>
#include "neuron.h"



using namespace std;
//--------------solvers------------------------------

template<class EqsTyp>
class Solvers
{
  public:
  EqsTyp *pObj;

  //solvers
  float RK4(float num, float (EqsTyp::*FEqs)(float)) //Runge Kutta order 4th
  {
    float k1,k2,k3,k4;
    k1 = (pObj->*FEqs)( num );
    k2 = (pObj->*FEqs)( num + DTIME*k1/2.0f );
    k3 = (pObj->*FEqs)( num + DTIME*k2/2.0f );
    k4 = (pObj->*FEqs)( num + DTIME*k3 );
    return num + DTIME*( k1 + 2.0f*k2 + 2.0f*k3 + k4 )/6.0f;
  };

  float IpvEul(float num, float (EqsTyp::*FEqs)(float)) //improved Euler or RK2
  {
    float k1 = (pObj->*FEqs)( num );
    float k2 = (pObj->*FEqs)( num + DTIME*k1 );
    return num + DTIME*(k1+k2)/2.0f;
  };

  float Eul(float num, float (EqsTyp::*FEqs)(float)) //Euler
  {
    float k1=(pObj->*FEqs)(num);
    return num + DTIME*k1;
  };
};
//---------------------------------------------------------------------------
EqsTmp
struct EqsItr
{
  RepP PItr; // receptor parameters
  RepV VItr; // receptors variables
  RepV Rep_begin,Rep_end; //for accumulation receptor current

  RepP2 PItr2; // receptor parameters
  RepV2 VItr2; // receptors variables
  RepV2 Rep2_begin,Rep2_end; //for accumulation receptor current

  Neu NItr; // neuron list
  Neu NItr_begin; // for spike
};


typedef EqsItr<RP1 *, vector<RV1>::iterator, RP2 *, vector<RV2>::iterator, vector<Neuron_LIF3>::iterator > ItrLIF3;
typedef EqsItr<RepPar1_gp2 *, vector<RepVar1_gp>::iterator, RepPar2_gp2 *, vector<RepVar1_gp>::iterator, vector<Neuron_gp4>::iterator > ItrLIF_gp4;



//----------------------------------------------
class EqsLIF3 : public ItrLIF3
{
  public:
  EqsLIF3(){};

  float FRK4Irep1Acc(); //accumulation receptor current for AMPA GABA in Rk4 solver
  float FRK4Irep2Acc(); //accumulation receptor current for NMDA in Rk4 solver
  float FIpvIrep1Acc(); //accumulation receptor current for AMPA GABA in improved eular solver
  float FIpvIrep2Acc(); //accumulation receptor current for NMDA in improved eular solver
  float FEulIrep1Acc(); //accumulation receptor current for AMPA GABA in eular solver
  float FEulIrep2Acc(); //accumulation receptor current for NMDA in eular solver

  vector<bool>::iterator SItr,Sti_begin;
  int SolType;
  void Activation(unsigned int);
  void firing();
  void Spike(unsigned int);

  void FRK4Vm();
  void FIpvVm();
  void FEulVm();

  //neumerical method function
  float Fxt(float);  // receptor rising
  float Fst(float);  // receptor decay
  float FFst(float);  // receptor decay
  float MemPot(float);  // receptor decay

  //slovers
  float RK4(float num, float (EqsLIF3::*FEqs)(float)) //Runge Kutta order 4th
  {
    float k1,k2,k3,k4;
    k1 = (this->*FEqs)( num );
    k2 = (this->*FEqs)( num + DTIME*k1/2.0f );
    k3 = (this->*FEqs)( num + DTIME*k2/2.0f );
    k4 = (this->*FEqs)( num + DTIME*k3 );
    return num + DTIME*( k1 + 2.0f*k2 + 2.0f*k3 + k4 )/6.0f;
  };

  float IpvEul(float num, float (EqsLIF3::*FEqs)(float)) //improved Euler or RK2
  {
    float k1 = (this->*FEqs)( num );
    float k2 = (this->*FEqs)( num + DTIME*k1 );
    return num + DTIME*(k1+k2)/2.0f;
  };

  float Eul(float num, float (EqsLIF3::*FEqs)(float)) //Euler
  {
    float k1=(this->*FEqs)(num);
    return num + DTIME*k1;
  };

};

//----------------------------------------------
class EqsLIF_gp4 : public ItrLIF_gp4
{
  public:
  EqsLIF_gp4(){};

  float SpiIn;

  void FRK4Irep1Acc(float,float); //accumulation receptor current for AMPA GABA in RK4 solver
  void FRK4Irep1AccG(float,float); //accumulation receptor current for AMPA GABA in RK4 solver
  void FRK4Irep2Acc(float,float); //accumulation receptor current for NMDA in RK4 solver
  void FIpvIrep1Acc(float,float); //accumulation receptor current for AMPA GABA in improved eular solver
  void FIpvIrep1AccG(float,float); //accumulation receptor current for AMPA GABA in improved eular solver
  void FIpvIrep2Acc(float,float); //accumulation receptor current for NMDA in eular improved solver
  void FEulIrep1Acc(float,float); //accumulation receptor current for AMPA GABA in eular solver
  void FEulIrep1AccG(float,float); //accumulation receptor current for AMPA GABA in eular solver
  void FEulIrep2Acc(float,float); //accumulation receptor current for NMDA in eular solver

  vector<NeuVar>::iterator VmItr;
  int SolType;
  void Activation(unsigned int);
  void firing();
  void Spike(unsigned int);
  void SpikeRK4();
  void SpikeIpv();
  void SpikeEul();

  void FRK4Vm();
  void FIpvVm();
  void FEulVm();

  //neumerical method function
  float Fst(float);  // receptor decay
  float FFst(float);  // receptor decay
  float MemPot(float);  // receptor decay

  float tst; // tau of st
  float txt; // tau of xt
  float Gxt; // xt
  float as; // alpha of st

  float FGS1(float);
  float FGS2(float);

  //slovers
  float RK4(float num, float (EqsLIF_gp4::*FEqs)(float)) //Runge Kutta order 4th
  {
    float k1,k2,k3,k4;
    k1 = (this->*FEqs)( num );
    k2 = (this->*FEqs)( num + DTIME*k1/2.0f );
    k3 = (this->*FEqs)( num + DTIME*k2/2.0f );
    k4 = (this->*FEqs)( num + DTIME*k3 );
    return num + DTIME*( k1 + 2.0f*k2 + 2.0f*k3 + k4 )/6.0f;
  };

  float IpvEul(float num, float (EqsLIF_gp4::*FEqs)(float)) //improved Euler or RK2
  {
    float k1 = (this->*FEqs)( num );
    float k2 = (this->*FEqs)( num + DTIME*k1 );
    return num + DTIME*(k1+k2)/2.0f;
  };

  float Eul(float num, float (EqsLIF_gp4::*FEqs)(float)) //Euler
  {
    float k1=(this->*FEqs)(num);
    return num + DTIME*k1;
  };
};



//----------------------------------------------

template<class EqsTyp,class NTyp>
class EqsThd
{
  public:
  EqsThd()
  {
    Nets_size=0;
    EqsTime=0.0f;
    SolType=ACCURATE;
  };

  ~EqsThd()
  {
    delete [] th;
    vector<EqsTyp> ().swap(Equations);
  };

  int SolType;
  float EqsTime; //simulation time
  thread *th;
  unsigned int Nets_size;
  vector<EqsTyp> Equations;
  unsigned int LoadCnt;
  bool LoadEmpty;
  mutex mtx;

  void iniEqsThd(vector<NTyp> &Nets,unsigned int TNum)
  {
    th=new thread[TNum];
    Nets_size=Nets.size();
    EqsTyp EqsNew;
    Equations.insert(Equations.end(),TNum,EqsNew);

    for(unsigned int i=0;i<Equations.size();i++)
    {
      Equations[i].SolType=SolType;
      Equations[i].NItr_begin=Nets.begin();
    }

  };

  void RunActSteps(unsigned int thid,unsigned int Nets_end)
  {
    for(unsigned int NID=0;NID<Nets_end;NID++)
    { Equations[thid].Activation(NID); }

    for(unsigned int NID=0;NID<Nets_end;NID++)
    { Equations[thid].Spike(NID); }
  };

  void RunAct(unsigned int thid)
  {
    unsigned int NID;
    bool Stop=false;
    while(!LoadEmpty)
    {
      mtx.lock();
      if(LoadCnt<Nets_size)
      { NID=LoadCnt++; }
      else
      {
        LoadEmpty=true;
        Stop=true;
      }
      mtx.unlock();

      if(!Stop)
      { Equations[thid].Activation(NID); }
    }
  };

  void RunSpi(unsigned int thid)
  {
    unsigned int NID;
    bool Stop=false;
    while(!LoadEmpty)
    {
      mtx.lock();
      if(LoadCnt<Nets_size)
      { NID=LoadCnt++; }
      else
      {
        LoadEmpty=true;
        Stop=true;
      }
      mtx.unlock();

      if(!Stop)
      { Equations[thid].Spike(NID); }
    }
  };

  void RunThd()
  {

    LoadEmpty=false;
    LoadCnt=0;
    for(unsigned int i=0;i<Equations.size();i++)
    {*(th+i)=thread(&EqsThd::RunAct,this,i); }

    for(unsigned int i=0;i<Equations.size();i++)
    {(th+i)->join();}

    LoadEmpty=false;
    LoadCnt=0;
    for(unsigned int i=0;i<Equations.size();i++)
    { *(th+i)=thread(&EqsThd::RunSpi,this,i); }

    for(unsigned int i=0;i<Equations.size();i++)
    {(th+i)->join();}
  };

};

//----------------------------------------------

class EqsGnl
{
  public:
  EqsGnl()
  {
    dt=DTIME;
    thsize=1;
  };

  NeuronGnl *NItr;
  NeuronGnl *NItr_begin;
  int SolType;
  void Activation(unsigned int);
  void firing();
  void Spike(unsigned int,unsigned int);
  float dt;
  unsigned int thsize;

  //hogkin huxily model
  float IHHRK4();
  float IHHIpv();
  float IHHEul();
  float An();
  float Bn();
  float Am();
  float Bm();
  float Ah();
  float Bh();
  float dx(float);
  float Vx;
  float Ax;
  float Bx;

  //LTP
  void BackProp();
  void ltpRK4();
  void ltpIpv();
  void ltpEul();

  //STP
  void stpRK4();
  void stpIpv();
  void stpEul();

  //membrane function
  void VmRK4();
  void VmIpv();
  void VmEul();

  //neumerical method function
  float MemPot(float);  // receptor decay
  float FCar(float);
  float FCap(float);
  float FF(float);
  float FD(float);

  float FFst(float,float);
  float Fst(float,float,float,float);

//slovers
//-------------------pass 1 Variable
  float RK4(float num, float (EqsGnl::*FEqs)(float)) //Runge Kutta order 4th
  {
    float k1,k2,k3,k4;
    k1 = (this->*FEqs)( num );
    k2 = (this->*FEqs)( num + dt*k1/2.0f );
    k3 = (this->*FEqs)( num + dt*k2/2.0f );
    k4 = (this->*FEqs)( num + dt*k3 );
    return num + dt*( k1 + 2.0f*k2 + 2.0f*k3 + k4 )/6.0f;
  };

  float IpvEul(float num, float (EqsGnl::*FEqs)(float)) //improved Euler or RK2
  {
    float k1 = (this->*FEqs)( num );
    float k2 = (this->*FEqs)( num + dt*k1 );
    return num + dt*(k1+k2)/2.0f;
  };

  float Eul(float num, float (EqsGnl::*FEqs)(float)) //Euler
  {
    float k1=(this->*FEqs)(num);
    return num + dt*k1;
  };

//-------------------pass 2 Variable
  float RK4Var2(float num1, float num2, float (EqsGnl::*FEqs)(float,float)) //Runge Kutta order 4th
  {
    float k1,k2,k3,k4;
    k1 = (this->*FEqs)( num1 , num2 );
    k2 = (this->*FEqs)( num1 + dt*k1/2.0f , num2 );
    k3 = (this->*FEqs)( num1 + dt*k2/2.0f , num2);
    k4 = (this->*FEqs)( num1 + dt*k3 , num2);
    return num1 + dt*( k1 + 2.0f*k2 + 2.0f*k3 + k4 )/6.0f;
  };

  float IpvEulVar2(float num1, float num2, float (EqsGnl::*FEqs)(float,float)) //improved Euler or RK2
  {
    float k1 = (this->*FEqs)( num1, num2 );
    float k2 = (this->*FEqs)( num1 + dt*k1, num2 );
    return num1 + dt*(k1+k2)/2.0f;
  };

  float Eul2Var(float num1, float num2,float (EqsGnl::*FEqs)(float,float)) //Euler
  {
    float k1=(this->*FEqs)(num1,num2);
    return num1 + dt*k1;
  };

//-------------------pass 4 Variable
  float RK4Var4(float num1, float num2, float num3, float num4, float (EqsGnl::*FEqs)(float,float,float,float)) //Runge Kutta order 4th
  {
    float k1,k2,k3,k4;
    k1 = (this->*FEqs)( num1 , num2 , num3 , num4);
    k2 = (this->*FEqs)( num1 + dt*k1/2.0f , num2 , num3 , num4 );
    k3 = (this->*FEqs)( num1 + dt*k2/2.0f , num2 , num3 , num4);
    k4 = (this->*FEqs)( num1 + dt*k3 , num2 , num3 , num4);
    return num1 + dt*( k1 + 2.0f*k2 + 2.0f*k3 + k4 )/6.0f;
  };

  float IpvEulVar4(float num1,float num2,float num3,float num4, float (EqsGnl::*FEqs)(float,float,float,float)) //improved Euler or RK2
  {
    float k1 = (this->*FEqs)( num1 , num2 , num3 , num4);
    float k2 = (this->*FEqs)( num1 + dt*k1 , num2 , num3 , num4);
    return num1 + dt*(k1+k2)/2.0f;
  };

  float Eul4Var(float num1, float num2, float num3, float num4,float (EqsGnl::*FEqs)(float,float,float,float)) //Euler
  {
    float k1=(this->*FEqs)(num1,num2,num3,num4);
    return num1 + dt*k1;
  };
};


class EqsGnlThd
{
  public:
  EqsGnlThd()
  {
    dt=DTIME;
    SolType=ACCURATE;
    thsize=1;
  };

  ~EqsGnlThd()
  {
    delete [] th;
    delete [] Equations;
    delete [] ActNum;
    delete [] SpiNum;
  };

  int SolType;
  float dt;
  unsigned int thsize;
  thread *th;
  EqsGnl *Equations;
  unsigned int *ActNum;
  unsigned int *SpiNum;

  //random engine
  mt19937 rdgen;

  void iniEqsThd(NeuronGnl *Nets, unsigned int TNum)
  {
    thsize=TNum;
    if(thsize>1)
    {
      th=new thread[thsize];
      Equations = new EqsGnl[thsize];
      for(unsigned int i=0;i<thsize;i++)
      {
        (Equations+i)->dt=dt;
        (Equations+i)->thsize=thsize;
        (Equations+i)->SolType=SolType;
        (Equations+i)->NItr_begin=Nets;
      }
    }
    else
    {
      Equations = new EqsGnl;
      Equations->dt=dt;
      Equations->thsize=thsize;
      Equations->SolType=SolType;
      Equations->NItr_begin=Nets;
    }
  };

  void RunAct(unsigned int thid,unsigned int NIDBeg,unsigned int NIDEnd)
  {
    for(unsigned int NID=NIDBeg;NID<NIDEnd;NID++)
    { (Equations+thid)->Activation(NID); }
  };

  void RunSpi(unsigned int thid,unsigned int NIDBeg,unsigned int NIDEnd)
  {
    for(unsigned int NID=NIDBeg;NID<NIDEnd;NID++)
    { (Equations+thid)->Spike(thid,NID); }
  };

  void RunThd()
  {
    if(thsize>1)
    {
      for(unsigned int i=0;i<thsize;i++)
      {
        if(i==0)
          *(th+i)=thread(&EqsGnlThd::RunAct,this,i,0,*ActNum);
        else
          *(th+i)=thread(&EqsGnlThd::RunAct,this,i,*(ActNum+i-1),*(ActNum+i));
      }

      for(unsigned int i=0;i<thsize;i++)
      {(th+i)->join();}

      for(unsigned int i=0;i<thsize;i++)
      {
        if(i==0)
          *(th+i)=thread(&EqsGnlThd::RunSpi,this,i,0,*SpiNum);
        else
          *(th+i)=thread(&EqsGnlThd::RunSpi,this,i,*(SpiNum+i-1),*(SpiNum+i));
      }

      for(unsigned int i=0;i<thsize;i++)
      {(th+i)->join();}
    }
    else
    {
      RunAct(0,0,*ActNum);
      RunSpi(0,0,*SpiNum);
    }
  };

};




class EqsGnl1_1
{
  public:
  EqsGnl1_1();
  int SolType;

  //void CntInit(void);
  void Activation(unsigned int,unsigned int);
  void firing(unsigned int);
  void firingHH(unsigned int);
  void Spike(unsigned int,unsigned int);
  //void BgInit(unsigned int);
  void BgInit(unsigned int,unsigned int);
  void VarInit(void);

  void ThdSynapse(unsigned int);
  void ThdSynapseSTP(unsigned int);
  void ThdSynapseSTD(unsigned int);
  void ThdPGatVar(unsigned int);
  void ThdPGatVarSTP(unsigned int);
  void ThdPGatVarSTD(unsigned int);

  void ThdSynapseIpvEul(unsigned int);
  void ThdSynapseSTPIpvEul(unsigned int);
  void ThdSynapseSTDIpvEul(unsigned int);
  void ThdPGatVarIpvEul(unsigned int);
  void ThdPGatVarSTPIpvEul(unsigned int);
  void ThdPGatVarSTDIpvEul(unsigned int);

  void ThdSynapseRK4(unsigned int);
  void ThdSynapseSTPRK4(unsigned int);
  void ThdSynapseSTDRK4(unsigned int);
  void ThdPGatVarRK4(unsigned int);
  void ThdPGatVarSTPRK4(unsigned int);
  void ThdPGatVarSTDRK4(unsigned int);

  void ThdGapSyn(unsigned int);
  void ThdPSynapse(unsigned int);

  float dst(float st,float tst); //delta st
  float ndst(float st,float tst,float xt,float as); //NMDA receptor delta st
  float dD(float D, float tD);
  float dt;

  thread *th;
  unsigned int thsize;

  //random engine
  mt19937 *rdgen;
  unsigned int *rseeds;


  //LTP
  void firLP(unsigned int);
  void ltpEul(unsigned int);
  void ltpIpvEul(unsigned int);
  void ltpRK4(unsigned int);
  void LTP(unsigned int);

  //STP
  void STP(unsigned int);
  void STPIpvEul(unsigned int);
  void STPRK4(unsigned int);
  void STD(unsigned int);
  void STDIpvEul(unsigned int);
  void STDRK4(unsigned int);


  //membrane function
  float deltaV(float Vm,unsigned int);
  void MemPot(unsigned int);


  void EqsOne(unsigned int);

//-----membrane variable
  NeuronGnl2 *Nets;
  NeuRepPar2<RP1,RP2> *RPs;
  NeuRepVar2 *RVs;

  unsigned int NetSize;
  unsigned int RPsSize;
  unsigned int RVsSize;
  unsigned int RGsSize;
  

  ParBase *NetPBse;
  VarBase *NetVBse;

  ParLIF *NetPLIF;
  VarLIF2 *NetVLIF;

  ParStp *NetPSTP;
  VarStp *NetVSTP;

  ParLP *NetPLTP;
  VarLP2 *NetVLTP;

  unsigned int NetPBseSize;
  unsigned int NetVBseSize;

  unsigned int NetPLIFSize;
  unsigned int NetVLIFSize;

  unsigned int NetPSTPSize;
  unsigned int NetVSTPSize;

  unsigned int NetPLTPSize;
  unsigned int NetVLTPSize;

// back ground noise
  BgSti1<RV1,RV2> *BgSt;
  ziggurat *zig;

  
//-----single neuron variable
  SynFast2 *SynFst1;
  SynSlow2 *SynSlo1;
  SynGap2  *SynCnd1;

  unsigned int SynAMPA1End;
  unsigned int SynGABA1End;
  unsigned int SynACH1End;
  unsigned int SynGCL1End;

  unsigned int SynNMDA1Size;
  unsigned int SynCnd1Size;


//-----thread level parallelism
  unsigned int AMPA1CntSize;
  unsigned int GABA1CntSize;
  unsigned int ACH1CntSize;
  unsigned int GCL1CntSize;
  unsigned int NMDA1CntSize;
  unsigned int Cnd1CntSize;

  NeuMap *NeuMapAMPA1;
  NeuMap *NeuMapGABA1;
  NeuMap *NeuMapACH1;
  NeuMap *NeuMapGCL1;
  NeuMap *NeuMapNMDA1;
  NeuMap *NeuMapCnd1;

  unsigned int PreAMPA1CntSize;
  unsigned int PreGABA1CntSize;
  unsigned int PreACH1CntSize;
  unsigned int PreGCL1CntSize;
  unsigned int PreNMDA1CntSize;
  unsigned int PreCnd1CntSize;

  NeuMap *NeuPreMapAMPA1;
  NeuMap *NeuPreMapGABA1;
  NeuMap *NeuPreMapACH1;
  NeuMap *NeuPreMapGCL1;
  NeuMap *NeuPreMapNMDA1;
  NeuMap *NeuPreMapCnd1;

  unsigned int *PreMapAMPA1;
  unsigned int *PreMapGABA1;
  unsigned int *PreMapACH1;
  unsigned int *PreMapGCL1;
  unsigned int *PreMapNMDA1;
  unsigned int *PreMapCnd1;



//-----group model

  unsigned int *PSynA;
  unsigned int *PSynG;
  unsigned int *PSynH;
  unsigned int *PSynL;
  unsigned int *PSynN;
  unsigned int *PSynCp; //to PreSynapseID
  unsigned int *PSynCo; //to PosNeuronID

  unsigned int PSynASize;
  unsigned int PSynGSize;
  unsigned int PSynHSize;
  unsigned int PSynLSize;
  unsigned int PSynNSize;
  unsigned int PSynCpSize;

  SynFast2 *PSynFstGatA;  //AMPA
  SynFast2 *PSynFstGatG;  //GABA
  SynFast2 *PSynFstGatH;  //ACH
  SynFast2 *PSynFstGatL;  //GCL
  SynSlow2 *PSynSloGat;
  PSynGap *PSynCndP; //integrate PreID current to PosID

  unsigned int PGatASize;
  unsigned int PGatGSize;
  unsigned int PGatHSize;
  unsigned int PGatLSize;
  unsigned int PGatNSize;
  unsigned int PGapPSize;

  NeuMap *NeuMapPAMPA;
  NeuMap *NeuMapPGABA;
  NeuMap *NeuMapPACH;
  NeuMap *NeuMapPGCL;
  NeuMap *NeuMapPNMDA;

  NeuMap *NeuMapPCndP; //integrate PreID current to PosID
  NeuMap *NeuMapPCndO; //integrate PosID current to PreID

  unsigned int PAMPACntSize;
  unsigned int PGABACntSize;
  unsigned int PACHCntSize;
  unsigned int PGCLCntSize;
  unsigned int PNMDACntSize;
  unsigned int PCndPCntSize;
  unsigned int PCndOCntSize;

  unsigned int *MapPAPosID;
  unsigned int *MapPGPosID;
  unsigned int *MapPHPosID;
  unsigned int *MapPLPosID;
  unsigned int *MapPNPosID;
  unsigned int *MapPCPosID;

  template<class C>
  NeuMap *BulBiasDNeuMap(C *c, unsigned int CBeg, unsigned int CEnd, unsigned int &MSize)
  {
    if(CBeg!=CEnd)
    {
      unsigned int now=c[CBeg].PosID;
      NeuMap *nArray;
      MSize=0;

      for(unsigned int i=CBeg;i<CEnd;)
      {
        while(now==c[i].PosID && i<CEnd)
        { i++; }

        if(i<CEnd)
        {
          now=c[i].PosID;
          MSize++;
        }
      }
      MSize++;
      nArray = new NeuMap[MSize];

      now=c[CBeg].PosID;
      nArray[0].Beg=CBeg;
      MSize=0;
      for(unsigned int i=CBeg;i<CEnd;)
      {
        while(now==c[i].PosID && i<CEnd)
        { i++; }

        if(i<CEnd)
        {
          nArray[MSize].End=i-1;
          now=c[i].PosID;
          MSize++;
          nArray[MSize].Beg=i;
        }
      }
      nArray[MSize].End=CEnd-1;
      MSize++;
      return nArray;
    }
    else
    {return nullptr;}
  };

  template<class C>
  NeuMap *BulDPreNeuMap(C *c, unsigned int CBeg, unsigned int CEnd,unsigned int &MSize,unsigned int NeuSize,unsigned int *PreArray)
  {
    if(CBeg!=CEnd)
    {
      vector<NeuMap> Array;
      for(unsigned int i=0,n=0,m=0;i<NeuSize;i++)
      {
        m=n;
        for(unsigned int j=CBeg;j<CEnd;j++)
        {
          if(c[j].PreID==i)
          {
            PreArray[n]=j;
            n++;
          }
        }
        if(m!=n)
        {
          NeuMap NewNeuMap;
          NewNeuMap.Beg=m;
          NewNeuMap.End=n-1;
          Array.push_back(NewNeuMap);
        }
      }
      MSize=Array.size();
      NeuMap *nArray=new NeuMap[Array.size()];
      for(unsigned int i=0;i<Array.size();i++)
      {  nArray[i]=Array[i]; }

      { vector<NeuMap> ().swap(Array);}

      return nArray;
    }
    else
    {return nullptr;}
  };


//slovers
//-------------------pass 2 Variable
  float RK4FF(float num1, float num2, float (EqsGnl1_1::*FEqs)(float,float)) //Runge Kutta order 4th
  {
    float k1,k2,k3,k4;
    k1 = (this->*FEqs)( num1 , num2 );
    k2 = (this->*FEqs)( num1 + dt*k1/2.0f , num2 );
    k3 = (this->*FEqs)( num1 + dt*k2/2.0f , num2);
    k4 = (this->*FEqs)( num1 + dt*k3 , num2);
    return dt*( k1 + 2.0f*k2 + 2.0f*k3 + k4 )/6.0f;
  };

  float IpvEulFF(float num1, float num2, float (EqsGnl1_1::*FEqs)(float,float)) //improved Euler or RK2
  {
    float k1 = (this->*FEqs)( num1, num2 );
    float k2 = (this->*FEqs)( num1 + dt*k1, num2 );
    return dt*(k1+k2)/2.0f;
  };

//-------------------pass 2 Variable
  float RK4FuI(float num1, unsigned int num2, float (EqsGnl1_1::*FEqs)(float,unsigned int)) //Runge Kutta order 4th, Float, unsigned int
  {
    float k1,k2,k3,k4;
    k1 = (this->*FEqs)( num1 , num2 );
    k2 = (this->*FEqs)( num1 + dt*k1/2.0f , num2 );
    k3 = (this->*FEqs)( num1 + dt*k2/2.0f , num2);
    k4 = (this->*FEqs)( num1 + dt*k3 , num2);
    return ( k1 + 2.0f*k2 + 2.0f*k3 + k4 )/6.0f;
  };

  float IpvEulFuI(float num1, unsigned int num2, float (EqsGnl1_1::*FEqs)(float,unsigned int)) //improved Euler or RK2, Float, unsigned int
  {
    float k1 = (this->*FEqs)( num1, num2 );
    float k2 = (this->*FEqs)( num1 + dt*k1, num2 );
    return (k1+k2)/2.0f;
  };

//-------------------pass 4 Variable
  float RK4FFFF(float num1, float num2, float num3, float num4, float (EqsGnl1_1::*FEqs)(float,float,float,float)) //Runge Kutta order 4th, Float, Float, Float, Float
  {
    float k1,k2,k3,k4;
    k1 = (this->*FEqs)( num1 , num2 , num3 , num4);
    k2 = (this->*FEqs)( num1 + dt*k1/2.0f , num2 , num3 , num4 );
    k3 = (this->*FEqs)( num1 + dt*k2/2.0f , num2 , num3 , num4);
    k4 = (this->*FEqs)( num1 + dt*k3 , num2 , num3 , num4);
    return dt*( k1 + 2.0f*k2 + 2.0f*k3 + k4 )/6.0f;
  };

  float IpvEulFFFF(float num1,float num2,float num3,float num4, float (EqsGnl1_1::*FEqs)(float,float,float,float)) //improved Euler or RK2, Float, Float, Float, Float
  {
    float k1 = (this->*FEqs)( num1 , num2 , num3 , num4);
    float k2 = (this->*FEqs)( num1 + dt*k1 , num2 , num3 , num4);
    return dt*(k1+k2)/2.0f;
  };

};


class EqsGnl2
{
  public:
  EqsGnl2();
  int SolType;

  //void CntInit(void);
  //void Activation(unsigned int);
  void Activation(unsigned int,unsigned int);
  void firing(unsigned int);
  void Spike(unsigned int,unsigned int);
  //void BgInit(unsigned int);
  void BgInit(unsigned int,unsigned int);
  void VarInit(void);

  void ThdSynapse(unsigned int);
  void ThdSynapseSTP(unsigned int);
  void ThdSynapseSTD(unsigned int);
  void ThdGapSyn(unsigned int);

  void ThdPGatVar(unsigned int);
  void ThdPGatVarSTP(unsigned int);
  void ThdPGatVarSTD(unsigned int);
  void ThdPSynapse(unsigned int);

  float dt;

  thread *th;
  unsigned int thsize;

  //random engine
  mt19937 *rdgen;
  unsigned int *rseeds;


  //LTP
  void firLP(unsigned int);
  void ltpEul(unsigned int);
  void LTP(unsigned int);

  //STP
  void STP(unsigned int);
  void STD(unsigned int);

  //membrane function
  void VmEul(unsigned int);

  //neumerical method function
  float MemPot(float,unsigned int); // Irep and Ileaky

  void EqsOne(unsigned int);

//-----membrane variable
  NeuronGnl2 *Nets;
  NeuRepPar2<RP3,RP4> *RPs;
  NeuRepVar2 *RVs;

  unsigned int NetSize;
  unsigned int RPsSize;
  unsigned int RVsSize;

  ParBase2 *NetPBse;
  VarBase *NetVBse;

  ParLIF *NetPLIF;
  VarLIF2 *NetVLIF;

  ParStp2 *NetPSTP;
  VarStp *NetVSTP;

  ParLP2 *NetPLTP;
  VarLP2 *NetVLTP;

  unsigned int NetPBseSize;
  unsigned int NetVBseSize;

  unsigned int NetPLIFSize;
  unsigned int NetVLIFSize;

  unsigned int NetPSTPSize;
  unsigned int NetVSTPSize;

  unsigned int NetPLTPSize;
  unsigned int NetVLTPSize;



// back ground noise
  //BgSti2<RV3,RV4> *BgSt;

  BgSti1<RV3,RV4> *BgSt;
  ziggurat *zig;

//-----single neuron variable
  SynFast2 *SynFst1;
  SynSlow2 *SynSlo1;
  SynGap2  *SynCnd1;

  unsigned int SynAMPA1End;
  unsigned int SynGABA1End;
  unsigned int SynACH1End;
  unsigned int SynGCL1End;

  unsigned int SynNMDA1Size;
  unsigned int SynCnd1Size;


//-----thread level parallelism
  unsigned int AMPA1CntSize;
  unsigned int GABA1CntSize;
  unsigned int ACH1CntSize;
  unsigned int GCL1CntSize;
  unsigned int NMDA1CntSize;
  unsigned int Cnd1CntSize;

  NeuMap *NeuMapAMPA1;
  NeuMap *NeuMapGABA1;
  NeuMap *NeuMapACH1;
  NeuMap *NeuMapGCL1;
  NeuMap *NeuMapNMDA1;
  NeuMap *NeuMapCnd1;

  unsigned int PreAMPA1CntSize;
  unsigned int PreGABA1CntSize;
  unsigned int PreACH1CntSize;
  unsigned int PreGCL1CntSize;
  unsigned int PreNMDA1CntSize;
  unsigned int PreCnd1CntSize;

  NeuMap *NeuPreMapAMPA1;
  NeuMap *NeuPreMapGABA1;
  NeuMap *NeuPreMapACH1;
  NeuMap *NeuPreMapGCL1;
  NeuMap *NeuPreMapNMDA1;
  NeuMap *NeuPreMapCnd1;

  unsigned int *PreMapAMPA1;
  unsigned int *PreMapGABA1;
  unsigned int *PreMapACH1;
  unsigned int *PreMapGCL1;
  unsigned int *PreMapNMDA1;
  unsigned int *PreMapCnd1;



//-----group model

  unsigned int *PSynA;
  unsigned int *PSynG;
  unsigned int *PSynH;
  unsigned int *PSynL;
  unsigned int *PSynN;
  unsigned int *PSynCp; //to PreSynapseID
  unsigned int *PSynCo; //to PosNeuronID

  unsigned int PSynASize;
  unsigned int PSynGSize;
  unsigned int PSynHSize;
  unsigned int PSynLSize;
  unsigned int PSynNSize;
  unsigned int PSynCpSize;

  SynFast2 *PSynFstGatA;  //AMPA
  SynFast2 *PSynFstGatG;  //GABA
  SynFast2 *PSynFstGatH;  //ACH
  SynFast2 *PSynFstGatL;  //GCL
  SynSlow2 *PSynSloGat;
  PSynGap *PSynCndP; //integrate PreID current to PosID

  unsigned int PGatASize;
  unsigned int PGatGSize;
  unsigned int PGatHSize;
  unsigned int PGatLSize;
  unsigned int PGatNSize;
  unsigned int PGapPSize;

  NeuMap *NeuMapPAMPA;
  NeuMap *NeuMapPGABA;
  NeuMap *NeuMapPACH;
  NeuMap *NeuMapPGCL;
  NeuMap *NeuMapPNMDA;

  NeuMap *NeuMapPCndP; //integrate PreID current to PosID
  NeuMap *NeuMapPCndO; //integrate PosID current to PreID

  unsigned int PAMPACntSize;
  unsigned int PGABACntSize;
  unsigned int PACHCntSize;
  unsigned int PGCLCntSize;
  unsigned int PNMDACntSize;
  unsigned int PCndPCntSize;
  unsigned int PCndOCntSize;

  unsigned int *MapPAPosID;
  unsigned int *MapPGPosID;
  unsigned int *MapPHPosID;
  unsigned int *MapPLPosID;
  unsigned int *MapPNPosID;
  unsigned int *MapPCPosID;

  template<class C>
  NeuMap *BulBiasDNeuMap(C *c, unsigned int CBeg, unsigned int CEnd, unsigned int &MSize)
  {
    if(CBeg!=CEnd)
    {
      unsigned int now=c[CBeg].PosID;
      NeuMap *nArray;
      MSize=0;

      for(unsigned int i=CBeg;i<CEnd;)
      {
        while(now==c[i].PosID && i<CEnd)
        { i++; }

        if(i<CEnd)
        {
          now=c[i].PosID;
          MSize++;
        }
      }
      MSize++;
      nArray = new NeuMap[MSize];

      now=c[CBeg].PosID;
      nArray[0].Beg=CBeg;
      MSize=0;
      for(unsigned int i=CBeg;i<CEnd;)
      {
        while(now==c[i].PosID && i<CEnd)
        { i++; }

        if(i<CEnd)
        {
          nArray[MSize].End=i-1;
          now=c[i].PosID;
          MSize++;
          nArray[MSize].Beg=i;
        }
      }
      nArray[MSize].End=CEnd-1;
      MSize++;
      return nArray;
    }
    else
    {return nullptr;}
  };

  template<class C>
  NeuMap *BulDPreNeuMap(C *c, unsigned int CBeg, unsigned int CEnd,unsigned int &MSize,unsigned int NeuSize,unsigned int *PreArray)
  {
    if(CBeg!=CEnd)
    {
      vector<NeuMap> Array;
      for(unsigned int i=0,n=0,m=0;i<NeuSize;i++)
      {
        m=n;
        for(unsigned int j=CBeg;j<CEnd;j++)
        {
          if(c[j].PreID==i)
          {
            PreArray[n]=j;
            n++;
          }
        }
        if(m!=n)
        {
          NeuMap NewNeuMap;
          NewNeuMap.Beg=m;
          NewNeuMap.End=n-1;
          Array.push_back(NewNeuMap);
        }
      }
      MSize=Array.size();
      NeuMap *nArray=new NeuMap[Array.size()];
      for(unsigned int i=0;i<Array.size();i++)
      {  nArray[i]=Array[i]; }

      { vector<NeuMap> ().swap(Array);}

      return nArray;
    }
    else
    {return nullptr;}
  };

};





class EqsGnl2_1
{
  public:
  EqsGnl2_1();
  int SolType;

  //void CntInit(void);
  void Activation(unsigned int,unsigned int);
  void firing(unsigned int);
  void firingHH(unsigned int);
  void Spike(unsigned int,unsigned int);
  void BgInit(unsigned int,unsigned int);
  void VarInit(void);

  void ThdSynapse(unsigned int);
  void ThdSynapseSTP(unsigned int);
  void ThdSynapseSTD(unsigned int);
  void ThdGapSyn(unsigned int);

  float dt;

  thread *th;
  unsigned int thsize;

  //random engine
  mt19937 *rdgen;
  unsigned int *rseeds;


  //LTP
  void firLP(unsigned int);
  void ltpEul(unsigned int);
  void LTP(unsigned int);

  //STP
  void STP(unsigned int);
  void STD(unsigned int);

  //membrane function
  void VmEul(unsigned int);
  void VmSOD(unsigned int);
  void VmHH(unsigned int);

  //neumerical method function
  float MemPot(float,unsigned int); // Irep and Ileaky

  void EqsOne(unsigned int);

//-----membrane variable
  NeuronGnl2 *Nets;
  NeuRepPar2<RP3,RP4> *RPs;
  NeuRepVar2 *RVs;
  NeuGatVar3 *RGs;

  unsigned int NetSize;
  unsigned int RPsSize;
  unsigned int RVsSize;
  unsigned int RGsSize;
  

  ParBase2 *NetPBse;
  VarBase *NetVBse;

  ParLIF *NetPLIF;
  VarLIF2 *NetVLIF;

  ParStp2 *NetPSTP;
  VarStp *NetVSTP;

  ParLP2 *NetPLTP;
  VarLP2 *NetVLTP;

  ParHH2 *NetPHH;
  VarHH *NetVHH;

  unsigned int NetPBseSize;
  unsigned int NetVBseSize;

  unsigned int NetPLIFSize;
  unsigned int NetVLIFSize;

  unsigned int NetPSTPSize;
  unsigned int NetVSTPSize;

  unsigned int NetPLTPSize;
  unsigned int NetVLTPSize;

// back ground noise
  BgSti1<float,RV4> *BgSt;
  ziggurat *zig;

  
//-----single neuron variable
  SynFast2 *SynFst1;
  SynSlow2 *SynSlo1;
  SynGap2  *SynCnd1;

  unsigned int SynAMPA1End;
  unsigned int SynGABA1End;
  unsigned int SynACH1End;
  unsigned int SynGCL1End;

  unsigned int SynNMDA1Size;
  unsigned int SynCnd1Size;


//-----thread level parallelism
  unsigned int AMPA1CntSize;
  unsigned int GABA1CntSize;
  unsigned int ACH1CntSize;
  unsigned int GCL1CntSize;
  unsigned int NMDA1CntSize;
  unsigned int Cnd1CntSize;

  NeuMap *NeuMapAMPA1;
  NeuMap *NeuMapGABA1;
  NeuMap *NeuMapACH1;
  NeuMap *NeuMapGCL1;
  NeuMap *NeuMapNMDA1;
  NeuMap *NeuMapCnd1;

  unsigned int PreAMPA1CntSize;
  unsigned int PreGABA1CntSize;
  unsigned int PreACH1CntSize;
  unsigned int PreGCL1CntSize;
  unsigned int PreNMDA1CntSize;
  unsigned int PreCnd1CntSize;

  NeuMap *NeuPreMapAMPA1;
  NeuMap *NeuPreMapGABA1;
  NeuMap *NeuPreMapACH1;
  NeuMap *NeuPreMapGCL1;
  NeuMap *NeuPreMapNMDA1;
  NeuMap *NeuPreMapCnd1;

  unsigned int *PreMapAMPA1;
  unsigned int *PreMapGABA1;
  unsigned int *PreMapACH1;
  unsigned int *PreMapGCL1;
  unsigned int *PreMapNMDA1;
  unsigned int *PreMapCnd1;


  template<class C>
  NeuMap *BulBiasDNeuMap(C *c, unsigned int CBeg, unsigned int CEnd, unsigned int &MSize)
  {
    if(CBeg!=CEnd)
    {
      unsigned int now=c[CBeg].PosID;
      NeuMap *nArray;
      MSize=0;

      for(unsigned int i=CBeg;i<CEnd;)
      {
        while(now==c[i].PosID && i<CEnd)
        { i++; }

        if(i<CEnd)
        {
          now=c[i].PosID;
          MSize++;
        }
      }


      MSize++;
      nArray = new NeuMap[MSize];

      now=c[CBeg].PosID;
      nArray[0].Beg=CBeg;
      MSize=0;
      for(unsigned int i=CBeg;i<CEnd;)
      {
        while(now==c[i].PosID && i<CEnd)
        { i++; }

        if(i<CEnd)
        {
          nArray[MSize].End=i-1;
          now=c[i].PosID;
          MSize++;
          nArray[MSize].Beg=i;
        }
      }
      nArray[MSize].End=CEnd-1;
      MSize++;
      return nArray;
    }
    else
    {return nullptr;}
  };

  template<class C>
  NeuMap *BulDPreNeuMap(C *c, unsigned int CBeg, unsigned int CEnd,unsigned int &MSize,unsigned int NeuSize,unsigned int *PreArray)
  {
    if(CBeg!=CEnd)
    {
      vector<NeuMap> Array;
      for(unsigned int i=0,n=0,m=0;i<NeuSize;i++)
      {
        m=n;
        for(unsigned int j=CBeg;j<CEnd;j++)
        {
          if(c[j].PreID==i)
          {
            PreArray[n]=j;
            n++;
          }
        }
        if(m!=n)
        {
          NeuMap NewNeuMap;
          NewNeuMap.Beg=m;
          NewNeuMap.End=n-1;
          Array.push_back(NewNeuMap);
        }
      }
      MSize=Array.size();
      NeuMap *nArray=new NeuMap[Array.size()];
      for(unsigned int i=0;i<Array.size();i++)
      {  nArray[i]=Array[i]; }

      { vector<NeuMap> ().swap(Array);}

      return nArray;
    }
    else
    {return nullptr;}
  };

};










