
#define _USE_MATH_DEFINES

// c++ v11 standard support
#define CPP11

//rabbitMQ
//#define RABBITMQ

//LEG worm SOCKET type
//#define _WINSOCKET
#define _LINSOCKET


// solver support
#define ACCURATE 0
#define MODERATE 1
#define ROUGH    2

// parameters
#define DTIME 0.0001f // time unit is 0.1 ms
#define OneMilliSecond 10 // 1/DTIME = steps in 1 ms
#define ReportClock 1000 // 100/DTIME = steps in every 100ms report on time

#define VRESTING -70.0f  // mV
#define VRESET -50.0f  // mV
#define CarBase 20.0f  // uM

#define MNoiseMean 0.0f
#define MNoiseSTD 0.0f
#define FExtSinAmp 100*1e-12f

#define MIN_TAU 2*dt
//parser keyword 0~255
#define KW_AMPA 0
#define KW_GABA 1
#define KW_NMDA 2
#define KW_GAP  3
#define KW_CUBA 4
#define KW_SINE 5
#define RSM_SYNCOND 6
#define RSM_GATVAR 7
#define RSM_MEMPOT 8
#define RSM_RNDSED 9
#define RSM_ALL 10
#define KW_ACH 11
#define KW_GCL 12

// main.cpp FSM
#define STATE_SimName 1
#define STATE_command 2
#define STATE_help 3
#define STATE_pro 4
#define STATE_conf 5
#define STATE_thread 6
#define STATE_OutMem 7
#define STATE_OutSpi 8
#define STATE_OutFR 9
#define STATE_SolverType 10
#define STATE_Repeat 11
#define STATE_NeuModel 12
#define STATE_Daemon 13
#define STATE_DTIME 14
#define STATE_UDFRSEED 15

//---------------------------network.conf FSM
#define STATE_Initial 0
#define STATE_NeuPopuName 1
#define STATE_NeuPopu 2
#define STATE_N 3

#define STATE_Taum 4
#define STATE_Vth 5
#define STATE_gl 7
#define STATE_Cm 8
#define STATE_Vresting 10
#define STATE_Vreset 11
#define STATE_RfPed 12
#define STATE_SpiDy 13
#define STATE_SelfCon 14

#define STATE_HH_AnA 15
#define STATE_HH_AnB 16
#define STATE_HH_AnC 17
#define STATE_HH_AnD 18
#define STATE_HH_BnA 19
#define STATE_HH_BnB 20
#define STATE_HH_BnC 21
#define STATE_HH_BnD 22

#define STATE_HH_AmA 23
#define STATE_HH_AmB 24
#define STATE_HH_AmC 25
#define STATE_HH_AmD 26
#define STATE_HH_BmA 27
#define STATE_HH_BmB 28
#define STATE_HH_BmC 29
#define STATE_HH_BmD 30

#define STATE_HH_AhA 31
#define STATE_HH_AhB 32
#define STATE_HH_AhC 33
#define STATE_HH_AhD 34
#define STATE_HH_BhA 35
#define STATE_HH_BhB 36
#define STATE_HH_BhC 37
#define STATE_HH_BhD 38

#define STATE_HH_Ek 39
#define STATE_HH_gk 40
#define STATE_HH_Ena 41
#define STATE_HH_gna 42

#define STATE_STP_pv 43
#define STATE_STP_tD 44
#define STATE_STP_aF 45
#define STATE_STP_tF 46
#define STATE_STP_aCap 47
#define STATE_STP_tCap 48
#define STATE_STP_aCar 49
#define STATE_STP_tCar 50

#define STATE_LTP_tLP 51
#define STATE_LTP_SpiLPDy 52
#define STATE_LTP_PosISI 53
#define STATE_LTP_NegISI 54

#define STATE_RepName 55
#define STATE_Receptor 56
#define STATE_RepTauST 57
#define STATE_RepTauXT 58
#define STATE_RepVrv 59
#define STATE_RepFExt 60
#define STATE_RepMnExtEff 61
#define STATE_RepMnExtCon 62

#define STATE_NeuConName 63
#define STATE_NeuCon 64
#define STATE_NeuCon_Conct 65
#define STATE_NeuCon_TarRep 66
#define STATE_NeuCon_MnEff 67
#define STATE_NeuCon_weight 68


#define STATE_VthSod 69


//------------------network.pro FSM

#define STATE_TimEv 65
#define STATE_TimEvDef 66
#define STATE_TimTyp 67
#define STATE_TimLbl 68
#define STATE_TimPp 69
#define STATE_TimSubTyp 70
#define STATE_TimVar1 71
#define STATE_TimVar2 72

#define STATE_OutCtr 75
#define STATE_OutName 76
#define STATE_OutType 77
#define STATE_OutRWin 78
#define STATE_OutFPnt 79
#define STATE_OutTimEv 80
#define STATE_OutAllPopu 81
#define STATE_OutPopu 82
#define STATE_OutSumAllPp 83
#define STATE_OutAlign 84

#define STATE_DefMac 85
#define STATE_MacGpName 86
#define STATE_MacGpMbs 87
#define STATE_MacRndSeds 88

// Events
#define EVTEXTFREQ 1
#define EVTCUNTINJ 2

//output file max tpye
#define FILEMAX 4

// output file type ID
#define ReadFromPro 0
#define MemPotID 1
#define SpikeID 2
#define FRateID 3
#define SpiMemPot 4
#define GatVarID 5
#define IsynSepSumID 6
#define SaveID 7
#define ResumeID 8
#define SynWgtID 9
#define dVOvfID 10

//output control flag
#define ALIGN_OF 1
#define SUMPARSALL_OF 2





#define MODE 2

//neuron models in main.cpp
#define SIMPLE 1
#define LIF 2
#define LIF2 3
#define LIF3 4
#define LIF_STP 5
#define LIF_BP 6
#define LIF_UNI 7
#define HH 8
#define LIF_GP 9
#define LIF_GP2 10
#define LIF_GP3 11
#define LIF_GP4 12
#define LIF_LTP 13
#define GNL 14
#define LIF4 15
#define GNL2 16
#define GNL3 17

//neuron states
#define FIRE_F 1 // firing flag, 1 for firing
#define SPIKE_F 2 // spike flag, 1 for spiking
#define LP_FIRE_F 4 // LTP firing flag, 1 for firing
#define LP_SPIKE_F 8 // LTP spike flag, 1 for spiking
#define MEMDV_OVF 16 // membrane current over limit flag, 1 for overflow
#define MEMDV_UDF 32 // membrane current under limit flag, 1 for underflow
#define INJ_EN 64 // injection current enable



//template abbreviation
#define EqsTmp template<class RepP,class RepV, class RepP2,class RepV2, class Neu>
//#define EqsGen Eqs< RepP, RepV, RepP2, RepV2, Neu >
#define EqsOMPGen EqsOMP< RepP, RepV, RepP2, RepV2, Neu >
#define ItrCUBA EqsItr<RepP, RepV, RepP2, RepV2, Neu>

#define EqsNewGen EqsOMP<EqsType>

#define NeuDatTmp template<class RepP,class RepV, class RepP2,class RepV2,class ConT>
#define NeuDatGen NeuDatBase<RepP, RepV, RepP2, RepV2, ConT>
#define NeuBaseGen NeuBase<RepP, RepV, RepP2, RepV2, ConT>

#define NeuGPBaseGen NeuGPBase<RepP, RepV, RepP2, RepV2, ConT>

#define NeuParTmp template<class RepP, class RepP2>

//Active flag mask
//receptor
#define AMPA_ACT 1
#define GABA_ACT 2
#define NMDA_ACT 4
#define GAP_ACT 8
#define CUBA_ACT 16
#define ACH_ACT 32
#define GCL_ACT 64


// injection current
#define INJ_ACT 128

//membrane equation
#define LIF_ACT 256
#define HH_ACT 512
#define FZNA_ACT 1024
#define IZH_ACT 2048
#define SOD_ACT 4096
#define EUL_ACT 131072

//synpse plasticity
#define LTP_ACT 8192
#define STD_ACT 16384
#define STP_ACT 32768
#define AXA_ACT 65536

//user define random seed
#define UDFRDSEED_ACT 16384



//index inset FSM
#define IdxIns 1
#define IdxFiPre 2
#define IdxRe 3

//Macros
#define SetSti1() VItr->st = VItr->st + 1.0f
#define SetSti2() VItr2->xt = VItr2->xt + 1.0f
#define SetSti1stp() VItr->st = VItr->st + VItr->ax
#define SetSti2stp() VItr2->xt = VItr2->xt + VItr2->ax
#define SetSti4() VItr->st = VItr->st + PItr2->as*(1.0f - VItr->st)

#define PreSynR (ConfSta[i].PpRep+k)
#define PoSyn ((Confs+i)->PpLink+j)
#define PoSynPp PoSyn->TarPpName
#define IdxBgn IdxNum[j]-ConfSta[PoSynPp].NeuNum

#define POSYN Confs[i].PpLink[j]
#define POSYNPP PostPp[j].TarPpName
#define POSYNREP PostPp[j].TarRep
#define POSYNCON PostPp[j].Connectivity
#define POSYNEFF PostPp[j].MeanEff
#define POSYNWGT PostPp[j].weight

#define PreNID (NItr_begin + SynItr->PreID)
#define PosNID (NItr_begin + SynItr->PosID)

// LEGO worm 
#define WEvts Worms[WID]->Evts[i]
#define WNets Worms[WID]->Nets[m]
#define WBGSTI Worms[WID]->BgSt[m]
#define WNVL Worms[WID]->NetVLIF[m]




//equations
#define HHA(A,B,C,D,V) A*(B-V)/(C+(float)exp((B-V)/D))
#define HHB(A,B,C,D,V) A*(C+(float)exp((B-V)/D))
#define HHC(A,B,C,D,V) A/(C+(float)exp((B-V)/D))
#define HHx(A,B,dt,x) (x*exp(-(A+B)*dt)+(1.0f-exp(-(A+B)*dt))*A/(A+B))
#define HHdx(A,B,x) (A*(1-x)+B*x)
#define HHx3(A,B,dt,x) ((x*(1.0f/dt-(A+B)/2.0f)+A)/(1.0f/dt+(A+B)/2.0f))
#define HHEul(A,B,dt,x) (x+dt*HHdx(A,B,x))
#define Ina(gna,m,h,Ena,V) gna*m*m*m*h*(Ena-V)
#define Ik(gk,n,Ek,V) gk*n*n*n*n*(Ek-V)

#define I_LIF(V) (NRV.gAMPA*(NRP.ParAMPA.Vrev-V)+\
NRV.gGABA*(NRP.ParGABA.Vrev-V)+\
NRV.gNMDA*(NRP.ParNMDA.Vrev-V)/( 1.0f + NRP.ParNMDA.mg * (float)exp( -0.062f * V ) / 3.57f )+\
NRV.gACH*(NRP.ParACH.Vrev-V)+\
NRV.gGCL*(NRP.ParGCL.Vrev-V))*1e-3f+\
NetVBse[NID].Isyn+NPB.gl*(NPB.El-V)

inline
float expApx(float x) {
  x = 1.0f + x / 1024.0f;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x; x *= x; x *= x;
  x *= x; x *= x;
  return x;
};



//receptor backward capability factor: GNL2 ==> GNL1
//#define REP_CAP_FACTOR 1.25714f


