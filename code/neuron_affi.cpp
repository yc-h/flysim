#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>
#include "neuron_affi.h"

using namespace std;

ConTab::ConTab()
{

  NeuID=0;

  AMPA_ID=0;
  GABA_ID=0;
  NMDA_ID=0;

  AMPA_Act=false;
  GABA_Act=false;
  NMDA_Act=false;

}

ConTab_gp::ConTab_gp()
{

  NeuID=0;

  AMPA_MeanG=0.0f;
  GABA_MeanG=0.0f;
  NMDA_MeanG=0.0f;

  AMPA_Act=false;
  GABA_Act=false;
  NMDA_Act=false;

}

Nlink::Nlink()
{
  AMPA_st=0.0f;
  GABA_st=0.0f;
  NMDA_xt=0.0f;
  NMDA_st=0.0f;

}


ConTab_gp2::ConTab_gp2()
{
  NeuID=0;

  AMPA_MeanG=0.0f;
  GABA_MeanG=0.0f;
  NMDA_MeanG=0.0f;

  //NMDA_as=1.0f;
  NMDA_as=0.6332f;

  AMPA_tst=2.0e-3f;
  GABA_tst=5.0e-3f;
  NMDA_txt=2.0e-3f;
  NMDA_tst=100.0e-3f;

  AMPA_Act=false;
  GABA_Act=false;
  NMDA_Act=false;

  GAP_MeanG=0.0f;  // pre-synapse to post-synapse connection mean efficacy
  GAP_Act=false; //GAP junction active flag, true for active
}
//-----------------------------------------------------------------------
NeuPar::NeuPar()
{
// membrance parameter
  Vth=-50.0e-3f; //-50mV
  El=-70.0e-3f; //-70mV
  gl=2.5e-9f;
  Cm=1.0e-9f;

  Vreset=VRESET;  // mV

  RfPed=19;  // 1.8ms = (19(total) - 1(firing))*0.1m
  SpiDy=0;  // spike is delay 18 DTIME = 1.8ms
}

NeuPar_stp::NeuPar_stp()
{
// membrance parameter
  Vth=-50.0e-3f; //-50mV
  El=-70.0e-3f; //-70mV
  gl=2.5e-9f;
  Cm=1.0e-9f;

  Vreset=VRESET;  // mV

  RfPed=19;  // 1.8ms = (19(total) - 1(firing))*0.1m
  SpiDy=0;  // spike is delay 18 DTIME = 1.8ms

// short-term plasticity
  pv=0.6f; // reduce factor of depression factor D
  tD=300.0e-3f; // decay time constant D

  aF=0.0035e-9f; // the alpha of facilitation F
  tF=7000.0e-3f; // the decay time constant of F

  aCap=250.0e-6f; // the alpha of Ca^2+ peak
  tCap=2.0e-3f; // decay time constant  Ca^2+ peak

  aCar=0.1e-6f; // alpha of residual Ca^2+
  tCar=2000.0e-3f; // decay time constant residual Ca^2+
}

NeuPar_HH::NeuPar_HH()
{
// membrance parameter
  Vth=-50.0e-3f; //-50mV
  El=-60.6e-3f; //-70mV
  gl=2.5e-9f;
  Cm=1.0e-9f;

  Vreset=VRESET;  // mV

  RfPed=19;  // 1.8ms = (19(total) - 1(firing))*0.1m
  SpiDy=0;  // spike is delay 18 DTIME = 1.8ms

// Hodgkin Huxley model
  Ek=-82.0e-3f;
  Ena=45.0e-3f;

  gk=36.0e-6f;
  gna=120.0e-6f;
}

NeuPar_HST::NeuPar_HST()
{
// membrance parameter
  Vth=-50.0e-3f; //-50mV
  El=-70.0e-3f; //-70mV
  gl=2.5e-9f;
  Cm=1.0e-9f;

  Vreset=VRESET;  // mV

  RfPed=19;  // 1.8ms = (19(total) - 1(firing))*0.1m
  SpiDy=0;  // spike is delay 18 DTIME = 1.8ms

// Hodgkin Huxley model
  Ek=-78.0e-3f; //-90mV
  Ena=60.0e-3f; //65mV

  gk=36.0e-9f;
  gna=120.0e-9f;

// short-term plasticity
  pv=0.6f; // reduce factor of depression factor D
  tD=300.0e-3f; // decay time constant D

  aF=0.0035e-9f; // the alpha of facilitation F
  tF=7000.0e-3f; // the decay time constant of F

  aCap=250.0e-6f; // the alpha of Ca^2+ peak
  tCap=2.0e-3f; // decay time constant  Ca^2+ peak

  aCar=0.1e-6f; // alpha of residual Ca^2+
  tCar=2000.0e-3f; // decay time constant residual Ca^2+
}


NeuPar_bp::NeuPar_bp()
{
// membrance parameter
  Vth=-50.0e-3f; //-50mV
  El=-70.0e-3f; //-70mV
  gl=2.5e-9f;
  Cm=1.0e-9f;

  Vreset=VRESET;  // mV

  RfPed=19;  // 1.8ms = (19(total) - 1(firing))*0.1m
  SpiDy=0;  // spike is delay 18 DTIME = 1.8ms


// dendritic part
  tx=0.3f; // taum factor

// these are about 0.2ms delay for dendrite to soma
  ax=1.0e3f;
  vx=1.0e-3f;
}

NeuPar_uni::NeuPar_uni()
{
// membrance parameter
  Vth=-50.0e-3f; //-50mV
  El=-70.0e-3f; //-70mV
  gl=2.5e-9f;
  Cm=1.0e-9f;

  Vreset=VRESET;  // mV

  RfPed=19;  // 1.8ms = (19(total) - 1(firing))*0.1m
  SpiDy=0;  // spike is delay 18 DTIME = 1.8ms

//dendrite part
  tx=0.3f; // taum factor
// these are about 0.2ms delay for dendrite to soma
  ax=1.0e3f;
  vx=1.0e-3f;

//axon part
  ta=0.3f; // taum factor
// these are about 0.2ms delay for dendrite to soma
  aa=1.0e3f;
  va=1.0e-3f;

}


NeuPar_gp::NeuPar_gp()
{
// membrance parameter
  Vth=-50.0e-3f; //-50mV
  El=-70.0e-3f; //-70mV
  gl=2.5e-9f;
  Cm=1.0e-9f;

  Vreset=VRESET;  // mV

  RfPed=19;  // 1.8ms = (19(total) - 1(firing))*0.1m
  SpiDy=0;  // spike is delay 18 DTIME = 1.8ms

  SelfSti=false;
}

//-----------------------------------------------------------------------
NeuVar::NeuVar()
{
  Isyn=0.0f;  // current of synapse
  V=VRESTING*1e-3f;  // resting menbrane potential

  Count=0;  // firing tabe index
  Firing=false;  // firing flag, true for firing
  SpiOut=false;  // spike of this neuron, true for spike
}

NeuVar_stp::NeuVar_stp()
{
  //NeuVar();
  Isyn=0.0f;  // current of synapse
  V=VRESTING*1e-3f;  // resting menbrane potential

  Count=0;  // firing tabe index
  Firing=false;  // firing flag, true for firing
  SpiOut=false;  // spike of this neuron, true for spike

  //VarStp();
  D=1.0f;  // depression factor D
  F=1.0f;  // facilitation factor F
  Cap=0.0e-6f; // Ca^2+ peak level
  Car=20.0e-6f;  // residual level of Ca^2+
  //Ca=0.0f;  // Ca^2+ level
}

NeuVar_HH::NeuVar_HH()
{
  //NeuVar();
  Isyn=0.0f;  // current of synapse
  V=VRESTING*1e-3f;  // resting menbrane potential

  Count=0;  // firing tabe index
  Firing=false;  // firing flag, true for firing
  SpiOut=false;  // spike of this neuron, true for spike

  //VarHH();
  n=0.320392546f;
  m=0.0534988507f;
  h=0.596120754f;
}

NeuVar_HST::NeuVar_HST()
{
  //NeuVar();
  Isyn=0;  // current of synapse
  V=VRESTING*1e-3f;  // resting menbrane potential

  Count=0;  // firing tabe index
  Firing=false;  // firing flag, true for firing
  SpiOut=false;  // spike of this neuron, true for spike

  //VarHH();
  n=0.320392546f;
  m=0.0534988507f;
  h=0.596120754f;

  //VarStp();
  D=1.0f;  // depression factor D
  F=1.0f;  // facilitation factor F
  Cap=0.0e-6f; // Ca^2+ peak level
  Car=20.0e-6f;  // residual level of Ca^2+
  //Ca=0.0f;  // Ca^2+ level
}

NeuVar_bp::NeuVar_bp()
{
  Isyn=0.0f;  // current of synapse
  V=VRESTING*1e-3f;  // resting menbrane potential

  Count=0;  // firing tabe index
  Firing=false;  // firing flag, true for firing
  SpiOut=false;  // spike of this neuron, true for spike
  
  Vd=VRESTING*1e-3f;
}

NeuVar_uni::NeuVar_uni()
{
  Isyn=0.0f;  // current of synapse
  V=VRESTING*1e-3f;  // resting menbrane potential

  Count=0;  // firing tabe index
  Firing=false;  // firing flag, true for firing
  SpiOut=false;  // spike of this neuron, true for spike
  
  Vd=VRESTING*1e-3f;
  Va=VRESTING*1e-3f;
}


