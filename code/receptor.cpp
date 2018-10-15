#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <omp.h>
#include "receptor.h"


using namespace std;

// Receptor1 parameter for AMPA and GABA

RepPar1::RepPar1()
{
  //as=10.0f; // the alpha st
  tst=2.0e-3f; //decay time constant for st
  Vrev=0.0f; // reversal potential
}

// Receptor2 parameter for NMDA

RepPar2::RepPar2()
{
  as=0.6332f; // the alpha st
  txt=2.0e-3f; // the decay time constant of xt
  tst=100.0e-3f; //decay time constant for st
  Vrev=0.0e-3f; // reversal potential

  mg=1.0f; // the concentration of magnesium
}


RepVar1::RepVar1()
{
  g=0.11e-9f; // synaptic conductance
  Sti=false;
  st=0.0f;
}

RepVar2::RepVar2()
{
  g=0.11e-9f; // synaptic conductance
  Sti=false;
  st=0.0f;
  xt=0.0f;
}

RepVar1_stp::RepVar1_stp()
{
  g=0.11e-9f; // synaptic conductance
  Sti=false;
  st=0.0f;
  ax=1.0f; // the alpha xt
}

RepVar2_stp::RepVar2_stp()
{
  g=0.11e-9f; // synaptic conductance
  Sti=false;
  st=0.0f;
  xt=0.0f;
  ax=1.0f; // the alpha xt
}

//---------group model------------------------------------

RepPar1_gp::RepPar1_gp()
{
  tst=2.0e-3f; //decay time constant for st
  Vrev=0.0f; // reversal potential
  GTotal=0.0f;
  g=0.0f;
  MeanG=0.0f;
}

RepPar2_gp::RepPar2_gp()
{
  as=0.6332f; // the alpha st
  txt=2.0e-3f; // the decay time constant of xt
  tst=100.0e-3f; //decay time constant for st
  Vrev=0.0e-3f; // reversal potential

  mg=1.0f; // the concentration of magnesium
  GTotal=0.0f;
  g=0.0f;
  MeanG=0.0f;
}

RepPar1_gp2::RepPar1_gp2()
{
  tst=2.0e-3f; //decay time constant for st
  Vrev=0.0f; // reversal potential
  g=0.0f;
  MeanG=0.0f;
}

RepPar2_gp2::RepPar2_gp2()
{
  as=0.6332f; // the alpha st
  txt=2.0e-3f; // the decay time constant of xt
  tst=100.0e-3f; //decay time constant for st
  Vrev=0.0e-3f; // reversal potential

  mg=1.0f; // the concentration of magnesium
  g=0.0f;
  MeanG=0.0f;
}

RepVar1_gp::RepVar1_gp()
{ st=0.0f; }

RepVar2_gp::RepVar2_gp()
{
  st=0.0f;
  xt=0.0f;
}



