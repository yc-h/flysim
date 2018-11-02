#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <omp.h>
#include <random>
#include <chrono>
#include "neuron.h"

using namespace std;

void NeuronGnl2::TexInfo(string& ostr)
{
/*
  ostringstream ostr_tmp;
  string str;
//----------------------------noise--------------------------------------
  ostr_tmp<<"noise parameters:"<<"\n";

  ostr_tmp
    <<"Injection current:"
    <<"mean="<<Inj.mean()
    <<" stddev="<<Inj.stddev()
    <<"\n";

  ostr_tmp<<"StiExt AMPA:"<<*pb<<"\n";
  ostr_tmp<<"StiExt GABA:"<<*(pb+1)<<"\n";
  ostr_tmp<<"StiExt NMDA:"<<*(pb+2)<<"\n";
  ostr_tmp<<"StiExt Ach:"<<*(pb+3)<<"\n";
  ostr_tmp<<"StiExt GluCl:"<<*(pb+4)<<"\n";

  ostr=ostr_tmp.str();
*/
}

