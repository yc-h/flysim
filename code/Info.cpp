#include <iostream>
#include <string>
#include <omp.h>
#include "Info.h"

using namespace std;

void Information::Info()
{
  string s1;
  TexInfo(s1);
  cout<<s1;
}

void Information::Info(string& s1)
  { TexInfo(s1); }




