
#include <string>
#include "Inc.h"


using namespace std;

class Information
{
  public:
  Information(){};
  virtual ~Information(){};

  void Info();
  void Info(string&);

  virtual void TexInfo(string&)=0;
};


