#include <string>
using namespace std;

class RandGen
{
  public:
  RandGen();
  virtual ~RandGen(){};

  void init();
  void Info();
  void Info(string&);

  void SetPar(float,float,int);
  void SetPar(float,float);
  void SetPar(float);
  void SetPar(int);
  float RVal(); // random value generating

  float Constant(float); //for mean
  float Uniform(float); //for mean
  float Exponential(float); //for mean
  int   Poisson(float); //for mean
  bool  PoiST(float,float); //for mean
  float Pareto(float); //for alpha
  float Normal(float,float); //for mean, standard devation
  int   Geometric(float);  //for p
  float Weibull(float,float); //for scale, for shape
  float Erlang(int,float); //for scale, for shape

  float Constant();
  float Uniform();
  float Exponential();
  int   Poisson();
  bool  PoiST();
  float Pareto();
  float Normal();
  int   Geometric();
  float Weibull();
  float Erlang();

  private:
  void TexInfo(string&);
  float a;  //for mean
  float b;  // for number of event space, standard devation
  int  c;  // for scale

};






