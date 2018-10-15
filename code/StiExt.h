
#include <string>
#include <random>
#include <thread>
#include "receptor.h"
#include "Info.h"
#include "RandGen.h"


using namespace std;
// global variable: random generator
extern random_device seed_gen;



class SinGen : public Information
{  // for extra stimulate with sine wave form

  public:
  SinGen();
  virtual ~SinGen(){};

  void SetPar(float,float,float); // set parameters
  void GetPar(float&,float&,float&); // get parameters

  float Irep(float);

  private:
  void TexInfo(string&);
  float amp;  // sine wave amplitude
  float freq; // sine wave frequency
  float phase;// sine wave phase

};

//----------------------------------------------------------------------------

class PulseGen : public Information
{  // for extra stimulate with Pulse

  public:
  PulseGen();
  virtual ~PulseGen(){};

  void SetPar(float,float,float,float); // set parameters
  void GetPar(float&,float&,float&,float&); // get parameters

  float Irep(float);

  private:
  void TexInfo(string&);
  float amp;  // Pulse amplitude
  float freq; // Pulse frequency
  float duty; // Pulse duty
  float phase;// Pulse phase
  float TBias;

};

//----------------------------------------------------------------------------
class StiExt
{
  public:
  StiExt();
  virtual ~StiExt(){};

  unsigned int rnd_seed;

  default_random_engine rdgen;

  normal_distribution<float> MemNoise;
  uniform_real_distribution<float> uni_dist;

  float p1;
  float p2;
  float p3;

  bool R1();
  bool R2();
  bool R3();

  SinGen RSin;  // sinusoidal stimulation
  PulseGen RPulse;  // pulse stimulation

  float ISti(float);
  float Igauss();
};

//----------------------------------------------------------------------------
class StiExt2
{
  public:
  StiExt2();
  virtual ~StiExt2(){delete pb;};

  unsigned int rnd_seed;

  default_random_engine rdgen;

  normal_distribution<float> MemNoise;
  uniform_real_distribution<float> uni_dist;

  float *pb;

  bool Rnd(float);

  SinGen RSin;  // sinusoidal stimulation
  PulseGen RPulse;  // pulse stimulation

  float ISti(float);
  float Igauss();
};

//----------------------------------------------------------------------------
template<class RepV, class RepV2>
class StiExt3
{
  public:
  StiExt3()
  {
    normal_distribution<float>::param_type par(MNoiseMean,MNoiseSTD);
    Inj.param(par);

    uniform_real_distribution<float>::param_type uni_par(0.0f,1.0f);
    uni_dist.param(uni_par);

    SetIgsEn();
  }


  ~StiExt3()
  {
    delete pb;
    delete ExtFst;
  };

  mt19937 *rdgen;
  normal_distribution<float> Inj;
  uniform_real_distribution<float> uni_dist;
  bool IgsEn;

  bool Rnd(float p)
  { return p>uni_dist(*rdgen)? true:false; }

  bool SetIgsEn(void)
  {
    if(Inj.mean()==0.0f && Inj.stddev()==0.0f)
      IgsEn=false;
    else
      IgsEn=true;

    return IgsEn;
  }


  float Igauss(void)
  { return Inj(*rdgen); }

  float *pb;
  RepV  *ExtFst; //AMPA,GABA,ACH,GCL
  RepV2 ExtSlo; //NMDA

};




struct BgSti // back ground stimulation
{
  float pb[5];  // probability
  uniform_real_distribution<float> uni_dist;
  normal_distribution<float> Inj;
};


template<class RepV, class RepV2>
struct BgSti1 // back ground stimulation
{
  RepV  SynFst1ES[4]; //AMPA,GABA,ACH,GCL
  RepV2 SynSlo1ES; //NMDA
  float pb[5];  // probability
  uniform_real_distribution<float> uni_dist;

  float mean,std;
};

template<class RepV, class RepV2>
struct BgSti2 : public BgSti1<RepV,RepV2>// back ground stimulation
{ normal_distribution<float> Inj; };

struct ziggurat // ziggurat method
{
  ziggurat()
  {
    kn=new uint32_t[128];
    fn=new float[128];
    wn=new float[128];
    r4_nor_setup(kn,fn,wn);
  };

  ~ziggurat()
  {
    // ziggurat gassuian distrubution
    delete kn;
    delete fn;
    delete wn;
  };

  uint32_t jsr,jcong;  //uni random seed and uni random
  uint32_t *kn;
  float *fn, *wn;


  float Igauss(float mean, float std)
  {
    jsr=cong_seeded(jcong);
    return mean + std*r4_nor(jsr, kn, fn, wn);
  };



  uint32_t cong_seeded( uint32_t &jcong )
  {
    uint32_t value;
    jcong = 69069 * ( jcong ) + 1234567;
    value = jcong;

    return value;
  };

  uint32_t shr3_seeded ( uint32_t &jsr )
  {
    uint32_t jsr_input;
    uint32_t value;

    jsr_input = jsr;
    jsr = ( jsr ^ ( jsr <<   13 ) );
    jsr = ( jsr ^ ( jsr >>   17 ) );
    jsr = ( jsr ^ ( jsr <<    5 ) );
    value = jsr_input + jsr;

    return value;
  };

  float r4_uni ( uint32_t &jsr )  //uniform distribution
  {
    uint32_t jsr_input;
    float value;

    jsr_input = jsr;
    jsr = ( jsr ^ ( jsr <<   13 ) );
    jsr = ( jsr ^ ( jsr >>   17 ) );
    jsr = ( jsr ^ ( jsr <<    5 ) );

    value = fmod(
    0.5 + (float)( jsr_input + jsr )/65536.0/65536.0, 1.0);

    return value;
  };


  float r4_nor ( uint32_t &jsr, uint32_t *kn, float *fn, float *wn) //normal distribution
  {
    int hz;
    uint32_t iz;
    const float r = 3.442620;
    float value;
    float x;
    float y;

    hz = ( int ) shr3_seeded ( jsr );
    iz = ( hz & 127 );

    if ( fabs ( hz ) < kn[iz] )
      value = ( float ) ( hz ) * wn[iz];
    else
    {
      for ( ; ; )
      {
        if ( iz == 0 )
        {
          for ( ; ; )
          {
            x = - 0.2904764 * log ( r4_uni ( jsr ) );
            y = - log ( r4_uni ( jsr ) );
            if ( x * x <= y + y )
              break;
          }

          if ( hz <= 0 )
            value = - r - x;
          else
            value = + r + x;

          break;
        }

        x = ( float ) ( hz ) * wn[iz];

        if ( fn[iz] + r4_uni ( jsr ) * ( fn[iz-1] - fn[iz] ) 
        < exp ( - 0.5 * x * x ) )
        {
          value = x;
          break;
        }

        hz = ( int ) shr3_seeded ( jsr );
        iz = ( hz & 127 );

        if ( fabs ( hz ) < kn[iz] )
        {
          value = ( float ) ( hz ) * wn[iz];
          break;
        }
      }
    }

    return value;
  };

  void r4_nor_setup ( uint32_t *kn, float *fn, float *wn) //setup normal distribution
  {
    double dn = 3.442619855899;
    int i;
    const double m1 = 2147483648.0;
    double q;
    double tn = 3.442619855899;
    const double vn = 9.91256303526217E-03;

    q = vn / exp ( - 0.5 * dn * dn );

    kn[0] = ( uint32_t ) ( ( dn / q ) * m1 );
    kn[1] = 0;

    wn[0] = ( float ) ( q / m1 );
    wn[127] = ( float ) ( dn / m1 );

    fn[0] = 1.0;
    fn[127] = ( float ) ( exp ( - 0.5 * dn * dn ) );

    for ( i = 126; 1 <= i; i-- )
    {
      dn = sqrt ( - 2.0 * log ( vn / dn + exp ( - 0.5 * dn * dn ) ) );
      kn[i+1] = ( uint32_t ) ( ( dn / tn ) * m1 );
      tn = dn;
      fn[i] = ( float ) ( exp ( - 0.5 * dn * dn ) );
      wn[i] = ( float ) ( dn / m1 );
    }

    return;
  };

};




