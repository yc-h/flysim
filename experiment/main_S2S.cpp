#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <thread> // thread
#include <mutex>  // mutex
#include <atomic>
#include <chrono>


// segment to segment distance calculation

#define SMALL_NUM   1e-8f // anything that avoids division overflow
#define dot(u,v)   ((u).x * (v).x + (u).y * (v).y + (u).z * (v).z)// dot product (3D) which allows vector operations in arguments
#define norm(v)    sqrt(dot(v,v))  // norm = length of  vector
#define abs(x)     ((x) >= 0.0f ? (x) : -(x))   //  absolute value

#define MAX_DISTANCE 1.0f

#define IDP s[i].idpar
#define IDPj s[j].idpar
#define PreR allswc[NID].s[i].R
#define PreRidp allswc[NID].s[allswc[NID].IDP].R
#define PoR allswc[n].s[j].R
#define PoRidp allswc[n].s[allswc[n].IDPj].R

#define UNITX 1.0f
#define UNITY 1.0f
#define UNITZ 1.0f

using namespace std;

mutex mtx;           // mutex for critical section
unsigned int ThdNum;
unsigned int op;
string FileConWgtHead;
vector<string> fstr;
atomic<unsigned int> LoadCnt;


struct RecPar
{
  float Beg;
  unsigned int Count;
  float Inc;
};

struct swcDef
{
  unsigned int R;
  float x,y,z;
  float D;
  unsigned int idpar;
};

struct swc
{
  swcDef *s;
  unsigned int swcSize;
};

swc *allswc;

struct point
{
  float x;
  float y;
  float z;
};

struct pointApr
{
  float x;
  float y;
  float z;
  unsigned int apr;
};

struct seg
{
  point P0;
  point P1;
};

void dpReplace(float f,string &s)
{
  size_t found;
  stringstream out;
  out<<f;

  s = out.str();

  found = s.find(".");
  if(found != string::npos)
  {  s.replace(found,1,"_");}
}

float S2S( point &uPri, point &vPri, point &wPri, float &sc_ret, float &tc_ret)
{
  point u;
  point v;
  point w;

  u.x=uPri.x*UNITX;
  u.y=uPri.y*UNITY;
  u.z=uPri.z*UNITZ;

  v.x=vPri.x*UNITX;
  v.y=vPri.y*UNITY;
  v.z=vPri.z*UNITZ;

  w.x=wPri.x*UNITX;
  w.y=wPri.y*UNITY;
  w.z=wPri.z*UNITZ;

  float  a = dot(u,u);         // always >= 0
  float  b = dot(u,v);
  float  c = dot(v,v);         // always >= 0
  float  d = dot(u,w);
  float  e = dot(v,w);
  float  D = a*c - b*b;        // always >= 0
  float  sc, sN, sD = D;     // sc = sN / sD, default sD = D >= 0
  float  tc, tN, tD = D;       // tc = tN / tD, default tD = D >= 0

    // compute the line parameters of the two closest points
    if (D < SMALL_NUM) { // the lines are almost parallel
        sN = 0.0f;         // force using point P0 on segment S1
        sD = 1.0f;         // to prevent possible division by 0.0f later
        tN = e;
        tD = c;
    }
    else {                 // get the closest points on the infinite lines
        sN = (b*e - c*d);
        tN = (a*e - b*d);
        if (sN < 0.0f) {        // sc < 0 => the s=0 edge is visible
            sN = 0.0f;
            tN = e;
            tD = c;
        }
        else if (sN > sD) {  // sc > 1  => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

    if (tN < 0.0f) {            // tc < 0 => the t=0 edge is visible
        tN = 0.0f;
        // recompute sc for this edge
        if (-d < 0.0f)
            sN = 0.0f;
        else if (-d > a)
            sN = sD;
        else {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD) {      // tc > 1  => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if ((-d + b) < 0.0f)
            sN = 0.0f;
        else if ((-d + b) > a)
            sN = sD;
        else {
            sN = (-d +  b);
            sD = a;
        }
    }
    // finally do the division to get sc and tc
    sc = (abs(sN) < SMALL_NUM ? 0.0f : sN / sD);
    tc = (abs(tN) < SMALL_NUM ? 0.0f : tN / tD);

    // get the difference of the two closest points
    point   dP;
    dP.x= w.x + (sc * u.x) - (tc * v.x);  // =  S1(sc) - S2(tc)
    dP.y= w.y + (sc * u.y) - (tc * v.y);  // =  S1(sc) - S2(tc)
    dP.z= w.z + (sc * u.z) - (tc * v.z);  // =  S1(sc) - S2(tc)
    sc_ret=sc;
    tc_ret=tc;
    return norm(dP);   // return the closest distance
};


//-------------------------------------------------------------------------1
void runPolarityConnection(unsigned int NID, ofstream *fo1, RecPar RecP)
{
  float d,sc_ret=0.0f,tc_ret=0.0f;
  unsigned int *Dist;
  point u;
  point v;
  point w;

  Dist = new unsigned int[RecP.Count];

  for(unsigned int n=0;n<fstr.size();n++)
  {
    if(n!=NID)
    {
      for(unsigned int p=0;p<RecP.Count;p++)
      { Dist[p]=0; }

      for(unsigned int i=0;i<allswc[NID].swcSize;i++)
      {
        if(PreR==2 && PreRidp==2)  // 2=Axon
        {
          u.x= allswc[NID].s[i].x - allswc[NID].s[allswc[NID].IDP].x;
          u.y= allswc[NID].s[i].y - allswc[NID].s[allswc[NID].IDP].y;
          u.z= allswc[NID].s[i].z - allswc[NID].s[allswc[NID].IDP].z;

          for(unsigned int j=0;j<allswc[n].swcSize;j++)
          {
            if(PoR==3 && PoRidp==3) // 3=dendrite
            {
              v.x = allswc[n].s[j].x - allswc[n].s[allswc[n].IDPj].x;
              v.y = allswc[n].s[j].y - allswc[n].s[allswc[n].IDPj].y;
              v.z = allswc[n].s[j].z - allswc[n].s[allswc[n].IDPj].z;

              w.x = allswc[NID].s[allswc[NID].IDP].x - allswc[n].s[allswc[n].IDPj].x;
              w.y = allswc[NID].s[allswc[NID].IDP].y - allswc[n].s[allswc[n].IDPj].y;
              w.z = allswc[NID].s[allswc[NID].IDP].z - allswc[n].s[allswc[n].IDPj].z;

              d=S2S(u,v,w,sc_ret,tc_ret);

              for(unsigned int p=0;p<RecP.Count;p++)
              { if(d<(p*RecP.Inc+RecP.Beg)) Dist[p]++; }

            }
          }
        }
      }

      for(unsigned int p=0;p<RecP.Count;p++)
      { if(Dist[p]!=0) fo1[p]<<NID<<" "<<n<<" "<<Dist[p]<<"\n"; }
    }
  }

}


//-------------------------------------------------------------------------2
void runNoPolarityConnection(unsigned int NID, ofstream *fo1, RecPar RecP)
{
  float d,sc_ret=0.0f,tc_ret=0.0f;
  unsigned int *Dist;
  point u;
  point v;
  point w;

  Dist = new unsigned int[RecP.Count];

  for(unsigned int n=0;n<fstr.size();n++)
  {
    if(n!=NID)
    {
      for(unsigned int p=0;p<RecP.Count;p++)
      { Dist[p]=0; }

      for(unsigned int i=0;i<allswc[NID].swcSize;i++)
      {
          u.x= allswc[NID].s[i].x - allswc[NID].s[allswc[NID].IDP].x;
          u.y= allswc[NID].s[i].y - allswc[NID].s[allswc[NID].IDP].y;
          u.z= allswc[NID].s[i].z - allswc[NID].s[allswc[NID].IDP].z;

          for(unsigned int j=0;j<allswc[n].swcSize;j++)
          {
              v.x = allswc[n].s[j].x - allswc[n].s[allswc[n].IDPj].x;
              v.y = allswc[n].s[j].y - allswc[n].s[allswc[n].IDPj].y;
              v.z = allswc[n].s[j].z - allswc[n].s[allswc[n].IDPj].z;

              w.x = allswc[NID].s[allswc[NID].IDP].x - allswc[n].s[allswc[n].IDPj].x;
              w.y = allswc[NID].s[allswc[NID].IDP].y - allswc[n].s[allswc[n].IDPj].y;
              w.z = allswc[NID].s[allswc[NID].IDP].z - allswc[n].s[allswc[n].IDPj].z;

              d=S2S(u,v,w,sc_ret,tc_ret);

              for(unsigned int p=0;p<RecP.Count;p++)
              { if(d<(p*RecP.Inc+RecP.Beg)) Dist[p]++; }

          }
      }

      for(unsigned int p=0;p<RecP.Count;p++)
      { if(Dist[p]!=0) fo1[p]<<NID<<" "<<n<<" "<<Dist[p]<<"\n"; }
    }
  }

}

//-------------------------------------------------------------------------3
void runPolarityPositionAxon(unsigned int NID, ofstream *fo1, ofstream *fo2, RecPar RecP)
{
  float d,sc_ret=0.0f,tc_ret=0.0f;
  unsigned int *Dist;
  point u;
  point v;
  point w;

  point NewPVec;
  vector<vector<point> > PVec;

  vector<pointApr> PVecRut;
  pointApr NewPVecRut;

  Dist = new unsigned int[RecP.Count];
  PVec.insert(PVec.begin(),RecP.Count,vector<point> ());

  for(unsigned int n=0;n<fstr.size();n++)
  {
    if(n!=NID)
    {
      for(unsigned int p=0;p<RecP.Count;p++)
      { Dist[p]=0; }

      for(unsigned int i=0;i<allswc[NID].swcSize;i++)
      {
        if(PreR==2 && PreRidp==2)  // 2=Axon
        {
          u.x= allswc[NID].s[i].x - allswc[NID].s[allswc[NID].IDP].x;
          u.y= allswc[NID].s[i].y - allswc[NID].s[allswc[NID].IDP].y;
          u.z= allswc[NID].s[i].z - allswc[NID].s[allswc[NID].IDP].z;

          for(unsigned int j=0;j<allswc[n].swcSize;j++)
          {
            if(PoR==3 && PoRidp==3) // 3=dendrite
            {
              v.x = allswc[n].s[j].x - allswc[n].s[allswc[n].IDPj].x;
              v.y = allswc[n].s[j].y - allswc[n].s[allswc[n].IDPj].y;
              v.z = allswc[n].s[j].z - allswc[n].s[allswc[n].IDPj].z;

              w.x = allswc[NID].s[allswc[NID].IDP].x - allswc[n].s[allswc[n].IDPj].x;
              w.y = allswc[NID].s[allswc[NID].IDP].y - allswc[n].s[allswc[n].IDPj].y;
              w.z = allswc[NID].s[allswc[NID].IDP].z - allswc[n].s[allswc[n].IDPj].z;

              d=S2S(u,v,w,sc_ret,tc_ret);

              for(unsigned int p=0;p<RecP.Count;p++)
              {
                if(d<(p*RecP.Inc+RecP.Beg))
                {
                  Dist[p]++;
                  NewPVec.x=allswc[NID].s[allswc[NID].IDP].x+sc_ret*u.x;
                  NewPVec.y=allswc[NID].s[allswc[NID].IDP].y+sc_ret*u.y;
                  NewPVec.z=allswc[NID].s[allswc[NID].IDP].z+sc_ret*u.z;
                  PVec[p].push_back(NewPVec);
                }
              }

            }
          }
        }
      }

      for(unsigned int p=0;p<RecP.Count;p++)
      {
        bool found;
        for(unsigned int ii=0;ii<PVec[p].size();ii++)
        {
          found=false;
          for(unsigned int jj=0;jj<PVecRut.size();jj++)
          {
            if((PVec[p][ii].x==PVecRut[jj].x)&&(PVec[p][ii].y==PVecRut[jj].y)&&(PVec[p][ii].z==PVecRut[jj].z))
            {
              found=true;
              PVecRut[jj].apr++;
            }
            
          }
          if(!found) // first appear in PVec
          {
            NewPVecRut.x=PVec[p][ii].x; 
            NewPVecRut.y=PVec[p][ii].y;
            NewPVecRut.z=PVec[p][ii].z;
            NewPVecRut.apr=1;

            PVecRut.push_back(NewPVecRut);
          }
        }

        for(unsigned int jj=0;jj<PVecRut.size();jj++)
        { fo2[p]<<NID<<" "<<n<<" "<<PVecRut[jj].x<<" "<<PVecRut[jj].y<<" "<<PVecRut[jj].z<<" "<<PVecRut[jj].apr<<"\n";  }

        if(Dist[p]!=0) fo1[p]<<NID<<" "<<n<<" "<<Dist[p]<<"\n";

        PVecRut.clear();
        PVec[p].clear();
      }

    }
  }

}

//-----------------------------------------------------------------------4
void runNoPolarityPositionAxon(unsigned int NID, ofstream *fo1, ofstream *fo2, RecPar RecP)
{
  float d,sc_ret=0.0f,tc_ret=0.0f;
  unsigned int *Dist;
  point u;
  point v;
  point w;

  point NewPVec;
  vector<vector<point> > PVec;

  vector<pointApr> PVecRut;
  pointApr NewPVecRut;

  Dist = new unsigned int[RecP.Count];
  PVec.insert(PVec.begin(),RecP.Count,vector<point> ());

  for(unsigned int n=0;n<fstr.size();n++)
  {
    if(n!=NID)
    {
      for(unsigned int p=0;p<RecP.Count;p++)
      { Dist[p]=0; }

      for(unsigned int i=0;i<allswc[NID].swcSize;i++)
      {
          u.x= allswc[NID].s[i].x - allswc[NID].s[allswc[NID].IDP].x;
          u.y= allswc[NID].s[i].y - allswc[NID].s[allswc[NID].IDP].y;
          u.z= allswc[NID].s[i].z - allswc[NID].s[allswc[NID].IDP].z;

          for(unsigned int j=0;j<allswc[n].swcSize;j++)
          {
              v.x = allswc[n].s[j].x - allswc[n].s[allswc[n].IDPj].x;
              v.y = allswc[n].s[j].y - allswc[n].s[allswc[n].IDPj].y;
              v.z = allswc[n].s[j].z - allswc[n].s[allswc[n].IDPj].z;

              w.x = allswc[NID].s[allswc[NID].IDP].x - allswc[n].s[allswc[n].IDPj].x;
              w.y = allswc[NID].s[allswc[NID].IDP].y - allswc[n].s[allswc[n].IDPj].y;
              w.z = allswc[NID].s[allswc[NID].IDP].z - allswc[n].s[allswc[n].IDPj].z;

              d=S2S(u,v,w,sc_ret,tc_ret);

              for(unsigned int p=0;p<RecP.Count;p++)
              {
                if(d<(p*RecP.Inc+RecP.Beg))
                {
                  Dist[p]++;
                  NewPVec.x=allswc[NID].s[allswc[NID].IDP].x+sc_ret*u.x;
                  NewPVec.y=allswc[NID].s[allswc[NID].IDP].y+sc_ret*u.y;
                  NewPVec.z=allswc[NID].s[allswc[NID].IDP].z+sc_ret*u.z;
                  PVec[p].push_back(NewPVec);
                }
              }

          }
      }

      for(unsigned int p=0;p<RecP.Count;p++)
      {
        bool found;
        for(unsigned int ii=0;ii<PVec[p].size();ii++)
        {
          found=false;
          for(unsigned int jj=0;jj<PVecRut.size();jj++)
          {
            if((PVec[p][ii].x==PVecRut[jj].x)&&(PVec[p][ii].y==PVecRut[jj].y)&&(PVec[p][ii].z==PVecRut[jj].z))
            {
              found=true;
              PVecRut[jj].apr++;
            }
            
          }
          if(!found) // first appear in PVec
          {
            NewPVecRut.x=PVec[p][ii].x; 
            NewPVecRut.y=PVec[p][ii].y;
            NewPVecRut.z=PVec[p][ii].z;
            NewPVecRut.apr=1;

            PVecRut.push_back(NewPVecRut);
          }
        }

        for(unsigned int jj=0;jj<PVecRut.size();jj++)
        { fo2[p]<<NID<<" "<<n<<" "<<PVecRut[jj].x<<" "<<PVecRut[jj].y<<" "<<PVecRut[jj].z<<" "<<PVecRut[jj].apr<<"\n";  }

        if(Dist[p]!=0) fo1[p]<<NID<<" "<<n<<" "<<Dist[p]<<"\n";

        PVecRut.clear();
        PVec[p].clear();
      }

    }
  }

}

//-------------------------------------------------------------------------5
void runPolarityPositionDendrite(unsigned int NID, ofstream *fo1, ofstream *fo2, RecPar RecP)
{
  float d,sc_ret=0.0f,tc_ret=0.0f;
  unsigned int *Dist;
  point u;
  point v;
  point w;

  point NewPVec;
  vector<vector<point> > PVec;

  point NewQVec;
  vector<vector<point> > QVec;

  Dist = new unsigned int[RecP.Count];
  PVec.insert(PVec.begin(),RecP.Count,vector<point> ());
  QVec.insert(QVec.begin(),RecP.Count,vector<point> ());

  for(unsigned int n=0;n<fstr.size();n++)
  {
    if(n!=NID)
    {
      for(unsigned int p=0;p<RecP.Count;p++)
      { Dist[p]=0; }

      for(unsigned int i=0;i<allswc[NID].swcSize;i++)
      {
        if(PreR==2 && PreRidp==2)  // 2=Axon
        {
          u.x= allswc[NID].s[i].x - allswc[NID].s[allswc[NID].IDP].x;
          u.y= allswc[NID].s[i].y - allswc[NID].s[allswc[NID].IDP].y;
          u.z= allswc[NID].s[i].z - allswc[NID].s[allswc[NID].IDP].z;

          for(unsigned int j=0;j<allswc[n].swcSize;j++)
          {
            if(PoR==3 && PoRidp==3) // 3=dendrite
            {
              v.x = allswc[n].s[j].x - allswc[n].s[allswc[n].IDPj].x;
              v.y = allswc[n].s[j].y - allswc[n].s[allswc[n].IDPj].y;
              v.z = allswc[n].s[j].z - allswc[n].s[allswc[n].IDPj].z;

              w.x = allswc[NID].s[allswc[NID].IDP].x - allswc[n].s[allswc[n].IDPj].x;
              w.y = allswc[NID].s[allswc[NID].IDP].y - allswc[n].s[allswc[n].IDPj].y;
              w.z = allswc[NID].s[allswc[NID].IDP].z - allswc[n].s[allswc[n].IDPj].z;

              d=S2S(u,v,w,sc_ret,tc_ret);

              for(unsigned int p=0;p<RecP.Count;p++)
              {
                if(d<(p*RecP.Inc+RecP.Beg))
                {
                  Dist[p]++;
                  NewPVec.x=allswc[NID].s[allswc[NID].IDP].x+sc_ret*u.x;
                  NewPVec.y=allswc[NID].s[allswc[NID].IDP].y+sc_ret*u.y;
                  NewPVec.z=allswc[NID].s[allswc[NID].IDP].z+sc_ret*u.z;
                  PVec[p].push_back(NewPVec);

                  NewQVec.x=allswc[n].s[allswc[n].IDPj].x+tc_ret*v.x;
                  NewQVec.y=allswc[n].s[allswc[n].IDPj].y+tc_ret*v.y;
                  NewQVec.z=allswc[n].s[allswc[n].IDPj].z+tc_ret*v.z;
                  QVec[p].push_back(NewQVec);

                }
              }

            }
          }
        }
      }

      for(unsigned int p=0;p<RecP.Count;p++)
      {
        for(unsigned int ii=0;ii<PVec[p].size();ii++)
        { fo2[p]<<NID<<" "<<n<<" "<<PVec[p][ii].x<<" "<<PVec[p][ii].y<<" "<<PVec[p][ii].z<<" "<<QVec[p][ii].x<<" "<<QVec[p][ii].y<<" "<<QVec[p][ii].z<<"\n"; }

        if(Dist[p]!=0) fo1[p]<<NID<<" "<<n<<" "<<Dist[p]<<"\n";

        PVec[p].clear();
        QVec[p].clear();
      }

    }
  }

}

//-----------------------------------------------------------------------6
void runNoPolarityPositionDendrite(unsigned int NID, ofstream *fo1, ofstream *fo2, RecPar RecP)
{
  float d,sc_ret=0.0f,tc_ret=0.0f;
  unsigned int *Dist;
  point u;
  point v;
  point w;

  point NewPVec;
  vector<vector<point> > PVec;

  point NewQVec;
  vector<vector<point> > QVec;

  Dist = new unsigned int[RecP.Count];
  PVec.insert(PVec.begin(),RecP.Count,vector<point> ());
  QVec.insert(QVec.begin(),RecP.Count,vector<point> ());

  for(unsigned int n=0;n<fstr.size();n++)
  {
    if(n!=NID)
    {
      for(unsigned int p=0;p<RecP.Count;p++)
      { Dist[p]=0; }

      for(unsigned int i=0;i<allswc[NID].swcSize;i++)
      {
          u.x= allswc[NID].s[i].x - allswc[NID].s[allswc[NID].IDP].x;
          u.y= allswc[NID].s[i].y - allswc[NID].s[allswc[NID].IDP].y;
          u.z= allswc[NID].s[i].z - allswc[NID].s[allswc[NID].IDP].z;

          for(unsigned int j=0;j<allswc[n].swcSize;j++)
          {
              v.x = allswc[n].s[j].x - allswc[n].s[allswc[n].IDPj].x;
              v.y = allswc[n].s[j].y - allswc[n].s[allswc[n].IDPj].y;
              v.z = allswc[n].s[j].z - allswc[n].s[allswc[n].IDPj].z;

              w.x = allswc[NID].s[allswc[NID].IDP].x - allswc[n].s[allswc[n].IDPj].x;
              w.y = allswc[NID].s[allswc[NID].IDP].y - allswc[n].s[allswc[n].IDPj].y;
              w.z = allswc[NID].s[allswc[NID].IDP].z - allswc[n].s[allswc[n].IDPj].z;

              d=S2S(u,v,w,sc_ret,tc_ret);

              for(unsigned int p=0;p<RecP.Count;p++)
              {
                if(d<(p*RecP.Inc+RecP.Beg))
                {
                  Dist[p]++;
                  NewPVec.x=allswc[NID].s[allswc[NID].IDP].x+sc_ret*u.x;
                  NewPVec.y=allswc[NID].s[allswc[NID].IDP].y+sc_ret*u.y;
                  NewPVec.z=allswc[NID].s[allswc[NID].IDP].z+sc_ret*u.z;
                  PVec[p].push_back(NewPVec);

                  NewQVec.x=allswc[n].s[allswc[n].IDPj].x+tc_ret*v.x;
                  NewQVec.y=allswc[n].s[allswc[n].IDPj].y+tc_ret*v.y;
                  NewQVec.z=allswc[n].s[allswc[n].IDPj].z+tc_ret*v.z;
                  QVec[p].push_back(NewQVec);
                }
              }

          }
      }

      for(unsigned int p=0;p<RecP.Count;p++)
      {
        for(unsigned int ii=0;ii<PVec[p].size();ii++)
        { fo2[p]<<NID<<" "<<n<<" "<<PVec[p][ii].x<<" "<<PVec[p][ii].y<<" "<<PVec[p][ii].z<<" "<<QVec[p][ii].x<<" "<<QVec[p][ii].y<<" "<<QVec[p][ii].z<<"\n"; }

        if(Dist[p]!=0) fo1[p]<<NID<<" "<<n<<" "<<Dist[p]<<"\n";

        PVec[p].clear();
        QVec[p].clear();
      }

    }
  }

}

//-----------------------------------------------------------------------------


void SWCTol(string &FileSWC)
{
  string str="";
  string s1="";
  ifstream fin;
  istringstream istream;

  fin.open(FileSWC.c_str());
  while(!fin.eof())
  {
    string::size_type startpos=0;

    getline(fin,s1);
    s1.find_first_of(" :;,'=\t\r\n",startpos);
    if(startpos!=string::npos && s1.size()!=0)
    { fstr.push_back(s1); }
  }
  fin.close();


  allswc = new swc[fstr.size()];

  for(unsigned int i=0,PreIdx=0;i<fstr.size();i++)
  {
    cout<<"read file:"<<100.0f*(float)(i+1)/fstr.size()<<"%        \r"<<flush;

    PreIdx=0;
    s1="";
    str="";

    fin.open(fstr[i].c_str());
    while(!fin.eof())
    {
      getline(fin,s1);
      if(s1.find_first_of("#", 0)!=0)
      {
        str = str + s1 +"\n";
        PreIdx++;
      }
    }
    PreIdx--;
    fin.close();

    allswc[i].swcSize=PreIdx;
    allswc[i].s = new swcDef[PreIdx];

    istream.str(str);

    int PreIdp=0;
    for(unsigned int m=0,SID=0,R=0;m<PreIdx;m++)
    {
      for(unsigned int n=0;n<7;n++)
      {
        switch(n)
        {
          case 0:
          istream>>SID;
          break;

          case 1:
          istream>>R;
          allswc[i].s[m].R=abs(R)/10;
          break;

          case 2:
          istream>>allswc[i].s[m].x;
          break;

          case 3:
          istream>>allswc[i].s[m].y;
          break;

          case 4:
          istream>>allswc[i].s[m].z;
          break;

          case 5:
          istream>>allswc[i].s[m].D;
          break;

          case 6:
          istream>>PreIdp;
          allswc[i].s[m].idpar=abs(PreIdp)-1;
          break;

        }
      }
    }
  }

}


void RunThd(unsigned int thid, RecPar r)
{
  typedef std::chrono::high_resolution_clock Clock;
  auto t1=Clock::now();
  auto t2=Clock::now();
  auto t3=Clock::now();
  auto t4=Clock::now();

  ofstream *fo1,*fo2;
  unsigned int NID;
  float f;
  string s;




  mtx.lock();
  fo1 = new ofstream[r.Count];

  if((op==3) || (op==4) || (op==5)|| (op==6))
  { fo2 = new ofstream[r.Count]; }

  for(unsigned int i=0;i<r.Count;i++)
  {
    f=r.Beg+r.Inc*i;
    dpReplace(f,s);
    string FileStr(FileConWgtHead+".d"+s+".th"+to_string(thid));
    fo1[i].open(FileStr.c_str());
    cout<<"thread:"<<thid<<",open file "<<i<<":"<<FileStr<<"\n";

  if((op==3) || (op==4) || (op==5)|| (op==6))
    {
        FileStr=FileConWgtHead+"_p.d"+s+".th"+to_string(thid);
        fo2[i].open(FileStr.c_str());
        cout<<"thread:"<<thid<<",open file "<<i<<":"<<FileStr<<"\n";
    }
  }
  mtx.unlock();

  t1=Clock::now();
  while((NID=LoadCnt++)<fstr.size())
  {
      t3=Clock::now();
      switch(op)
      {
        case 2:
        runNoPolarityConnection(NID,fo1,r);
        break;

        case 3:
        runPolarityPositionAxon(NID,fo1,fo2,r);
        break;

        case 4:
        runNoPolarityPositionAxon(NID,fo1,fo2,r);
        break;

        case 5:
        runPolarityPositionDendrite(NID,fo1,fo2,r); //just fix runNoPolarityPositionDendrite to runPolarityPositionDendrite here from main31.cpp
        break;

        case 6:
        runNoPolarityPositionDendrite(NID,fo1,fo2,r);
        break;

        case 1:
        default:
        runPolarityConnection(NID,fo1,r);
        break;
      }

      t4=Clock::now();
      cout<<"Neuron ID:"<<NID<<", t="<<(chrono::duration_cast<chrono::nanoseconds>(t4-t3).count())*1e-9<<", thread:"<<thid<<" t="<<(chrono::duration_cast<chrono::nanoseconds>(t4-t1).count())*1e-9<<endl;
  }

  t2=Clock::now();
  cout<<"In thread:"<<thid<<", total time="<<(chrono::duration_cast<chrono::nanoseconds>(t2-t1).count())*1e-9<<endl;

  for(unsigned int i=0;i<r.Count;i++)
  { fo1[i].close(); }

  if((op==3) || (op==4) || (op==5)|| (op==6))
  {
    for(unsigned int i=0;i<r.Count;i++)
    { fo2[i].close(); }
  }

}






int main (int argc, char * const argv[])
{
  string FileSWC;
  thread *th;
  LoadCnt=0;
  op=(unsigned int)atoi(argv[1]);

  if(argc==1)
  {
    ThdNum=1;
    FileSWC="AllSWC.txt";
    FileConWgtHead="ConWgt.txt";
  }
  else
  {
    ThdNum=(unsigned int)atoi(argv[2]);
    FileSWC=argv[3];
    FileConWgtHead=argv[4];
  }

  cout<<"example: s2s_v32.out option threads InputIndex.txt OutputHeader begin count inc\n"
      <<"  option:\n"
      <<"  1=Polarity + Connection Number only\n"
      <<"  2=No Polarity + Connection Number only\n"
      <<"  3=Polarity + Connection Number Position + AxonPosition\n"
      <<"  4=No Polarity + Connection + Position + AxonPosition\n"
      <<"  5=Polarity + Connection + Position + Axon Dendrite Pair\n"
      <<"  6=No Polarity + Connection + Position + Axon Dendrite Pair\n"
      <<"\n"
      <<" threads=1\n"
      <<" InputIndex=AllSWC.txt\n"
      <<" OutputHeader:ConWgt.txt\n"
      <<"\n"
      <<" distance criterion:\n"
      <<" begin:1.0\n"
      <<" count:3\n"
      <<" inc:0.32\n";

  SWCTol(FileSWC);
  cout<<"\noption="<<op<<",threads="<<ThdNum<<", Inputindex.txt:"<<FileSWC<<", OutputHead:"<<FileConWgtHead<<", Neurons="<<fstr.size()<<endl;


  RecPar r;
  r.Beg = (float)atof(argv[5]);
  r.Count = (unsigned int)atoi(argv[6]);
  r.Inc = (float)atof(argv[7]);
  cout<<"begin:"<<r.Beg<<", Count:"<<r.Count<<", Inc:"<<r.Inc<<", End:"<<r.Beg+r.Inc*(r.Count-1);
  cout<<endl;

  if(ThdNum>fstr.size())
  { ThdNum=fstr.size(); }

  th=new thread[ThdNum];

  for(unsigned int i=0;i<ThdNum;i++)
  { *(th+i)=thread(RunThd,i,r); }

  for(unsigned int i=0;i<ThdNum;i++)
  {(th+i)->join();}

  return 0;
}


