#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <omp.h>
#include <random>
#include <chrono>
#include <thread>
#include <mutex>
#include <string.h>
#include <sys/stat.h>
#include "NetSim.h"

//windows socket
#ifdef _WINSOCKET
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <winsock2.h>
#include <ws2tcpip.h>
#pragma comment (lib, "Ws2_32.lib")
#pragma comment (lib, "Mswsock.lib")
#pragma comment (lib, "AdvApi32.lib")
WSADATA wsaData;  //winsocket
#endif


//linux socket
#ifdef _LINSOCKET
#include <unistd.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/tcp.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
void closesocket(int socket) {close(socket);}
#endif

using namespace std;
//-----------------------------------------------------------------------------
//------------------NetSimLEGWorm

NetSimLEGWorm::NetSimLEGWorm()
{
  TSteps=0;
  CID=0;
}

NetSimLEGWorm::~NetSimLEGWorm()
{
  for(unsigned int i=0;i<Worms.size();i++)
  { delete Worms[i]; }

  {vector<NetSimGnl2 *> ().swap(Worms);}
}

void NetSimLEGWorm::LdEst(int TNum)  //load estimate
{}


void NetSimLEGWorm::TexInfo(string& ostr)
{
  ostringstream ostr_tmp;
  string str;

  ostr=ostr_tmp.str();
}


unsigned int NetSimLEGWorm::AddWorm(string &conf,int thsize)
{
  string s1;
  NetSimGnl2 *NewWorm=new NetSimGnl2;
  NewWorm->ActCtr=AMPA_ACT | GABA_ACT | NMDA_ACT | ACH_ACT | GCL_ACT | GAP_ACT | LIF_ACT;
  NewWorm->dtime=0.1f;
  NewWorm->ActCtr &= ~UDFRDSEED_ACT;

  NewWorm->CSta(conf,NewWorm->dtime);  // statistic of netwrok configuration file
  NewWorm->parse(conf,NewWorm->dtime);  // read netwrok configuration file
  NewWorm->SetPp(thsize,ACCURATE);  // setup population

  s1="WormID"+to_string(Worms.size())+".log";
  NewWorm->logDataWrite(s1.c_str());
  NewWorm->VarInit();

  Worms.push_back(NewWorm);

  return (Worms.size()-1);
}



void NetSimLEGWorm::DelWorm(unsigned int WID)
{
  //only release memory but not delete address
  delete Worms[WID];
  Worms[WID]=nullptr;
}

void NetSimLEGWorm::DoEvt(string &pro,unsigned int WID)
{
  {vector<OutCtr> ().swap(Worms[WID]->OFCtr);}
  {vector<TmEvt> ().swap(Worms[WID]->Evts);}

  Worms[WID]->DefMacro(pro);
  Worms[WID]->parse(pro,Worms[WID]->Evts,Worms[WID]->OFCtr);  // read time events, and output file events
  for(unsigned int i=0;i<Worms[WID]->Evts.size();i++)
  {
    if(WEvts.Type==EVTEXTFREQ)
    {
      for(unsigned int m=Worms[WID]->NIdx[WEvts.Pp];m<(Worms[WID]->NIdx[WEvts.Pp]+Worms[WID]->ConfSta[WEvts.Pp].NeuNum);m++)
      {
        switch(WEvts.SubTyp)
        {
          case KW_AMPA:
          if(WEvts.Var1==0.0f)
            WNets.ActCtr &= ~AMPA_ACT;
          else
          {
            WNets.ActCtr |= AMPA_ACT;
            WBGSTI.pb[0]=(float)(WEvts.Var1*Worms[WID]->dt*1e-3f);
          }
          break;

          case KW_GABA:
          if(WEvts.Var1==0.0f)
            WNets.ActCtr &= ~GABA_ACT;
          else
          {
            WNets.ActCtr |= GABA_ACT;
            WBGSTI.pb[1]=(float)(WEvts.Var1*Worms[WID]->dt*1e-3f);
          }
          break;

          case KW_ACH:
          if(WEvts.Var1==0.0f)
            WNets.ActCtr &= ~ACH_ACT;
          else
          {
            WNets.ActCtr |= ACH_ACT;
            WBGSTI.pb[2]=(float)(WEvts.Var1*Worms[WID]->dt*1e-3f);
          }
          break;

          case KW_GCL:
          if(WEvts.Var1==0.0f)
            WNets.ActCtr &= ~GCL_ACT;
          else
          {
            WNets.ActCtr |= GCL_ACT;
            WBGSTI.pb[3]=(float)(WEvts.Var1*Worms[WID]->dt*1e-3f);
          }
          break;

          case KW_NMDA:
          if(WEvts.Var1==0.0f)
            WNets.ActCtr &= ~NMDA_ACT;
          else
          {
            WNets.ActCtr |= NMDA_ACT;
            WBGSTI.pb[4]=(float)(WEvts.Var1*Worms[WID]->dt*1e-3f);
          }
          break;

          default:;
        }
      }
    }
    else if(WEvts.Type==EVTCUNTINJ)
    {
      for(unsigned int m=Worms[WID]->NIdx[WEvts.Pp];m<(Worms[WID]->NIdx[WEvts.Pp]+Worms[WID]->ConfSta[WEvts.Pp].NeuNum);m++)
     {
        WBGSTI.mean=WEvts.Var1;
        WBGSTI.std=WEvts.Var2;
        if(WEvts.Var1==0.0f && WEvts.Var2==0.0f)
          WNets.ActCtr &= ~INJ_ACT;
        else
          WNets.ActCtr |= INJ_ACT;
     }
    }
  }

}


void NetSimLEGWorm::RunWorm(int setps)
{
  bool RunFlag=false;
  for(unsigned int i=0;i<Worms.size();i++)
  {
    if(Worms[i]!=nullptr)
    {

      for(int i=0;i<setps;i++)
      {  Worms[i]->EqsOne(Worms[i]->ActCtr);} // only 1 thread
      RunFlag=true;
    }
  }

  if(RunFlag) TSteps++;
}

/*
void NetSimLEGWorm::Ack(int sockfd)
{
  char AckInfo[] = "Flysim acknowledge!";
  int clientfd;
  int rval=0;

  //struct sockaddr_in client_addr;
  //socklen_t addrlen = sizeof(client_addr);
  //clientfd = accept(sockfd, (struct sockaddr *) &client_addr, &addrlen);// Wait and Accept connection
  clientfd = accept(sockfd, NULL,NULL);// Wait and Accept connection

  int flag = 1; 
  setsockopt(clientfd, IPPROTO_TCP, TCP_NODELAY, (char *) &flag, sizeof(int));
  // Send acknowledge message to client
  rval=send(clientfd, AckInfo, sizeof(AckInfo),0);

  while(rval >=0 )
  { rval=CltReq(clientfd);}

  closesocket(clientfd);
}
*/


int NetSimLEGWorm::CltReq(int clientfd)
{
  char buf32[32];
  int rev;
  int ret=0;

  // receive command
  memset(buf32,0,32);
  rev=recv(clientfd,buf32,sizeof(buf32),0);
  cout<<"command "<<CID++<<":"<<buf32<<", rev="<<rev<<endl;

  if(strcmp(buf32,"FILE_WRITE")==0)
  {
    rev=recv(clientfd,buf32,sizeof(buf32),0);
    cout<<"write file name:"<<buf32<<", rev="<<rev<<endl;

    /* receive file*/
    char ch[16];
    ofstream fo;
    fo.open(buf32);

    memset(ch,0,16);
    while( (rev=recv(clientfd,ch,sizeof(ch),0)) > 0  )
    {
      if(strcmp(ch,"END_FILE_WRITE")==0) { break; }
      fo<<ch;
      memset(ch,0,16);
    }

    fo.close();
    cout<<"file write\n";
  }
  if(strcmp(buf32,"FILE_SEND")==0)
  {
    rev=recv(clientfd,buf32,sizeof(buf32),0);
    cout<<"write file name:"<<buf32<<", rev="<<rev<<endl;

    /* receive file*/
    ofstream fo;
    fo.open(buf32);

    memset(buf32,0,32);
    while( (rev=recv(clientfd,buf32,sizeof(buf32),0)) > 0  )
    {
      if(strcmp(buf32,"END_FILE_SEND")==0)
      { break; }

      fo<<buf32;
      memset(buf32,0,32);
    }

    fo.close();
    cout<<"file write\n";
  }
  else if(strcmp(buf32,"FILE_READ")==0)
  {
    rev=recv(clientfd,buf32,sizeof(buf32),0);
    cout<<"read file name:"<<buf32<<", rev="<<rev<<endl;

    /* Send file*/
    char ch[16];
    ifstream fin;
    fin.open(buf32);

    if(fin)
    {
      while(!fin.eof())
      {
        fin.get(ch,15,EOF);
        rev=send(clientfd, ch, sizeof(ch),0);
      }
      rev=send(clientfd, "END_FILE_READ", sizeof("END_FILE_READ"),0);
    }
    else
    { cout<<".conf open file error!"<<endl; }

    fin.close();
    cout<<"file read\n";
  }
  else if(strcmp(buf32,"ADD_WORM")==0)
  {
    rev=recv(clientfd,buf32,sizeof(buf32),0);
    string conf=buf32;

    AddWorm(conf,1);
    cout<<"conf file name:"<<conf<<", rev="<<rev<<", worms:"<<Worms.size()<<endl;

    memset(buf32,0,32);
    sprintf(buf32,"WID=%d END_ADD_WORM",(int)(Worms.size()-1));
    rev=send(clientfd, buf32, sizeof(buf32),0);
  }
  else if(strcmp(buf32,"DEL_WORM")==0)
  {
    rev=recv(clientfd,buf32,sizeof(buf32),0);
    unsigned int WID=(unsigned int)atoi(buf32);
    cout<<"DelWorm WID:"<<WID<<", rev="<<rev<<endl;

    DelWorm(WID);

    memset(buf32,0,32);
    sprintf(buf32,"WID=%d END_DEL_WORM",WID);
    rev=send(clientfd, buf32, sizeof(buf32),0);
  }
  else if(strcmp(buf32,"DO_EVENTS")==0)
  {
    rev=recv(clientfd,buf32,sizeof(buf32),0);
    string pro=buf32;
    cout<<"pro file name:"<<pro<<", rev="<<rev<<endl;

    memset(buf32,0,32);
    rev=recv(clientfd,buf32,sizeof(buf32),0);
    unsigned int WID=(unsigned int)atoi(buf32); // receive WID
    cout<<"Do events WID:"<<WID<<", rev="<<rev<<endl;

    DoEvt(pro,WID);
    rev=send(clientfd, "END_DO_EVENTS", sizeof("END_DO_EVENTS"),0);
  }
  else if(strcmp(buf32,"DATA_READ")==0)
  {
    rev=recv(clientfd,buf32,sizeof(buf32),0);
    cout<<"Data type:"<<buf32<<", rev="<<rev<<endl;

    if(strcmp(buf32,"Membrane_Potential_Header")==0)
    {
      int jj=0;
      memset(buf32,0,32);
      sprintf(buf32,"TotalSteps:%d ",jj++);
      rev=send(clientfd, buf32, sizeof(buf32),0);

      for(unsigned int i=0;i<Worms.size();i++) // Send nets data
      {
        if(Worms[i]!=nullptr)
        {
          for(unsigned int j=0;j<Worms[i]->NetSize;j++)
          {
            memset(buf32,0,32);
            sprintf(buf32,"WormID %d:%d ",i,jj++);
            rev=send(clientfd, buf32, sizeof(buf32),0);

            memset(buf32,0,32);
            sprintf(buf32,"NeuronID %d:%d ",j,jj++);
            rev=send(clientfd, buf32, sizeof(buf32),0);

            memset(buf32,0,32);
            sprintf(buf32,"Vm:%d ",jj++);
            rev=send(clientfd, buf32, sizeof(buf32),0);
          }
        }
      }
      rev=send(clientfd, "END_DATA_HEADER", sizeof("END_DATA_HEADER"),0);
    }
    else if(strcmp(buf32,"Membrane_Potential")==0)
    {
      memset(buf32,0,32);
      sprintf(buf32,"%d ",TSteps);
      rev=send(clientfd, buf32, sizeof(buf32),0);

      for(unsigned int i=0;i<Worms.size();i++) // Send nets data
      {
        if(Worms[i]!=nullptr)
        {
          for(unsigned int j=0;j<Worms[i]->NetSize;j++)
          {
            memset(buf32,0,32);
            sprintf(buf32,"%d %d %6f ",i,j,Worms[i]->NetVBse[j].V*1e-3f);
            rev=send(clientfd, buf32, sizeof(buf32),0);
          }
        }
      }
      rev=send(clientfd, "END_DATA_READ", sizeof("END_DATA_READ"),0);
    }
    else if(strcmp(buf32,"Membrane_Potential_Cont")==0)
    {
      memset(buf32,0,32);
      rev=recv(clientfd,buf32,sizeof(buf32),0);
      cout<<"repeat times:"<<buf32<<", rev="<<rev<<endl;
      int End=atoi(buf32);

      int Start=TSteps;
      for(int ii=Start;ii<Start+End;ii++)
      {
        memset(buf32,0,32);
        sprintf(buf32,"%d ",TSteps);
        rev=send(clientfd, buf32, sizeof(buf32),0);

        for(unsigned int i=0;i<Worms.size();i++) // Send nets data
        {
          if(Worms[i]!=nullptr)
          {
            for(unsigned int j=0;j<Worms[i]->NetSize;j++)
            {
              memset(buf32,0,32);
              sprintf(buf32,"%d %d %f ",i,j,Worms[i]->NetVBse[j].V*1e-3f);
              rev=send(clientfd, buf32, sizeof(buf32),0);
            }
          }
        }

        RunWorm(1);
      }
      rev=send(clientfd, "END_DATA_READ", sizeof("END_DATA_READ"),0);
    }
    else if(strcmp(buf32,"Spikes_Header")==0)
    {
      int jj=0;
      memset(buf32,0,32);
      sprintf(buf32,"TotalSteps:%d ",jj++);
      rev=send(clientfd, buf32, sizeof(buf32),0);

      memset(buf32,0,32);
      sprintf(buf32,"WormID:%d ",jj++);
      rev=send(clientfd, buf32, sizeof(buf32),0);

      memset(buf32,0,32);
      sprintf(buf32,"NeuronID:%d ",jj++);
      rev=send(clientfd, buf32, sizeof(buf32),0);

      rev=send(clientfd, "END_DATA_HEADER", sizeof("END_DATA_HEADER"),0);
    }
    else if(strcmp(buf32,"Spikes")==0)
    {
      for(unsigned int i=0;i<Worms.size();i++) // Send nets data
      {
        if(Worms[i]!=nullptr)
        {
          for(unsigned int j=0;j<Worms[i]->NetSize;j++)
          {
            //if(Worms[i]->Nets[j].SpiOut)
            if(Worms[i]->NetVLIF[j].NState & SPIKE_F)
            {
              memset(buf32,0,32);
              sprintf(buf32,"%d %d %d ",TSteps,i,j);
              rev=send(clientfd, buf32, sizeof(buf32),0);
            }
          }
        }
      }
      rev=send(clientfd, "END_DATA_READ", sizeof("END_DATA_READ"),0);
    }
    else if(strcmp(buf32,"Spikes_Cont")==0)
    {
      memset(buf32,0,32);
      rev=recv(clientfd,buf32,sizeof(buf32),0);
      cout<<"repeat times:"<<buf32<<", rev="<<rev<<endl;
      int End=atoi(buf32);

      int Start=TSteps;
      for(int ii=Start;ii<Start+End;ii++)
      {
        for(unsigned int i=0;i<Worms.size();i++) // Send nets data
        {
          if(Worms[i]!=nullptr)
          {
            for(unsigned int j=0;j<Worms[i]->NetSize;j++)
            {
              //if(Worms[i]->Nets[j].SpiOut)
              if(Worms[i]->NetVLIF[j].NState & SPIKE_F)
              {
                memset(buf32,0,32);
                sprintf(buf32,"%d %d %d ",TSteps,i,j);
                rev=send(clientfd, buf32, sizeof(buf32),0);
              }
            }
          }
        }

        RunWorm(1);
      }
      rev=send(clientfd, "END_DATA_READ", sizeof("END_DATA_READ"),0);
    }
  }
  else if(strcmp(buf32,"END_CLIENT")==0)
  { ret=-1;}

cout<<"\n";

  if(rev<0)
  {
    cout<<"client has some errors!"<<endl;
    ret=-1;
  }

  return ret;
}

void NetSimLEGWorm::DaemonMode(int port)
{
  struct addrinfo *servinfo=NULL;
  struct addrinfo serv;
  int status;
  int sockfd;

  char AckInfo[] = "Flysim acknowledge!";
  int clientfd;
  int rval=0;

#ifdef _WINSOCKET
  status = WSAStartup(0x0202, &wsaData);//winsocket
  if (status != 0)
  { printf("WSAStartup failed with error: %d\n", status); }
#endif
  
  // initialize value in dest
  memset(&serv,0,sizeof(serv));
  serv.ai_family = AF_INET;
  serv.ai_socktype = SOCK_STREAM;
  serv.ai_protocol = IPPROTO_TCP;
  serv.ai_flags = AI_PASSIVE;
  char str[8];
  sprintf(str,"%d",port);

  if ((status = getaddrinfo(NULL, str, &serv, &servinfo)) != 0)
  {
    fprintf(stderr, "getaddrinfo error: %s\n", gai_strerror(status));
    exit(1);
    }

  sockfd = socket(servinfo->ai_family, servinfo->ai_socktype, servinfo->ai_protocol);
  bind(sockfd, servinfo->ai_addr, (int)servinfo->ai_addrlen);
  freeaddrinfo(servinfo);
  listen(sockfd, 20);



  while(1)
  {
    clientfd = accept(sockfd, NULL,NULL);// Wait and Accept connection
    if(clientfd>=0)
    {
    int flag = 1; 
    setsockopt(clientfd, IPPROTO_TCP, TCP_NODELAY, (char *) &flag, sizeof(int));
    // Send acknowledge message to client
    rval=send(clientfd, AckInfo, sizeof(AckInfo),0);

    while(rval >=0 )
    { rval=CltReq(clientfd);}

    closesocket(clientfd);
    }
  }

  closesocket(sockfd);
#ifdef _WINSOCKET
  WSACleanup(); //winsocket
#endif
}






