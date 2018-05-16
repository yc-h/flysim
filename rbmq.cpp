
#include <iostream>
#include <stdio.h>
#include "NetSim.h"



#ifdef RABBITMQ
#include <amqp_tcp_socket.h>
#include <amqp.h>
#include <amqp_framing.h>


using namespace std;

void rbmq::MQini()
{
  MQport=5672;
  int MQstatus;

  MQconn = amqp_new_connection();

  MQsocket = amqp_tcp_socket_new(MQconn);
  if(!MQsocket) { cout<<"can't creating MQ socket\n"; }

  MQstatus = amqp_socket_open(MQsocket, MQhost, MQport);
  if(MQstatus)  { cout<<"can't opening MQ socket\n"; }

  amqp_login(MQconn, "/", 0, 131072, 0, AMQP_SASL_METHOD_PLAIN, "guest", "guest");

  amqp_channel_open(MQconn, 1);
  amqp_exchange_declare(MQconn, 1, amqp_cstring_bytes(EXCHANGE_NAME), amqp_cstring_bytes(exchangetype),
                        0, 0, amqp_empty_table);//equal to JAVA:channel.exchangeDeclare(EXCHANGE_NAME, "direct");

  props._flags = AMQP_BASIC_CONTENT_TYPE_FLAG | AMQP_BASIC_DELIVERY_MODE_FLAG;
  props.content_type = amqp_cstring_bytes("text/plain");
  props.delivery_mode = 2;
}

void rbmq::CMQini()
{
  CMQport=5672;
  int CMQstatus;
  amqp_bytes_t queuename;
  MQListop=false;

  CMQconn = amqp_new_connection();

  CMQsocket = amqp_tcp_socket_new(CMQconn);
  if(!CMQsocket) { cout<<"can't creating CMQ socket\n"; }

  CMQstatus = amqp_socket_open(CMQsocket, CMQhost, CMQport);
  if(CMQstatus)  { cout<<"can't opening CMQ socket\n"; }

  amqp_login(CMQconn, "/", 0, 131072, 0, AMQP_SASL_METHOD_PLAIN, "guest", "guest");

  amqp_channel_open(CMQconn, 1);

  amqp_queue_declare_ok_t *r = amqp_queue_declare(CMQconn, 1, amqp_empty_bytes, 0, 0, 0, 1, amqp_empty_table);
  queuename = amqp_bytes_malloc_dup(r->queue);
  amqp_exchange_declare(CMQconn, 1, amqp_cstring_bytes(EXCHANGE_NAME2), amqp_cstring_bytes(exchangetype),
                        0, 0, amqp_empty_table);//equal to JAVA:channel.exchangeDeclare(EXCHANGE_NAME, "direct");
  amqp_queue_bind(CMQconn, 1, queuename, amqp_cstring_bytes(EXCHANGE_NAME2), amqp_cstring_bytes(bindingkey),amqp_empty_table);
  amqp_basic_consume(CMQconn, 1, queuename, amqp_empty_bytes, 0, 1, 0, amqp_empty_table);
}

void rbmq::MQSend(string &s1)
{
  amqp_basic_publish(MQconn, 1, amqp_cstring_bytes(EXCHANGE_NAME), amqp_cstring_bytes(routingkey), 0, 0, &props, amqp_cstring_bytes(s1.c_str()) );
}

void rbmq::MQListen()
{
  while(!MQListop)
  {
    amqp_rpc_reply_t res;
    amqp_envelope_t envelope;
    amqp_maybe_release_buffers(CMQconn);
    res = amqp_consume_message(CMQconn, &envelope, NULL, 0);

    if (AMQP_RESPONSE_NORMAL != res.reply_type) { break;}

    mtx.lock();
    sprintf(chbuff,"%.*s",(int) envelope.message.body.len, (char *) envelope.message.body.bytes);
    mtx.unlock();

    amqp_destroy_envelope(&envelope);
  }

}


void rbmq::MQclose(amqp_connection_state_t &conn)
{
  amqp_channel_close(conn, 1, AMQP_REPLY_SUCCESS);
  amqp_connection_close(conn, AMQP_REPLY_SUCCESS);
  amqp_destroy_connection(conn);
}





#endif

