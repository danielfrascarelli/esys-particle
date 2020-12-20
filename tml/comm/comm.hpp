/////////////////////////////////////////////////////////////
//                                                         //
// Copyright (c) 2003-2017 by The University of Queensland //
// Centre for Geoscience Computing                         //
// http://earth.uq.edu.au/centre-geoscience-computing      //
//                                                         //
// Primary Business: Brisbane, Queensland, Australia       //
// Licensed under the Open Software License version 3.0    //
// http://www.apache.org/licenses/LICENSE-2.0              //
//                                                         //
/////////////////////////////////////////////////////////////


//--- point to point communication primitives for TML_Comm ---
#include "tml/message/packed_message.h"

//--- STL ---
#include <map>
using std::map;

/*!
  send a C-array of data with known dimensions
  
  \param data the data to be sent
  \param ndata the size of the array (nr. of elements)
  \param dest the rank of the destination process the data is sent to
  \param tag the message tag
  \warning no checks
*/
template <typename T> 
void TML_Comm::send_array(T* data,int ndata,int dest,int tag)
{
  MPI_Send(data,ndata,GetType(*data),dest,tag,m_comm);
}

/*!
  receive a C-array of data with known dimensions
  
  \param data the data to be received
  \param ndata the number of integers to be received
  \param source the rank of the destination process the data comes from
  \param tag the message tag
  \warning no checks
*/
template <typename T> 
void TML_Comm::receive_array(T* data,int ndata,int source,int tag)
{
  MPI_Recv(data,ndata,GetType(*data),source,tag,m_comm,&m_status);
}

/*!
  send and receive a C-array of data with known dimensions
  
  \param send_data the data to be sent
  \param send_count the number of integers to be sent
  \param recv_data the data to be received
  \param recv_count the number of integers to be received
  \param dest the rank of the destination process the data is sent to
  \param source the rank of the destination process the data comes from
  \param tag the message tag
  \warning no checks
*/
template <typename T,typename P> 
void TML_Comm::sendrecv_array(T* send_data,int send_count,P* recv_data,int recv_count,int dest,int source,int tag)
{
  MPI_Sendrecv(send_data,send_count,GetType(*send_data),dest,tag,recv_data,recv_count,GetType(*recv_data),source,tag,m_comm,&m_status);
}

/*!
  send single data
  
  \param data the data to be sent
  \param dest the rank of the destination process the data is sent to
  \param tag the message tag
*/
template <typename T> 
void TML_Comm::send(T data,int dest,int tag)
{
  MPI_Send(&data,1,GetType(data),dest,tag,m_comm);
}


/*!
  receive single data
  
  \param data the data to be received
  \param source the rank of the destination process the data comes from
  \param tag the message tag
  \warning no checks
*/
template <typename T> 
void TML_Comm::receive(T& data,int source ,int tag)
{
  MPI_Recv(&data,1,GetType(data),source,tag,m_comm,&m_status);
}

/*!
  send and receive single data
  
  \param send_data the data to be sent
  \param recv_data the data to be received
  \param dest the rank of the destination process the data is sent to
  \param source the rank of the destination process the data comes from
  \param tag the message tag
  \warning no checks
*/template <typename T,typename P> 
void TML_Comm::sendrecv(T send_data,P& recv_data,int dest,int source,int tag)
{
  MPI_Sendrecv(&send_data,1,GetType(send_data),dest,tag,&recv_data,1,GetType(recv_data),source,tag,m_comm,&m_status);
}

/*!
  Send an STL container or anything that has iterators, begin() and end(). Uses 2 send operations for size and data. 
  
  \param data the data to be sent
  \param dest the rank of the destination process the data is sent to
  \param tag the message tag
*/
template <typename T> 
void  TML_Comm::send_cont(const T& data, int dest, int tag)
{
  int data_size=data.size();

  // setup buffer
  typename T::value_type *buffer=new typename T::value_type[data_size];

  // put data into buffer
  int count=0;
  
  for(typename T::const_iterator iter=data.begin();
      iter!=data.end();
      iter++){
    void* buf=reinterpret_cast<void*>(&(buffer[count])); // get memory adress for buffer element
    new(buf)(typename T::value_type)(*iter); // initialize object at this adress
    // the placement new stuff (see Stroustrup, p.255ff) is necessary
    // because assignment buffer[count]=(*iter) doesn't work if the 
    // T::value_type contains constant member data. Initialization by 
    // calling the copy constructor directly doesn't work for builtin
    // types
    count++;
  }

  //send size
  send(data_size,dest,tag);

  //send data
  send_array(buffer,data_size,dest,tag+1024);

  // clean up
  delete [] buffer;
}

/*!
  Receive an STL container or anything that has iterators, begin() and end(). Single item
  insert (a.insert(p,t)) is used instead of range insert to be compatible with both 
  sequence and associative containers;
  
  \param data the data to be received
  \param source the rank of the destination process the data is coming from
  \param tag the message tag
*/
template <typename T> 
void TML_Comm::receive_cont(T& data,int source,int tag)
{
  int data_size;

  //get size
  receive(data_size,source,tag);
  // setup recv buffer
  typename T::value_type *buffer=new typename T::value_type[data_size];

  //get data
  receive_array(buffer,data_size,source,tag+1024);
  // insert into container
  for(int i=0;i<data_size;i++){
    data.insert(data.end(),buffer[i]);
  }
  delete [] buffer;
}


/*!
  send and receive an STL container or anything that has iterators, begin() and end(). Uses 2 MPI_Sendrecv operations for size and data
  
  \param send_data the data to be sent
  \param recv_data the data to be received
  \param dest the rank of the destination process the data is sent to
  \param source the rank of the destination process the data comes from
  \param tag the message tag
*/
template <typename T,typename P> 
void TML_Comm::sendrecv_cont(T send_data,P& recv_data,int dest,int source,int tag)
{
  int send_count=send_data.size();
  int recv_count;

  // setup send buffer
  typename T::value_type *send_buffer=new typename T::value_type[send_count];

  // put data into send buffer
  int count=0;
  for(typename T::const_iterator iter=send_data.begin();
      iter!=send_data.end();
      iter++){
    void* buf=reinterpret_cast<void*>(&(send_buffer[count])); // get memory adress for buffer element
    new(buf)(typename T::value_type)(*iter); // initialize object at this adress (see send_cont)
    count++;
  }

  //send/receive size
  sendrecv(send_count,recv_count,dest,source,tag);
  // check for recv from Null process -> set recv_count to 0 if so
  if(source==MPI_PROC_NULL){
    recv_count=0;
  }

  // setup recv buffer
  typename T::value_type *recv_buffer=new typename T::value_type[recv_count];

  //send/receive data
  sendrecv_array(send_buffer,send_count,recv_buffer,recv_count,dest,source,tag+1024);

  // insert into container
  for(int i=0;i<recv_count;i++){
    recv_data.insert(recv_data.end(),recv_buffer[i]);
  }

  delete [] send_buffer;
  delete [] recv_buffer;
}

/*!
  send and receive an STL container or anything that has iterators, begin() and end(), erase() and insert(). The input data is replaced (via erase/insert) by the output data. Uses separate send/receive buffers and MPI_Sendrecv instead of a single buffer and MPI_Sendrecv_replace to accomodate different size send and received messages. Uses 2 MPI_Sendrecv operations for size and data.
  
  \param data the data to be sent and received
  \param dest the rank of the destination process the data is sent to
  \param source the rank of the destination process the data comes from
  \param tag the message tag
*/
template <typename T> 
void TML_Comm::sendrecv_cont_replace(T& data,int dest,int source,int tag)
{
  int send_count=data.size();
  int recv_count;

  // setup send buffer
  typename T::value_type *send_buffer=new typename T::value_type[send_count];

  // put data into send buffer
  int count=0;
  for(typename T::const_iterator iter=data.begin();
      iter!=data.end();
      iter++){
    void* buf=reinterpret_cast<void*>(&(send_buffer[count])); // get memory adress for buffer element
    new(buf)(typename T::value_type)(*iter); // initialize object at this adress
    count++;
  }

  //send/receive size
  sendrecv(send_count,recv_count,dest,source,tag);
  // check for recv from Null process -> set recv_count to 0 if so
  if(source==MPI_PROC_NULL){
    recv_count=0;
  }

  // setup recv buffer
  typename T::value_type *recv_buffer=new typename T::value_type[recv_count];

  //send/receive data
  sendrecv_array(send_buffer,send_count,recv_buffer,recv_count,dest,source,tag+1024);

  // replace data
  data.erase(data.begin(),data.end());
    // insert into container
  for(int i=0;i<recv_count;i++){
    data.insert(data.end(),recv_buffer[i]);
  }


  delete [] send_buffer;
  delete [] recv_buffer;
}

/*!
  Send an STL container or anything that has iterators, begin() and end(). Uses 2 send operations for size and data. If the "checked" option is set, a checked message buffer is used to detect type mismatches between pack and unpack operations. The TML_Packed_Message has to grow dynamically. It is not possible to predetermine the size of the message because while T always has a size() function, it is not guaranteed that sizeof() or size() work for the valuetype of T.
  
  \param data the data to be sent
  \param dest the rank of the destination process the data is sent to
  \param checked use checked message buffer if true
  \param tag the message tag
  \warning checked stuff not yet implemented
*/
template <typename T> 
void TML_Comm::send_cont_packed(T data, int dest ,bool checked ,int tag)
{
  TML_Packed_Message* msg=new TML_Packed_Message(m_comm);
  int nb_data=data.size();
  
  msg->pack(nb_data); // pack number of items first
  // pack data
  for(typename T::const_iterator iter=data.begin();
      iter!=data.end();
      iter++){
    msg->pack(*iter);
  }

  //send size
  send(msg->size(),dest,tag);

  //send data
  send_array(msg->buffer(),msg->size(),dest,tag+1024);
  
  delete msg;
}

/*!
  Receive an STL container or anything that has iterators and insert(). Uses 2 receive operations for size and data. If the "checked" option is set, a checked message buffer is used to detect type mismatches between pack and unpack operations. 
  
  \param data the data to be received
  \param source the rank of the destination process the data is coming from
  \param checked use checked message buffer if true
  \param tag the message tag
  \warning checked stuff not yet implemented
*/
template <typename T> 
void TML_Comm::receive_cont_packed(T& data,int source,bool checked,int tag)
{
  int msg_size; // total size of the message
  int nb_data;  // number of data (i.e. the expected data.size())

  //get size
  receive(msg_size,source,tag);

  TML_Packed_Message* msg=new TML_Packed_Message(m_comm,msg_size);

  //get data
  receive_array(msg->buffer(),msg_size,source,tag+1024);

  // extract nuber of items
  nb_data=msg->pop_int();

  // unpack data
  for(int i=0;i<nb_data;i++){
    typename T::value_type tv;
    msg->unpack(tv);
    data.insert(data.end(),tv);
  }
  delete msg;
}

/*!
  send and receive an STL container or anything that has iterators, begin() and end() and insert(). Uses 2 MPI_Sendrecv operations for size and data. If the "checked" option is set, a checked message buffer is used to detect type mismatches between pack and unpack operations.
  
  \param send_data the data to be sent
  \param recv_data the data to be received
  \param dest the rank of the destination process the data is sent to
  \param source the rank of the destination process the data comes from
  \param tag the message tag
  \warning checked stuff not yet implemented
*/
template <typename T,typename P> 
void TML_Comm::sendrecv_cont_packed(T send_data, P& recv_data, int dest, int source, bool checked,int tag)
{
  TML_Packed_Message* send_msg=new TML_Packed_Message(m_comm);

  int send_nb_data=send_data.size();
  int send_msg_size;
  int recv_nb_data;
  int recv_msg_size; 

  send_msg->pack(send_nb_data); // pack number of items first
  // pack data
  for(typename T::const_iterator iter=send_data.begin();
      iter!=send_data.end();
      iter++){
    send_msg->pack(*iter);
  }
  send_msg_size=send_msg->size();
   
  //send/receive size
  sendrecv(send_msg_size,recv_msg_size,dest,source,tag);
  // check for recv from Null process -> set recv_count to 0 if so
  if(source==MPI_PROC_NULL){
    recv_msg_size=0;
  }

  // setup receive message 
  TML_Packed_Message* recv_msg=new TML_Packed_Message(m_comm,recv_msg_size);

  //send/receive data
  sendrecv_array(send_msg->buffer(),send_msg_size,recv_msg->buffer(),recv_msg_size,dest,source,tag+1024);

  if(source!=MPI_PROC_NULL){ // if the source exists
    // extract nuber of items
    recv_nb_data=recv_msg->pop_int();

    // unpack data
    typename T::value_type tv;
    for(int i=0;i<recv_nb_data;i++){
      recv_msg->unpack(tv);
      recv_data.insert(recv_data.end(),tv);
    }
  }
  delete send_msg;
  delete recv_msg;
}

/*!
  send and receive an STL container or anything that has iterators, begin() and end() and insert(). In input data is replaced by the received data. Uses 2 MPI_Sendrecv operations for size and data. If the "checked" option is set, a checked message buffer is used to detect type mismatches between pack and unpack operations.
  
  \param data the data to be sent and received
  \param dest the rank of the destination process the data is sent to
  \param source the rank of the destination process the data comes from
  \param tag the message tag
  \warning checked stuff not yet implemented
*/
template <typename T> 
void TML_Comm::sendrecv_cont_packed_replace(T& data, int dest, int source, bool checked,int tag)
{
  TML_Packed_Message* send_msg=new TML_Packed_Message(m_comm);

  int send_nb_data=data.size();
  int send_msg_size;
  int recv_nb_data;
  int recv_msg_size; 

  send_msg->pack(send_nb_data); // pack number of items first
  // pack data
  for(typename T::const_iterator iter=data.begin();
      iter!=data.end();
      iter++){
    send_msg->pack(*iter);
  }
  send_msg_size=send_msg->size();
   
  //send/receive size
  sendrecv(send_msg_size,recv_msg_size,dest,source,tag);
  // check for recv from Null process -> set recv_count to 0 if so
  if(source==MPI_PROC_NULL){
    recv_msg_size=0;
  }

  // setup receive message 
  TML_Packed_Message* recv_msg=new TML_Packed_Message(m_comm,recv_msg_size);

  //send/receive data
  sendrecv_array(send_msg->buffer(),send_msg_size,recv_msg->buffer(),recv_msg_size,dest,source,tag+1024);

  // extract nuber of items
  recv_nb_data=recv_msg->pop_int();

  // replace data
  data.erase(data.begin(),data.end());
  for(int i=0;i<recv_nb_data;i++){
    typename T::value_type tv;
    recv_msg->unpack(tv);
    data.insert(data.end(),tv);
  }

  delete send_msg;
  delete recv_msg;
}

