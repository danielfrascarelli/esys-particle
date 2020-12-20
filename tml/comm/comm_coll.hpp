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


// --- collective communication primitives for TML_Comm ---
#include "tml/message/packed_multi_message.h"

/*!
  broadcast single data to all other nodes

  \param data the data to be broadcast
*/
template <typename T> 
void TML_Comm::broadcast(T data)
{
  MPI_Bcast(&data,1,GetType(data),rank(),m_comm);
}

/*!
  broadcast an array of known size

  \param data the array to be broadcast
  \param ndata the size of the array (nr. of elements)
*/
template <typename T> 
void TML_Comm::broadcast_array(T* data,int ndata)
{
  MPI_Bcast(data,ndata,GetType(*data),rank(),m_comm);
}
/*!
  broadcast the content of a STL container of simple types. Uses
  2 MPI broadcasts for size and data.

  \param data the data to be broadcast
*/
template <typename T> 
void TML_Comm::broadcast_cont(const T& data)
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

  // broadcast size 
  broadcast(data_size);

  // broadcast buffer
  broadcast_array(buffer,data_size);

  // clean up
  delete [] buffer;
}

/*!
  broadcast the content of a STL container of packable objects. Uses
  2 MPI broadcasts for size and data.

  \param data the data to be broadcast
*/ 
template <typename T> 
void TML_Comm::broadcast_cont_packed(const T &data)
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
  // broadcast size
  broadcast(msg->size());

  // broadcast data
  broadcast_array(msg->buffer(),msg->size());

  delete msg;
}

/*!
  Receive broadcast of single data from a given node
  
  \param data the data to be broadcast
  \param root the node which sent the broadcast
*/
template <typename T> 
void TML_Comm::recv_broadcast(T& data,int root)
{
  MPI_Bcast(&data,1,GetType(data),root,m_comm);
}

/*!
  Receive broadcast of an array of known size

  \param data the array to be broadcast
  \param ndata the size of the array (nr. of elements)
  \param root the node which sent the broadcast
*/
template <typename T> 
void TML_Comm::recv_broadcast_array(T* data,int ndata,int root)
{
  MPI_Bcast(data,ndata,GetType(*data),root,m_comm);
}

/*!
  Receive broadcast of the content of a STL container of simple types
  from a given node. Uses 2 MPI broadcasts for size and data.
  
  \param data the data to be broadcast
  \param root the node which sent the broadcast
*/
template <typename T> 
void TML_Comm::recv_broadcast_cont(T& data,int root)
{
  int data_size;

  //get size
  recv_broadcast(data_size,root);
  // setup recv buffer
  typename T::value_type *buffer=new typename T::value_type[data_size];

  //get data
  recv_broadcast_array(buffer,data_size,root);
  // insert into container
  for(int i=0;i<data_size;i++){
    data.insert(data.end(),buffer[i]);
  }

  //clean up
  delete [] buffer;
}
/*!
  Receive broadcast of the content of a STL container of packable objects
  from a given node. Uses 2 MPI broadcasts for size and data.
  
  \param data the data to be broadcast
  \param root the node which sent the broadcast
*/
template <typename T> 
void TML_Comm::recv_broadcast_cont_packed(T& data,int root)
{
  int msg_size; // total size of the message
  int nb_data;  // number of data (i.e. the expected data.size())

  //get size
  recv_broadcast(msg_size,root);

  //setup message
  TML_Packed_Message* msg=new TML_Packed_Message(m_comm,msg_size);
  
  //get data
  recv_broadcast_array(msg->buffer(),msg_size,root);
   // extract nuber of items
  nb_data=msg->pop_int();

  // unpack data
  for(int i=0;i<nb_data;i++){
    typename T::value_type tv;
    msg->unpack(tv);
    data.insert(data.end(),tv);
  }

  // clean up
  delete msg;
}

/*!
  scatter the content of a multimap. The key of the multimap is used to decide where to
  send the data. Uses one MPI_Scatter (sizes) and one MPI_Scatterv (data) call.

  \param data the multimap containing the data to be scattered
  \todo checks, handling of out of range indices
*/
template <typename T> 
void TML_Comm::scatter(const multimap<int,T> data)
{
  // put data into buffer
  int total_size=data.size();
  T* buffer=new T[total_size];
  int comm_size=size();
  int *size_buffer=new int[comm_size];
  int *offs_buffer=new int[comm_size];
  int count=0;
  int count_p=0;
  for(int i=0;i<comm_size;i++){
    typename multimap<int,T>::const_iterator begin=data.find(i);
    typename multimap<int,T>::const_iterator end=data.upper_bound(i);
    if(begin!=data.end()){
      for(typename multimap<int,T>::const_iterator iter=begin;
	  iter!=end;
	  iter++){
	buffer[count]=iter->second;
	count++;
      }
    }
    size_buffer[i]=count-count_p;
    count_p=count;
  }
  // send size info
  int dummy;
  MPI_Scatter(size_buffer,1,MPI_INT,&dummy,1,MPI_INT,rank(),m_comm);
  // construct offsets from sizes
  offs_buffer[0]=0;
  for(int i=1;i<comm_size;i++){
    offs_buffer[i]=offs_buffer[i-1]+size_buffer[i-1];
  }
  // send data
  T dummy2;
  MPI_Scatterv(buffer,size_buffer,offs_buffer,GetType(*buffer),&dummy2,0,GetType(*buffer),rank(),m_comm);
  // clean up
  delete [] size_buffer;
  delete [] offs_buffer;
  delete [] buffer;
}

/*!
  receive scattered data

  \param data the received data
  \param root the process which scattered the data
*/
template <typename T> 
void TML_Comm::recv_scatter(T& data,int root)
{
  // receive size
  int size;
  MPI_Scatter(NULL,0,MPI_INT,&size,1,MPI_INT,root,m_comm);
  // allocate buffer buffer
  typename T::value_type *buffer=new typename T::value_type[size];
  // receive data
  MPI_Scatterv(NULL,NULL,NULL,MPI_INT,buffer,size,GetType(*buffer),root,m_comm);
  // put data in container
  for(int i=0;i<size;i++){
    data.insert(data.end(),buffer[i]);
  }
  // clean up
  delete [] buffer;
}

/*!
  Gather data from all other processes in the communicator into a multimap. The multimap-key 
  will be set according to the rank of the process where the data came from

  \param data the multimap
*/
template <typename T> 
void TML_Comm::gather(multimap<int,T>& data)
{
  int dummy=0;
  int comm_size=size();
  int *size_buffer=new int[comm_size];
  int *offs_buffer=new int[comm_size];
  // receive sizes
  MPI_Gather(&dummy,1,MPI_INT,size_buffer,1,MPI_INT,rank(),m_comm);
  int totalsize=0;
  for(int i=0;i<comm_size;i++){
    totalsize+=size_buffer[i];
  }
  // setup receive buffer
  T *buffer=new T[totalsize];
  // construct offsets from sizes
  offs_buffer[0]=0;
  for(int i=1;i<comm_size;i++){
    offs_buffer[i]=offs_buffer[i-1]+size_buffer[i-1];
  }
  // receive data
  T dummy2;
  MPI_Gatherv(&dummy2,0,GetType(dummy),buffer,size_buffer,offs_buffer,GetType(*buffer),rank(),m_comm);
  // put data into multimap
  for(int i=0;i<comm_size;i++){
    for(int j=offs_buffer[i];j<offs_buffer[i]+size_buffer[i];j++){
      data.insert(pair<int,T>(i,buffer[j]));
    }
  }
  // clean up
  delete [] size_buffer;	
  delete [] offs_buffer;
  delete [] buffer;
}

/*!
  Gather data from all other processes in the communicator into a multimap. The multimap-key 
  will be set according to the rank of the process where the data came from. Debug version 
  (output data)

  \param data the multimap
*/
template <typename T> 
void TML_Comm::gather_debug(multimap<int,T>& data)
{
  int dummy=0;
  int comm_size=size();
  int *size_buffer=new int[comm_size];
  int *offs_buffer=new int[comm_size];
  // receive sizes
  MPI_Gather(&dummy,1,MPI_INT,size_buffer,1,MPI_INT,rank(),m_comm);
  int totalsize=0;
  for(int i=0;i<comm_size;i++){
    console.Debug() << "buffer size " << i << " - " << size_buffer[i] << "\n";
    totalsize+=size_buffer[i];
  }
  // setup receive buffer
  T *buffer=new T[totalsize];
  // construct offsets from sizes
  offs_buffer[0]=0;
  for(int i=1;i<comm_size;i++){
    offs_buffer[i]=offs_buffer[i-1]+size_buffer[i-1];
  }
  // receive data
  T dummy2;
  MPI_Gatherv(&dummy2,0,GetType(dummy),buffer,size_buffer,offs_buffer,GetType(*buffer),rank(),m_comm);
  // put data into multimap
  for(int i=0;i<comm_size;i++){
    for(int j=offs_buffer[i];j<offs_buffer[i]+size_buffer[i];j++){
      data.insert(pair<int,T>(i,buffer[j]));
    }
  }
  // clean up
  delete [] size_buffer;	
  delete [] offs_buffer;
  delete [] buffer;
}

/*!
  Send data to be gathered by a root process.

  \param data the data (a container)
  \param root the root process
*/
template <typename T> 
void TML_Comm::send_gather(T& data,int root)
{
  // send size
  int size=data.size();
  MPI_Gather(&size,1,MPI_INT,NULL,0,MPI_INT,root,m_comm);
  // setup send buffer
  typename T::value_type *buffer=new typename T::value_type[size];
  // put data into send buffer
  int count=0;
  for(typename T::const_iterator iter=data.begin();
      iter!=data.end();
      iter++){
    buffer[count]=*iter;
    count++;
  }
  // send data
  MPI_Gatherv(buffer,size,GetType(*buffer),NULL,NULL,NULL,MPI_INT,root,m_comm);
  // clean up
  delete [] buffer;
}

 
/*!
  Send data to be gathered by a root process. Debug version.

  \param data the data (a container)
  \param root the root process
*/
template <typename T> 
void TML_Comm::send_gather_debug(T& data,int root)
{
  // send size
  int size=data.size();
  console.Debug() << "send size :" << size << "\n";
  MPI_Gather(&size,1,MPI_INT,NULL,0,MPI_INT,root,m_comm);
  // setup send buffer
  typename T::value_type *buffer=new typename T::value_type[size];
  // put data into send buffer
  int count=0;
  for(typename T::const_iterator iter=data.begin();
      iter!=data.end();
      iter++){
    buffer[count]=*iter;
    count++;
  }
  console.Debug() << "send count :" << count << "\n";
  // send data
  MPI_Gatherv(buffer,size,GetType(*buffer),NULL,NULL,NULL,MPI_INT,root,m_comm);
  // clean up
  delete [] buffer;
}

 
/*!
  Scatter the content of a multimap of packable objects. The key of the multimap is used to decide where to
  send the data. Uses one MPI_Scatter (sizes) and one MPI_Scatterv (data) call.

  \param data the multimap containing the data to be scattered
  \todo checks, handling of out of range indices
*/
template <typename T> 
void TML_Comm::scatter_packed(const multimap<int,T> data)
{
  int comm_size=size();
  // construct new packed multimessage
  TML_PackedMultiMessage *msg=new TML_PackedMultiMessage(m_comm);
  // pack data
  for(int i=1;i<comm_size;i++){
    typename multimap<int,T>::const_iterator begin=data.find(i);
    typename multimap<int,T>::const_iterator end=data.upper_bound(i);
    (*msg)[i].pack(int(data.count(i)));
    if(begin!=data.end()){
      for(typename multimap<int,T>::const_iterator iter=begin;
	  iter!=end;
	  iter++){
	(*msg)[i].pack(iter->second);
      }
    }
  }
  // send size
  int dummy;
  MPI_Scatter(msg->sizes(),1,MPI_INT,&dummy,1,MPI_INT,rank(),m_comm);  
  // send data
  T dummy2;
  MPI_Datatype data_type=GetType(*(msg->buffer()));
  MPI_Scatterv(msg->buffer(),msg->sizes(),msg->offsets(),data_type,
  	       &dummy2,0,data_type,rank(),m_comm);
  // clean up
  delete msg;
}

/*!
  receive scattered packed data

  \param data the received data
  \param root the process which scattered the data
*/
template <typename T> 
void TML_Comm::recv_scatter_packed(T& data,int root)
{
  // receive size
  int msg_size;
  MPI_Scatter(NULL,0,MPI_INT,&msg_size,1,MPI_INT,root,m_comm);
  // construct packed message of sufficient size
  TML_Packed_Message* msg=new TML_Packed_Message(m_comm,msg_size);

  // recceive data
  MPI_Datatype data_type=GetType(*(msg->buffer()));
  MPI_Scatterv(NULL,NULL,NULL,MPI_INT,msg->buffer(),msg_size,data_type,root,m_comm);
  // extract number of items
  msg->begin_unpack();
  int nitems=msg->pop_int();
  // unpack data
  typename T::value_type item;
  for(int i=0;i<nitems;i++){
    msg->unpack(item);
    data.insert(data.end(),item);
  }
  // clean up
  delete msg;
}

template <typename T> 
void TML_Comm::gather_packed(multimap<int,T> &data)
{
  console.Debug() << "TML_Comm::gather_packed: enter\n";
  int dummy=0;
  int comm_size=size();
  int *size_buffer=new int[comm_size];
  int *offs_buffer=new int[comm_size];
  // receive sizes
  console.Debug()
     << "TML_Comm::gather_packed: gathering sizes\n";
  MPI_Gather(&dummy,1,MPI_INT,size_buffer,1,MPI_INT,rank(),m_comm);
  int totalsize=0;
  for(int i=0;i<comm_size;i++){
    //console.Debug()
    //  << "TML_Comm::gather_packed:"
    //  << " size from rank " << i << " = " << size_buffer[i] << "\n";
    totalsize+=size_buffer[i];
  }
  console.Debug()
     << "TML_Comm::gather_packed: total msg size = " << totalsize << "\n";
  // setup receive buffer
    //setup message
  TML_Packed_Message* msg=new TML_Packed_Message(m_comm,totalsize);
  
  // construct offsets from sizes
  offs_buffer[0]=0;
  for(int i=1;i<comm_size;i++){
    offs_buffer[i]=offs_buffer[i-1]+size_buffer[i-1];
  }
  // receive data
  T dummy2;
  console.Debug()
     << "TML_Comm::gather_packed: gathering data\n";
  MPI_Gatherv(
    &dummy2,0, GetType(dummy), msg->buffer(), size_buffer, offs_buffer,
    GetType(*(msg->buffer())),rank(),m_comm
  );
  // put data into multimap
  console.Debug()
     << "TML_Comm::gather_packed: unpacking into multi-map\n";
  for(int i=0;i<comm_size;i++){
    if (size_buffer[i] > 0)
    {
      const int numElems = msg->pop_int();
      for(int j=0; j < numElems; j++){
        //console.Debug()
        //  << "TML_Comm::gather_packed:"
        //  << " unpacking object (" << i << "," << j << ")\n";
        T t;
        msg->unpack(t);
        data.insert(pair<int,T>(i, t));
      }
    }
  }
  console.Debug() << "TML_Comm::gather_packed: done unpacking into multi-map\n";
  // clean up
  delete [] size_buffer;
  delete [] offs_buffer;
  delete msg;
  console.Debug() << "TML_Comm::gather_packed: exit\n";
}

template <typename T> 
void TML_Comm::send_gather_packed(const T &data, int root)
{
  console.Debug() << "TML_Comm::send_gather_packed: enter\n";
  // setup send buffer
  TML_Packed_Message* msg =
    new TML_Packed_Message(
      m_comm,
      std::max(static_cast<size_t>(64), data.size()*sizeof(typename T::value_type))
  );
  // put data into send buffer
  const int numElems = data.size();
  msg->pack(numElems);
  for(typename T::const_iterator iter=data.begin();
      iter!=data.end();
      iter++){
    msg->pack(*iter);
  }

  // send size
  int size = msg->size();
  console.Debug() << "TML_Comm::send_gather_packed: sending data size...\n";
  MPI_Gather(&size,1,MPI_INT,NULL,0,MPI_INT,root,m_comm);
  
  // send data
  console.Debug() << "TML_Comm::send_gather_packed: sending data...\n";
  MPI_Gatherv(
    msg->buffer(), msg->size(),
    GetType(*(msg->buffer())),
    NULL, NULL, NULL, MPI_INT, root, m_comm
  );
  // clean up
  delete msg;
  console.Debug() << "TML_Comm::send_gather_packed: exit\n";
}

template <typename T> 
T TML_Comm::sum_all(const T& data)
{
  T res;

  MPI_Datatype data_type=GetType(data);
  MPI_Allreduce((void*)(&data),(void*)(&res),1,data_type,MPI_SUM,m_comm);

  return res;
}
