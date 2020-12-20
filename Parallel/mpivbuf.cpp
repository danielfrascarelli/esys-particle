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

//--- Project includes ---
#include "mpivbuf.h"
#include "console.h"

//--- system includes ---
#include <cstdlib> // for exit()
#include <string.h>

/*!
  Constructor. Allocates the buffer and sets the MPI communicator to be used for send/receive operations. If the initial buffer size is not given a buffer of initial size 16 is allocated. 

  \param comm the MPI communicator
  \param s the initial size of the buffer, defaults to 16
*/ 
CVarMPIBuffer::CVarMPIBuffer(MPI_Comm comm,int size):AMPIBufferPP(comm)
{
  m_buffersize=size;
  m_buffer=new char[m_buffersize];
  for(int i=0;i<m_buffersize;i++) m_buffer[i]=char(0);
  // get increments for int,double
  MPI_Pack_size(1,MPI_INT,m_comm,&m_int_increment);
  MPI_Pack_size(1,MPI_DOUBLE,m_comm,&m_dbl_increment);
  m_position=0;
  m_lock=false;
}

CVarMPIBuffer::~CVarMPIBuffer()
{
  delete [] m_buffer;
}

/*! 
  Sends the contents of the buffer to a given destination. There are actually two messages sent, the first one announces the size of the buffer so the buffer on the receiving end can be grown if necessary, the second one does the transfer of the data. It uses buffered sends (MPI_Bsend) and it thus deadlock-save 

  \param dest the rank of the destination process in the current communicator 
  \param tag the message tag

  \warning It is not checked if the destination actually exists. 
  \warning Overlapping requests can deadlock.
*/
void CVarMPIBuffer::sendTo(int dest,int tag)
{
  // send message size (m_position), not buffer size !
  MPI_Send((void*)&m_position,1,MPI_INT,dest,tag,m_comm);
  MPI_Send(m_buffer,m_position,MPI_PACKED,dest,tag+2048,m_comm);
}

/*! 
  Nonblocking version of CVarMPIBuffer::sendTo. Uses MPI_Isend and is (should be) thus deadlock-safe.

  \param dest the rank of the destination process in the current communicator 
  \param tag the message tag

  \warning It is not checked if the destination actually exists.
*/
void CVarMPIBuffer::NBsendTo(int dest,int tag)
{
  MPI_Request req[2];
  MPI_Status stat[2];
  // send message size (m_position), not buffer size !
  MPI_Isend((void*)&m_position,1,MPI_INT,dest,tag,m_comm,&(req[0]));
  MPI_Isend(m_buffer,m_position,MPI_PACKED,dest,tag+2048,m_comm,&(req[1]));
  MPI_Waitall(2,req,stat);
}

/*!
  Initate send,lock buffer and immediately return (equivalent to MPI_Isend)
 
  \param dest the rank of the destination process in the current communicator 
  \param tag the message tag
 */
void CVarMPIBuffer::initSendTo(int dest,int tag)
{
  if(m_lock){ // change to exeption
    console.Critical() << "trying to use locked buffer, aborting !\n";
    exit(1);
  } else {
    MPI_Isend((void*)&m_position,1,MPI_INT,dest,tag,m_comm,&(m_req[0]));
    MPI_Isend(m_buffer,m_position,MPI_PACKED,dest,tag+2048,m_comm,&(m_req[1]));
    m_lock=true;
  }
}

/*!
  Wait for completion of transaction on this buffer. If completed, unlock buffer and return
 */
void CVarMPIBuffer::wait()
{
  MPI_Waitall(2,m_req,m_stat);
  m_lock=false;
}

/*!
  Receives a message from a given source and stores it in the buffer.The size of the buffer is automatically adjusted so it will be big enough to fit the message. For this reason 2 messages are received, the first one for the size of the data, the second one for the data. If no source and no tag are given, any message from any source is accepted.

  \param src rank of the sender in the current communicator, defaults to MPI_ANY_SOURCE
  \param tag the message tag, defaults to MPI_ANY_TAG
  
*/
void CVarMPIBuffer::receiveFrom(int src,int tag)
{
  int size;
  MPI_Recv((void*)&size,1,MPI_INT,src,tag,m_comm,&m_status);
  if(size>m_buffersize){
    growTo(size);
  }
  int n_src=m_status.MPI_SOURCE; // necessary in case of overlapping MPI_ANY_SOURCE transactions
  MPI_Recv(m_buffer,size,MPI_PACKED,n_src,tag+2048,m_comm,&m_status) ;
}

/*!
  Grows the buffer to a given size. If the buffer is already larger that the given size, nothing is done. Used by receiveFrom.
  
  \param size size to which the buffer is grown
*/
void CVarMPIBuffer::growTo(int size)
{
  if(size>m_buffersize){
    char *temp = new char[size];
    memcpy((void*)temp,m_buffer,size_t(m_position));
    delete [] m_buffer;
    m_buffer = temp;    
    m_buffersize = size;
  }
}

/*!
  Grows the buffer to twice its current size, thus guaranteeing that append works in amortized constant time.  
*/
void CVarMPIBuffer::grow()
{
  growTo(m_buffersize*2);
}

/*!
  Append an integer to the buffer. If necessary, the buffer is enlarged.

  \warning currently does not check if there is enough free space to allocate larger buffer
*/
void CVarMPIBuffer::append(int i)
{
  if(m_position+m_int_increment>m_buffersize){ // to small, grow
    grow();
  }
  MPI_Pack(&i,1,MPI_INT,m_buffer,m_buffersize,&m_position,m_comm);
}

/*!
  Append a double to the buffer. If necessary, the buffer is enlarged.

  \warning currently does not check if there is enough free space to allocate larger buffer
  \sa grow()
*/
void CVarMPIBuffer::append(double d)
{
  if(m_position+m_dbl_increment>m_buffersize){ // to small, grow
    grow();
  }
  MPI_Pack(&d,1,MPI_DOUBLE,m_buffer,m_buffersize,&m_position,m_comm);
}

/*!
  Append a string to the buffer. If necessary, the buffer is enlarged.

  \warning currently does not check if there is enough free space to allocate larger buffer
  \sa grow()
*/
void CVarMPIBuffer::append(const char* str)
{
  int len=strlen(str);
  while (m_position+m_int_increment+len>m_buffersize) {
    grow();
  }
  MPI_Pack(&len,1,MPI_INT,m_buffer,m_buffersize,&m_position,m_comm);
  MPI_Pack((void *)str,len,MPI_CHAR,m_buffer,m_buffersize,&m_position,m_comm);
}

/*!
  Pops an integer from the buffer, i.e. it pops the last sizeof(MPI_INT) bytes of the buffer, interpreting them ans an int. 

  \warning No check for underflow
  \return The int.
 */
int CVarMPIBuffer::pop_int()
{
  int res;
  MPI_Unpack(m_buffer,m_buffersize,&m_position,&res,1,MPI_INT,m_comm);

  return res; 
}

/*! 
  Pops a double from the buffer.
  \warning No check for underflow
  \return the double.
  \sa CMPIBuffer::pop_int()
 */
double CVarMPIBuffer::pop_double()
{
  double res;
  MPI_Unpack(m_buffer,m_buffersize,&m_position,&res,1,MPI_DOUBLE,m_comm);

  return res;
}

/*!
  pop an array of doubles from a buffer
*/
void CVarMPIBuffer::pop_doubles(double *dbl, int ndb)
{
  MPI_Unpack(m_buffer,m_buffersize,&m_position,dbl,ndb,MPI_DOUBLE,m_comm);
}

/*! 
  Pops a string from the buffer. The first for bytes are interpreted as int, giving the length of the string (without terminating '\0'), the rest as the characters.

  \warning no consistency check, i.e. it is not checked if the length is smaller than the buffersize.
  \return the double.
  \sa CVarMPISingle::pop_int()
*/
std::string CVarMPIBuffer::pop_string()
{
  int len = 0;
  MPI_Unpack(m_buffer,m_buffersize,&m_position,&len,1,MPI_INT,m_comm);
  char *res=new char[len+1]; // +1 for terminating '\0'
  MPI_Unpack(m_buffer,m_buffersize,&m_position,res,len,MPI_CHAR,m_comm);
  res[len]='\0';

  std::string poppedString = res;
  delete [] res;
  return poppedString;
}

/*!
  Broadcast a message to all members of the communicator. 

  \param root the root of the broadcast
*/
void CVarMPIBuffer::broadcast(int root)
{
  MPI_Bcast(&m_position,1,MPI_INT,root,m_comm);
  MPI_Bcast(m_buffer,m_position,MPI_PACKED,root,m_comm);
}

/*!
  receive broadcast

  \param root the root of the broadcast
*/
void CVarMPIBuffer::receiveBroadcast(int root)
{
  int size;
  MPI_Bcast(&size,1,MPI_INT,root,m_comm);
  if(size>m_buffersize){
    growTo(size);
  }
  MPI_Bcast(m_buffer,size,MPI_PACKED,root,m_comm);
}
