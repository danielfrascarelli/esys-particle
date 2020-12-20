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

#include "mpibuf.h"
#include <string.h>

AMPIBufferPP::AMPIBufferPP(MPI_Comm comm):AMPIBuffer(comm)
{}

void AMPIBuffer::append(const Vec3 & V)
{
  append(V.X()) ;
  append(V.Y()) ;
  append(V.Z()) ;
}

Vec3 AMPIBuffer::pop_vector()
{
  double x=pop_double();
  double y=pop_double();
  double z=pop_double();
  return Vec3(x,y,z);
}  

/*!
  Constructor. Allocates the buffer and sets the MPI communicator to be used for send/receive operations.

  \param comm the MPI communicator
  \param s the size of the buffer
 */ 
CMPIBuffer::CMPIBuffer(MPI_Comm comm,int s):AMPIBufferPP(comm)
{
  m_buffersize=s;
  m_buffer=new char[m_buffersize];
  // get increments for int,double
  MPI_Pack_size(1,MPI_INT,m_comm,&m_int_increment);
  MPI_Pack_size(1,MPI_DOUBLE,m_comm,&m_dbl_increment);
  m_position=0;
}

CMPIBuffer::~CMPIBuffer()
{
  delete m_buffer;
}

/*! 
  Sends the contents of the buffer to a given destination.
  \param dest the rank of the destination process in the current communicator 
  \param tag the message tag

  \warning It is not checked if the destination actually exists.
*/
void CMPIBuffer::sendTo(int dest,int tag)
{
  MPI_Send(m_buffer,m_position,MPI_PACKED,dest,tag,m_comm);
}

/*!
  Recieves a message from a given source and stores it in the buffer. It is assumed that the buffer is large enough the take the message. If no source and no tag are given, any message from any source is accepted.

  \warning No check if the buffer is big enough
  \param src rank of the sender in the current communicator, defaults to MPI_ANY_SOURCE
  \param tag the message tag, defaults to MPI_ANY_TAG
  
*/
void CMPIBuffer::receiveFrom(int src,int tag)
{
  MPI_Recv(m_buffer,m_buffersize,MPI_PACKED,src,tag,m_comm,&m_status) ;
  m_position=0;
}

/*!
  Append an integer to the buffer.

  \warning No check for overflow
*/
void CMPIBuffer::append(int i)
{
  MPI_Pack(&i,1,MPI_INT,m_buffer,m_buffersize,&m_position,m_comm);
}

/*!
  Append a double to the buffer.

  \warning No check for overflow
*/
void CMPIBuffer::append(double d)
{
  MPI_Pack(&d,1,MPI_DOUBLE,m_buffer,m_buffersize,&m_position,m_comm);
}

/*!
  Append a string to the buffer. The string appended is a normal (zero-terminated) C-string, but is internally handeled by packing the length frist and then the string.

  \warning No check for overflow
*/
void CMPIBuffer::append(const char* str)
{
  int len=strlen(str);
  MPI_Pack(&len,1,MPI_INT,m_buffer,m_buffersize,&m_position,m_comm);
  MPI_Pack((void *)str,len,MPI_CHAR,m_buffer,m_buffersize,&m_position,m_comm);
}
  

/*!
  Append an integer to the buffer with overflow check. If the buffer is big enough the integer is appended, if not nothing is done.

  \return true if the append succeded, false otherwise
*/
bool CMPIBuffer::append_checked(int i)
{
  bool res;

  if(m_position+m_int_increment<m_buffersize){
    MPI_Pack(&i,1,MPI_INT,m_buffer,m_buffersize,&m_position,m_comm);
    res=true;
  } else {
    res=false;
  }
  return res;
}

/*!
  Append a double to the buffer with overflow check.

  \sa CMPIBuffer::append_checked(int i)
*/
bool CMPIBuffer::append_checked(double d)
{
  bool res;

  if(m_position+m_dbl_increment<m_buffersize){
    MPI_Pack(&d,1,MPI_DOUBLE,m_buffer,m_buffersize,&m_position,m_comm);
    res=true;
  } else {
    res=false;
  }
  return res;
}

/*!
  Pops an integer from the buffer, i.e. it pops the last sizeof(MPI_INT) bytes of the buffer, interpreting them as an int. 

  \warning No check for underflow
  \return the int.
 */
int CMPIBuffer::pop_int()
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
double CMPIBuffer::pop_double()
{
  double res;
  MPI_Unpack(m_buffer,m_buffersize,&m_position,&res,1,MPI_DOUBLE,m_comm);

  return res;
}

/*!
  pop an array of doubles from a buffer
*/
void CMPIBuffer::pop_doubles(double *dbl, int ndb)
{
  MPI_Unpack(m_buffer,m_buffersize,&m_position,dbl,ndb,MPI_DOUBLE,m_comm);
}


/*! 
  Pops a string from the buffer. The first for bytes are interpreted as int, giving the length of the string (without terminating '\0'), the rest as the characters.

  \warning no consistency check, i.e. it is not checked if the length is smaller than the buffersize.
  \return the double.
  \sa CMPISingle::pop_int()
*/
std::string CMPIBuffer::pop_string()
{
  int len = 0;
  MPI_Unpack(m_buffer,m_buffersize,&m_position,&len,1,MPI_INT,m_comm);
  char *res = new char[len+1]; // +1 for terminating '\0'
  MPI_Unpack(m_buffer,m_buffersize,&m_position,res,len,MPI_CHAR,m_comm);
  res[len]='\0';

  std::string poppedString = res;
  delete [] res;
  return poppedString;
}
