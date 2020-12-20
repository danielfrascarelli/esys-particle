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


#include "tml/message/packed_multi_message.h"

//--- I/O ---
#include <iostream>
#include <cstring>
using std::cout;
using std::endl;
using std::flush;

/*!
  Constructor for TML_PackedMultiMessage

  \param comm the MPI communicator
  \param isize initial buffer size per slice, default 64 byte
 */
TML_PackedMultiMessage::TML_PackedMultiMessage(MPI_Comm comm,int isize)
{
  m_comm=comm;
  MPI_Comm_size(m_comm,&m_size);
  m_vbuffersize=isize;
  m_vbuffer=new char[m_vbuffersize*m_size];
  m_position=new int[m_size];
  m_rpos=new int[m_size];
  m_recvcount=new int[m_size];
  // setup initial displacements in buffer
  m_displ=new int[m_size];
  for(int i=0;i<m_size;i++){
    m_displ[i]=i*m_vbuffersize;
    m_position[i]=i*m_vbuffersize;
    m_rpos[i]=0;
  }
  m_position[0]=0;
  m_recvcount[0]=0;
  MPI_Pack_size(1,MPI_INT,m_comm,&m_int_increment);
  MPI_Pack_size(1,MPI_DOUBLE,m_comm,&m_dbl_increment);
}

TML_PackedMultiMessage::~TML_PackedMultiMessage()
{
  delete m_vbuffer;
  delete m_position;
  delete m_recvcount;
  delete m_displ;
  delete m_rpos;
}

/*!
  return a slab
*/
TML_PackedMultiMessageSlab TML_PackedMultiMessage::operator[](int i)
{
  return TML_PackedMultiMessageSlab(this,i);
}


/*!
  Grows the buffer to twice its current size, thus guaranteeing that append works in amortized constant time. Currently grows the buffer homogeneously, i.e. all slices have the same size. 

  \warning no check if there is enough space for the new buffer
*/
void TML_PackedMultiMessage::grow()
{
  char *temp=m_vbuffer;
  int vbs_old=m_vbuffersize;
  m_vbuffersize+=m_vbuffersize;
  m_vbuffer=new char[m_vbuffersize*m_size];
  for(int i=0;i<m_size;i++){
    memcpy((void*)(&m_vbuffer[i*m_vbuffersize]),(void*)(&temp[i*vbs_old]),size_t(m_position[i]-m_displ[i]));
    m_position[i]=m_position[i]+i*vbs_old;;
    m_displ[i]+=m_displ[i];
    m_rpos[i]=m_position[i]-m_displ[i];
  }
  delete temp;
}

/*!
  Grow buffer to a specified size

  \param size the size to grow to
*/
void TML_PackedMultiMessage::growTo(int size)
{
  if(size>m_vbuffersize){
    char *temp=m_vbuffer;
    int vbs_old=m_vbuffersize;
    m_vbuffersize=size;
    m_vbuffer=new char[m_vbuffersize*m_size];
    for(int i=0;i<m_size;i++){
      memcpy((void*)(&m_vbuffer[i*m_vbuffersize]),(void*)(&temp[m_displ[i]]),size_t(m_position[i]-m_displ[i]));
      m_position[i]=m_position[i]+i*(m_vbuffersize-vbs_old);
      m_displ[i]=i*m_vbuffersize;
      m_rpos[i]=m_position[i]-m_displ[i];
    }
    delete temp;
  }  
}

/*!
  clear message buffer, i.e. reset all positions to 0
*/
void TML_PackedMultiMessage::clear()
{
  for(int i=1;i<m_size;i++) {
    m_position[i]=m_displ[i];
    m_rpos[i]=0;
  }
}
  
/*!
  reset single packing posn to 0
*/
void TML_PackedMultiMessage::begin_pack(int i)
{
  m_position[i]=m_displ[i];
}

/*!
  reset single unpacking posn to 0
*/
void TML_PackedMultiMessage::begin_unpack(int i)
{
  m_position[i]=m_displ[i];
}



/*!
  Append an integer to a given slice of the buffer.

  \param i the integer
  \param nslice the nr. of the slice

  \warning No check for overflow
*/
void TML_PackedMultiMessage::append(int i,int nslice)
{
  if((m_position[nslice]-m_displ[nslice])+m_int_increment>m_vbuffersize){ // to small, grow
    grow();
  }
  MPI_Pack(&i,1,MPI_INT,m_vbuffer,m_vbuffersize*m_size,&(m_position[nslice]),m_comm);
  m_rpos[nslice]=m_position[nslice]-m_displ[nslice];
}

/*!
  Append a double to a given slice of the buffer.

  \param d the double
  \param nslice the nr. of the slice

  \warning No check for overflow
*/
void TML_PackedMultiMessage::append(double d,int nslice)
{
  if((m_position[nslice]-m_displ[nslice])+m_dbl_increment>m_vbuffersize){ // to small, grow
    grow();
  }
  MPI_Pack(&d,1,MPI_DOUBLE,m_vbuffer,m_vbuffersize*m_size,&(m_position[nslice]),m_comm);
  m_rpos[nslice]=m_position[nslice]-m_displ[nslice];
}

/*!
  Append a STL-string to a given slice of the buffer.

  \param str the string
  \param nslice the nr. of the slice

  \warning No check for overflow
*/
void TML_PackedMultiMessage::append(const string& str,int nslice)
{
  int len=str.size();
  if((m_position[nslice]-m_displ[nslice])+m_int_increment+len>m_vbuffersize){ // to small, grow
    grow();
  }
  MPI_Pack(&len,1,MPI_INT,m_vbuffer,m_vbuffersize*m_size,&(m_position[nslice]),m_comm);
  MPI_Pack((void *)str.c_str(),len,MPI_CHAR,m_vbuffer,m_vbuffersize*m_size,&(m_position[nslice]),m_comm);
  m_rpos[nslice]=m_position[nslice]-m_displ[nslice];
}

/*!
  Append boolean value to a given slice of the buffer.
  
  \param b the boolean
  \param nslice the nr. of the slice
*/
void TML_PackedMultiMessage::append(bool b, int nslice)
{
  int i;

  if(b) {i=1;}
  else {i=0;}
  append(i,nslice);
}

/*!
  Pops an integer from a given slice of the the buffer, i.e. it pops the last sizeof(MPI_INT) bytes of the buffer, interpreting them as an int. 

  \param nslice the nr. of the slice
  \return the int.

  \warning No check for underflow
*/
int TML_PackedMultiMessage::pop_int(int nslice)
{
  int res;
  MPI_Unpack(m_vbuffer,m_vbuffersize*m_size,&(m_position[nslice]),&res,1,MPI_INT,m_comm);
  return res;  
}

/*!
  Pops a double from a given slice of the the buffer. 

  \param nslice the nr. of the slice
  \return the double.

  \warning No check for underflow
*/
double TML_PackedMultiMessage::pop_double(int nslice)
{
  double res;
  MPI_Unpack(m_vbuffer,m_vbuffersize*m_size,&(m_position[nslice]),&res,1,MPI_DOUBLE,m_comm);

  return res;
}

/*!
  Pops a string from a given slice of the the buffer. 

  \param nslice the nr. of the slice
  \return the string.

  \warning Not implemented
  \todo implement
*/
//string TML_PackedMultiMessage::pop_string(int nslice)
string TML_PackedMultiMessage::pop_string()
{
  return "";
}

/*!
  Pops a boolean from a given slice of the the buffer, i.e. it pops the last sizeof(MPI_INT) bytes of the buffer, interpreting them as a boolean (via pop_int()). 

  \param nslice the nr. of the slice
  \return the boolean value

  \warning No check for underflow
*/
bool TML_PackedMultiMessage::pop_bool(int nslice)
{
  int i=pop_int(nslice);

  return (i==1);  
}
