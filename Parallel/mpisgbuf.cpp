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

#include "mpisgbuf.h"
#include <cstdlib> // NULL
#include <string.h>

//-- AMPISGBufferRoot member functions----

/*!
  Constructor for AMPISGBufferRoot

  \param comm the MPI communicator
 */
AMPISGBufferRoot::AMPISGBufferRoot(MPI_Comm comm)
{
  m_comm=comm;
  MPI_Comm_rank(comm,&m_rank);
  MPI_Comm_size(comm,&m_size);
  MPI_Pack_size(1,MPI_INT,m_comm,&m_int_increment);
  MPI_Pack_size(1,MPI_DOUBLE,m_comm,&m_dbl_increment);
}

void AMPISGBufferRoot::append(const Vec3 &v,int nslice)
{
  append(v.X(),nslice);
  append(v.Y(),nslice);
  append(v.Z(),nslice);
}

Vec3 AMPISGBufferRoot::pop_vector(int nslice)
{
  double x[3] ;
  pop_doubles(nslice,x,3);
  return Vec3(x[0],x[1],x[2]);
}

//-- CMPISGBufferRoot member functions----
/*!
  Constructor for CMPISGBufferRoot

  \param comm the MPI communicator
  \param buffersize buffer size per slice
 */
CMPISGBufferRoot::CMPISGBufferRoot(MPI_Comm comm,int buffersize):AMPISGBufferRoot(comm)
{
  m_buffersize=buffersize;
  m_buffer=new char[m_buffersize*m_size];
  m_dummy_buffer=new char[m_buffersize];
  m_position=new int[m_size];
  for(int i=0;i<m_size;i++){
    m_position[i]=i*m_buffersize;
  }
}

CMPISGBufferRoot::~CMPISGBufferRoot()
{
  delete [] m_buffer;
  delete [] m_dummy_buffer;
  delete [] m_position;
}

void CMPISGBufferRoot::clear()
{
  for(int i=0;i<m_size;i++){
    m_position[i]=i*m_buffersize;
  }
}

/*!
  Get data from all other members of the communicator, using MPI_Gather 
*/
void CMPISGBufferRoot::gather()
{
  MPI_Gather(m_dummy_buffer,m_buffersize,MPI_PACKED,m_buffer,m_buffersize,MPI_PACKED,m_rank,m_comm);
}

/*!
  Send data to all other  members of the communicator, using MPI_Scatter
*/
void CMPISGBufferRoot::scatter()
{
  MPI_Scatter(m_buffer,m_buffersize,MPI_PACKED,m_dummy_buffer,m_buffersize,MPI_PACKED,m_rank,m_comm);
}

/*!
  Append an integer to a given slice of the buffer.

  \param i the integer
  \param nslice the nr. of the slice

  \warning No check for overflow
*/
void CMPISGBufferRoot::append(int i,int nslice)
{
  MPI_Pack(&i,1,MPI_INT,m_buffer,m_buffersize*m_size,&(m_position[nslice]),m_comm);
}

/*!
  Append an double to a given slice of the buffer.

  \param d the double
  \param nslice the nr. of the slice

  \warning No check for overflow
*/
void CMPISGBufferRoot::append(double d,int nslice)
{
  MPI_Pack(&d,1,MPI_DOUBLE,m_buffer,m_buffersize*m_size,&(m_position[nslice]),m_comm);
}
/*!
  Append an C string (char*) to a given slice of the buffer.

  \param str the string
  \param nslice the nr. of the slice

  \warning No check for overflow
*/
void CMPISGBufferRoot::append(const char* str,int nslice)
{
  int len=strlen(str);
  MPI_Pack(&len,1,MPI_INT,m_buffer,m_buffersize*m_size,&(m_position[nslice]),m_comm);
  MPI_Pack((void *)str,len,MPI_CHAR,m_buffer,m_buffersize*m_size,&(m_position[nslice]),m_comm);
}

/*!
  Pops an integer from a given slice of the the buffer, i.e. it pops the last sizeof(MPI_INT) bytes of the buffer, interpreting them as an int. 

  \param nslice the nr. of the slice
  \return the int.

  \warning No check for underflow
*/
int CMPISGBufferRoot::pop_int(int nslice)
{
  int res;
  MPI_Unpack(m_buffer,m_buffersize*m_size,&(m_position[nslice]),&res,1,MPI_INT,m_comm);
  
  return res;  
}


/*!
  Pops an double from a given slice of the the buffer. 

  \param nslice the nr. of the slice
  \return the double.

  \warning No check for underflow
*/
double CMPISGBufferRoot::pop_double(int nslice)
{
  double res;
  MPI_Unpack(m_buffer,m_buffersize,&(m_position[nslice]),&res,1,MPI_DOUBLE,m_comm);

  return res;
}

void CMPISGBufferRoot::pop_doubles(int nslice, double *dbl, int ndb)
{
  MPI_Unpack(m_buffer,m_buffersize,&(m_position[nslice]),dbl,ndb,MPI_DOUBLE,m_comm);
}

//-- AMPISGBufferLeaf member functions----

/*!
  Constuctor for AMPISGBufferLeaf

  \param comm the MPI communicator
  \param root rank of the root process
*/
AMPISGBufferLeaf::AMPISGBufferLeaf(MPI_Comm comm, int root):AMPIBuffer(comm)
{
  m_root=root;
  MPI_Pack_size(1,MPI_INT,m_comm,&m_int_increment);
  MPI_Pack_size(1,MPI_DOUBLE,m_comm,&m_dbl_increment);
}

//-- CMPISGBufferLeaf member functions----

/*!
  Constuctor for CMPISGBufferLeaf

  \param comm the MPI communicator
  \param root rank of the root process
  \param buffersize size of the communication buffer
*/
CMPISGBufferLeaf::CMPISGBufferLeaf(MPI_Comm comm,int root,int buffersize):AMPISGBufferLeaf(comm,root)
{
  m_buffersize=buffersize;
  m_buffer = new char[m_buffersize];
  for (int i = 0; i < m_buffersize; i++)
  {
    m_buffer[i] = '\0';
  }
  m_position=0;
}

CMPISGBufferLeaf::~CMPISGBufferLeaf()
{
  delete [] m_buffer;
}

void CMPISGBufferLeaf::clear()
{
  m_position=0;
}

/*!
  Send data to the root process, using MPI_Gather
*/
void CMPISGBufferLeaf::send()
{
  MPI_Gather(m_buffer,m_buffersize,MPI_PACKED,NULL,0,MPI_PACKED,m_root,m_comm);
}

/*!
  Receive data from root process, using MPI_Scatter
*/
void CMPISGBufferLeaf::receive()
{
  MPI_Scatter(NULL,0,MPI_PACKED,m_buffer,m_buffersize,MPI_PACKED,m_root,m_comm); 
}

/*!
  Append an integer to the buffer.

  \warning No check for overflow
*/
void CMPISGBufferLeaf::append(int i)
{
  MPI_Pack(&i,1,MPI_INT,m_buffer,m_buffersize,&m_position,m_comm);
}

/*!
  Append a double to the buffer.

  \warning No check for overflow
*/
void CMPISGBufferLeaf::append(double d)
{
  MPI_Pack(&d,1,MPI_DOUBLE,m_buffer,m_buffersize,&m_position,m_comm);
}

/*!
  Append a string to the buffer. The string appended is a normal (zero-terminated) C-string, but is internally handeled by packing the length frist and then the string.

  \warning No check for overflow
*/
void CMPISGBufferLeaf::append(const char* str)
{
  int len=strlen(str);
  MPI_Pack(&len,1,MPI_INT,m_buffer,m_buffersize,&m_position,m_comm);
  MPI_Pack((void *)str,len,MPI_CHAR,m_buffer,m_buffersize,&m_position,m_comm);
}

/*!
  Pops an integer from the buffer, i.e. it pops the last sizeof(MPI_INT) bytes of the buffer, interpreting them as an int. 

  \warning No check for underflow
  \return the int.
 */
int CMPISGBufferLeaf::pop_int()
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
double CMPISGBufferLeaf::pop_double()
{
  double res;
  MPI_Unpack(m_buffer,m_buffersize,&m_position,&res,1,MPI_DOUBLE,m_comm);

  return res;
}

void CMPISGBufferLeaf::pop_doubles(double *dbl, int ndb)
{
  MPI_Unpack(m_buffer,m_buffersize,&m_position,dbl,ndb,MPI_DOUBLE,m_comm);
}
/*! 
  Pops a string from the buffer. The first for bytes are interpreted as int, giving the length of the string (without terminating '\0'), the rest as the characters.

  \warning no consistency check, i.e. it is not checked if the length is smaller than the buffersize.
  \return the double.
  \sa CMPISingle::pop_int()
*/
std::string CMPISGBufferLeaf::pop_string()
{
  int len;
  char* res;

  MPI_Unpack(m_buffer,m_buffersize,&m_position,&len,1,MPI_INT,m_comm);
  res=new char[len+1]; // +1 for terminating '\0'
  MPI_Unpack(m_buffer,m_buffersize,&m_position,res,len,MPI_CHAR,m_comm);
  res[len]='\0';

  std::string poppedString = res;
  delete [] res;
  return poppedString;
}
