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

#include "mpisgvbuf.h"
//#include "console.h"

#include <string.h>

//-- CMPIVarSGBufferRoot member functions ----
/*!
  Constructor for CMPISGBufferRoot

  \param comm the MPI communicator
  \param isize initial buffer size per slice, default 16 byte
 */
CMPIVarSGBufferRoot::CMPIVarSGBufferRoot(MPI_Comm comm,int isize):AMPISGBufferRoot(comm)
{
  m_vbuffersize=isize;
  m_vbuffer=new char[m_vbuffersize*m_size];
  m_dummy_vbuffer=new char[m_vbuffersize];
  m_position=new int[m_size];
  m_rpos=new int[m_size];
  m_recvcount=new int[m_size];
  // setup initial displacements in buffer
  m_displ=new int[m_size];
  for(int i=0;i<m_size;i++){
    m_displ[i]=i*m_vbuffersize;
    m_position[i]=i*m_vbuffersize;
  }
  m_position[0]=0;
  m_recvcount[0]=0;
}


CMPIVarSGBufferRoot::~CMPIVarSGBufferRoot()
{
  delete m_vbuffer;
  delete m_dummy_vbuffer;
  delete m_position;
  delete m_recvcount;
  delete m_displ;
  delete m_rpos;
}

/*!
  Grows the buffer to twice its current size, thus guaranteeing that append
  works in amortized constant time. Currently grows the buffer homogeneously,
  i.e. all slices have the same size. 

  \warning no check if there is enough space for the new buffer
*/
void CMPIVarSGBufferRoot::grow()
{
  char *temp=m_vbuffer;
  int vbs_old=m_vbuffersize;
  m_vbuffersize+=m_vbuffersize;
  m_vbuffer=new char[m_vbuffersize*m_size];
  for(int i=0;i<m_size;i++){
    memcpy((void*)(&m_vbuffer[i*m_vbuffersize]),(void*)(&temp[i*vbs_old]),size_t(m_position[i]-m_displ[i]));
    m_position[i]=m_position[i]+i*vbs_old;;
    m_displ[i]+=m_displ[i];
  }
  delete temp;
}

void CMPIVarSGBufferRoot::growTo(int size)
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
    }
    delete temp;
  }  
}

void CMPIVarSGBufferRoot::clear()
{
  for(int i=1;i<m_size;i++) m_position[i]=m_displ[i];
}

/*!
  Send data to the root process, using MPI_Gather and MPI_Gatherv.
  The receive buffer grows to fit the data if neccessary.
*/
void CMPIVarSGBufferRoot::gather()
{
  int i;

  // zero positions (rel.)
  for(i=1;i<m_size;i++) m_position[i]=m_displ[i];
  // get data sizes
  MPI_Gather(m_dummy_vbuffer,1,MPI_INT,m_recvcount,1,MPI_INT,m_rank,m_comm);
  // if buffer to small, grow
  int max_recv=0;
  for(i=1;i<m_size;i++){ // from 1, i.e. excluding the dummy !
    max_recv=(m_recvcount[i]>max_recv) ? m_recvcount[i] : max_recv;
  }
  if(max_recv>m_vbuffersize){
    growTo(max_recv);
  }
  m_recvcount[0]=1;
  // get data
  MPI_Gatherv(m_dummy_vbuffer,1,MPI_PACKED,m_vbuffer,m_recvcount,m_displ,MPI_PACKED,m_rank,m_comm);
}

/*!
  Send data to all other  members of the communicator, using MPI_Scatter/MPI_Scatterv
*/
void CMPIVarSGBufferRoot::scatter()
{
  //console.XDebug() << "CMPIVarSGBufferRoot::scatter()\n";
  //console.XDebug() << "m_vbuffersize " << m_vbuffersize << "\n";
  // distribute buffer size
  for(int j=1;j<m_size;j++){
    m_rpos[j]=m_position[j]-m_displ[j];
    //console.XDebug() << "posn: " << j << " , " << m_rpos[j] << "\n";
  }
  m_rpos[0]=1;
  MPI_Scatter((void *)m_rpos,1,MPI_INT,&m_ndummy,1,MPI_INT,m_rank,m_comm);
  // send the actual data
  MPI_Scatterv(m_vbuffer,m_rpos,m_displ,MPI_PACKED,m_dummy_vbuffer,1,MPI_PACKED,m_rank,m_comm); 
}

/*!
  Append an integer to a given slice of the buffer.

  \param i the integer
  \param nslice the nr. of the slice

  \warning No check for overflow
*/
void CMPIVarSGBufferRoot::append(int i,int nslice)
{
  if((m_position[nslice]-m_displ[nslice])+m_int_increment>m_vbuffersize){ // to small, grow
    grow();
  }
  MPI_Pack(&i,1,MPI_INT,m_vbuffer,m_vbuffersize*m_size,&(m_position[nslice]),m_comm);
}


/*!
  Append a double to a given slice of the buffer.

  \param d the double
  \param nslice the nr. of the slice

  \warning No check for overflow
*/
void CMPIVarSGBufferRoot::append(double d,int nslice)
{
  if((m_position[nslice]-m_displ[nslice])+m_dbl_increment>m_vbuffersize){ // to small, grow
    grow();
  }
  MPI_Pack(&d,1,MPI_DOUBLE,m_vbuffer,m_vbuffersize*m_size,&(m_position[nslice]),m_comm);
}

/*!
  Append a C-string to a given slice of the buffer.

  \param str the string
  \param nslice the nr. of the slice

  \warning No check for overflow
*/
void CMPIVarSGBufferRoot::append(const char* str,int nslice)
{
  int len=strlen(str);
  //console.XDebug()<< "append string,len " << str  << " , " <<  len << "\n";
  if((m_position[nslice]-m_displ[nslice])+m_int_increment+len>m_vbuffersize){ // to small, grow
    grow();
  }
  MPI_Pack(&len,1,MPI_INT,m_vbuffer,m_vbuffersize*m_size,&(m_position[nslice]),m_comm);
  MPI_Pack((void *)str,len,MPI_CHAR,m_vbuffer,m_vbuffersize*m_size,&(m_position[nslice]),m_comm);
}

/*!
  Pops an integer from a given slice of the the buffer, i.e. 
  it pops the last sizeof(MPI_INT) bytes of the buffer, interpreting
  them as an int. 

  \param nslice the nr. of the slice
  \return the int.

  \warning No check for underflow
*/
int CMPIVarSGBufferRoot::pop_int(int nslice)
{
  int res;
  MPI_Unpack(m_vbuffer,m_vbuffersize*m_size,&(m_position[nslice]),&res,1,MPI_INT,m_comm);
  return res;  
}

/*!
  Pops an double from a given slice of the the buffer. 

  \param nslice the nr. of the slice
  \return the double.

  \warning No check for underflow
*/
double CMPIVarSGBufferRoot::pop_double(int nslice)
{
  double res;
  MPI_Unpack(m_vbuffer,m_vbuffersize*m_size,&(m_position[nslice]),&res,1,MPI_DOUBLE,m_comm);

  return res;
}

void CMPIVarSGBufferRoot::pop_doubles(int nslice, double *dbl, int ndb)
{
  MPI_Unpack(m_vbuffer,m_vbuffersize*m_size,&(m_position[nslice]),dbl,ndb,MPI_DOUBLE,m_comm);
}

//-- CMPIVarSGBufferLeaf member functions ----  

/*!
  Constuctor for CMPISGBufferLeaf

  \param comm the MPI communicator
  \param root rank of the root process
  \param isize initial size of the communication buffer, default 16
*/
CMPIVarSGBufferLeaf::CMPIVarSGBufferLeaf(MPI_Comm comm,int root,int isize):AMPISGBufferLeaf(comm,root)
{
  m_vbuffersize=isize;
  m_vbuffer=new char[m_vbuffersize];
  m_position=0;
}

CMPIVarSGBufferLeaf::~CMPIVarSGBufferLeaf()
{
  delete m_vbuffer;
}

/*!
  Grows the buffer to twice its current size, thus guaranteeing that append works in amortized constant time.  
*/
void CMPIVarSGBufferLeaf::grow()
{
  char *temp=m_vbuffer;
  m_vbuffersize+=m_vbuffersize;
  m_vbuffer=new char[m_vbuffersize];
  memcpy((void*)m_vbuffer,(void*)temp,size_t(m_position));
  delete temp;
}

/*!
  Grows the buffer to a given size. If the buffer is already larger
  that the given size, nothing is done. Used by receiveFrom.
  
  \param size size to which the buffer is grown
*/
void CMPIVarSGBufferLeaf::growTo(int size)
{
  if(size>m_vbuffersize){
    char *temp=m_vbuffer;
    m_vbuffersize=size;
    m_vbuffer=new char[m_vbuffersize];
    memcpy((void*)m_vbuffer,(void*)temp,size_t(m_position));
    delete temp;
  }
}

void CMPIVarSGBufferLeaf::clear()
{
  m_position=0;  
}

/*!
  Send data to the root process, using MPI_Gather/MPI_Gatherv
*/
void CMPIVarSGBufferLeaf::send()
{
  // send size of data
  MPI_Gather(&m_position,1,MPI_INT,NULL,0,MPI_INT,m_root,m_comm);
  // send data
  MPI_Gatherv(m_vbuffer,m_position,MPI_PACKED,NULL,NULL,NULL,MPI_PACKED,m_root,m_comm);
}

/*!
  Receive data from root process, using MPI_Scatter/MPI_Scatterv.
  The buffer grows to fit the data, if neccesary. 
*/
void CMPIVarSGBufferLeaf::receive()
{
  //console.XDebug() <<  "CMPIVarSGBufferLeaf::receive()\n";
  // get size of data
  MPI_Scatter(NULL,0,MPI_INT,&m_data_size,1,MPI_INT,m_root,m_comm);
  //console.XDebug() <<  "recv. data size : " << m_data_size << "\n";
  // if buffer to small, grow
  if(m_data_size>m_vbuffersize){
    growTo(m_data_size);
  }
  //console.XDebug() << "buffersize : " << m_vbuffersize << "\n";
  // get data
  MPI_Scatterv(NULL,NULL,NULL,MPI_PACKED,m_vbuffer,m_data_size,MPI_PACKED,m_root,m_comm);
}

/*!
  Append an integer to the buffer. If necessary, the buffer is enlarged.

  \warning currently does not check if there is enough free space to allocate larger buffer
*/
void CMPIVarSGBufferLeaf::append(int i)
{
  if(m_position+m_int_increment>m_vbuffersize){ // to small, grow
    grow();
  }
  MPI_Pack(&i,1,MPI_INT,m_vbuffer,m_vbuffersize,&m_position,m_comm);
}

/*!
  Append a double to the buffer. If necessary, the buffer is enlarged.

  \warning currently does not check if there is enough free space to allocate larger buffer
  \sa grow()
*/
void CMPIVarSGBufferLeaf::append(double d)
{
  if(m_position+m_dbl_increment>m_vbuffersize){ // to small, grow
    grow();
  }
  MPI_Pack(&d,1,MPI_DOUBLE,m_vbuffer,m_vbuffersize,&m_position,m_comm);
}

/*!
  Append a string to the buffer. If necessary, the buffer is enlarged.

  \warning currently does not check if there is enough free space to allocate larger buffer
  \sa grow()
*/
void CMPIVarSGBufferLeaf::append(const char* str)
{
  int len=strlen(str);
  if(m_position+m_int_increment+len<m_vbuffersize){
    grow();
  }
  MPI_Pack(&len,1,MPI_INT,m_vbuffer,m_vbuffersize,&m_position,m_comm);
  MPI_Pack((void *)str,len,MPI_CHAR,m_vbuffer,m_vbuffersize,&m_position,m_comm);
}

/*!
  Pops an integer from the buffer, i.e. it pops the last sizeof(MPI_INT)
  bytes of the buffer, interpreting them as an int. 

  \warning No check for underflow
  \return the int.
 */
int CMPIVarSGBufferLeaf::pop_int()
{
  int res;
  MPI_Unpack(m_vbuffer,m_vbuffersize,&m_position,&res,1,MPI_INT,m_comm);

  return res;
}


/*! 
  Pops a double from the buffer.
  \warning No check for underflow
  \return the double.
  \sa CMPIBuffer::pop_int()
*/
double CMPIVarSGBufferLeaf::pop_double()
{
  double res;
  MPI_Unpack(m_vbuffer,m_vbuffersize,&m_position,&res,1,MPI_DOUBLE,m_comm);

  return res;
}

void CMPIVarSGBufferLeaf::pop_doubles(double *dbl, int ndb)
{
  MPI_Unpack(m_vbuffer,m_vbuffersize,&m_position,dbl,ndb,MPI_DOUBLE,m_comm);
}

/*! 
  Pops a string from the buffer. The first for bytes are interpreted
  as int, giving the length of the string (without terminating '\0'),
  the rest as the characters.

  \warning no consistency check, i.e. it is not checked if the length
  is smaller than the buffersize.
  \return the string.
  \sa CVarMPISingle::pop_int()
*/
std::string CMPIVarSGBufferLeaf::pop_string()
{
  int len = 0;
  MPI_Unpack(m_vbuffer,m_vbuffersize,&m_position,&len,1,MPI_INT,m_comm);
  char *res=new char[len+1]; // +1 for terminating '\0'
  //console.XDebug()<< "pop string,len " << len  << "\n";
  MPI_Unpack(m_vbuffer,m_vbuffersize,&m_position,res,len,MPI_CHAR,m_comm);
  res[len]='\0';

  std::string poppedString = res;
  delete [] res;
  return poppedString;
}
