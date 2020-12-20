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


#include "tml/message/packed_message.h"
#include <cstring>

/*!
  Constructor. Allocates the buffer and sets the MPI communicator to be used for pack/unpack operations. If the initial buffer size is not given a buffer of initial size 64 is allocated.

  \param comm the MPI communicator
  \param s the initial size of the buffer
*/ 
TML_Packed_Message::TML_Packed_Message(MPI_Comm comm,unsigned int s)
{
  m_comm=comm;
  m_buffersize=s;
  m_buffer=new char[m_buffersize];
  m_pack_pos=0;
  m_unpack_pos=0;
  MPI_Pack_size(1,MPI_INT,m_comm,&m_int_increment);
  MPI_Pack_size(1,MPI_DOUBLE,m_comm,&m_dbl_increment);
}

/*!
  Destructor. Free buffer.
*/
TML_Packed_Message::~TML_Packed_Message()
{
  delete [] m_buffer;
}

/*!
  Grows the buffer to a given size. If the buffer is already larger that the given size, nothing is done. Used by receiveFrom.
  
  \param size size to which the buffer is grown
*/
void TML_Packed_Message::growTo(int size)
{
  if(size>m_buffersize){
    char *temp = new char[size];
    memcpy((void*)temp,m_buffer,size_t(m_pack_pos));
    delete [] m_buffer;
    m_buffer = temp;
    m_buffersize = size;
  }
}

/*!
  Grows the buffer to twice its current size, thus guaranteeing that append works in amortized constant time.  
*/
void TML_Packed_Message::grow()
{
  growTo(2*m_buffersize);
}

/*!
  Append an integer to the buffer. If necessary, the buffer is enlarged.

  \warning currently does not check if there is enough free space to allocate larger buffer
*/
void TML_Packed_Message::append(int i)
{
  while(m_pack_pos+m_int_increment>m_buffersize){ // to small, grow
    grow();
  }
  MPI_Pack(&i,1,MPI_INT,m_buffer,m_buffersize,&m_pack_pos,m_comm);
}

/*!
  Append a double to the buffer. If necessary, the buffer is enlarged.

  \warning currently does not check if there is enough free space to allocate larger buffer
  \sa grow()
*/
void TML_Packed_Message::append(double d)
{
  while(m_pack_pos+m_dbl_increment>m_buffersize){ // to small, grow
    grow();
  }
  MPI_Pack(&d,1,MPI_DOUBLE,m_buffer,m_buffersize,&m_pack_pos,m_comm);
}

/*!
  Append a STL string to the buffer. The string is internally handeled by packing the length frist and then the string. If necessary, the buffer is enlarged.

  \warning currently does not check if there is enough free space to allocate larger buffer
  \sa grow()
*/
void TML_Packed_Message::append(const string& str)
{
  int len=str.size();
  while (m_pack_pos+m_int_increment+len>m_buffersize){
    grow();
  }
  MPI_Pack(&len,1,MPI_INT,m_buffer,m_buffersize,&m_pack_pos,m_comm);
  MPI_Pack((void *)str.c_str(),len,MPI_CHAR,m_buffer,m_buffersize,&m_pack_pos,m_comm);
}

/*!
  Append a Vec3 to the message buffer. Calls append(double) per element
*/
void TML_Packed_Message::append(const Vec3& v)
{
  append(v[0]);
  append(v[1]);
  append(v[2]);
}

/*!
  Append a Matrix3 to the message buffer. Calls append(double) per element
*/
void TML_Packed_Message::append(const Matrix3& m)
{
    append(m(0,0));
    append(m(0,1));
    append(m(0,2));
    append(m(1,0));
    append(m(1,1));
    append(m(1,2));
    append(m(2,0));
    append(m(2,1));
    append(m(2,2));
}


/*!
  Append a boolean to the message buffer. The bool gest transported as an int (1/0) because
  MPI doesn't have a native boolean type. Therefore calls append(int).
*/
void TML_Packed_Message::append(bool b)
{
  int i;
  if(b){i=1;}
  else {i=0;}
  append(i);
}

/*!
  Pops an integer from the buffer, i.e. it pops the last sizeof(MPI_INT) bytes of the buffer, interpreting them as an int. 

  \warning No check for underflow
  \return the int.
 */
int TML_Packed_Message::pop_int()
{
  int res;
  MPI_Unpack(m_buffer,m_buffersize,&m_unpack_pos,&res,1,MPI_INT,m_comm);

  return res;
}

/*! 
  Pops a double from the buffer.
  \warning No check for underflow
  \return the double.
  \sa CMPIBuffer::pop_int()
 */
double TML_Packed_Message::pop_double()

{
  double res;
  MPI_Unpack(m_buffer,m_buffersize,&m_unpack_pos,&res,1,MPI_DOUBLE,m_comm);

  return res;
}

/*!
  pop a C-array of doubles from a buffer. Faster than doing multiple pop_double operations

  \param dbl the array
  \param ndb the number of doubles to be popped
  \warning No check for underflow
*/
void TML_Packed_Message::pop_doubles(double *dbl,int ndb)
{
  MPI_Unpack(m_buffer,m_buffersize,&m_unpack_pos,dbl,ndb,MPI_DOUBLE,m_comm);
}

/*! 
  Pops a string from the buffer. The first for bytes are interpreted as int, giving the length of the string (without terminating '\0'), the rest as the characters.

  \warning no consistency check, i.e. it is not checked if the length is smaller than the buffersize.
  \return the double.
  \sa CMPISingle::pop_int()
*/
string TML_Packed_Message::pop_string()
{
  int len;
  char* res;

  MPI_Unpack(m_buffer,m_buffersize,&m_unpack_pos,&len,1,MPI_INT,m_comm);
  res=new char[len+1]; // +1 for terminating '\0'
  MPI_Unpack(m_buffer,m_buffersize,&m_unpack_pos,res,len,MPI_CHAR,m_comm);
  res[len]='\0';
  string str=string(res); 
  delete [] res;

  return str;
}

/*!
  Pop a Vec3 of the buffer. Calls pop_double per element
*/
Vec3 TML_Packed_Message::pop_vec3()
{
  Vec3 res;

  res[0]=pop_double();
  res[1]=pop_double();
  res[2]=pop_double();

  return res;
}

/*!
  Pop a Matrix3 of the buffer. Calls pop_doubles[] 
*/
Matrix3 TML_Packed_Message::pop_matrix3()
{
  double db[9];
  pop_doubles(db, 9); // stress tensor
  Matrix3 res(db);

  return res;
}

/*!
  Pop a boolean value of the buffer. Booleans are transported as int (0/1)
*/
bool TML_Packed_Message::pop_bool()
{
  bool res;
  int i;

  i=pop_int();
  res=(i==1);

  return res;
}
