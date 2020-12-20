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

/*!
  Shift a single value along a cartesian direction. Gets the source/destination with MPI_Cart_shift and moves the data with sendrecv. If the direction is out of range, nothing is done.

  \param send_data data to be sent
  \param recv data data to be received
  \param dir direction
  \param dist the shift distance along dir
  \param tag the message tag
*/
template <typename T,typename P> 
void TML_CartComm::shift(T send_data,P& recv_data,int dir,int dist,int tag)
{
  int source,dest;

  if(dir<m_ndims){
    MPI_Cart_shift(m_comm,dir,dist,&source,&dest);
    sendrecv(send_data,recv_data,dest,source,tag);
  }
}

/*!
  Shift C-arrays of known size value along a cartesian direction. Gets the source/destination with MPI_Cart_shift and moves the data with sendrecv. If the direction is out of range, nothing is done.
  
  \param send_data data to be sent
  \param size of arry to be sent
  \param recv data data to be received
  \param size of arry to be received
  \param dir direction
  \param dist the shift distance along dir
  \param tag the message tag
*/
template <typename T,typename P> 
void TML_CartComm::shift_array(T* send_data,int send_count,P* recv_data,int recv_count,int dir,int dist,int tag)
{
  int source,dest;

  if(dir<m_ndims){
    MPI_Cart_shift(m_comm,dir,dist,&source,&dest);
    sendrecv(send_data,send_count,recv_data,recv_count,dest,source,tag);
  }
}


/*!
  Shift STL containers along a cartesian direction. Gets the source/destination with MPI_Cart_shift and moves the data with sendrecv. If the direction is out of range, nothing is done.

  \param send_data data to be sent
  \param recv data data to be received
  \param dir direction
  \param dist the shift distance along dir
  \param tag the message tag
*/
template <typename T,typename P> 
void TML_CartComm::shift_cont(T send_data,P& recv_data,int dir,int dist,int tag)
{
  int source,dest;

  if(dir<m_ndims){
    MPI_Cart_shift(m_comm,dir,dist,&source,&dest);
    sendrecv_cont(send_data,recv_data,dest,source,tag);
  }
}

/*!
  Shift a single packable object along a cartesian direction. Gets the source/destination with MPI_Cart_shift and moves the data with sendrecv.  If the direction is out of range, nothing is done.

  \param send_data data to be sent
  \param recv data data to be received
  \param dir direction
  \param dist the shift distance along dir
  \param tag the message tag
*/
template <typename T,typename P> 
void TML_CartComm::shift_packed(T send_data,P& recv_data,int dir,int dist,int tag)
{
  int source,dest;

  if(dir<m_ndims){
    MPI_Cart_shift(m_comm,dir,dist,&source,&dest);
    sendrecv_packed(send_data,recv_data,dest,source,tag);
  }
}

/*!
  Shift C-arrays of packable objects of known size value along a cartesian direction. Gets the source/destination with MPI_Cart_shift and moves the data with sendrecv.  If the direction is out of range, nothing is done.
  
  \param send_data data to be sent
  \param size of arry to be sent
  \param recv data data to be received
  \param size of arry to be received
  \param dir direction
  \param dist the shift distance along dir
  \param tag the message tag
*/
template <typename T,typename P> 
void TML_CartComm::shift_array_packed(T* send_data,int send_count,P* recv_data,int recv_count,int dir,int dist,int tag)
{
  int source,dest;

  if(dir<m_ndims){
    MPI_Cart_shift(m_comm,dir,dist,&source,&dest);
    sendrecv_array_packed(send_data,send_count,recv_data,recv_count,dest,source,tag);
  }
}


/*!
  Shift STL containers of packable objects along a cartesian direction. Gets the source/destination with MPI_Cart_shift and moves the data with sendrecv. If the direction is out of range, nothing is done.

  \param send_data data to be sent
  \param recv data data to be received
  \param dir direction
  \param dist the shift distance along dir
  \param tag the message tag
*/
template <typename T,typename P> 
void TML_CartComm::shift_cont_packed(T send_data,P& recv_data,int dir,int dist,int tag)
{
  int source,dest;

  if(dir<m_ndims){
    MPI_Cart_shift(m_comm,dir,dist,&source,&dest);
    sendrecv_cont_packed(send_data,recv_data,dest,source,false,tag);
  }
}
