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

#ifndef __COMM_H
#define __COMM_H

//--- MPI ---
#include <mpi.h>

//--- project includes ---
#include "Foundation/console.h"

//--- TML includes ---
#include "tml/type/gettype.h"

//--- STL includes ---
#include <vector>
#include <map>
#include <utility>
#include <string>

#undef barrier

using std::vector;
using std::multimap;
using std::pair;
using std::string;

/*!
  \class TML_Comm
  \brief abstract base class for communicator

  \author Steffen Abe
  $Revision$
  $Date$
*/
class TML_Comm
{
protected:
  MPI_Status m_status;
  MPI_Comm m_comm;

public:
  bool isNull() const {return m_comm==MPI_COMM_NULL;};
  int rank() const;
  int size();
  MPI_Comm comm() const {return m_comm;};

  // assignment
  TML_Comm& operator=(const TML_Comm&);
  void setComm(MPI_Comm);

  // construction of new communicators
  TML_Comm();
  TML_Comm(MPI_Comm);
  TML_Comm include(const vector<int>&);
  TML_Comm exclude(const vector<int>&);

  // send/recv for single data
  template <typename T> void send(T,int,int=0);
  template <typename T> void receive(T&,int,int=MPI_ANY_TAG);

  // send/recv for C-arrays with known dimension
  template <typename T> void send_array(T*,int,int,int=0);
  template <typename T> void receive_array(T*,int,int,int=MPI_ANY_TAG);

  // send/recv for STL containers
  template <typename T> void send_cont(const T&,int,int=0);
  template <typename T> void receive_cont(T&,int,int=MPI_ANY_TAG);

  // send/recv for STL containers of packable objects
  template <typename T> void send_cont_packed(T,int,bool,int=0);
  template <typename T> void receive_cont_packed(T&,int,bool,int=MPI_ANY_TAG);

  // sendrecv for single data
  template <typename T,typename P> void sendrecv(T,P&,int,int,int=0);
  
  // sendrecv for C-arrays with known dimension
  template <typename T,typename P> void sendrecv_array(T*,int,P*,int,int,int,int=0);

  // sendrecv for STL containers 
  template <typename T,typename P> void sendrecv_cont(T,P&,int,int,int=0);
  template <typename T> void sendrecv_cont_replace(T&,int,int,int=0);

  // sendrecv for STL containers of packable objects
  template <typename T,typename P> void sendrecv_cont_packed(T,P&,int,int,bool,int=0);
  template <typename T> void sendrecv_cont_packed_replace(T&,int,int,bool,int=0);

  // broadcast
  template <typename T> void broadcast(T);
  template <typename T> void broadcast_array(T*,int);
  template <typename T> void broadcast_cont(const T&);
  template <typename T> void broadcast_cont_packed(const T &);

  // receive_broadcast
  template <typename T> void recv_broadcast(T&,int);
  template <typename T> void recv_broadcast_array(T*,int,int);
  template <typename T> void recv_broadcast_cont(T&,int);
  template <typename T> void recv_broadcast_cont_packed(T&,int);

  // scatter/gather 
  template <typename T> void scatter(const multimap<int,T>);
  template <typename T> void recv_scatter(T&,int);
  template <typename T> void gather(multimap<int,T>&);
  template <typename T> void send_gather(T&,int);

  // debug versions
  template <typename T> void gather_debug(multimap<int,T>&);
  template <typename T> void send_gather_debug(T&,int);

  // scatter/gather packed
  template <typename T> void scatter_packed(const multimap<int,T>);
  template <typename T> void recv_scatter_packed(T&,int);
  template <typename T> void gather_packed(multimap<int,T>&);
  template <typename T> void send_gather_packed(const T &,int);

  // reduce ops
  template <typename T> T sum_all(const T&);

  // syncronisation
  void barrier();
  void barrier(const string&);
};

#include "tml/comm/comm.hpp"
#include "tml/comm/comm_coll.hpp"

#endif //__COMM_H
