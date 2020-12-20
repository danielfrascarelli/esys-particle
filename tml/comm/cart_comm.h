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

#ifndef __CARTCOMM_H
#define __CARTCOMM_H

//--- MPI ---
#include <mpi.h>

//--- project includes ---
#include "tml/comm/comm.h"

//--- STL ---
#include <vector>
using std::vector;
/*!
  \class TML_CartComm
  \brief class for a cartesian communicator

  \author Steffen Abe
  $Revision$
  $Date$
*/
class TML_CartComm : public TML_Comm
{
private:
  int m_ndims;
  vector<int> m_dims;
 
 public:
  // constructor
  TML_CartComm(TML_Comm *comm, vector<unsigned int> dims, vector<bool> circular);
  TML_CartComm(TML_Comm *comm, unsigned int ndims, const vector<int> &dims, const vector<bool> &circular);
  TML_CartComm(TML_Comm*,unsigned int,int*,int*);

  vector<int> get_coords(int); //!< get coords of a process
  vector<int> get_coords() const; //!< get own coords
  vector<int> get_all_dims() const; //!< get size of communicator in all directions 
  int get_dim(int); //!< get size of communicator in one direction 
  int get_ndim() const {return m_ndims;}; 

  //! shift ops
  template <typename T,typename P> void shift(T,P&,int,int,int=0);
  template <typename T,typename P> void shift_array(T*,int,P*,int,int,int,int=0);
  template <typename T,typename P> void shift_cont(T,P&,int,int,int=0);

  //! packed shift ops
  template <typename T,typename P> void shift_packed(T,P&,int,int,int=0);
  template <typename T,typename P> void shift_array_packed(T*,int,P*,int,int,int,int=0);
  template <typename T,typename P> void shift_cont_packed(T,P&,int,int,int=0);
  // shift_replace ops?
};

#include "tml/comm/cart_comm.hpp"

#endif // __CARTCOMM_H
