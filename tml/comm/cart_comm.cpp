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


#include "tml/comm/cart_comm.h"
#include "Foundation/console.h"

//--- I/0 ---
#include <iostream>
#include <sstream>

/*!
  Constructor, using an STL vector for boundary conditions and optionally prescribing
  dimensions. Dimensions are choosen according to the size of the communicator
  via MPI_Create_dims.

  \param old_comm the old communicator
  \param dims the dimensions
  \param circular circular boundaries
*/
TML_CartComm::TML_CartComm(TML_Comm* old_comm, vector<unsigned int> dims, vector<bool> circular)
{
  const unsigned int ndims = 3;
  while (dims.size() < ndims)
  {
    dims.push_back((dims.size() > 0 ? 1 : 0));
  }
  int *mpi_dims=new int[ndims];
  for(unsigned int i=0;i<ndims;i++){
    mpi_dims[i]=dims[i];
  }

  // get dims
  MPI_Dims_create(old_comm->size(),ndims,mpi_dims);
  
  console.Debug() << "ndims = " << ndims << "\n";
  std::stringstream str;
  str << "dims = [";
  if (ndims > 0) {
    str << mpi_dims[0];
    for (unsigned int i = 1; i < ndims; i++) {
      str << ", " << mpi_dims[i];
    }
  }
  str << "]";
  console.Debug() << str.str() << "\n";

  for(unsigned int i=0;i<ndims;i++){ 
    m_dims.push_back(mpi_dims[i]);
  }
  // C-array boundary conditions
  int *l_circ=new int[ndims];

  // fill in
  for(unsigned int i=0;i<ndims;i++){
    l_circ[i]=circular[i];
  }

  // generate new communicator
  MPI_Cart_create(old_comm->comm(), ndims, mpi_dims, l_circ, 0, &m_comm);
  m_ndims=ndims;

  // clean up
  delete [] mpi_dims;
  delete [] l_circ;
}


/*!
  Constructor, using STL vectors for dimensions and boundary conditions 

  \param old_comm the old communicator
  \param ndims the number of dimensions
  \param dims the dimensions
  \param circular circular boundaries
*/
TML_CartComm::TML_CartComm(TML_Comm* old_comm, unsigned int ndims, const vector<int> &dims, const vector<bool> &circular)
{
  // do checks
  if((dims.size()!=ndims) || (circular.size()!=ndims)){
    std::cerr << "wrong nr. of dims in TML_CartComm constructor" << std::endl;
    exit(1);
  } 
  m_dims=dims;
  int nelem=1;
  for(std::vector<int>::const_iterator iter=dims.begin();
      iter!=dims.end();
      iter++){
    nelem*=*iter;
  }
  if(nelem!=old_comm->size()){
    std::cerr << "wrong nr. of processes TML_CartComm constructor" << std::endl;
    exit(1);
  }

  // C-arrays for dims and boundary conditions
  int *l_dims=new int[ndims];
  int *l_circ=new int[ndims];

  // fill in
  for(unsigned int i=0;i<ndims;i++){
    l_dims[i]=dims[i];
    l_circ[i]=circular[i];
  }

  // construct the new communicator
  MPI_Cart_create(old_comm->comm(),ndims,l_dims,l_circ,0,&m_comm);
  m_ndims=ndims;

  delete [] l_dims;
  delete [] l_circ;
}

/*!
  Constructor, using C arrays for dimensions and boundary conditions 

  \param old_comm the old communicator
  \param ndims the number of dimensions
  \param dims the dimensions
  \param circular circular boundaries
  \warning no checking of nr. of dims and boundary cond.
*/
TML_CartComm::TML_CartComm(TML_Comm* old_comm,unsigned int ndims,int* dims,int* circular)
{
  int nelem=1;

  // check number of processes
  for(unsigned int i=0;i<ndims;i++){
    nelem*=dims[i];
  }
  for(unsigned int i=0;i<ndims;i++){ 
    m_dims.push_back(dims[i]);
  }

  if(nelem!=old_comm->size()){
    std::cerr << "wrong nr. of processes TML_CartComm constructor" << std::endl;
    exit(1);
  }

  // generate new communicator
  MPI_Cart_create(old_comm->comm(),ndims,dims,circular,0,&m_comm);
  m_ndims=ndims;
}

/*!
  Get cartesian coordinates of a given process in the communicator

  \param rank the rank of the process
*/
vector<int> TML_CartComm::get_coords(int rank) 
{
  int *coords=new int[m_ndims];
  
  MPI_Cart_coords(m_comm,rank,m_ndims,coords);
  
  vector<int> res=vector<int>(coords,&(coords[m_ndims]));
  delete [] coords;

  return res;
}

/*!
  Get cartesian coordinates of local process
*/
vector<int> TML_CartComm::get_coords() const
{
  int *coords=new int[m_ndims];
  
  MPI_Cart_coords(m_comm,rank(),m_ndims,coords);
  
  vector<int> res=vector<int>(coords,&(coords[m_ndims]));
  delete [] coords;
  
  return res;
} 
 
/*!
  Get  size of the communicator in all directions
*/
vector<int> TML_CartComm::get_all_dims() const
{
  vector<int> res;

  for(int i=0;i<m_ndims;i++){
    res.push_back(m_dims[i]);
  }

  return res; 
}

/*!
  get size of communicator in direction i

  \param i the number of the direction 
*/
int TML_CartComm::get_dim(int i)
{
  int res;
  if((0<=i)&&(i<m_ndims)){
    res=m_dims[i];
  }else{
    res=0;
  }
  return res;
}
