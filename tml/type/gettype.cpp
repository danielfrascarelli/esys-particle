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


#include "tml/type/gettype.h"

// get MPI type for buildin types
// should eventually cover all MPI data types 
template <> 
MPI_Datatype SGetType::operator()<signed char>(const signed char&)
{
  return MPI_CHAR;
}

template <> 
MPI_Datatype SGetType::operator()<char>(const char&)
{
  return MPI_CHAR;
}

template <> 
MPI_Datatype SGetType::operator()<int>(const int&)
{
  return MPI_INT;
}

template <> 
MPI_Datatype SGetType::operator()<float>(const float&)
{
  return MPI_FLOAT;
}

template <> 
MPI_Datatype SGetType::operator()<double>(const double&)
{
  return MPI_DOUBLE;
}
