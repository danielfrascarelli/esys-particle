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

////////////////////////////////////////////////////////////////////////////////

#include "pis/pi_storage.h"

/*!
  generate new scalar history field saver from the PIS. Dummy implementation: 
	- only generates a error message
	- returns a NULL pointer
  Proper implementation only in bonded PIS.

  \param comm 
  \param fieldname
  \param is_tagged
  \param tag
  \param mask
*/
AFieldSlave* AParallelInteractionStorage::generateNewScalarHistoryFieldSlave(TML_Comm* comm,const string& fieldname,int is_tagged,int tag,int mask)
{
	std::cerr << "generateNewScalarFieldSlave not implemented in this PIS" << std::endl;
	
	return NULL;
}