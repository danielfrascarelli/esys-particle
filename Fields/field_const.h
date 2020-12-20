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

#ifndef __FIELD_CONST_H
#define __FIELD_CONST_H

enum {
    WRITE_TYPE_DX = 0,
    WRITE_TYPE_POV,
    WRITE_TYPE_SUM,
    WRITE_TYPE_MAX,
    WRITE_TYPE_RAW_SERIES,
    WRITE_TYPE_RAW2,
    WRITE_TYPE_RAW,
    WRITE_TYPE_RAW_WITH_ID,
    WRITE_TYPE_RAW_WITH_POS_ID,
    WRITE_TYPE_SILO,
    /****fluid contents: begin****/
    WRITE_TYPE_RAW_WITH_POS,
    WRITE_TYPE_RAW_WITH_ID_POS,
    WRITE_TYPE_RAW_WITH_PARTICLE,
    WRITE_TYPE_VTI,
    WRITE_TYPE_VTU
    /****fluid contents: end****/
};

enum {
    COLL_TYPE_FULL = 1,
    COLL_TYPE_SUM,
    COLL_TYPE_MAX,
    COLL_TYPE_FULL2 = 5,
    COLL_TYPE_FULL_DX,
    COLL_TYPE_FULL_WITH_ID,
    COLL_TYPE_FULL_WITH_POS_ID,
    /****fluid contents: begin****/
    COLL_TYPE_FULL_WITH_POS,
    COLL_TYPE_FULL_WITH_ID_POS,
    COLL_TYPE_FULL_WITH_PARTICLE
    /****fluid contents: end****/
};

#endif // __FIELD_CONST_H
