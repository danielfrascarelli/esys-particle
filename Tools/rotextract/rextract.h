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

#ifndef __REXTRACT_H
#define __REXTRACT_H


//--- STL includes ---
#include <map>
#include <set>
#include <string>
#include <vector>

using std::set;
using std::vector;
using std::map;
using std::string;

/*!
  \class Rextract
  \brief class for the extraction of RMS and Stddev of angvel from snapshots
*/
class Rextract
{
 private:
  vector<double> m_angvel_rms;
  vector<double> m_rad;
  set<int> m_tag_set;
  string m_infilebase;
  string m_outfilename;
  int m_tag;
  int m_count;
  int m_sumpart;
  bool m_initialized;

 public:
  Rextract(const string&,const string&,int);
  void read_frame(int);
  void write_data();
  void write_data_bin(double,double,int);
};

#endif // __REXTRACT_H
