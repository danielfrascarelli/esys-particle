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

#include "SnapFileHelp.h"

// --- Project includes ---
#include "Foundation/vec3.h"

// --- STL includes ---
#include <iterator>
#include <fstream>
#include <set>

using std::istream_iterator;
using std::back_inserter;
using std::ifstream;
using std::make_pair;
using std::set;

int get_version(const string& infilename)
{
  string dummystring;
  int version;
  ifstream headerfile(infilename.c_str()); 
  // read token  
  headerfile >> dummystring;

  if(dummystring=="V"){ // if V -> new version 
    headerfile >> version ;
    cout << "version : " << version << endl;
  } else {
    cout << "pre- V.1 version" << endl;
    version=0;
  }
  headerfile.close();

  return version;
}

vector<string> get_filenames(const string& infilename, int version)
{
  cout << "infilename : " << infilename << endl;
  ifstream headerfile(infilename.c_str()); 
  float dummy,xmax,ymax,zmax,xmin,ymin,zmin;
  vector<string> filenames;
  string dummystring;

  if(version==0){
    headerfile >> dummy >> dummy >> dummy;
    headerfile >> dummystring >> dummy;
  } else if ((version==1) || (version==2) || (version==3)){
    headerfile >> dummystring >> dummy;
    headerfile >> dummy >> dummy >> dummy;
    headerfile >> dummystring >> dummy;
  } else {
    cerr << "unknown checkpoint version " << version << endl;
  }
  // get bounding box
  headerfile >> dummystring;
  headerfile >> xmin >> ymin >> zmin >> xmax >> ymax >> zmax ;

  // ignore periodic bdry
  headerfile >> dummystring >> dummy >> dummy >> dummy;

  // ignore dimension
  headerfile >> dummystring >> dummystring;

  // get file names
  copy(istream_iterator<string>(headerfile),istream_iterator<string>(),back_inserter(filenames));
  
  headerfile.close();

  cout << "nr. of filenames: " << filenames.size() << endl;
  return filenames;
}
