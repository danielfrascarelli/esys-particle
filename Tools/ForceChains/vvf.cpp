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

#include "vvf.h"

//--- IO includes ---
#include <fstream>
#include <iostream>

using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;
using std::flush;

//--- project includes ---
#include "vec3.h"

/*!
  Convert a file containing interaction forces in RAW2 format into a .vvf file (PFC dump)

  \param ifname name of the input file
  \param ofname name of the output file
*/
void convert_to_vvf(const string& ifname,const string& ofname)
{
  Vec3 ppos1; // particle 1 position
  Vec3 ppos2; // particle 2 position
  Vec3 ipos; // interaction position
  Vec3 force; // interaction force
  double r1,r2; // particle radii;
  int cnt=0; // counter

  // open input file
  ifstream infile(ifname.c_str());
  // open output file
  ofstream outfile(ofname.c_str());

  // write vvf header
  outfile << "1025,1000000" << endl;
  outfile << "I7,S11,S9,F19.11,F19.11,F19.11,F19.11,F19.11,F19.11,F19.11,I3" << endl << flush;
  outfile << "EVENT_LABEL,EVENT_DATE,EVENT_TIME,PFC_CFORCE_LOCX,PFC_CFORCE_LOCY,PFC_CFORCE_LOCZ,PFC_CFORCE_LENGTH,PFC_CFORCE_FX,PFC_CFORCE_FY,PFC_CFORCE_FZ,PFC_CFORCE_TENSION" << endl;
  while(!infile.eof()){// until end of input file
    // read line
    infile >> ppos1 >> r1 >> ppos2 >> r2 >> ipos >> force;
    // if force non-zero -> write line
    if(force.norm()>0.0) {
      cnt++;
      // event label (count)
      outfile.flags(ios_base::fixed|ios_base::right);
      outfile.width(7);
      outfile << cnt;
      // dummy event date
      outfile.flags(ios_base::fixed|ios_base::right);
      outfile.width(11);
      outfile << "01-04-2000";
      // dummy event time
      outfile.flags(ios_base::fixed|ios_base::right);
      outfile.width(9);
      outfile << "01:23:45";
      outfile.flags(ios_base::scientific|ios_base::right);
      outfile.precision(11);
      outfile.width(19);
      outfile << ipos.X() ;
      outfile.width(19);
      outfile << ipos.Y(); 
      outfile.width(19);
      outfile << ipos.Z();
      outfile.width(19);
      outfile << r1+r2;
      outfile.width(19);
      outfile << force.X();
      outfile.width(19);
      outfile << force.Y();
      outfile.width(19);
      outfile << force.Z();
      // dummy "tension"
      outfile.flags(ios_base::fixed|ios_base::right);
      outfile.width(3);
      outfile << "1" << endl;
    }
  }
  cout << cnt << " lines written" << endl; 
  // close files
  infile.close();
  outfile.close();
}
