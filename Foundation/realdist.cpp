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

#include "realdist.h"

// -- system includes --
#include <cmath>
using std::floor;

// -- IO includes --
#include <fstream>
#include <sstream>
using std::ostringstream;
using std::ofstream;


RealDist::RealDist (double MinSize, double MaxSize, int Nbins)
{
   minsize = MinSize;
   maxsize = MaxSize;
   nbins = Nbins;
   Nevents = 0;
  
   Create ();
}

void RealDist::Create ()
{
   int i; 

   binsize = (maxsize-minsize)/nbins;

   Edist = new long[nbins];

   for (i=0;i<nbins;i++) Edist[i] = 0;
}

void RealDist::Destroy ()
{
   delete [] (Edist);
}

RealDist::~RealDist ()
{
   Destroy ();
}

void RealDist::AddSample (double evsize)
{
   AddEvSize (evsize);
   Nevents++;
}

void RealDist::Write (const string& filename)
{
   ostringstream evfilename;
   int i;

   evfilename << filename << ".r";

   ofstream evfile(evfilename.str().c_str());

   for (i=0;i<nbins;i++) {
     if (Edist[i] > 0) evfile << double(i)*binsize+minsize+binsize*0.5 << " " << double(Edist[i]) << std::endl;
   }

   evfile.close();
}

void RealDist::AddEvSize (double evsize)
{
   int index;
   double Esize;

   Esize = (evsize - minsize);
   index = (int)(floor(Esize/binsize));
   
   if ((index < nbins)&&(index >= 0)) Edist[index]++;
}

void RealDist::Clear()
{
  for (int i=0;i<nbins;i++) Edist[i]=0;
  Nevents=0;
}
