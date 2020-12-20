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

#ifndef PROBDIST_H

#define PROBDIST_H

#include <stdio.h>
#include <stdlib.h>

class ProbDist {
   public:
      ProbDist (double MinSize, double MaxSize, double BaseConst, int DistType);
      ~ProbDist ();
      void AddSample (double evsize);
      void AddSample (double evsize,double weight);
      void Write (const char *filename, double EvRate);
   private:
      long nbins; 
      int disttype;
      double maxsize, minsize, binsize, base;
      long Nevents;
      void Create ();
      void Destroy ();
      long *Edist;
      double *Vdist;
      void AddEvSize (double evsize);
      void AddEvSize (double evsize,double weight);
};

#endif
