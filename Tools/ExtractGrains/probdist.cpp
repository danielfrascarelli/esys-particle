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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "probdist.h"

ProbDist::ProbDist (double MinSize, double MaxSize, double BaseConst, int DistType)
{
   minsize = MinSize;
   maxsize = MaxSize;
   base = BaseConst;

   nbins = (long)(log(maxsize/minsize)/log(base));
   nbins++;
   //fprintf (stderr, "Info: PDF nbins = %7d\n", nbins);

   Nevents = 0;

   disttype = DistType;

   Create ();
}

void ProbDist::Create ()
{
   Edist = new long[nbins];
   Vdist = new double[nbins];

   for (long i=0;i<nbins;i++) Edist[i] = 0;
   for (long i=0;i<nbins;i++) Vdist[i] = 0.0;
}

void ProbDist::Destroy ()
{
   delete [] Edist;
   delete [] Vdist;
}

ProbDist::~ProbDist ()
{
   Destroy ();
}

void ProbDist::AddSample (double evsize)
{
//   fprintf (stderr, "Here Too\n");
   AddEvSize (evsize);
//   fprintf (stderr, "Here Too b\n");
}

void ProbDist::AddSample (double size, double weight)
{
//   fprintf (stderr, "Here Too\n");
  AddEvSize (size,weight);
//   fprintf (stderr, "Here Too b\n");
}

void ProbDist::Write (const char *filename, double EvRate)
{
   FILE *evfile;
   char evfilename[50];
   long i;
   double sizebin, midbin, probdensity;

//   fprintf (stderr, "Info: PDF Write\n");
   sprintf(evfilename,"%s.p", filename);
   //fprintf (stderr, "Info: PDF Write %s\n", evfilename);

   evfile = fopen(evfilename, "w");
   //fprintf (stderr, "Info: PDF Open \n");

   for (i=0;i<nbins;i++) {
   //fprintf (stderr, "Info: i=%ld Edist[i]=%ld\n", i, Edist[i]);
      if (Edist[i] > 0) {
         //fprintf (stderr, "Info: i=%ld ; Edist[i] = %ld\n", i, Edist[i]);
         sizebin = minsize*powf(base, i)*(base - 1.0);
         midbin = minsize*(powf(base, i) + powf(base, i+1))/2.0;
         //fprintf (stderr, "Info: sizebin=%8.6e ; midbin = %8.6e\n", sizebin, midbin);
	 if (disttype==0) {
  	    probdensity = (double)(Edist[i])/Nevents/sizebin;
         //fprintf (stderr, "Info: probdensity=%8.6e ; Nevents = %ld\n", probdensity, Nevents);
           fprintf (evfile, "%8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %8.6e %7ld\n", midbin, probdensity, midbin*EvRate, probdensity/EvRate, log(midbin),log(probdensity),Vdist[i], Edist[i]);
	 }
	 else if (disttype==1) {
  	    probdensity = (double)(Edist[i]);
fprintf (evfile, "%8.6e %8.6e %8.6e %8.6e %8.6e %8.6e\n", midbin, probdensity, midbin*EvRate, probdensity/EvRate, log(midbin),log(probdensity));
	 }
	 else if (disttype==2) {
	   fprintf (evfile, "%8.6e %7ld %8.6e %8.6e %8.6e %8.6e \n", midbin,Edist[i],Vdist[i], log(midbin),log(double(Edist[i])),log(Vdist[i]));
	 }
      }
   }

   fclose (evfile);
}

void ProbDist::AddEvSize (double evsize)
{
   long index, i;

   if (disttype==0) {
      index = (long)(log(evsize/minsize)/log(base));
      if ((index>=0)&&(index<nbins)) {
         Edist[index]++;
	 Vdist[index]+=evsize;
         Nevents++;
      }
      //fprintf(stderr, "Info: AddSample (%8.6e) index=%ld Edist[index]=%ld\n", evsize, index, Edist[index]);
   }
   else if (disttype==1) {
      index = (long)(log(evsize/minsize)/log(base));
      for (i=0;i<=index;i++) {
         Edist[i]++;
         Nevents++;
      }
   }
}

void ProbDist::AddEvSize (double size,double weight)
{
   long index, i;

   if (disttype==0) {
      index = (long)(log(size/minsize)/log(base));
      if ((index>=0)&&(index<nbins)) {
         Edist[index]+=weight;
	 Vdist[index]+=weight;
         Nevents++;
      }
   }
   else if (disttype==1) {
      index = (long)(log(size/minsize)/log(base));
      for (i=0;i<=index;i++) {
         Edist[i]+=weight;
         Nevents++;
      }
   }
   else if (disttype==2) {
     
      index = (long)(log(size/minsize)/log(base));
      for (i=index;i<nbins;i++) {
	Edist[i]++;
	Vdist[i]+=weight;
	Nevents++;
      }
   }
}


