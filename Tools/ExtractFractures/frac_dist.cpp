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

#include "frac_dist.h"
//--- IO includes ---
#include <fstream>
#include <iostream>

using std::ofstream;
using std::endl;

/*!
  \param data the fracture data (pos,dir,size)
  \param dmin distribution minimum
  \param dmax distribution maximum
  \param nbin nr. of bins for distribution
*/
FracDist::FracDist(const vector<FracFrame::fdata> data ,double dmin,double dmax,int nbin)
{
  m_dist=vector<int>(nbin,0);
  m_dmin=dmin;
  m_dmax=dmax;
  m_nbin=nbin;
  m_ntotal=0;
  m_dx=(m_dmax-m_dmin)/double(m_nbin);
  for(vector<FracFrame::fdata>::const_iterator iter1=data.begin();
      iter1!=data.end();
      iter1++){
    Vec3 pos1=iter1->pos;
    for(vector<FracFrame::fdata>::const_iterator iter2=iter1+1;
	iter2!=data.end();
	iter2++){
      Vec3 pos2=iter2->pos;
      double dist=(pos2-pos1).norm();
      int idx=int(floor((dist-m_dmin)/m_dx));
      m_ntotal++;
      if((idx>=0) && (idx<m_nbin)){
	m_dist[idx]++;
      }
    }
  }
}

/*!
  \param filename
*/
void FracDist::write(const string& filename)
{
  std::cout << "write distribution (" <<   filename << " , " << m_dmin << " , " << m_dmax << " , " << m_nbin << std::endl;
  ofstream outfile(filename.c_str());
  for(int i=0;i<m_nbin;i++){
    double bin_mid=m_dmin+(double(i)+0.5)*m_dx;
    outfile << bin_mid << " " << double(m_dist[i])/double(m_ntotal) << endl;
  }
  outfile.close();
}
