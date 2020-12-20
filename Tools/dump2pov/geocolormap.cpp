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

#include "geocolormap.h"

//-- system includes --
#include <cstdlib>
using std::rand;

GeoColorMap::GeoColorMap(const Vec3& c0,const Vec3& cm,double x0,double xm,int nlayer,double rd)
  :ColorMap(c0,cm,x0,xm)
{
  // calc initial layer thickness
  double sum_l=0.0;
  for(int i=0;i<nlayer;i++){
    double l=1.0-rd+(2.0*rd*((double)(rand())/(double)(RAND_MAX)));
    m_bdry.push_back(l);
    sum_l+=l;
    std::cout << "layer: " << l << std::endl;
  }
  std::cout << "layer sum: " << sum_l << std::endl;
  // scale layer thickness
  double scale=(x_max-x_min)/sum_l;
  for(vector<double>::iterator iter=m_bdry.begin();
      iter!=m_bdry.end();
      iter++){
    (*iter)=(*iter)*scale;
    std::cout << "new layer: " << *iter << std::endl;
  }
}

Vec3 GeoColorMap::getColor(double x) const
{
  Vec3 res; 
  if(x<x_min){
    res=c_min;
  } else if (x>x_max){
    res=c_max;
  } else {
    double l=x_min;
    int count=0;
    vector<double>::const_iterator iter=m_bdry.begin();
    while((iter!=m_bdry.end()) && (l<x)){
      l+=(*iter);
      count++;
      iter++;
    }
    if((count%2)==0){
      res=c_min;
    } else {
      res=c_max;
    }
  }
  

  return res;
} 
