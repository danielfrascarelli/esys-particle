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

#include "Foundation/console.h"

#include <utility>
using std::make_pair;

/*!
  Call a constant member function of Triangle taking no argument and returning a 
  value for all Triangles and collect the return values in a container. The 
  container has to be an STL sequence container (vector,list...) or something 
  with the same interface. The template parameter P is a type of container of 
  the return type of the particle member function, not the return type itself.

  \param cont the container
  \param rdf the particle member function
*/
template <typename P> 
void TriMesh::forAllTrianglesGet(P& cont,typename P::value_type (Triangle::*rdf)() const)
{
  for(vector<Triangle>::iterator iter=m_triangles.begin();
      iter!=m_triangles.end();
      iter++){
    cont.push_back(((*iter).*rdf)());
  }
}

/*!
  \param rdf the particle member function
*/
template <typename P> 
vector<pair<int,P> > TriMesh::forAllTrianglesGetIndexed(P (Triangle::*rdf)() const)
{
  console.XDebug() << "TriMesh::forAllTrianglesGetIndexed\n";  
  vector<pair<int,P> > res;

  for(vector<Triangle>::iterator iter=m_triangles.begin();
      iter!=m_triangles.end();
      iter++){
	// console.XDebug() << "ID: " <<  iter->getID() << " data: " <<  ((*iter).*rdf)() << "\n";
    res.push_back(make_pair(iter->getID(),((*iter).*rdf)()));
  }

  return res;
}
