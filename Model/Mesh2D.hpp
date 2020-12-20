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

#include <utility>
using std::make_pair;

/*!
  Call a constant member function of Edge2D taking no argument and returning a 
  value for all edges and collect the return values in a container. The 
  container has to be an STL sequence container (vector,list...) or something 
  with the same interface. The template parameter P is a type of container of 
  the return type of the particle member function, not the return type itself.

  \param cont the container
  \param rdf the particle member function
*/
template <typename P> 
void Mesh2D::forAllEdgesGet(P& cont,typename P::value_type (Edge2D::*rdf)() const)
{
  for(vector<Edge2D>::iterator iter=m_edges.begin();
      iter!=m_edges.end();
      iter++){
    cont.push_back(((*iter).*rdf)());
  }
}

/*!
  \param rdf the particle member function
*/
template <typename P> 
vector<pair<int,P> > Mesh2D::forAllEdgesGetIndexed(P (Edge2D::*rdf)() const)
{
  vector<pair<int,P> > res;

  for(vector<Edge2D>::iterator iter=m_edges.begin();
      iter!=m_edges.end();
      iter++){
    res.push_back(make_pair(iter->getID(),((*iter).*rdf)()));
  }

  return res;
}
