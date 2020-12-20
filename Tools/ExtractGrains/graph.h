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

#ifndef __GRAPH_H
#define __GRAPH_H

// --- project includes ---
#include "Foundation/vec3.h"

// --- STL includes ---
#include <list>
#include <vector>
#include <utility>
#include <map>
#include <iostream>
#include <string>

using std::list;
using std::vector;
using std::pair;
using std::map;
using std::ostream;
using std::endl;
using std::string;

struct pdata{
  int tag;
  double mass;
  Vec3 pos;
  double rad;
  Vec3 vel;
  Vec3 angvel;
  
  pdata(int t=-1,
	double m=0.0,
	Vec3 p=Vec3(0.0,0.0,0.0),
	double r=0.0, 
	Vec3 v=Vec3(0.0,0.0,0.0),
	Vec3 av=Vec3(0.0,0.0,0.0))
    :tag(t),mass(m),pos(p),rad(r),vel(v),angvel(av){};
};

struct pdata2d{
  int gid; // grain id
  Vec3 pos2d; // projected position
  double rad; // radius
};

/*!
  \class Graph
  \brief Graph class, partially based on Sedgewick, "Alg. in C++", progs. 17.1, 17.9 and 17.10

*/
class Graph
{
 public: // types
  struct Edge{
    int i,j;
    Edge(int v=-1,int w=-1):i(v),j(w){}
  };
  
  struct Node{
    int id;
    Node(int i=-1):id(i){};
  };
  typedef list<int>::iterator adjIterator;

 private:
  map<int,list<int> > m_data;
  map<int,int> cid; // component id for conn. component
  map<int,pdata> m_vertex_data;
  map<int,double> m_grain_mass;
  map<int,Vec3> m_grain_rot;
  int ccnt;

  void ccR(int);
  void dfsIter(int);

 public: // methods
  Graph();
  ~Graph();
  int numV() const;
  int numE() const;
  int getGrainID(int) const;
  double getParticleMass(int) const;
  void insert(const Edge&);
  void insert(const pair<int,int>&);
  void setVertexData(int,const pdata&);
  pdata getVertexData(int i) const;
  void remove(const Edge&);
  bool isEdge(int,int);
  void removeDoubles();
  adjIterator IterBegin(int);
  adjIterator IterEnd(int);
  map<int,int>::const_iterator cid_begin() const {return cid.begin();};
  map<int,int>::const_iterator cid_end() const {return cid.end();};
  void makeConnComp();
  void printGrainPCount(ostream&);
  void printGrainMass(ostream&);
  void printIdList(const string&);
  void printRotList(const string&);
  void printGrainCountDist(const string&);
  void printGrainDiamDist(const string&,double,int);
  void printGrainMassDist(const string&,double,int);
  void printSieveDist(const string&,double);
  double getPercentile(double);
  void writeAvgGrainSizeProfile(const string&,double,double,int) ;
  void writeAvgGrainSizeGrid(const string&,double,double,double,double,double,double,double) ;
  void writeMatrixFractionProfile(const string&,double,double,int,double) ;
  void writeMatrixFractionGrid(const string&,double,double,double,double,double,double,double,double) ;
  void printGrainsAsVtk(const string&,double);
  void printAllAsVtk(const string&);
  void printCrossSection(const Vec3&,const Vec3&,const Vec3&,const string&,int,int,double,double,double,double,bool,bool);
  void printGrainCenterPosition(const string&);
};

#endif //__GRAPH_H
