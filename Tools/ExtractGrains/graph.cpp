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

#include "graph.h"
#include "Matrix3.h"
#include "Triangle2d.h"
#include <iostream>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <cmath>

using std::pow;
using std::cout;
using std::cerr;
using std::endl;
using std::sort;
using std::multimap;
using std::make_pair;
using std::ostringstream;
using std::ofstream;

// --- 
#include "probdist.h"

Graph::Graph()
{}

Graph::~Graph()
{}

/*!
  get number of vertices
*/
int Graph::numV() const
{
  return m_data.size();
}

/*!
  get number of egdes
*/
int Graph::numE() const
{
  int nedges=0;
  for(map<int,list<int> >::const_iterator iter=m_data.begin();
      iter!=m_data.end();
      iter++){
    nedges+=iter->second.size();
  }
  return nedges/2;
}

/*!
  insert a new edge

  \param E the edge
*/
void Graph::insert(const Edge& E)
{
  // check if the id's are inside
  //  if((E.i<m_data.size()) && (E.j<m_data.size())){
    m_data[E.i].push_back(E.j);
    m_data[E.j].push_back(E.i);
    //  }
}

/*!
  insert a new edge as pair<int,int>

  \param p the edge
*/
void Graph::insert(const pair<int,int>&  p)
{
  // check if the id's are inside
  //  if((p.first<m_data.size()) && (p.second<m_data.size())){
    m_data[p.first].push_back(p.second);
    m_data[p.second].push_back(p.first);
    //  }
}

void Graph::removeDoubles()
{
  for(map<int,list<int> >::iterator iter=m_data.begin();
      iter!=m_data.end();
      iter++){
    iter->second.sort();
    iter->second.unique();
  }
}

void Graph::remove(const Edge& E)
{
}

bool Graph::isEdge(int i,int j)
{
  list<int>::iterator iter=find(m_data[i].begin(),m_data[i].end(),j);
  return (iter!=m_data[i].end());
}

int Graph::getGrainID(int i) const 
{
  int res;

  map<int,int>::const_iterator it=cid.find(i);
  if(it!=cid.end()){
    res=it->second;
  } else {
    res=-1;
  }

  return res;
}

double Graph::getParticleMass(int i) const
{
  double res;

  map<int,pdata>::const_iterator it=m_vertex_data.find(i);
  if(it!=m_vertex_data.end()){
    res=it->second.mass;
  } else {
    res=0.0;
  }

  return res;
}

Graph::adjIterator Graph::IterBegin(int i)
{
  return m_data[i].begin();
}

Graph::adjIterator Graph::IterEnd(int i)
{
  return m_data[i].end();
}

void Graph::makeConnComp()
{
  ccnt=1;
  for(map<int,list<int> >::iterator iter=m_data.begin();
      iter!=m_data.end();
      iter++){
    if(cid[iter->first]==0){
      //ccR(iter->first);
      dfsIter(iter->first);
      ccnt++;
    }
  }
  cout << ccnt-1 << " components found" << endl;
}

void Graph::ccR(int i)
{
  cid[i]=ccnt;
  for(list<int>::iterator iter=m_data[i].begin();
      iter!=m_data[i].end();
      iter++){
    if(cid[*iter]==0){
      ccR(*iter);
    }
  }
}

  
/*!
  Iterative depth first search. Avoids possible issues with recursion depth
  in very large components

  \param u starting node
*/
void Graph::dfsIter(int u)
{
  list<int> to_visit;
  //  int vcnt=0; // count visited
  if(cid[u]==0){ // node is "white"
    to_visit.push_back(u);
    while(!to_visit.empty()){
      int v=to_visit.back();
      to_visit.pop_back();
      if(cid[v]==0){ // v has not been visited yet
	cid[v]=ccnt; // mark v as visited
	//	vcnt++;
	//	std::cerr << " v " << vcnt << std::endl;
	// add unvisited neigbours to "to visit" list
	for(list<int>::iterator iter=m_data[v].begin();
	    iter!=m_data[v].end();
	    iter++){
	  if(cid[*iter]==0){
	    to_visit.push_back(*iter);
	  }
	}
      }
    }
  }
}

void Graph::setVertexData(int n,const pdata& pd)
{
  m_vertex_data[n]=pd;
  cid[n]=0;
  m_data[n].sort(); // dummy to make sure the node is inserted
}

// fix undefined if not found
pdata Graph::getVertexData(int i) const
{
  pdata res;

  map<int,pdata>::const_iterator iter=m_vertex_data.find(i);
  if(iter!=m_vertex_data.end()){
    res=iter->second;
  } 
  return res;
}


void Graph::printGrainPCount(ostream& ost)
{
  vector<int> pcv(ccnt,0);
  
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    pcv[(iter->second)-1]++;
  }

  for(int i=0;i<ccnt-1;i++){
    ost << i << " " << pcv[i] << endl;
  }
}

void Graph::printGrainMass(ostream& ost)
{
  vector<double> pmv(ccnt,0);
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    int id=(iter->second)-1;
    pdata pd=m_vertex_data[iter->first];
    pmv[id]+=pd.mass;
  }

  for(int i=0;i<ccnt-1;i++){
    ost << i << " " << pmv[i] << endl;
  }
}

void Graph::printIdList(const string& filename)
{
  ofstream outfile(filename.c_str());

  // calculate grain masses
  vector<double> pmv(ccnt,0);
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    int id=(iter->second)-1;
    pdata pd=m_vertex_data[iter->first];
    pmv[id]+=pd.mass;
  }
  
  // write pid,gid,gmass
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    outfile << iter->first << " " << iter->second << " " << pmv[(iter->second)-1] << endl;
  }

  outfile.close();
}

/*!
  print grain list with per-grain rotation vectors

  \param filename the name of the file written
*/
void Graph::printRotList(const string& filename)
{
  // setup grain - particle map
  multimap<int,int> grain_particle_map;
  for(map<int,int>::const_iterator iter=cid_begin();
      iter!=cid_end();
      iter++){
    int part_idx=iter->first;
    int grain_idx=(iter->second)-1;
    grain_particle_map.insert(make_pair(grain_idx,part_idx));
  }

  // calculate center of mass for each grain
  vector<Vec3> grain_center(ccnt,Vec3(0.0,0.0,0.0));
  vector<Vec3> grain_vel(ccnt,Vec3(0.0,0.0,0.0));
  vector<double> grain_mass(ccnt,0);
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    int id=(iter->second)-1;
    pdata pd=m_vertex_data[iter->first];
    grain_mass[id]+=pd.mass;
    grain_center[id]+=pd.pos*pd.mass;
    grain_vel[id]+=pd.vel*pd.mass;
  }
  for(int i=0;i<ccnt;i++){
    if(grain_mass[i]>0){
      grain_center[i]=grain_center[i]/grain_mass[i];
      grain_vel[i]=grain_vel[i]/grain_mass[i];
    } else {
      grain_center[i]=Vec3(0.0,0.0,0.0);
      grain_vel[i]=Vec3(0.0,0.0,0.0);
    }
  }
  

  // calculate rotations
  vector<Vec3> grain_angvel(ccnt,Vec3(0.0,0.0,0.0));
  for(int i=0;i<ccnt;i++){
    multimap<int,int>::iterator grain_begin=grain_particle_map.lower_bound(i);
    multimap<int,int>::iterator grain_end=grain_particle_map.upper_bound(i);
    int npart_g=grain_particle_map.count(i); // nr. of particles in grain
    npart_g=npart_g>5 ? 5 : npart_g;
    if(npart_g>1){
      // calc axis direction
      Vec3 axis=Vec3(0.0,0.0,0.0);
      multimap<int,int>::iterator iter1=grain_begin;
      multimap<int,int>::iterator iter2=grain_begin;
      iter2++;
      int cnt=0;
      std::cout << "grain " << i << " ---- " << std::endl;
      for(int j=0;j<npart_g-1;j++){
	Vec3 v1=(m_vertex_data[iter1->second]).vel-grain_vel[j];
	for(int k=1;k<npart_g;k++){
	  Vec3 v2=(m_vertex_data[iter2->second]).vel-grain_vel[k];
	  Vec3 temp_axis=cross(v1,v2).unit();
	  // force positive z
	  if(temp_axis.Z()<0.0) temp_axis=-1.0*temp_axis;
	  axis+=temp_axis;
	  cnt++;
	  //	std::cout << temp_axis << std::endl;
	  iter2++;
	}
	iter1++;
      }
      axis=(axis/double(cnt)).unit();
      std::cout << "axis,mass " << axis << " " << grain_mass[i] << std::endl; 
      // convert per-particle velocity to grain rotation  
      for(multimap<int,int>::iterator iter=grain_begin;
	  iter!=grain_end;
	  iter++){
	// vector from axis to position
	Vec3 d=(m_vertex_data[iter->second]).pos-grain_center[iter->first];
	Vec3 r=d-axis*(axis*d);
	// rotation should be velocity x r / (r.r)
	Vec3 vcm=(m_vertex_data[iter->second]).vel-grain_vel[iter->first]; // velocity rel. to center of mass 
	Vec3 rot=cross(r,vcm)/(r*r);
	// std::cout << "rot, d_axis " << rot << " " << rot.unit()-axis << std::endl;
	grain_angvel[i]+=rot*(m_vertex_data[iter->second]).mass;
      } 
    } else {
      multimap<int,int>::iterator grain_begin=grain_particle_map.find(i);
      grain_angvel[i]=(m_vertex_data[grain_begin->second]).angvel.unit();
    }
  }

  
  // write gid,mass,center, angvel
  ofstream outfile(filename.c_str());
  for(int i=0;i<ccnt;i++){
    if(grain_mass[i]>0){
      outfile << i << " " << grain_mass[i] << " " << (grain_angvel[i]/grain_mass[i]).norm() << " "  << grain_center[i] <<  " " << grain_angvel[i]/grain_mass[i]  << endl;
    }
  }
  outfile.close();
}

void Graph::printGrainCountDist(const string& filename)
{
  // get counts and min/max 
  vector<int> pcv(ccnt-1,0);
  
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    pcv[(iter->second)-1]++;
  }
  sort(pcv.begin(),pcv.end());
  cout << "min, max counts: " << pcv.front() << " , " << pcv.back() << endl;

  // setup distribution
  ProbDist PD(double(pcv.front()),double(pcv.back()),2.5,0);
  // fill in data
  for(vector<int>::iterator iter=pcv.begin();
      iter!=pcv.end();
      iter++){
    PD.AddSample(double(*iter));
  }
  // write to file
  PD.Write(filename.c_str(),1.0);
}

/*!
 output grain mass distribution

 \param filename the output filename
 \param base the base for the exponential bin size
 \param cum cumulative or not
*/
void Graph::printGrainMassDist(const string& filename,double base,int cum)
{
  cout << "Graph::printGrainMassDist( " << filename << " " << base << " " << cum << " )" << std::endl;
  // get masses and min/max
  vector<double> pmv(ccnt-1,0);
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    int id=(iter->second)-1;
    pdata pd=m_vertex_data[iter->first];
    pmv[id]+=pd.mass;
  }

  sort(pmv.begin(),pmv.end());
  cout << "min, max mass : " << pmv.front() << " , " << pmv.back() << endl;
  
  // setup distribution
  ProbDist PD(pmv.front(),pmv.back(),base,cum);
  // fill in data
  for(vector<double>::iterator iter=pmv.begin();
      iter!=pmv.end();
      iter++){
    PD.AddSample(*iter);
  }
  
  PD.Write(filename.c_str(),1.0);
}

/*!
  print grain size distribution - cumulative mass over diameter

 \param filename the output filename
 \param base the base for the exponential bin size
 \param cum cumulative or not
*/
void Graph::printGrainDiamDist(const string& filename,double base,int cum)
{
  cout << "Graph::printGrainDiamDist( " << filename << " " << base << " " << cum << " )" << std::endl;
  // get masses and min/max
  vector<double> pmv(ccnt-1,0);
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    int id=(iter->second)-1;
    pdata pd=m_vertex_data[iter->first];
    pmv[id]+=pd.mass;
  }

  sort(pmv.begin(),pmv.end());

  // setup distribution
  double min_dia=cbrt(pmv.front());
  double max_dia=cbrt(pmv.back());
  cout << "min, max diameter : " << min_dia << " , " << max_dia << endl;
  ProbDist PD(min_dia,max_dia,base,cum);
  // fill in data
  for(vector<double>::iterator iter=pmv.begin();
      iter!=pmv.end();
      iter++){
    double grain_mass=*iter;
    double diameter=cbrt(grain_mass);
    PD.AddSample(diameter);
  }
  
  PD.Write(filename.c_str(),1.0);
}

/*!
  print grain size distribution - cumulative mass over diameter

 \param filename the output filename
 \param base the base for the exponential bin size
 \param cum cumulative or not
*/
void Graph::printSieveDist(const string& filename,double base)
{
  cout << " Sieving: " << filename << " " << base  << std::endl;
  // get masses and min/max
  vector<double> pmv(ccnt-1,0);
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    int id=(iter->second)-1;
    pdata pd=m_vertex_data[iter->first];
    pmv[id]+=pd.mass;
  }

  sort(pmv.begin(),pmv.end());

  // setup distribution
  double min_dia=cbrt(pmv.front());
  double max_dia=cbrt(pmv.back());
  cout << "min, max diameter : " << min_dia << " , " << max_dia << endl;
  ProbDist PD(min_dia,max_dia,base,2);
  // fill in data
  for(vector<double>::iterator iter=pmv.begin();
      iter!=pmv.end();
      iter++){
    double grain_mass=*iter;
    double diameter=cbrt(grain_mass);
    PD.AddSample(diameter,grain_mass);
  }
  
  PD.Write(filename.c_str(),1.0);
}

/*!
  get the sieving size of the x-th percentile

  \param p the percentile (in percent!)
*/
double Graph::getPercentile(double p)
{
  double res=0.0;

  // get masses and min/max
  double tmass=0.0; // total mass
  vector<double> pmv(ccnt-1,0);
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    int id=(iter->second)-1;
    pdata pd=m_vertex_data[iter->first];
    pmv[id]+=pd.mass;
    tmass+=pd.mass;
  }

  double pmass=(p/100.0)*tmass; // percentile mass

  sort(pmv.begin(),pmv.end());

  double smass=0.0; // sum mass
  vector<double>::iterator iter=pmv.begin();
  while((iter!=pmv.end()) && (smass<pmass)){
    smass+=*iter;
    res=cbrt(*iter);
    iter++;
  }
  
  return res;
}

/*!
  write y-profile of average grain size (mass, equiv. diameter)

  \param filename
  \param ymin
  \param ymax
  \param nbin
*/
void Graph::writeAvgGrainSizeProfile(const string& filename,double ymin, double ymax, int nbin) 
{
  cout << "writeAvgGrainSizeProfile ( " << filename << " , " << ymin << " , " << ymax << " , " << nbin << " ) " << endl;

  // get masses 
  vector<double> pmv(ccnt-1,0);
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    int id=(iter->second)-1;
    pdata pd=m_vertex_data[iter->first];
    pmv[id]+=pd.mass;
  }

  // for each particle, get posn and grain mass -> update profile
  vector<double> mass_vec=vector<double>(nbin,0.0);
  vector<double> gd_vec=vector<double>(nbin,0.0);
  vector<double> gm_vec=vector<double>(nbin,0.0);
  double binsize=(ymax-ymin)/double(nbin);
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    int id=(iter->second)-1;
    pdata pd=m_vertex_data[iter->first];
    double ypos=pd.pos.Y(); // y-position of particle
    double gmass=pmv[id]; // mass of grain the particle belongs to 
    double pmass=pd.mass; // particle mass
    // calc bin index
    int bin_idx=int(floor(ypos-ymin)/binsize);
    if((bin_idx>=0) && (bin_idx<nbin)) {
      mass_vec[bin_idx]+=pmass;
      gm_vec[bin_idx]+=pmass*gmass;
      gd_vec[bin_idx]+=pmass*pow(gmass,1.0/3.0);
   } 
  }

  // write profile
  ofstream outfile(filename.c_str());
  for(int i=0;i<nbin;i++){
    double gmv,gdv;
    if(mass_vec[i]!=0.0){
      gmv=gm_vec[i]/mass_vec[i];
      gdv=gd_vec[i]/mass_vec[i];
    } else {
      gmv=0.0;
      gdv=0.0;
    }
    outfile << ymin+(double(i)+0.5)*binsize << " " << gmv << "  "  << gdv << endl;
  }
  outfile.close();
}

/*!
  write VTK header

  \param outfile the output file
  \param nx
  \param ny
  \param nz
  \param x0
  \param dx
  \param y0
  \param dy
  \param z0
  \param dz
*/
void write_vtk_header(ofstream&  outfile, int nx, int ny, int nz, double x0,double dx, double y0,double dy, double z0,double dz)
{ 
  outfile << "# vtk DataFile Version 2.0"  << std::endl;
  outfile << "displacement data"  << std::endl;
  outfile << "ASCII" << std::endl;
  outfile << "DATASET RECTILINEAR_GRID" << std::endl;
  outfile << "DIMENSIONS " << nx+1  << " " << ny+1 << "  " << nz+1 << std::endl;
  outfile << "X_COORDINATES " << nx+1 << " float"  << std::endl;
  for(int i=0;i<nx+1;i++){
    outfile << x0+double(i)*dx << " "; 
  }
  outfile << std::endl;
  outfile << "Y_COORDINATES " << ny+1 << " float" << std::endl;
  for(int i=0;i<ny+1;i++){
    outfile << y0+double(i)*dy << " ";
  }
  outfile << std::endl;
  outfile << "Z_COORDINATES " << nz+1 << " float" << std::endl;
  for(int i=0;i<nz+1;i++){
    outfile << z0+double(i)*dz << " ";
  }

  outfile << std::endl;

  outfile << "CELL_DATA " <<  nz*ny*nx << std::endl;
  outfile << "SCALARS displacement float" << std::endl;
  outfile << "LOOKUP_TABLE default" << std::endl;
}


/*!
  write distribution of average grain size (mass, equiv. diameter) as VTK regular grid file (3D)

  \param filename the output file name
  \param xmin minimum of the x-range
  \param xmax maximum of the x-range
  \param ymin minimum of the y-range
  \param ymax maximum of the y-range
  \param zmin minimum of the z-range
  \param zmax maximum of the z-range
  \param cellsize the size of a grid cell (1 dimension, cells are (roughly) cubes)
*/
void Graph::writeAvgGrainSizeGrid(const string& filename,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,double cellsize) 
{
  //cout << "writeAvgGrainSizeGrid ( " << filename << " , " << ymin << " , " << ymax << " , " << nbin << " ) " << endl;

  // get masses 
  vector<double> pmv(ccnt-1,0);
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    int id=(iter->second)-1;
    pdata pd=m_vertex_data[iter->first];
    pmv[id]+=pd.mass;
  }

  int xcell=int(ceil((xmax-xmin)/cellsize));
  double xcelldim=(xmax-xmin)/double(xcell);
  int ycell=int(ceil((ymax-ymin)/cellsize));
  double ycelldim=(ymax-ymin)/double(ycell);
  int zcell=int(ceil((zmax-zmin)/cellsize));
  double zcelldim=(zmax-zmin)/double(zcell);
  int ncell=xcell*ycell*zcell;

  // for each particle, get posn and grain mass -> update profile
  vector<double> mass_vec=vector<double>(ncell,0.0);
  vector<double> gd_vec=vector<double>(ncell,0.0);
  vector<double> gm_vec=vector<double>(ncell,0.0);
  
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    int id=(iter->second)-1;
    pdata pd=m_vertex_data[iter->first];
    Vec3 pos=pd.pos; // position of particle
    double gmass=pmv[id]; // mass of grain the particle belongs to 
    double pmass=pd.mass; // particle mass
    // calc bin index
    int xbin_idx;
    int ybin_idx;
    int zbin_idx;
  
    xbin_idx=int(floor(pos.X()-xmin)/xcelldim); 
    ybin_idx=int(floor(pos.Y()-ymin)/ycelldim); 
    zbin_idx=int(floor(pos.Z()-zmin)/zcelldim); 
       
    // calc index
    if((xbin_idx>=0)&&(xbin_idx<xcell) && (ybin_idx>=0)&&(ybin_idx<ycell) && (zbin_idx>=0)&&(zbin_idx<zcell) ){
      int grid_idx=xbin_idx*ycell*zcell+ybin_idx*zcell+zbin_idx;
      mass_vec[grid_idx]+=pmass;
      gm_vec[grid_idx]+=pmass*gmass;
      gd_vec[grid_idx]+=pmass*pow(gmass,1.0/3.0);
    }  
  }

  // write profile
  ofstream outfile(filename.c_str());
  // VTK header
  write_vtk_header(outfile,xcell,ycell,zcell,xmin,xcelldim,ymin,ycelldim,zmin,zcelldim);

  for(int iz=0;iz<zcell;iz++){
    for(int iy=0;iy<ycell;iy++){
      for(int ix=0;ix<xcell;ix++){
	int grid_idx=ix*ycell*zcell+iy*zcell+iz;
	//double gmv;
	double gdv;
	if(mass_vec[grid_idx]!=0.0){
	  //gmv=gm_vec[grid_idx]/mass_vec[grid_idx];
	  gdv=gd_vec[grid_idx]/mass_vec[grid_idx];
	} else {
	  //gmv=0.0;
	  gdv=0.0;
	}
	outfile << gdv << " ";
      }
    }
    outfile << std::endl;
  }
  outfile << std::endl;

  outfile.close();
}

/*!
  write y-profile of average grain size (mass, equiv. diameter)

  \param filename
  \param ymin
  \param ymax
  \param nbin
  \param size_limit matrix cutoff (mass)
*/
void Graph::writeMatrixFractionProfile(const string& filename,double ymin, double ymax, int nbin, double size_limit) 
{
  cout << "writeMatrixFractionProfile ( " << filename << " , " << ymin << " , " << ymax << " , " << nbin << " , " << size_limit << " ) " << endl;

  // get masses 
  vector<double> pmv(ccnt-1,0);
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    int id=(iter->second)-1;
    pdata pd=m_vertex_data[iter->first];
    pmv[id]+=pd.mass;
  }

  // for each particle, get posn and grain mass -> update profile
  vector<double> mass_vec=vector<double>(nbin,0.0);
  vector<double> matrix_vec=vector<double>(nbin,0.0);
  double binsize=(ymax-ymin)/double(nbin);
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    int id=(iter->second)-1;
    pdata pd=m_vertex_data[iter->first];
    double ypos=pd.pos.Y(); // y-position of particle
    double gmass=pmv[id]; // mass of grain the particle belongs to 
    double pmass=pd.mass; // particle mass
    // calc bin index
    int bin_idx=int(floor(ypos-ymin)/binsize);
    if((bin_idx>=0) && (bin_idx<nbin)) {
      mass_vec[bin_idx]+=pmass;
      if(gmass < size_limit){
	matrix_vec[bin_idx]+=pmass;
      }
   } 
  }

  // write profile
  ofstream outfile(filename.c_str());
  for(int i=0;i<nbin;i++){
    double matrix_fraction;
    if(mass_vec[i]==0.0){
      matrix_fraction=-1.0;
    } else {
      matrix_fraction=matrix_vec[i]/mass_vec[i];
    }
    double bin_mid=ymin+(double(i)+0.5)*binsize; 
    outfile << bin_mid << "  " << matrix_fraction << std::endl;
  }
  outfile.close();
}

/*!
  write spatial distribution of matrix percentage (mass fraction of grains below threshold) as VTK regular grid file (3D)

  \param filename the output file name
  \param xmin minimum of the x-range
  \param xmax maximum of the x-range
  \param ymin minimum of the y-range
  \param ymax maximum of the y-range
  \param zmin minimum of the z-range
  \param zmax maximum of the z-range
  \param cellsize the size of a grid cell (1 dimension, cells are (roughly) cubes)
  \param size_limit the size limit below which a grain counts as "matrix"
*/
void Graph::writeMatrixFractionGrid(const string& filename,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,double cellsize,double size_limit) 
{
  //cout << "writeAvgGrainSizeGrid ( " << filename << " , " << ymin << " , " << ymax << " , " << nbin << " ) " << endl;

  // get masses 
  vector<double> pmv(ccnt-1,0);
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    int id=(iter->second)-1;
    pdata pd=m_vertex_data[iter->first];
    pmv[id]+=pd.mass;
  }

  int xcell=int(ceil((xmax-xmin)/cellsize));
  double xcelldim=(xmax-xmin)/double(xcell);
  int ycell=int(ceil((ymax-ymin)/cellsize));
  double ycelldim=(ymax-ymin)/double(ycell);
  int zcell=int(ceil((zmax-zmin)/cellsize));
  double zcelldim=(zmax-zmin)/double(zcell);
  int ncell=xcell*ycell*zcell;

  // for each particle, get posn and grain mass -> update profile
  vector<double> mass_vec=vector<double>(ncell,0.0);
  vector<double> matrix_vec=vector<double>(ncell,0.0);
  
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    int id=(iter->second)-1;
    pdata pd=m_vertex_data[iter->first];
    Vec3 pos=pd.pos; // position of particle
    double gmass=pmv[id]; // mass of grain the particle belongs to 
    double pmass=pd.mass; // particle mass
    // calc bin index
    int xbin_idx;
    int ybin_idx;
    int zbin_idx;
  
    xbin_idx=int(floor(pos.X()-xmin)/xcelldim); 
    ybin_idx=int(floor(pos.Y()-ymin)/ycelldim); 
    zbin_idx=int(floor(pos.Z()-zmin)/zcelldim); 
       
    // calc index
    if((xbin_idx>=0)&&(xbin_idx<xcell) && (ybin_idx>=0)&&(ybin_idx<ycell) && (zbin_idx>=0)&&(zbin_idx<zcell) ){
      int grid_idx=xbin_idx*ycell*zcell+ybin_idx*zcell+zbin_idx;
      mass_vec[grid_idx]+=pmass;
      if(gmass < size_limit){
	matrix_vec[grid_idx]+=pmass;
      }
    }  
  }

  // write profile
  ofstream outfile(filename.c_str());
  // VTK header
  write_vtk_header(outfile,xcell,ycell,zcell,xmin,xcelldim,ymin,ycelldim,zmin,zcelldim);

  for(int iz=0;iz<zcell;iz++){
    for(int iy=0;iy<ycell;iy++){
      for(int ix=0;ix<xcell;ix++){
	int grid_idx=ix*ycell*zcell+iy*zcell+iz;
	double matrix_fraction;
	if(mass_vec[grid_idx]==0.0){
	  matrix_fraction=-1.0;
	} else {
	  matrix_fraction=matrix_vec[grid_idx]/mass_vec[grid_idx];
	}
	outfile << matrix_fraction << " ";
      }
    }
    outfile << std::endl;
  }
  outfile << std::endl;

  outfile.close();
}



/*!
  write out all grains larger than a certain size as VTK-XML files

  \param basefilename base file name, files written are basefilename.N.xml
  \param min_mass minimum mass above which grains are saved
*/
void Graph::printGrainsAsVtk(const string& basefilename,double min_mass)
{
  // get masses
  vector<double> pmv(ccnt-1,0);
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    int id=(iter->second)-1;
    pdata pd=m_vertex_data[iter->first];
    pmv[id]+=pd.mass;
  }

  // setup inverse map
  multimap<int,int> inv_map;
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    inv_map.insert(make_pair(iter->second,iter->first));
  }
  
  // write out data
  int count=0;
  for(int i=0;i<ccnt-1;i++){
    if(pmv[i]>min_mass){
      // setup file name
      ostringstream filename;
      filename << basefilename << "." << count << ".xml";
      ofstream vtkfile(filename.str().c_str());
      // write the file
      // write header 
      vtkfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
      vtkfile << "<UnstructuredGrid>\n";
      vtkfile << "<Piece NumberOfPoints=\"" << inv_map.count(i+1) << "\" NumberOfCells=\"0\">\n";
      // get bounds
      pair<multimap<int,int>::iterator,multimap<int,int>::iterator> bounds=inv_map.equal_range(i+1);
      // write particle pos
      vtkfile << "<Points>\n";
      vtkfile << "<DataArray NumberOfComponents=\"3\" type=\"Float64\" format=\"ascii\">\n";
      for(multimap<int,int>::iterator iter=bounds.first;
	  iter!=bounds.second;
	  iter++){  
	pdata pd=m_vertex_data[iter->second];
	vtkfile << pd.pos<< endl;;
      }  
      vtkfile << "</DataArray>\n";
      vtkfile << "</Points>\n";
	
      // --- write particle data ---
      // radius
      vtkfile << "<PointData Scalars=\"radius\">\n";
      vtkfile << "<DataArray type=\"Float64\" Name=\"radius\" NumberOfComponents=\"1\" format=\"ascii\">\n";
      for(multimap<int,int>::iterator iter=bounds.first;
	  iter!=bounds.second;
	  iter++){  
	pdata pd=m_vertex_data[iter->second];
	vtkfile << pd.rad<< endl;;
      }  
      
      vtkfile << "</DataArray>\n";
      vtkfile << "</PointData>\n";

      // write empty cell block
      vtkfile << "<Cells>\n";
      vtkfile << "<DataArray type=\"Int32\" NumberOfComponents=\"1\" Name=\"connectivity\" format=\"ascii\">\n";
      vtkfile << "</DataArray>\n";
      vtkfile << "<DataArray type=\"UInt8\" NumberOfComponents=\"1\" Name=\"types\" format=\"ascii\">\n";
      vtkfile << "</DataArray>\n";
      vtkfile << "</Cells>\n";
      
      // write footer
      vtkfile << "</Piece>\n";
      vtkfile << "</UnstructuredGrid>\n";
      vtkfile << "</VTKFile>\n";
      
      // close file
      vtkfile.close();

      // update counter
      count++;
    }
  }
  
}

/*!
  write out all grains larger than a certain size as VTK-XML files

  \param filename ile name
*/
void Graph::printAllAsVtk(const string& filename)
{
  // get masses
  vector<double> pmv(ccnt-1,0);
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    int id=(iter->second)-1;
    pdata pd=m_vertex_data[iter->first];
    pmv[id]+=pd.mass;
  }

  // setup inverse map
  multimap<int,int> inv_map;
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    inv_map.insert(make_pair(iter->second,iter->first));
  }
  
  // write out data
  // open file
  ofstream vtkfile(filename.c_str());
  int count=0;
  // write the file
  // write header 
  vtkfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
  vtkfile << "<UnstructuredGrid>\n";
  vtkfile << "<Piece NumberOfPoints=\"" << numV() << "\" NumberOfCells=\"" << "0" << "\">\n";
  
  // write particle pos
  vtkfile << "<Points>\n";
  vtkfile << "<DataArray NumberOfComponents=\"3\" type=\"Float64\" format=\"ascii\">\n";
  for(map<int,pdata>::iterator iter=m_vertex_data.begin();
      iter!=m_vertex_data.end();
      iter++){  
    vtkfile << iter->second.pos<< endl;;
  }  
  vtkfile << "</DataArray>\n";
  vtkfile << "</Points>\n";
  
  // --- write particle data ---
  // radius
  vtkfile << "<PointData Scalars=\"radius\">\n";
  vtkfile << "<DataArray type=\"Float64\" Name=\"radius\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(map<int,pdata>::iterator iter=m_vertex_data.begin();
      iter!=m_vertex_data.end();
      iter++){  
    vtkfile << iter->second.rad<< endl;;
  }  
  
  vtkfile << "</DataArray>\n";
  
  // grain mass
  vtkfile << "<DataArray type=\"Float64\" Name=\"grainmass\" NumberOfComponents=\"1\" format=\"ascii\">\n";  
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    vtkfile << pmv[(iter->second)-1] << endl;
  }
 
  vtkfile << "</DataArray>\n";
  vtkfile << "</PointData>\n";
  
  // write empty cell block
  vtkfile << "<Cells>\n";
  vtkfile << "<DataArray type=\"Int32\" NumberOfComponents=\"1\" Name=\"connectivity\" format=\"ascii\">\n";
  vtkfile << "</DataArray>\n";
  vtkfile << "<DataArray type=\"UInt8\" NumberOfComponents=\"1\" Name=\"types\" format=\"ascii\">\n";
  vtkfile << "</DataArray>\n";
  vtkfile << "</Cells>\n";
      
  // write footer
  vtkfile << "</Piece>\n";
  vtkfile << "</UnstructuredGrid>\n";
  vtkfile << "</VTKFile>\n";
      
  // close file
  vtkfile.close();

  // update counter
  count++;
}

/*!
  generate cross section

  \param Pos the position of the plane
  \param W spanning vector in-plane
  \param V spanning vector in-plane
  \param filename the filename 
  \param xdim image width
  \param ydim image height
  \param write_ppm if true, output PPM file
  \param filter_singles if true, filter out single particle grains
*/
void Graph::printCrossSection(const Vec3& Pos,const Vec3& W,const Vec3& V,const string& filename, int xdim, int ydim,double xmin, double xmax, double ymin, double ymax,bool write_ppm,bool filter_singles)
{
  // define plane normal and local coordinate transform
  Vec3 normal=cross(W,V).unit();  
  Matrix3 invtrans=Matrix3(W,V,normal);
  invtrans.invert();

  map<int,pdata2d> pmap;
  vector<Triangle2D> trivec;
  // get particles intersecting plane
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    int id=iter->second;
    pdata pd=m_vertex_data[iter->first];
    double dist=(pd.pos-Pos)*normal;
    if(fabs(dist)<pd.rad){
      double rad=sqrt(pd.rad*pd.rad-dist*dist);
      Vec3 pos2d=invtrans*((pd.pos-Pos)-dist*normal);
      // check if z-comp is 0
      if(pos2d.Z() > 1e-5){
	std::cout << "z>0, particle pos: " << pd.pos << std::endl;
      }  
      pdata2d pd2;
      int pid=iter->first;
      pd2.gid=id+1;
      pd2.pos2d=pos2d;
      pd2.rad=rad;
      pmap[pid]=pd2;
      std::cout << "found particle " << pid << " grain id : " << id << "  -  orig:" << pd.pos << " | " << pd.rad << "  section: " << pos2d << " | " << rad << std::endl; 
    }
  }

  // generate triangles
  for(map<int,pdata2d>::iterator iter=pmap.begin(); // for all particles in cross section
      iter!=pmap.end();
      iter++){
    int id=iter->first;
    //    std::cout << "id : " << id << "  conn: " ;
    vector<int> cpv;
    // get connected particles
    for(list<int>::iterator liter=m_data[id].begin();
	liter!=m_data[id].end();
	liter++){
      int id2=*liter;
      // check if connected particle is in cross section
      map<int,pdata2d>::iterator citer=pmap.find(id2);
      if((citer!=pmap.end()) && (id2>id)){
	cpv.push_back(id2);
	//	std:: cout << id2 << " ";
      }
    }
    //    std::cout << std::endl;
    if(cpv.size()>2){
      for(unsigned int i=0;i<cpv.size()-1;i++){
	for(unsigned int j=i+1;j<cpv.size();j++){
	  //	  std::cout << "trying : " << id << " - " << cpv[i] << " - " << cpv[j] << std::endl;    
	  if(isEdge(cpv[i] ,cpv[j] )){
	    //	    std::cout << "Triangle found" << std::endl;
	    int tid=iter->second.gid;
	    Vec3 P0=iter->second.pos2d;
	    Vec3 P1=pmap[cpv[i]].pos2d;
	    Vec3 P2=pmap[cpv[j]].pos2d;
	    trivec.push_back(Triangle2D(P0,P1,P2,tid));
	  }
	}
      }
    }
   
  }

  // open file
  ofstream outfile(filename.c_str());
  if(write_ppm){
    outfile << "P3" << " " << xdim << " " << ydim << " 255" << endl;
  } 
  // write image data
  int npc=int(ceil(cbrt(double(ccnt))));
  int cf=127/(npc+1);
  std::cout << "npc : " << npc << " cf: " << cf << endl;
  for (int iy=0;iy<ydim;iy++){
    if(iy%100==0) cout << "iy=" << iy << endl;
    for(int ix=0;ix<xdim;ix++){
      // calc position of pixel in cross-section plane
      double px=xmin+double(ix)*((xmax-xmin)/double(xdim));
      double py=ymin+double(iy)*((ymax-ymin)/double(ydim));
      //      std::cout << "pos: " << px << " , " << py << endl;
      int id=0;
      // check if point is inside particle
      bool inside=false;
      for(map<int,pdata2d>::iterator iter=pmap.begin();
	  iter!=pmap.end();
	  iter++){
	double dx=iter->second.pos2d.X()-px;
	double dy=iter->second.pos2d.Y()-py;
	double dist=sqrt(dx*dx+dy*dy);
	if(dist<iter->second.rad){
	  id=iter->second.gid;
	  inside=true;
	}
      }
      // check if inside triangle
      // for(vector<Triangle2D>::iterator iter=trivec.begin();
// 	  iter!=trivec.end();
// 	  iter++){
// 	if(iter->isIn(Vec3(px,py,0.0))){
// 	     id=iter->Id();
// 	     //	     std::cout << " in Triangle" << std::endl;
// 	}
//       }
      // check if between 3 particles
      if(!inside){
	map<double,pdata2d> distmap;
	for(map<int,pdata2d>::iterator iter=pmap.begin();
	    iter!=pmap.end();
	    iter++){
	  double dx=iter->second.pos2d.X()-px;
	  double dy=iter->second.pos2d.Y()-py;
	  double dist=sqrt(dx*dx+dy*dy)-iter->second.rad;
	  if(distmap.size()<3) {
	    distmap[dist]=iter->second;
	  } else {
	    double dist3=distmap.rbegin()->first;
	    if (dist<dist3){		     
	      distmap[dist]=iter->second;
	      if(distmap.size()>3)  distmap.erase(dist3);
	    }
	  }
	}
	if(distmap.size()==3){
	  map<double,pdata2d>::iterator iter=distmap.begin();
	  pdata2d pd1=iter->second;
	  iter++;
	  pdata2d pd2=iter->second;
	  iter++;
	  pdata2d pd3=iter->second;
	  if((pd1.gid==pd2.gid) && (pd1.gid==pd3.gid) && (iter->first<0.5)){
	    Triangle2D Tr(pd1.pos2d,pd2.pos2d,pd3.pos2d,0);
	    if(Tr.isIn(Vec3(px,py,0.0))){
	      id=pd1.gid;
	    }
	  }
	}
      }
      if(write_ppm){
	if(id>0){
	  // convert id to color
	  int red=128+(id/(npc*npc))*cf;
	  int green=128+((id/npc)%npc)*cf;
	  int blue=128+(id%npc)*cf;
	  //std::cout << "id, rgb: " << id << "[ " << red << " , " << green << " , " << blue << " ]" << endl;
	  outfile << red << " " << green << " " << blue << endl;
	} else {
	  outfile << "0 0 0" << endl;
	}
      } else {
	outfile << id << " " ;
      }    
    }
    if(!write_ppm) outfile << endl;
  }

  // close file
  outfile.close();
}

/*!
  Write a list of grain id / grain center position (x,y,z) pairs to a file

  \param filename the name of the output file
*/
void Graph::printGrainCenterPosition(const string& filename)
{
  // calculate center of mass for each grain
  vector<Vec3> grain_center(ccnt,Vec3(0.0,0.0,0.0));
  vector<double> grain_mass(ccnt,0);
  for(map<int,int>::iterator iter=cid.begin();
      iter!=cid.end();
      iter++){
    int id=(iter->second)-1;
    pdata pd=m_vertex_data[iter->first];
    grain_mass[id]+=pd.mass;
    grain_center[id]+=pd.pos*pd.mass;
  }
  for(int i=0;i<ccnt;i++){
    if(grain_mass[i]>0){
      grain_center[i]=grain_center[i]/grain_mass[i];
    } else {
      grain_center[i]=Vec3(0.0,0.0,0.0);
    }
  }
  
  // write gid,center
  ofstream outfile(filename.c_str());
  for(int i=0;i<ccnt;i++){
    if(grain_mass[i]>0){
      outfile << i << " " << grain_center[i]  << endl;
    }
  }
  outfile.close();
}
