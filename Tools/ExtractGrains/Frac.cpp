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

#include "Frac.h"

#include <utility>
#include <sstream>
#include <fstream>

using std::make_pair;
using std::cout;
using std::endl;
using std::cerr;
using std::ostringstream;
using std::ofstream;

Frac::Frac(const Graph& g_old,const Graph& g_new)
{
  // setup old grain map from old graph
  for(map<int,int>::const_iterator iter=g_old.cid_begin();
      iter!=g_old.cid_end();
      iter++){
    int part_idx=iter->first;
    int grain_idx=iter->second;
    m_old_grain_map.insert(make_pair(grain_idx,part_idx));
//     cout << "particle " << part_idx << " in old grain " << grain_idx << endl;
  }
  // get inverse map for new graph
  for(map<int,int>::const_iterator iter=g_new.cid_begin();
      iter!=g_new.cid_end();
      iter++){
    int part_idx=iter->first;
    int grain_idx=iter->second;
    m_new_grain_map.insert(make_pair(part_idx,grain_idx));
    m_new_grain_to_particle_map.insert(make_pair(grain_idx,part_idx));
//     cout << "particle " << part_idx << " in new grain " << grain_idx << endl;
  }
  // setup old->new map
  for(multimap<int,int>::const_iterator iter=m_old_grain_map.begin();
      iter!=m_old_grain_map.end();
      iter++){
    int new_grain_id=m_new_grain_map[iter->second];
    m_grain_to_grain_map[iter->first].insert(new_grain_id);
  }
  // setup grain mass maps
  // old 
  for(map<int,int>::const_iterator iter=g_old.cid_begin();
      iter!=g_old.cid_end();
      iter++){
    int part_idx=iter->first;
    int grain_idx=iter->second;
    m_old_grain_mass_map[grain_idx]+=g_old.getParticleMass(part_idx);
  }
  // new 
  for(map<int,int>::const_iterator iter=g_new.cid_begin();
      iter!=g_new.cid_end();
      iter++){
    int part_idx=iter->first;
    int grain_idx=iter->second;
    m_new_grain_mass_map[grain_idx]+=g_new.getParticleMass(part_idx);
  }
  // setup mapping from old grain idx to particle tag (assuming all particles of a grain have the same tag)
  for(map<int,int>::const_iterator iter=g_old.cid_begin();
      iter!=g_old.cid_end();
      iter++){
    int part_idx=iter->first;
    int grain_idx=iter->second;
    int tag=g_old.getVertexData(part_idx).tag;
    m_tag_map.insert(make_pair(grain_idx,tag));
//     cout << "particle " << part_idx << " in old grain " << grain_idx << endl;
  }
}

void Frac::writeMassRatio(ostream& ost,double minmass)
{
  for(map<int,set<int> >::iterator m_iter=m_grain_to_grain_map.begin();
      m_iter!=m_grain_to_grain_map.end();
      m_iter++){
    list<double> mass_list;
    if(m_iter->second.size()>1){
      for(set<int>::iterator s_iter=m_iter->second.begin();
	  s_iter!=m_iter->second.end();
	  s_iter++){
	mass_list.push_back(m_new_grain_mass_map[*s_iter]);
// 	cout << "old [ mass ], new [ mass ] : " << m_iter->first << " [ " << m_old_grain_mass_map[m_iter->first] << 
// 	  " ] , " << *s_iter << " [ " << m_new_grain_mass_map[*s_iter] << " ] " << endl;
      }
      mass_list.sort();
      double ratio=mass_list.back()/m_old_grain_mass_map[m_iter->first];
      if(ratio>1){
	cerr << "error: ratio > 1" << endl;
      } else if (m_old_grain_mass_map[m_iter->first]>minmass) {
	ost << ratio << " " << m_old_grain_mass_map[m_iter->first] << endl;
      } 
    }
  }
}

/*!
  write out all fragment masses
*/
int Frac::writeAllMass(ostream& ost,double minmass,int count,bool with_tag)
{
  int local_count=0;
  for(map<int,set<int> >::iterator m_iter=m_grain_to_grain_map.begin();
      m_iter!=m_grain_to_grain_map.end();
      m_iter++){
    list<double> mass_list;
    if(m_iter->second.size()>1){
      for(set<int>::iterator s_iter=m_iter->second.begin();
	  s_iter!=m_iter->second.end();
	  s_iter++){
	mass_list.push_back(m_new_grain_mass_map[*s_iter]);
	cout << "old [ mass ], new [ mass ] : " << m_iter->first << " [ " << m_old_grain_mass_map[m_iter->first] << 
	  " ] , " << *s_iter << " [ " << m_new_grain_mass_map[*s_iter] << " ] " << endl;
      }
      mass_list.sort();
      double ratio=mass_list.back()/m_old_grain_mass_map[m_iter->first];
      if(ratio>1){
	cerr << "error: ratio > 1" << endl;
      } else if (m_old_grain_mass_map[m_iter->first]>minmass){
	// ost << m_iter->first << " " << m_old_grain_mass_map[m_iter->first] << " ";
	for(list<double>::reverse_iterator iter=mass_list.rbegin();
	    iter!=mass_list.rend();
	    iter++){
	  ost << *iter/m_old_grain_mass_map[m_iter->first] << " " << count+local_count;
	  if(with_tag) {
	    ost << " " << m_tag_map[m_iter->first] << " ";
	  }
	  ost << endl;
	}
	local_count++;
      } 
    }
  }
  return local_count;
}


/*!
  output fractured particles as VTK-XML files. The particles are "exploded", i.e. the
  fragments are moved out from the center for better visibility.

  \param basefilename the base filename, i.e. the files will be called basfilename.0.xml, basefilename.1.xml...
  \param exp the relative amount by which the new fragments are moved
  \param g_new
  \param fakevec
*/ 
void Frac::writeAsVtk(const string& basefilename, float exp, const Graph& g_new, bool fakevec)
{
  int old_grain_count=0;

  // get vectors for "exploding" grains
  map<int,Vec3> exp_vec=get_move_vectors(g_new);
  // get fragment masses
  map<int,double> mass_map=get_grain_mass(g_new);

  // run trough grain to grain map
  for(map<int,set<int> >::iterator m_iter=m_grain_to_grain_map.begin();
      m_iter!=m_grain_to_grain_map.end();
      m_iter++){
    // only use grains which have actually broken
    if(m_iter->second.size()>1){
      // for each "old" grain, generate a new filename 
      ostringstream filename;
      filename << basefilename << "." << old_grain_count << ".xml";
      cout << "filename: " << filename.str() << endl;
      // open file
      ofstream vtkfile(filename.str().c_str());
      // --  write out the grain (using new grain ID as a additional field)
      // for each particle in the old grain
      multimap<int,int>::iterator grain_begin=m_old_grain_map.lower_bound(m_iter->first);
      multimap<int,int>::iterator grain_end=m_old_grain_map.upper_bound(m_iter->first);
      int np=m_old_grain_map.count(m_iter->first);
      cout << "--" << m_iter->first << "--[" << np << "]--" << endl;
      // write file header 
      vtkfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
      vtkfile << "<UnstructuredGrid>\n";
      vtkfile << "<Piece NumberOfPoints=\"" << np << "\" NumberOfCells=\"0\">\n";

      // write particle pos
      vtkfile << "<Points>\n";
      vtkfile << "<DataArray NumberOfComponents=\"3\" type=\"Float64\" format=\"ascii\">\n";
      for(multimap<int,int>::iterator iter=grain_begin;
	  iter!=grain_end;
	  iter++){
	// get the data
	pdata pd=g_new.getVertexData(iter->second);
	int ngid=m_new_grain_map[iter->second];
	vtkfile << pd.pos+exp_vec[ngid]*exp << endl;
      }
      vtkfile << "</DataArray>\n";
      vtkfile << "</Points>\n";
	
      // --- write particle data ---
      // radius
      vtkfile << "<PointData Scalars=\"radius\">\n";
      vtkfile << "<DataArray type=\"Float64\" Name=\"radius\" NumberOfComponents=\"1\" format=\"ascii\">\n";

      for(multimap<int,int>::iterator iter=grain_begin;
	  iter!=grain_end;
	  iter++){
	// get the data
	pdata pd=g_new.getVertexData(iter->second);
	vtkfile << pd.rad << endl;
      }
      vtkfile << "</DataArray>\n";
	
      // new grain id
      // setup per-grain lookup table
      map<int,int> look_up;
      int idcount=0;
      for(multimap<int,int>::iterator iter=grain_begin;
	    iter!=grain_end;
	    iter++){
	int ngid=m_new_grain_map[iter->second];
	if(look_up.find(ngid)==look_up.end()){ // ngid not yet in table
	  look_up.insert(make_pair(ngid,idcount));
	  idcount++;
	}
      }
      // write data
      if(fakevec){
	vtkfile << "<DataArray type=\"Int32\" Name=\"new grain id\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	for(multimap<int,int>::iterator iter=grain_begin;
	    iter!=grain_end;
	    iter++){
	  vtkfile << look_up[m_new_grain_map[iter->second]] << " 0 0" << endl;
	}
      } else {
	vtkfile << "<DataArray type=\"Int32\" Name=\"new grain id\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	for(multimap<int,int>::iterator iter=grain_begin;
	    iter!=grain_end;
	    iter++){
	  vtkfile << look_up[m_new_grain_map[iter->second]] << endl;
	}
      } 
      // clear lookup table
      look_up.erase(look_up.begin(),look_up.end());
      vtkfile << "</DataArray>\n";
      // fragment mass
      if(fakevec){
	vtkfile << "<DataArray type=\"Float32\" Name=\"fragment mass\" NumberOfComponents=\"3\" format=\"ascii\">\n";
	for(multimap<int,int>::iterator iter=grain_begin;
	    iter!=grain_end;
	    iter++){
	  vtkfile << mass_map[m_new_grain_map[iter->second]] << " 0 0" << endl;
	}
      } else {
	vtkfile << "<DataArray type=\"Float32\" Name=\"fragment\" NumberOfComponents=\"1\" format=\"ascii\">\n";
	for(multimap<int,int>::iterator iter=grain_begin;
	    iter!=grain_end;
	    iter++){
	  vtkfile << mass_map[m_new_grain_map[iter->second]] << endl;
	}
      } 
      vtkfile << "</DataArray>\n";
      vtkfile << "</PointData>\n";
	
      // write empty cell block (no bonds)
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
            
      // close the file
      vtkfile.close();
      old_grain_count++;
    }
  }
}

/*!
  calculate the vectors for "exploding" the grain -> get centers of the new fragment

  \param g the graph from which the positions are taken
*/
map<int,Vec3> Frac::get_move_vectors(const Graph& g)
{
  map<int,Vec3> res;
  map<int,double> mass_map;
  Vec3 center;
  double total_mass=0.0;

  for(multimap<int,int>::iterator iter=m_old_grain_map.begin();
      iter!=m_old_grain_map.end();
      iter++){
    int ngid=m_new_grain_map[iter->second];
    Vec3 pos=g.getVertexData(iter->second).pos;
    double mass=g.getVertexData(iter->second).mass;
    
    res[ngid]+=mass*pos;
    mass_map[ngid]+=mass;
  }
  // calculate center*mass
  for(map<int,Vec3>::iterator iter=res.begin();
      iter!=res.end();
      iter++){
    center+=iter->second;
  }
  // calc total mass
  for(map<int,double>::iterator iter=mass_map.begin();
      iter!=mass_map.end();
      iter++){
    total_mass+=iter->second;
  }
  center=center/total_mass;

  // calc fragment centers
  for(map<int,Vec3>::iterator iter=res.begin();
      iter!=res.end();
      iter++){
    iter->second=(iter->second/mass_map[iter->first])-center;
  }

  return res;
}

/*!
  get mass of new grain fragments

  \param g the graph from which the masses are taken
*/
map<int,double> Frac::get_grain_mass(const Graph& g)
{
  map<int,double> mass_map;

  for(multimap<int,int>::iterator iter=m_old_grain_map.begin();
      iter!=m_old_grain_map.end();
      iter++){
    int ngid=m_new_grain_map[iter->second];
    double mass=g.getVertexData(iter->second).mass;
    
    mass_map[ngid]+=mass;
  }

  return mass_map;
}

/*
  get the amount of new surface generated, split between 
  grain splitting and abrasion

  \param orig the original set of bonds
  \param ratio the size ratio threshold above which a fracture event counts as splitting
  \param minmass minimum mass of the original grain
*/ 
pair<double,double> Frac::getSplitAbrasion(Graph &orig, double ratio, double minmass)
{
  double split_area=0.0;
  double ab_area=0.0;

  // iterate through grain-to-grain map
  for(map<int,set<int> >::iterator m_iter=m_grain_to_grain_map.begin();
      m_iter!=m_grain_to_grain_map.end();
      m_iter++){
    // only use grains which have actually broken and are larger than minmass
    if((m_iter->second.size()>1)&& (m_old_grain_mass_map[m_iter->first]>minmass)){
      // std::cout << "old grain id: " << m_iter->first << std::endl;
      // get largest new grain
      double maxmass=0.0;
      for(set<int>::iterator s_iter1=m_iter->second.begin();
	  s_iter1!=m_iter->second.end();
	  s_iter1++){
	double nmass=m_new_grain_mass_map[*s_iter1];
	maxmass=(nmass>maxmass) ? nmass : maxmass; 
	//std::cout << "new grain id, mass: " << *s_iter1 << " , " << nmass << std::endl; 
      }
      // mass ratio
      double mrat=maxmass/m_old_grain_mass_map[m_iter->first];
      //      std::cout << "old grain mass, ratio" << m_old_grain_mass_map[m_iter->first]<< " , " << mrat  << std::endl;
      // for all pairs of new grains
      for(set<int>::iterator s_iter1=m_iter->second.begin();
	  s_iter1!=m_iter->second.end();
	  s_iter1++){
	for(set<int>::iterator s_iter2=m_iter->second.begin();
	  s_iter2!=m_iter->second.end();
	  s_iter2++){
	  if(*s_iter1<*s_iter2){
	    // std::cout << "  id1,id2: " << *s_iter1 << " , " << *s_iter2 << std::endl;
	    // for all pairs of particles within the pair of grains
	    for(multimap<int,int>::const_iterator p_iter1=m_new_grain_to_particle_map.lower_bound(*s_iter1);
		p_iter1!=m_new_grain_to_particle_map.upper_bound(*s_iter1);
		p_iter1++){
	      for(multimap<int,int>::const_iterator p_iter2=m_new_grain_to_particle_map.lower_bound(*s_iter2);
		p_iter2!=m_new_grain_to_particle_map.upper_bound(*s_iter2);
		p_iter2++){
		if(orig.isEdge(p_iter1->second,p_iter2->second)) {
		  double r1=orig.getVertexData(p_iter1->second).rad;
		  double r2=orig.getVertexData(p_iter2->second).rad;
		  double rb=0.5*(r1+r2);
		  double area=0.785398*rb*rb;
		  // std::cout << p_iter1->second << " - " << p_iter2->second << " A: " << area << std::endl;
		  if(mrat>ratio) {
		    ab_area+=area;
		  } else {
		    split_area+=area;
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  return make_pair(split_area,ab_area);
}
