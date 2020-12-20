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

#include "frame_vtk.h"

// --- project includes ---
#include "vec3.h"


// --- STL includes ---
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <map>
#include <utility>
#include <sstream>
#include <boost/filesystem.hpp>

using std::vector;
using std::string;
using std::ifstream;
using std::istream_iterator;
using std::back_inserter;
using std::map;
using std::pair;
using std::make_pair;
using std::ostringstream;

struct nr_part{
  Vec3 pos;
  Vec3 init_pos;
  Vec3 vel;
  Vec3 force;
  Vec3 circ_shift;
  double rad;
  double mass;
  int tag;
  int proc_id;
};

struct r_part{
  Vec3 pos;
  Vec3 init_pos;
  Vec3 vel;
  Vec3 force;
  Vec3 circ_shift;
  double rad;
  double mass;
  int tag;
  double q1,q2,q3,q4;
  Vec3 angvel;
  int proc_id;
};

struct bond
{
  int id1;
  int id2;
  int tag;

  bond(int,int,int);
};

bond::bond(int i1,int i2,int t)
{
  id1=i1;
  id2=i2;
  tag=t;
}

/*!
  Write bonds to vtk file.

  \param vtkfile the output file stream
  \param data_map STL map containing particle data, particle ID as key. Used for bond strains.
  \param bonds bond data: particle ids and bond tag
  \param id2idx mapping between particle ID and index 
*/
template <typename T>
void write_bonds(ostream& vtkfile,map<int,T>& data_map,vector<bond>& bonds, map<int,int>& id2idx)
{
  vtkfile << "<Cells>\n";
  vtkfile << "<DataArray type=\"Int32\" NumberOfComponents=\"1\" Name=\"connectivity\" format=\"ascii\">\n";
  for(vector<bond>::iterator iter=bonds.begin();
      iter!=bonds.end();
      iter++){
    vtkfile << id2idx[iter->id1] << " " << id2idx[iter->id2] << endl;
  }
  vtkfile << "</DataArray>";
  // offsets
  vtkfile << "<DataArray type=\"Int32\" NumberOfComponents=\"1\" Name=\"offsets\" format=\"ascii\">\n";
  for (size_t i = 1; i < bonds.size()*2; i+=2) vtkfile << i+1 << "\n";
  vtkfile << "</DataArray>\n";
  // element type
  vtkfile << "<DataArray type=\"UInt8\" NumberOfComponents=\"1\" Name=\"types\" format=\"ascii\">\n";
  const int CELL_LINE_TYPE = 3;
  for (size_t i = 0; i < bonds.size(); i++) vtkfile << CELL_LINE_TYPE << "\n";
  vtkfile << "</DataArray>\n";  
  vtkfile << "</Cells>\n";
  // bond data
  vtkfile << "<CellData>\n";
  // bond strain 
  vtkfile << "<DataArray type=\"Float64\" Name=\"bondStrain\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(vector<bond>::iterator iter=bonds.begin();
      iter!=bonds.end();
      iter++){
    Vec3 pos1=data_map[iter->id1].pos;
    Vec3 pos2=data_map[iter->id2].pos;
    double r1=data_map[iter->id1].rad;
    double r2=data_map[iter->id2].rad;
    double strain=((pos1-pos2).norm()-(r1+r2))/(r1+r2);
    vtkfile << strain << endl;
    if ((pos1-pos2).norm() > 2.1) {
      // cerr << "excessive bond length " << (pos1-pos2).norm() << " in bond between " << iter->first << " at [" << pos1 << "] tag " << data_map[iter->first].tag <<
      //	" and " << iter->second << " at [" << pos2 << "] tag " << data_map[iter->second].tag <<endl;
    }
  }
  vtkfile << "</DataArray>\n"; 
  // bond tag
  vtkfile << "<DataArray type=\"Int32\" Name=\"bondTag\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(vector<bond>::iterator iter=bonds.begin();
      iter!=bonds.end();
      iter++){
    vtkfile << iter->tag << "\n";
  }
  vtkfile << "</DataArray>\n";  
  vtkfile << "</CellData>\n";
}


/*!
  Write bonds to vtk file. Version with bonded mesh IG.

  \param vtkfile the output file stream
  \param data_map STL map containing particle data, particle ID as key. Used for bond strains.
  \param bonds bond data: particle ids and bond tag
  \param mesh_bonds mesh bond data: particle ids and triangle id -> bond tag
  \param id2idx mapping between particle ID and index 
*/
template <typename T>
void write_bonds_with_mesh(ostream& vtkfile,
			   map<int,T>& data_map,
			   vector<bond>& bonds, 
			   vector<bond>& mesh_bonds, 
			   map<int,int>& id2idx)
{
  vtkfile << "<Cells>\n";
  vtkfile << "<DataArray type=\"Int32\" NumberOfComponents=\"1\" Name=\"connectivity\" format=\"ascii\">\n";
  // bonds
  for(vector<bond>::iterator iter=bonds.begin();
      iter!=bonds.end();
      iter++){
    vtkfile << id2idx[iter->id1] << " " << id2idx[iter->id2] << endl;
  }
  // mesh bonds
  for(vector<bond>::iterator iter=mesh_bonds.begin();
      iter!=mesh_bonds.end();
      iter++){
    vtkfile << id2idx[iter->id1] << " " << id2idx[iter->id2] << endl;
  }
  vtkfile << "</DataArray>";
  // offsets
  vtkfile << "<DataArray type=\"Int32\" NumberOfComponents=\"1\" Name=\"offsets\" format=\"ascii\">\n";
  for (size_t i = 1; i < (bonds.size()+mesh_bonds.size())*2; i+=2) vtkfile << i+1 << "\n";
  vtkfile << "</DataArray>\n";
  // element type
  vtkfile << "<DataArray type=\"UInt8\" NumberOfComponents=\"1\" Name=\"types\" format=\"ascii\">\n";
  const int CELL_LINE_TYPE = 3;
  for (size_t i = 0; i < (bonds.size()+mesh_bonds.size()); i++) vtkfile << CELL_LINE_TYPE << "\n";
  vtkfile << "</DataArray>\n";  
  vtkfile << "</Cells>\n";
  // bond data
  vtkfile << "<CellData>\n";
  // bond strain 
  vtkfile << "<DataArray type=\"Float64\" Name=\"bondStrain\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  // bonds
  for(vector<bond>::iterator iter=bonds.begin();
      iter!=bonds.end();
      iter++){
    Vec3 pos1=data_map[iter->id1].pos;
    Vec3 pos2=data_map[iter->id2].pos;
    double r1=data_map[iter->id1].rad;
    double r2=data_map[iter->id2].rad;
    double strain=((pos1-pos2).norm()-(r1+r2))/(r1+r2);
    vtkfile << strain << endl;
  }
  // mesh bonds -> currently set to 0
  for(vector<bond>::iterator iter=mesh_bonds.begin();
      iter!=mesh_bonds.end();
      iter++){
    vtkfile << "0.0" << endl;
  }
  vtkfile << "</DataArray>\n"; 
  // bond tag
  vtkfile << "<DataArray type=\"Int32\" Name=\"bondTag\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  // bonds
  for(vector<bond>::iterator iter=bonds.begin();
      iter!=bonds.end();
      iter++){
    vtkfile << iter->tag << "\n";
  }
  // mesh bonds
  for(vector<bond>::iterator iter=mesh_bonds.begin();
      iter!=mesh_bonds.end();
      iter++){
    vtkfile << iter->tag << "\n";
  }
  vtkfile << "</DataArray>\n";  
  // "mesh" flag - only written if there is a bonded mesh interaction
  // used to distinguish between normal and mesh bonded interactions
  vtkfile << "<DataArray type=\"Int32\" Name=\"meshFlag\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(vector<bond>::iterator iter=bonds.begin();
      iter!=bonds.end();
      iter++){
    vtkfile << "0 \n";
  }
  for(vector<bond>::iterator iter=mesh_bonds.begin();
      iter!=mesh_bonds.end();
      iter++){
    vtkfile << "1 \n";
  }
  vtkfile << "</DataArray>\n";
  vtkfile << "</CellData>\n";
}

/*!
  write vtu file header

  \param vtkfile the ofstream connected to the output file
  \param np nr. of particles
  \param nb nr. of bonds
*/
void write_vtu_header(ofstream& vtkfile,int np, int nb)
{
  vtkfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
  vtkfile << "<UnstructuredGrid>\n";
  vtkfile << "<Piece NumberOfPoints=\"" << np << "\" NumberOfCells=\"" << nb << "\">\n";
}

/*!
  write base vtu file footer

  \param vtkfile the ofstream connected to the output file
*/
void write_vtu_footer(ofstream& vtkfile)
{
  vtkfile << "</Piece>\n";
  vtkfile << "</UnstructuredGrid>\n";
  vtkfile << "</VTKFile>\n";

}

/*!
  finalize particle data section in vtu file
  
  \param vtkfile the ofstream connected to the output file
*/
void write_vtu_particle_data_end(ofstream& vtkfile)
{
  vtkfile << "</PointData>\n";
}

/*!
  write base particle data to vtu file

  \param vtkfile the ofstream connected to the output file
  \param datamap 
  \param unwrap
*/
void write_vtu_rot_base_data(ofstream& vtkfile,map<int,r_part>& datamap, bool unwrap)
{
  // write particle pos
  vtkfile << "<Points>\n";
  vtkfile << "<DataArray NumberOfComponents=\"3\" type=\"Float64\" format=\"ascii\">\n";
  for(map<int,r_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    if(unwrap){
      vtkfile << (iter->second).pos-(iter->second).circ_shift << endl;
    } else {
      vtkfile << (iter->second).pos << endl;
    }
  }  
  vtkfile << "</DataArray>\n";
  vtkfile << "</Points>\n";

  // --- write particle data ---
  // radius
  vtkfile << "<PointData Scalars=\"radius\">\n";
  vtkfile << "<DataArray type=\"Float64\" Name=\"radius\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(map<int,r_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    vtkfile << (iter->second).rad << endl;;
  } 
  vtkfile << "</DataArray>\n";
  // tag
  vtkfile << "<DataArray type=\"Int32\" Name=\"particleTag\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(map<int,r_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    vtkfile << (iter->second).tag << endl;;
  } 
  vtkfile << "</DataArray>\n";
  // tag
  vtkfile << "<DataArray type=\"Float64\" Name=\"particleMass\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(map<int,r_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    vtkfile << (iter->second).mass << endl;;
  } 
  vtkfile << "</DataArray>\n";
  // id 
  vtkfile << "<DataArray type=\"Int32\" Name=\"Id\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(map<int,r_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    vtkfile << iter->first << endl;;
  } 
  vtkfile << "</DataArray>\n";
  // velocity
  vtkfile << "<DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(map<int,r_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    vtkfile << (iter->second).vel << endl;;
  } 
  vtkfile << "</DataArray>\n";
  // angular velocity
  vtkfile << "<DataArray type=\"Float64\" Name=\"angular velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(map<int,r_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    vtkfile << (iter->second).angvel << endl;;
  } 
  vtkfile << "</DataArray>\n";
  // displacement
  vtkfile << "<DataArray type=\"Float64\" Name=\"displacement\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(map<int,r_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    vtkfile << ((iter->second).pos-(iter->second).circ_shift)-(iter->second).init_pos << endl;
    
  } 
  vtkfile << "</DataArray>\n";
  // initial position
  vtkfile << "<DataArray type=\"Float64\" Name=\"initial position\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(map<int,r_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    vtkfile << (iter->second).init_pos << endl;
    
  } 
  vtkfile << "</DataArray>\n";
  // process id
  vtkfile << "<DataArray type=\"Int32\" Name=\"proc_id\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(map<int,r_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    vtkfile << (iter->second).proc_id << endl;
    
  } 
  vtkfile << "</DataArray>\n";
}

/*!
  write base particle data to vtu file

  \param vtkfile the ofstream connected to the output file
  \param datamap 
  \param unwrap
*/
void write_vtu_nr_base_data(ofstream& vtkfile,map<int,nr_part>& datamap, bool unwrap)
{
  // write particle pos
  vtkfile << "<Points>\n";
  vtkfile << "<DataArray NumberOfComponents=\"3\" type=\"Float64\" format=\"ascii\">\n";
  for(map<int,nr_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    if(unwrap){
      vtkfile << (iter->second).pos-(iter->second).circ_shift << endl;
    } else {
      vtkfile << (iter->second).pos << endl;
    }
  }  
  vtkfile << "</DataArray>\n";
  vtkfile << "</Points>\n";

  // --- write particle data ---
  // radius
  vtkfile << "<PointData Scalars=\"radius\">\n";
  vtkfile << "<DataArray type=\"Float64\" Name=\"radius\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(map<int,nr_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    vtkfile << (iter->second).rad << endl;;
  } 
  vtkfile << "</DataArray>\n";
  // tag
  vtkfile << "<DataArray type=\"Int32\" Name=\"particleTag\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(map<int,nr_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    vtkfile << (iter->second).tag << endl;;
  } 
  vtkfile << "</DataArray>\n";
  // tag
  vtkfile << "<DataArray type=\"Float64\" Name=\"particleMass\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(map<int,nr_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    vtkfile << (iter->second).mass << endl;;
  } 
  vtkfile << "</DataArray>\n";
  // id 
  vtkfile << "<DataArray type=\"Int32\" Name=\"Id\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(map<int,nr_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    vtkfile << iter->first << endl;;
  } 
  vtkfile << "</DataArray>\n";
  // velocity
  vtkfile << "<DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(map<int,nr_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    vtkfile << (iter->second).vel << endl;;
  } 
  vtkfile << "</DataArray>\n";
  // displacement
  vtkfile << "<DataArray type=\"Float64\" Name=\"displacement\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(map<int,nr_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    vtkfile << ((iter->second).pos-(iter->second).circ_shift)-(iter->second).init_pos << endl;
    
  } 
  vtkfile << "</DataArray>\n";
  // initial position
  vtkfile << "<DataArray type=\"Float64\" Name=\"initial position\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(map<int,nr_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    vtkfile << (iter->second).init_pos << endl;
    
  } 
  vtkfile << "</DataArray>\n";
  // process id
  vtkfile << "<DataArray type=\"Int32\" Name=\"proc_id\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(map<int,nr_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    vtkfile << (iter->second).proc_id << endl;
    
  } 
  vtkfile << "</DataArray>\n";
}

/*!
  Get snapshot version.

  \param infilename the file name of the header file (incl. the _0.txt part) 
*/
int get_version(const string& infilename)
{
  string dummystring;
  int version;
  ifstream headerfile(infilename.c_str()); 
  // read token  
  headerfile >> dummystring;

  if(dummystring=="V"){ // if V -> new version 
    headerfile >> version ;
    cout << "version : " << version << endl;
  } else {
    cout << "pre- V.1 version" << endl;
    version=0;
  }
  headerfile.close();

  return version;
}

/*!
  Get data file names from the header file

  \param infilename the file name of the header file (incl. the _0.txt part)
  \param the snapshot version
  \param relpath
*/
vector<string> get_filenames(const string& infilename, int version, bool relpath=false)
{
  cout << "infilename : " << infilename << endl;
  ifstream headerfile(infilename.c_str()); 
  if(!headerfile){
    std::cerr << "header file " << infilename << " doesn't exist - exiting " << std::endl;
    std::exit(1);
  }

  float dummy,xmax,ymax,zmax,xmin,ymin,zmin,geo_version;
  vector<string> filenames;
  string dummystring;

  if(version==0){
    headerfile >> dummy >> dummy >> dummy;
    headerfile >> dummystring >> dummy;
  } else if ((version==1)||(version==2) || (version==3)) {
    headerfile >> dummystring >> dummy;
    headerfile >> dummy >> dummy >> dummy;
    std::getline(headerfile,dummystring);
  } else {
    cerr << "unknown checkpoint version " << version << endl;
  }
  // get bounding box
  std::getline(headerfile,dummystring);
  headerfile >> dummystring >> xmin >> ymin >> zmin >> xmax >> ymax >> zmax ;

  // ignore periodic bdry
  headerfile >> dummystring >> dummy >> dummy >> dummy;

  // v. 1.1 geometry files didn't have dimension info
  if((geo_version>1.15) || (version>2) ){
    headerfile >> dummystring >> dummystring;
  }

  // get file names
  copy(istream_iterator<string>(headerfile),istream_iterator<string>(),back_inserter(filenames));
  
  headerfile.close();

  cout << "nr. of filenames: " << filenames.size() << endl;

  if(relpath){
    // get path to header file
    boost::filesystem::path if_path(infilename);
    if_path.remove_filename();

    // change paths in date file names
    vector<string> mod_fn;
    for(vector<string>::iterator iter=filenames.begin();
	iter!=filenames.end();
	iter++)
      {
	boost::filesystem::path np(*iter);
	boost::filesystem::path p=if_path;
	p/=(np.filename());
	mod_fn.push_back(p.string());
      }
    filenames=mod_fn;
  }

  return filenames;
}

/*!
  Read snapshot with non-rotational particles and write it to vtk unstructured grid (.vtu) file

  \param infilename the name of the header file
  \param outfilename the "base" name of the output file
  \param iframe the frame number
  \param with_list
  \param listfilename
*/
void do_single_frame_vtk(const string& infilename,const string& outfilename,int iframe,bool with_list,const string& listfilename,bool remove_xbonds,double bond_remove_dist)
{
  int version=get_version(infilename);
  vector<string> filenames=get_filenames(infilename,version);
  map<int,nr_part> datamap;
  vector<bond> bonds;
  map<int,int> id2idx;
  bool hasMeshBondedInteractions=false;

  int proc_cnt=0;
  for(vector<string>::iterator iter=filenames.begin();
      iter!=filenames.end();
      iter++){
    proc_cnt++;
    cout << *iter << endl;
    ifstream datafile(iter->c_str());
    vector<double> pdata;
    // get particles
    Vec3 circ_shift;
    int npart;
    nr_part data;
    int id;
    Vec3 pos;
    Vec3 oldpos;
    Vec3 initpos;
    datafile >> npart;
    if(version<2){
      for(int i=0;i<npart;i++){
	datafile >> data.pos  >> data.rad >> id >> data.tag >> data.mass >> data.init_pos >> oldpos >> data.vel >> data.force;
	data.proc_id=proc_cnt;
	datamap[id]=data;
      }
    } else {
      for(int i=0;i<npart;i++){
	datafile >> data.pos  >> data.rad >> id >> data.tag >> data.mass >> data.init_pos >> oldpos >> data.vel >> data.force >> data.circ_shift;
	data.proc_id=proc_cnt;
	datamap[id]=data;
      }
    }

    // get bonds
    int ngrp;
    int nbond;
    string type;

    datafile >> ngrp;
    for(int ni=0;ni<ngrp;ni++){
      if(version>0){
	datafile >> type;
	if(type=="Bonded"){
	  datafile >> nbond;
	  if(remove_xbonds){
	    for(int i=0;i<nbond;i++){
	      int b1,b2,tag;
	      
	      datafile >> b1 >> b2 >> tag;
	      double dist=(datamap[b1].pos-datamap[b2].pos).norm();
	      if(dist < bond_remove_dist){
		if(version > 2){
		  bonds.push_back(bond(b1,b2,tag));
		} else {
		  bonds.push_back(bond(b1,b2,ni));
		}
	      }
	    } 
	  }  else {
	    for(int i=0;i<nbond;i++){
	      int b1,b2,tag;
	      
	      datafile >> b1 >> b2 >> tag;
	      if(version > 2){
		bonds.push_back(bond(b1,b2,tag));
	      } else {
		bonds.push_back(bond(b1,b2,ni));
	      }
	    }
	  }
	}
      } else { // pre - V1 snapshot -> assume bondend pair IG
	datafile >> nbond;
	for(int i=0;i<nbond;i++){
	  int b1,b2,tag;
	  
	  datafile >> b1 >> b2 >> tag;
	  if(version > 2){
	    bonds.push_back(bond(b1,b2,tag));
	  } else {
	    bonds.push_back(bond(b1,b2,ni));
	  }
	}
      }
    }


    datafile.close();
  }

  // generate reverse mapping between particle id and point idx
  int count=0;
  for(map<int,nr_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    id2idx.insert(make_pair(iter->first,count));
    count++;
  }
  // open output file
  ostringstream vtkfilename;
  vtkfilename << outfilename << iframe << ".vtu";
  ofstream vtkfile(vtkfilename.str().c_str());
  // write header 
  write_vtu_header(vtkfile, datamap.size(), bonds.size());
  
  // write particle base data
  write_vtu_nr_base_data(vtkfile,datamap,false);

  // end particle data
  write_vtu_particle_data_end(vtkfile);

 
  // write bonds 
  if(hasMeshBondedInteractions) {
    // write_bonds_with_mesh(vtkfile,datamap,bonds,id2idx);
  } else {
    write_bonds(vtkfile,datamap,bonds,id2idx);
  }

  // write footer
  write_vtu_footer(vtkfile);

  // close file
  vtkfile.close();
}

/*!
  Read snapshot with rotational particles and write it to vtk unstructured grid (.vtu) file

  \param infilename the name of the header file
  \param outfilename the "base" name of the output file
  \param iframe the frame number
  \param with_list
  \param listfilename
  \param remove_xbonds if "true", bonds crossing a circular boundary are not written to the file
  \param bond_remove_dist distance threshold to detect if a  bond is crossing a circular boundary
  \param with_brklist
  \param brklistname
  \param unwrap if "true", the model is unwrapped, i.e. movements across a circular boundary are taken into account
  \param relpath if "true" the filenames in the header file are taken to be relative to the path of the header
*/

void do_single_frame_vtk_r(const string& infilename,const string& outfilename,int iframe,bool with_list,const string& listfilename,bool remove_xbonds,double bond_remove_dist,bool with_brklist,const string& brklistname, bool unwrap, bool relpath)
{  
  int version=get_version(infilename);
  vector<string> filenames=get_filenames(infilename,version,relpath);
  map<int,r_part> datamap;
  map<int,int> ng_id_map;
  map<int,float> ng_mass_map;
  vector<bond> bonds;
  vector<bond> mesh_bonds;
  map<int,int> id2idx;
  map<int,int> id_brk_map;
  bool hasMeshBondedInteractions=false;

  int proc_cnt=0;
  int meshbond_count=1;
  for(vector<string>::iterator iter=filenames.begin();
      iter!=filenames.end();
      iter++){
    proc_cnt++;
    cout << *iter << " , " << proc_cnt << endl;
    ifstream datafile(iter->c_str());
    vector<double> pdata;
    // get particles
    int npart;
    r_part data;
    int id;
    Vec3 pos;
    Vec3 oldpos;
    datafile >> npart;
    if(version < 2){
      for(int i=0;i<npart;i++){
	datafile >> data.pos  >> data.rad >> id >> data.tag >> data.mass >> data.init_pos >> oldpos >> data.vel >> data.force
		 >> data.q1 >> data.q2 >> data.q3 >> data.q4 >> data.angvel;
	data.proc_id=proc_cnt;
	datamap[id]=data;
      } 
    } else {
      for(int i=0;i<npart;i++){
	datafile >> data.pos  >> data.rad >> id >> data.tag >> data.mass >> data.init_pos >> oldpos >> data.vel >> data.force >> data.circ_shift
		 >> data.q1 >> data.q2 >> data.q3 >> data.q4 >> data.angvel;
	data.proc_id=proc_cnt;
	if(unwrap){
	  
	  data.pos=data.pos;
	}
	datamap[id]=data;
      } 
    }  

    // get bonds
    int ngrp;
    int nbond;
    string type;

    
    datafile >> ngrp;
    for(int ni=0;ni<ngrp;ni++){
      if(version>0){
	datafile >> type;
	if((type=="Bonded") || (type=="RotBonded")){
	  datafile >> nbond;
	  if(remove_xbonds){
	    for(int i=0;i<nbond;i++){
	      int b1,b2,tag;
	    
	      datafile >> b1 >> b2 >> tag;
	      double dist=(datamap[b1].pos-datamap[b2].pos).norm();
	      if(dist < bond_remove_dist){
		if(version > 2){
		  bonds.push_back(bond(b1,b2,tag));
		} else {
		  bonds.push_back(bond(b1,b2,ni));
		}
	      }
	    }
	  } else {
	    for(int i=0;i<nbond;i++){
	      int b1,b2,tag;
	      
	      datafile >> b1 >> b2 >> tag;
	      if(version > 2){
		bonds.push_back(bond(b1,b2,tag));
	      } else {
		bonds.push_back(bond(b1,b2,ni));
	      }
	    }
	  }
	} else if ((type=="BondedMesh2D") || (type=="BondedTriMesh")){
	  hasMeshBondedInteractions=true;
	  datafile >> nbond;
	  for(int i=0;i<nbond;i++){
	    int b1,b2;
	    Vec3 ap_pos;
    
	    datafile >> b1 >> b2 >> ap_pos;	    
	    // add anchor point to point list
	    r_part ap;
	    ap.pos=ap_pos;
	    ap.init_pos=ap_pos;
	    ap.vel=Vec3(0,0,0);
	    ap.force=Vec3(0,0,0);
	    ap.circ_shift=Vec3(0,0,0);
	    ap.rad=0.0;
	    ap.mass=0.0;
	    ap.tag=-1*b2;
	    ap.q1=0;
	    ap.q2=0;
	    ap.q3=0;
	    ap.q4=0;
	    ap.angvel=Vec3(0,0,0);
	    ap.proc_id=proc_cnt;
	    int ap_id=-1*meshbond_count;
	    datamap[ap_id]=ap;
	    mesh_bonds.push_back(bond(b1,ap_id,b2));
	    meshbond_count++;
	  }
	} 
      } else { // pre - V1 snapshot -> assume bondend pair IG
	datafile >> nbond;
	for(int i=0;i<nbond;i++){
	  int b1,b2,tag;
	
	  datafile >> b1 >> b2 >> tag;
	  if(version > 2){
	    bonds.push_back(bond(b1,b2,tag));
	  } else {
	    bonds.push_back(bond(b1,b2,ni));
	  }
	}
      }
    }
    datafile.close();
  }

  // generate reverse mapping between particle id and point idx
  int count=0;
  for(map<int,r_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    id2idx.insert(make_pair(iter->first,count));
    count++;
  }

  // if option set, read list file into ng_id_map;
  if(with_list){
    ifstream listfile(listfilename.c_str());
    int id, nid;
    float gmass;
    while(!listfile.eof()){
      listfile >> id >> nid >> gmass;
      ng_id_map[id]=nid;
      ng_mass_map[id]=gmass;
    }
    listfile.close();
  }
  // if option set, read fracture list file int id_brk_map;
  if(with_brklist){
    std::cerr << "opening " << brklistname << std::endl;
    ifstream brklistfile(brklistname.c_str());
    int id, nbrk;
    while(!brklistfile.eof()){
      brklistfile >> id >> nbrk;
      id_brk_map[id]=nbrk;
    }
    brklistfile.close();
  }
  
  // open output file
  ostringstream vtkfilename;
  vtkfilename << outfilename << iframe << ".vtu";
  ofstream vtkfile(vtkfilename.str().c_str());
  // write header 
  write_vtu_header(vtkfile, datamap.size(), bonds.size()+mesh_bonds.size());
  // write particle base data
  write_vtu_rot_base_data(vtkfile,datamap,unwrap);

  if(with_list){
    
    // id
    vtkfile << "<DataArray type=\"Int32\" Name=\"Current GrainId\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for(map<int,r_part>::iterator iter=datamap.begin();
	iter!=datamap.end();
	iter++){
      vtkfile << ng_id_map[iter->first] << endl;;
    }
    vtkfile << "</DataArray>\n";
    
    // mass
    vtkfile << "<DataArray type=\"Float32\" Name=\"GrainMass\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for(map<int,r_part>::iterator iter=datamap.begin();
	iter!=datamap.end();
	iter++){
      vtkfile << ng_mass_map[iter->first] << endl;;
    }
    vtkfile << "</DataArray>\n";
  } 

  // if option set , number of broken bonds per particle 
  if(with_brklist){
    vtkfile << "<DataArray type=\"Int32\" Name=\"broken bonds\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    for(map<int,r_part>::iterator iter=datamap.begin();
	iter!=datamap.end();
	iter++){
      vtkfile << id_brk_map[iter->first] << endl;;
    } 
    vtkfile << "</DataArray>\n";
  }

  // end particle data
  write_vtu_particle_data_end(vtkfile);

  // write bonds 
  if(hasMeshBondedInteractions) {
    write_bonds_with_mesh(vtkfile,datamap,bonds,mesh_bonds,id2idx);
  } else {
    write_bonds(vtkfile,datamap,bonds,id2idx);
  }

  // write footer
  write_vtu_footer(vtkfile);

  // close file
  vtkfile.close();
}


/*!
  write particles of single tag from frame

  \param infilename the input file name
  \param outfilename the output file name base
  \param iframe frame nr.
  \param ptag particle tag to write out
  \param unwrap unwrap x-circular boundaries
*/
void do_single_frame_vtk_single_r(const string& infilename,const string& outfilename,int iframe,int ptag,bool unwrap,bool remove_xbonds,double bond_remove_dist)
{
  int version=get_version(infilename);
  vector<string> filenames=get_filenames(infilename,version);
  map<int,r_part> datamap;
  map<int,int> ng_id_map;
  vector<bond> bonds;
  map<int,int> id2idx;
  bool hasMeshBondedInteractions=false;

  for(vector<string>::iterator iter=filenames.begin();
      iter!=filenames.end();
      iter++){
    cout << *iter << endl;
    ifstream datafile(iter->c_str());
    vector<double> pdata;
    // get particles
    Vec3 circ_shift;
    int npart;
    r_part data;
    int id;
    Vec3 pos;
    Vec3 oldpos;
    Vec3 initpos;
    datafile >> npart;
    if(version < 2){
      for(int i=0;i<npart;i++){
	datafile >> data.pos  >> data.rad >> id >> data.tag >> data.mass >> initpos >> oldpos >> data.vel >> data.force
		 >> data.q1 >> data.q2 >> data.q3 >> data.q4 >> data.angvel;
	if(data.tag==ptag){
	  datamap[id]=data;
	}
      } 
    } else {
     for(int i=0;i<npart;i++){
       datafile >> data.pos  >> data.rad >> id >> data.tag >> data.mass >> initpos >> oldpos >> data.vel >> data.force >> circ_shift
		 >> data.q1 >> data.q2 >> data.q3 >> data.q4 >> data.angvel;
       if(data.tag==ptag){
	 datamap[id]=data;
       }
      } 
    }  

    // get bonds
// get bonds
    int ngrp;
    int nbond;
    string type;

    
    datafile >> ngrp;
    for(int ni=0;ni<ngrp;ni++){
      if(version>0){
	datafile >> type;
	if((type=="Bonded") || (type=="RotBonded")){
	  datafile >> nbond;
	  if(remove_xbonds){
	    for(int i=0;i<nbond;i++){
	      int b1,b2,tag;
	    
	      datafile >> b1 >> b2 >> tag;
	      double dist=(datamap[b1].pos-datamap[b2].pos).norm();
	      if((b1==ptag) && (b2==ptag)){
		if(dist < bond_remove_dist){
		  if(version > 2){
		    bonds.push_back(bond(b1,b2,tag));
		  } else {
		    bonds.push_back(bond(b1,b2,ni));
		  }
		}
	      }
	    }
	  } else {
	    for(int i=0;i<nbond;i++){
	      int b1,b2,tag;
	      
	      datafile >> b1 >> b2 >> tag;
	      if((datamap[b1].tag==ptag) && (datamap[b2].tag==ptag)){
		if(version > 2){
		  bonds.push_back(bond(b1,b2,tag));
		} else {
		  bonds.push_back(bond(b1,b2,ni));
		}
	      }
	    }
	  }
	} else if ((type=="BondedMesh2D") || (type=="BondedTriMesh")){
	  hasMeshBondedInteractions=true;
	  datafile >> nbond;
	  for(int i=0;i<nbond;i++){
	    int b1,b2;
	    Vec3 ap;
    
	    datafile >> b1 >> b2 >> ap;	    
	  }
	} 
      } else { // pre - V1 snapshot -> assume bondend pair IG
	datafile >> nbond;
	for(int i=0;i<nbond;i++){
	  int b1,b2,tag;
	
	  datafile >> b1 >> b2 >> tag;
	  if((datamap[b1].tag==ptag) && (datamap[b2].tag==ptag)){
	    bonds.push_back(bond(b1,b2,ni));
	  }
	}
      }
    }   
    datafile.close();
  }

  // generate reverse mapping between particle id and point idx
  int count=0;
  for(map<int,r_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    id2idx.insert(make_pair(iter->first,count));
    count++;
  }

  // open output file
  ostringstream vtkfilename;
  vtkfilename << outfilename << iframe << ".vtu";
  ofstream vtkfile(vtkfilename.str().c_str());
  // write header 
  write_vtu_header(vtkfile, datamap.size(), bonds.size());
  
  // write particle base data
  write_vtu_rot_base_data(vtkfile,datamap,unwrap);

  // end particle data
  write_vtu_particle_data_end(vtkfile);


  // write bonds 
  if(hasMeshBondedInteractions) {
    //write_bonds_with_mesh(vtkfile,datamap,bonds,id2idx);
  } else {
    write_bonds(vtkfile,datamap,bonds,id2idx);
  }

  // write footer
  write_vtu_footer(vtkfile);

  // close file
  vtkfile.close();
}

void do_single_frame_sliced_vtk(const string&,const string&,int,bool,const string&,double,double)
{}

/*!
  write a slice , doesn't write bonds
*/
void do_single_frame_sliced_vtk_r(const string& infilename,const string& outfilename,int iframe,bool with_list,const string& listfilename, double slz_min,double slz_max)
{  
  cout << "do_single_frame_sliced_vtk_r(.. " << slz_min << "," << slz_max << ")" << endl;

  int version=get_version(infilename);
  vector<string> filenames=get_filenames(infilename,version);
  map<int,r_part> datamap;
  map<int,int> ng_id_map;
  vector<pair<int,int> > bonds;
  map<int,int> id2idx;
  bool hasMeshBondedInteractions=false;
 
  for(vector<string>::iterator iter=filenames.begin();
      iter!=filenames.end();
      iter++){
    cout << *iter << endl;
    ifstream datafile(iter->c_str());
    vector<double> pdata;
    // get particles
    //float dummy=0.0;
    Vec3 circ_shift;
    int npart=0;
    r_part data;
    int id;
    Vec3 pos;
    Vec3 oldpos;
    Vec3 initpos;
    datafile >> npart;
    if(version < 2){
      for(int i=0;i<npart;i++){
	datafile >> data.pos  >> data.rad >> id >> data.tag >> data.mass >> initpos >> oldpos >> data.vel >> data.force
		 >> data.q1 >> data.q2 >> data.q3 >> data.q4 >> data.angvel;
	if((data.pos.Z()>slz_min) && (data.pos.Z()<slz_max)){ // if in slice, add to map
	  datamap[id]=data;
	}
      } 
    } else {
     for(int i=0;i<npart;i++){
       datafile >> data.pos  >> data.rad >> id >> data.tag >> data.mass >> initpos >> oldpos >> data.vel >> data.force >> circ_shift
		 >> data.q1 >> data.q2 >> data.q3 >> data.q4 >> data.angvel;
       if((data.pos.Z()>slz_min) && (data.pos.Z()<slz_max)){ // if in slice, add to map
	 datamap[id]=data;
       }
      } 
    }  

    datafile.close();
  }

  // generate reverse mapping between particle id and point idx
  int count=0;
  for(map<int,r_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    id2idx.insert(make_pair(iter->first,count));
    count++;
  }

  // if option set, read list file int ng_id_map;
  if(with_list){
    ifstream listfile(listfilename.c_str());
    int id, nid;
    while(!listfile.eof()){
      listfile >> id >> nid;
      ng_id_map[id]=nid;
    }
    listfile.close();
  }
  // open output file
  ostringstream vtkfilename;
  vtkfilename << outfilename << iframe << ".xml";
  ofstream vtkfile(vtkfilename.str().c_str());
  // write header 
  vtkfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";
  vtkfile << "<UnstructuredGrid>\n";
  vtkfile << "<Piece NumberOfPoints=\"" << datamap.size()<< "\" NumberOfCells=\"" << bonds.size() << "\">\n";
  // write particle pos (flatten in Z to middle of slice area)
  vtkfile << "<Points>\n";
  vtkfile << "<DataArray NumberOfComponents=\"3\" type=\"Float64\" format=\"ascii\">\n";
  for(map<int,r_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    vtkfile << (iter->second).pos.X() << " " << (iter->second).pos.Y() << " " << (slz_min+slz_max)*0.5 << endl;;
  }  
  vtkfile << "</DataArray>\n";
  vtkfile << "</Points>\n";

  // --- write particle data ---
  // radius
  vtkfile << "<PointData Scalars=\"radius\">\n";
  vtkfile << "<DataArray type=\"Float64\" Name=\"radius\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(map<int,r_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    vtkfile << (iter->second).rad << endl;;
  } 
  vtkfile << "</DataArray>\n";
  // tag
  vtkfile << "<DataArray type=\"Int32\" Name=\"particleTag\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(map<int,r_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    vtkfile << (iter->second).tag << endl;;
  } 
  vtkfile << "</DataArray>\n";
  // id 
  vtkfile << "<DataArray type=\"Int32\" Name=\"Id\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  for(map<int,r_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    vtkfile << iter->first << endl;;
  } 
  vtkfile << "</DataArray>\n";
  // velocity
  vtkfile << "<DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(map<int,r_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    vtkfile << (iter->second).vel << endl;;
  } 
  vtkfile << "</DataArray>\n";
  // angular velocity
  vtkfile << "<DataArray type=\"Float64\" Name=\"angular velocity\" NumberOfComponents=\"3\" format=\"ascii\">\n";
  for(map<int,r_part>::iterator iter=datamap.begin();
      iter!=datamap.end();
      iter++){
    vtkfile << (iter->second).angvel << endl;;
  } 
  vtkfile << "</DataArray>\n";
  // if option set , current grain id
  if(with_list){
    vtkfile << "<DataArray type=\"Int32\" Name=\"Current GrainId\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for(map<int,r_part>::iterator iter=datamap.begin();
	iter!=datamap.end();
	iter++){
      vtkfile << ng_id_map[iter->first] << " 0 0" << endl;;
    }
    vtkfile << "</DataArray>\n";
  }
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

  // close file
  vtkfile.close();
}

void writeMeshFile(const string& infilename,const string& outfilename, int snapNum){
  ifstream datafile(infilename.c_str());

  string header,skip;
  int numMeshIG;

  while (datafile >> header){
    if (header == "TMIG"){
      datafile >> numMeshIG;
      for(int ni=0;ni<numMeshIG;ni++){
        int numNodes;
        int numTri;
        string meshName;
        datafile >> meshName;
        datafile >> skip;
        datafile >> numNodes;

        ostringstream vtkfilename;
        vtkfilename << outfilename << "_Mesh" << ni << "_" << snapNum << ".vtk";
        ofstream outputfile(vtkfilename.str().c_str());

        outputfile << "# vtk DataFile Version 1.0\n";
        outputfile << meshName << " \n";
        outputfile << "ASCII\n\n";
        outputfile << "DATASET POLYDATA\n";
        outputfile << "POINTS " << numNodes << " float\n";

        std::map<int,int> node_map; // we need an id->position map for the nodes in case node numbering doesn't start at 0
        for(int nn=0;nn<numNodes;nn++){
           double X, Y, Z; 
	   int nid; // node id
           datafile >> nid >> skip >> skip >> X >> Y >> Z;
           outputfile << X << " " << Y << " " << Z << "\n";
	   node_map[nid]=nn;
        }

        datafile >> skip;
        datafile >> numTri;

        outputfile << "\nPOLYGONS " << numTri << " " << numTri*4 << "\n";

        for(int nn=0;nn<numTri;nn++){
           int P1, P2, P3; // corner ids in the snapshot
           datafile >> skip >> skip >> P1 >> P2 >> P3;
	   // convert node ids to node positions
	   auto cp1=node_map.find(P1);
	   auto cp2=node_map.find(P2);
	   auto cp3=node_map.find(P3);
	   
	   if(cp1==node_map.end() || cp2==node_map.end() || cp3==node_map.end()){
		std::cerr << "wrong node id found while processing mesh " << meshName << std::endl;  
	   } else {
		outputfile << "3 " << cp1->second << " " << cp2->second << " " << cp3->second << "\n";
	   }
        }
      }
    }
  }
}

