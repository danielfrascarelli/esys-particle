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

#include "read.h"

// --- STL includes ---
#include <vector>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <map>

using std::istream_iterator;
using std::back_inserter;
using std::ifstream;
using std::vector;
using std::map;
using std::cerr;
using std::cout;
using std::endl;
using std::floor;

// --- project includes ---
#include "Foundation/vec3.h"

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

vector<string> get_filenames(const string& infilename, int version)
{
  cout << "infilename : " << infilename << endl;
  ifstream headerfile(infilename.c_str()); 
  float dummy,xmax,ymax,zmax,xmin,ymin,zmin;
  vector<string> filenames;
  string dummystring;

  if(version==0){
    headerfile >> dummy >> dummy >> dummy;
    headerfile >> dummystring >> dummy;
  } else if ((version==1) || (version==2) || (version==3)){
    headerfile >> dummystring >> dummy;
    headerfile >> dummy >> dummy >> dummy;
    headerfile >> dummystring >> dummy;
  } else {
    cerr << "unknown checkpoint version " << version << endl;
  }
  // get bounding box
  headerfile >> dummystring;
  headerfile >> xmin >> ymin >> zmin >> xmax >> ymax >> zmax ;

  // ignore periodic bdry
  headerfile >> dummystring >> dummy >> dummy >> dummy;

  // ignore dimension
  headerfile >> dummystring >> dummystring;

  // get file names
  copy(istream_iterator<string>(headerfile),istream_iterator<string>(),back_inserter(filenames));
  
  headerfile.close();

  cout << "nr. of filenames: " << filenames.size() << endl;
  return filenames;
}

/*!
  read snapshot and write y-profile of x-displacement

  \param infilename name of the input file (the *_0.txt)
  \param outfilename name of the output file
  \param ymin minimum of the y-range
  \param ymax maximum of the y-range
  \param nbin number of bins
  \param dir direction of binning : 0->x, 1->y, 2->z
  \param mintag minimum particle tag considered 
*/
void read_and_write_profile_r(const string& infilename,const string& outfilename,double min,double max,int nbin, bool debug_on,bool grad,int dir,int mintag)
{
  int version=get_version(infilename);
  if(version<2) {
    cerr << "snapshot version (" << version << ") < 2 - can't calculate displacements" << endl;
  } else {

    vector<string> filenames=get_filenames(infilename,version);
    double binsize=(max-min)/double(nbin);
    vector<double> mass_vec=vector<double>(nbin,0.0);
    vector<double> disp_vec=vector<double>(nbin,0.0);
    
    for(vector<string>::iterator iter=filenames.begin();
	iter!=filenames.end();
	iter++){
      cout << *iter << endl;
      ifstream datafile(iter->c_str());
      Vec3 pos,initpos,oldpos,vel,force,angvel,circ_shift;
      double rad,mass,q1,q2,q3,q4;
      int id,tag,npart;
      
      datafile >> npart;
      for(int i=0;i<npart;i++){
	int bin_idx=-1;
	datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force >> circ_shift
		 >> q1 >> q2 >> q3 >> q4 >> angvel;
	if(tag>=mintag){
	  double xdisp=pos.X()-(circ_shift.X()+initpos.X());
	  switch(dir){
	  case 0 : bin_idx=int(floor(pos.X()-min)/binsize); break;
	  case 1 : bin_idx=int(floor(pos.Y()-min)/binsize); break;
	  case 2 : bin_idx=int(floor(pos.Z()-min)/binsize); 
	  }
	  // check validity of bin
	  if((bin_idx>=0) && (bin_idx<nbin)) {
	    if(debug_on){
	      cout << pos.Y() << " - " << bin_idx << " - " << xdisp << "[" << circ_shift.X() << "]" << endl;
	    }
	    mass_vec[bin_idx]+=mass;
	    disp_vec[bin_idx]+=mass*xdisp;
	  }
	} else {
	  if(debug_on){
	    cout << "outside: " << pos.Y() << " - " << bin_idx << endl;
	  }
	}
      }
      datafile.close();
    }
    
    // generate profile
    vector<double> disp=vector<double>(nbin,0.0);;
    for(int i=0;i<nbin;i++){
      if(mass_vec[i]!=0.0){
	disp[i]=disp_vec[i]/mass_vec[i];
      } else {
	disp[i]=0.0;
      }
    }
    // write profile
    ofstream outfile(outfilename.c_str());
    if(grad){
       for(int i=1;i<nbin;i++){
	 outfile << min+(double(i))*binsize << " " << (disp[i-1]-disp[i])/binsize << endl;
      }
    } else {
      for(int i=0;i<nbin;i++){
	outfile << min+(double(i)+0.5)*binsize << " " << disp[i] << endl;
      }
    }

    outfile.close();
  }
}
/*!
  read 2 snapshot and write y-profile of x-displacement between snapshots

  \param infilename1 name of the 1st input file (the *_0.txt)
  \param infilename2 name of the 2nd input file (the *_0.txt)
  \param outfilename name of the output file
  \param ymin minimum of the y-range
  \param ymax maximum of the y-range
  \param nbin number of bins
  \param dir direction of binning : 0->x, 1->y, 2->z
*/
void read_and_write_profile_rel(const string& infilename1,
				const string& infilename2,
				const string& outfilename,
				double min,
				double max,
				int nbin, 
				bool debug_on,
				bool grad,
				int dir)
{
  int version=get_version(infilename1);
  if(version<2) {
    cerr << "snapshot version (" << version << ") < 2 - can't calculate displacements" << endl;
  } else {

    double binsize=(max-min)/double(nbin);
    // read 1st file
    vector<string> filenames=get_filenames(infilename1,version);
    map<int,double> posmap_1; // positions 
    map<int,double> dispmap_1; // displacements 
    map<int,double> massmap_1; // mass

    for(vector<string>::iterator iter=filenames.begin();
	iter!=filenames.end();
	iter++){
      cout << *iter << endl;
      ifstream datafile(iter->c_str());
      Vec3 pos,initpos,oldpos,vel,force,angvel,circ_shift;
      double rad,mass,q1,q2,q3,q4;
      int id,tag,npart;
      
      datafile >> npart;
      for(int i=0;i<npart;i++){
	datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force >> circ_shift
		 >> q1 >> q2 >> q3 >> q4 >> angvel;
	double xdisp=pos.X()-(circ_shift.X()+initpos.X());
	switch(dir){
	case 0 : posmap_1[id]=pos.X();break;
	case 1 : posmap_1[id]=pos.Y();break;
	case 2 : posmap_1[id]=pos.Z();
	}
	dispmap_1[id]=xdisp;
	massmap_1[id]=mass;
      }
      datafile.close();
    }
    
    // read 2nd file
    filenames=get_filenames(infilename2,version);
     map<int,double> posmap_2; // positions 
    map<int,double> dispmap_2; // displacements 
    map<int,double> massmap_2; // mass
    map<int,double> rmap; // radius
 
    for(vector<string>::iterator iter=filenames.begin();
	iter!=filenames.end();
	iter++){
      cout << *iter << endl;
      ifstream datafile(iter->c_str());
      Vec3 pos,initpos,oldpos,vel,force,angvel,circ_shift;
      double rad,mass,q1,q2,q3,q4;
      int id,tag,npart;
      
      datafile >> npart;
      for(int i=0;i<npart;i++){
	datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force >> circ_shift
		 >> q1 >> q2 >> q3 >> q4 >> angvel;
	double xdisp=pos.X()-(circ_shift.X()+initpos.X());
	switch(dir){
	case 0 : posmap_2[id]=pos.X();break;
	case 1 : posmap_2[id]=pos.Y();break;
	case 2 : posmap_2[id]=pos.Z();
	}
	dispmap_2[id]=xdisp;
	massmap_2[id]=mass;
	rmap[id]=rad;
      }
      datafile.close();
    }
    

    // binning
    vector<double> disp_vec=vector<double>(nbin,0.0);
    vector<double> mass_vec=vector<double>(nbin,0.0);
    for(map<int,double>::iterator iter=posmap_2.begin();
	iter!=posmap_2.end();
	iter++){
      int idx=iter->first;
      // calc min/max dims of sphere
      double xmin=iter->second-rmap[idx];
      double xmax=iter->second+rmap[idx];
      
      // get min/max bins
      int bin_idx_min=int(floor((xmin-min)/binsize));
      bin_idx_min = bin_idx_min < 0 ? 0 : bin_idx_min;
      int bin_idx_max=int(floor((xmax-min)/binsize));
       bin_idx_max= bin_idx_max > nbin-1 ? nbin-1 : bin_idx_max;

       //     std::cerr << iter->second << " " << rmap[idx] << " " << xmin << " " << xmax << " - " <<  bin_idx_min << " " << bin_idx_max << std::endl;

       // calculate how much goes into which bin
       for(int i=bin_idx_min;i<=bin_idx_max;i++){
	 double bin_floor=double(i)*binsize+min;
	 double bin_ceil=double(i+1)*binsize+min;
	 if((xmin>bin_floor) && (xmax<bin_ceil)){ // fully inside bin
	   mass_vec[i]+=massmap_2[idx];
	   disp_vec[i]+=massmap_2[idx]*(dispmap_2[idx]-dispmap_1[idx]);
	 } else if ((xmin>bin_floor) && (xmax>bin_ceil)){ // lower cap
	   // calc cap volume
	   double h=xmin-bin_floor;
	   double R=rmap[idx];
	   double vcap=1.0471976*h*h*(3*R-h); // 1/3*pi*h^2*(3R-h)
	   double vsph=4.1887902*R*R*R; // 4/3*pi*R^3
	   // calc cap mass
	   double mass=massmap_2[idx]*(vcap/vsph);
	   mass_vec[i]+=mass;
	   disp_vec[i]+=mass*(dispmap_2[idx]-dispmap_1[idx]);
	 } else if ((xmin<bin_floor) && (xmax<bin_ceil)){ // upper cap
	    // calc cap volume
	   double h=bin_ceil-xmax;
	   double R=rmap[idx];
	   double vcap=1.0471976*h*h*(3*R-h); // 1/3*pi*h^2*(3R-h)
	   double vsph=4.1887902*R*R*R; // 4/3*pi*R^3
	   // calc cap mass
	   double mass=massmap_2[idx]*(vcap/vsph);
	   mass_vec[i]+=mass;
	   disp_vec[i]+=mass*(dispmap_2[idx]-dispmap_1[idx]);
	 } else if ((xmin<bin_floor) && (xmax>bin_ceil)){ // segment
	   //calc segment volume
	   double h=binsize;
	   double d=iter->second-bin_floor;
	   double R=rmap[idx];
	   double vseg=3.1415927*h*(R*R-d*d-h*d-h*h/3.0); // pi*h*(R^2-h^2-hd-1/3h^2)
	   double vsph=4.1887902*R*R*R; // 4/3*pi*R^3
	   // calc segment mass
	   double mass=massmap_2[idx]*(vseg/vsph);
	   mass_vec[i]+=mass;
	   disp_vec[i]+=mass*(dispmap_2[idx]-dispmap_1[idx]);
	 } else { // can't happen
	   std::cerr << "impossible case, xmin=" << xmin << " xmax= " << xmax << "  bin:" << bin_floor << " - " << bin_ceil << std::endl;
	 }

      } 
    }

    // generate profile
    vector<double> disp=vector<double>(nbin,0.0);
    
    for(int i=0;i<nbin;i++){
      if(mass_vec[i]!=0.0){
	disp[i]=(disp_vec[i]/mass_vec[i]);
      } else {
	disp[i]=0.0;
      }
    }
    // write profile
    ofstream outfile(outfilename.c_str());
    if(grad){
       for(int i=1;i<nbin;i++){
	 outfile << min+(double(i))*binsize << " " << (disp[i-1]-disp[i])/binsize << endl;
      }
    } else {
      for(int i=0;i<nbin;i++){
	outfile << min+(double(i)+0.5)*binsize << " " << disp[i] << endl;
      }
    }

    outfile.close();
  }
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
  read snapshot and write 3D grid VTK file of x-displacement

  \param infilename name of the input file (the *_0.txt)
  \param outfilename name of the output file
  \param xmin minimum of the x-range
  \param xmax maximum of the x-range
  \param ymin minimum of the y-range
  \param ymax maximum of the y-range
  \param zmin minimum of the z-range
  \param zmax maximum of the z-range
  \param cellsize the size of a grid cell (1 dimension, cells are cubes)
  \param grad displacement (gdim=false) or strain (grad=true)
  \param udim component of displacement
  \param gdim direction of derivative
*/
void read_and_write_disp_grid(const string& infilename,const string& outfilename,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,double cellsize,bool grad,int udim,int gdim)
{
  const double negvalue=-1e5; // large negative value marking "no data" boxes 

  int version=get_version(infilename);
  if(version<2) {
    cerr << "snapshot version (" << version << ") < 2 - can't calculate displacements" << endl;
  } else {

    vector<string> filenames=get_filenames(infilename,version);
    int xcell=int(ceil((xmax-xmin)/cellsize));
    double xcelldim=(xmax-xmin)/double(xcell);
    int ycell=int(ceil((ymax-ymin)/cellsize));
    double ycelldim=(ymax-ymin)/double(ycell);
    int zcell=int(ceil((zmax-zmin)/cellsize));
    double zcelldim=(zmax-zmin)/double(zcell);
    int ncell=xcell*ycell*zcell;

    std::cout << "grid dim: " << xcell << " , " << ycell << " , " << zcell << "   total cells : " << ncell << std::endl; 

    vector<double> mass_vec=vector<double>(ncell,0.0);
    vector<double> disp_vec=vector<double>(ncell,0.0);
    
    for(vector<string>::iterator iter=filenames.begin();
	iter!=filenames.end();
	iter++){
      cout << *iter << endl;
      ifstream datafile(iter->c_str());
      Vec3 pos,initpos,oldpos,vel,force,angvel,circ_shift;
      double rad,mass,q1,q2,q3,q4;
      int id,tag,npart;
      
      datafile >> npart;
      for(int i=0;i<npart;i++){
	int xbin_idx;
	int ybin_idx;
	int zbin_idx;

	datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force >> circ_shift
		 >> q1 >> q2 >> q3 >> q4 >> angvel;
	double xdisp=pos.X()-(circ_shift.X()+initpos.X());
	xbin_idx=int(floor(pos.X()-xmin)/xcelldim); 
	ybin_idx=int(floor(pos.Y()-ymin)/ycelldim); 
	zbin_idx=int(floor(pos.Z()-zmin)/zcelldim); 
       
	// calc index
	if((xbin_idx>=0)&&(xbin_idx<xcell) && (ybin_idx>=0)&&(ybin_idx<ycell) && (zbin_idx>=0)&&(zbin_idx<zcell) ){
	  int grid_idx=xbin_idx*ycell*zcell+ybin_idx*zcell+zbin_idx;
	  mass_vec[grid_idx]+=mass;
	  disp_vec[grid_idx]+=mass*xdisp;
	}
      }
      datafile.close();
    }
    
    // generate profile
    vector<double> disp=vector<double>(ncell,0.0);;
    for(int i=0;i<ncell;i++){
      if(mass_vec[i]!=0.0){
	disp[i]=disp_vec[i]/mass_vec[i];
      } else {
	disp[i]=negvalue;
      }
    }

    // write data
    ofstream outfile(outfilename.c_str());
 
    // VTK header
    if(grad){
      switch(gdim){
	
      case 1:  // dx
	write_vtk_header(outfile,xcell-1,ycell,zcell,xmin+0.5*xcelldim,xcelldim,ymin,ycelldim,zmin,zcelldim);

	for(int iz=0;iz<zcell;iz++){
	  for(int iy=0;iy<ycell;iy++){
	    for(int ix=1;ix<xcell;ix++){
	      int grid_idx=ix*ycell*zcell+iy*zcell+iz;
	      int grid_idx2=(ix-1)*ycell*zcell+iy*zcell+iz;
	      if((disp[grid_idx]==negvalue) || (disp[grid_idx2]==negvalue)) {
		outfile << negvalue << " " ;
	      } else {
		outfile << (disp[grid_idx]- disp[grid_idx2])/xcelldim<< " ";
	      }
	    }
	  }
	  outfile << std::endl;
	} break;
      case 2:  // dy
	write_vtk_header(outfile,xcell,ycell-1,zcell,xmin,xcelldim,(ymin+ycelldim*0.5),ycelldim,zmin,zcelldim);

	for(int iz=0;iz<zcell;iz++){
	  for(int iy=1;iy<ycell;iy++){
	    for(int ix=0;ix<xcell;ix++){
	      int grid_idx=ix*ycell*zcell+iy*zcell+iz;
	      int grid_idx2=ix*ycell*zcell+(iy-1)*zcell+iz;
	         if((disp[grid_idx]==negvalue) || (disp[grid_idx2]==negvalue)) {
		outfile << negvalue << " " ;
	      } else {
		   outfile << (disp[grid_idx]- disp[grid_idx2])/ycelldim<< " ";
		 }
	    }
	  }
	  outfile << std::endl;
	} break;
      case 3:  // dz
	write_vtk_header(outfile,xcell,ycell,zcell-1,xmin,xcelldim,ymin,ycelldim,zmin+zcelldim*0.5,zcelldim);

	for(int iz=1;iz<zcell;iz++){
	  for(int iy=0;iy<ycell;iy++){
	    for(int ix=0;ix<xcell;ix++){
	      int grid_idx=ix*ycell*zcell+iy*zcell+iz;
	      int grid_idx2=ix*ycell*zcell+iy*zcell+iz-1;
	         if((disp[grid_idx]==negvalue) || (disp[grid_idx2]==negvalue)) {
		   outfile << negvalue << " " ;
		 } else {
		   outfile << (disp[grid_idx]- disp[grid_idx2])/zcelldim<< " ";
		 }
	    }
	  }
	  outfile << std::endl;
	} break;
	
      }
    } else {
      write_vtk_header(outfile,xcell,ycell,zcell,xmin,xcelldim,ymin,ycelldim,zmin,zcelldim);

      for(int iz=0;iz<zcell;iz++){
	for(int iy=0;iy<ycell;iy++){
	  for(int ix=0;ix<xcell;ix++){
	    int grid_idx=ix*ycell*zcell+iy*zcell+iz;
	    outfile << disp[grid_idx] << " ";
	  }
	}
	outfile << std::endl;
      }
    }
    outfile << std::endl;
    outfile.close();
  }
}


/*!
  read snapshot and write 3D grid VTK file of porosity

  \param infilename name of the input file (the *_0.txt)
  \param outfilename name of the output file
  \param xmin minimum of the x-range
  \param xmax maximum of the x-range
  \param ymin minimum of the y-range
  \param ymax maximum of the y-range
  \param zmin minimum of the z-range
  \param zmax maximum of the z-range
  \param cellsize the size of a grid cell (1 dimension, cells are cubes)
*/
void read_and_write_poros_grid(const string& infilename,const string& outfilename,double xmin,double xmax,double ymin,double ymax,double zmin,double zmax,double cellsize)
{
  int version=get_version(infilename);
  if(version<2) {
    cerr << "snapshot version (" << version << ") < 2 - can't calculate displacements" << endl;
  } else {

    vector<string> filenames=get_filenames(infilename,version);
    int xcell=int(ceil((xmax-xmin)/cellsize));
    double xcelldim=(xmax-xmin)/double(xcell);
    int ycell=int(ceil((ymax-ymin)/cellsize));
    double ycelldim=(ymax-ymin)/double(ycell);
    int zcell=int(ceil((zmax-zmin)/cellsize));
    double zcelldim=(zmax-zmin)/double(zcell);
    int ncell=xcell*ycell*zcell;

    std::cout << "grid dim: " << xcell << " , " << ycell << " , " << zcell << "   total cells : " << ncell << std::endl; 

    vector<double> mass_vec=vector<double>(ncell,0.0);
    
    for(vector<string>::iterator iter=filenames.begin();
	iter!=filenames.end();
	iter++){
      cout << *iter << endl;
      ifstream datafile(iter->c_str());
      Vec3 pos,initpos,oldpos,vel,force,angvel,circ_shift;
      double rad,mass,q1,q2,q3,q4;
      int id,tag,npart;
      
      datafile >> npart;
      for(int i=0;i<npart;i++){
	int xbin_idx;
	int ybin_idx;
	int zbin_idx;

	datafile >> pos  >> rad >> id >> tag >> mass >> initpos >> oldpos >> vel >> force >> circ_shift
		 >> q1 >> q2 >> q3 >> q4 >> angvel;
	xbin_idx=int(floor(pos.X()-xmin)/xcelldim); 
	ybin_idx=int(floor(pos.Y()-ymin)/ycelldim); 
	zbin_idx=int(floor(pos.Z()-zmin)/zcelldim); 
       
	// calc index
	if((xbin_idx>=0)&&(xbin_idx<xcell) && (ybin_idx>=0)&&(ybin_idx<ycell) && (zbin_idx>=0)&&(zbin_idx<zcell) ){
	  int grid_idx=xbin_idx*ycell*zcell+ybin_idx*zcell+zbin_idx;
	  mass_vec[grid_idx]+=rad*rad*rad*4.1887; // (4/3*pi)
	}
      }
      datafile.close();
    }
    
    // write data
    ofstream outfile(outfilename.c_str());
 
    // VTK header
    write_vtk_header(outfile,xcell,ycell,zcell,xmin,xcelldim,ymin,ycelldim,zmin,zcelldim);

    double cellvol=xcelldim*ycelldim*zcelldim;
    for(int iz=0;iz<zcell;iz++){
      for(int iy=0;iy<ycell;iy++){
	for(int ix=0;ix<xcell;ix++){
	  int grid_idx=ix*ycell*zcell+iy*zcell+iz;
	  outfile << 1.0-(mass_vec[grid_idx]/cellvol) << " ";
	}
      }
      outfile << std::endl;
      }
    outfile << std::endl;
    outfile.close();
  }
}
