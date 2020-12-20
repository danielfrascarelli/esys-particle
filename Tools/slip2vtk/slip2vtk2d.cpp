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

#include "slip2vtk2d.h"

// --- STL includes ---
#include <sstream>
#include <map>
#include <utility>
#include <vector>

using std::ostringstream;
using std::map;
using std::make_pair;
using std::vector;

// --- I/O includes ---
#include <fstream>

using std::ofstream;
using std::ifstream;

// --- project includes ---
#include "Foundation/vec3.h"
 
/*!
  write VTK header

  \param outfile the output file
  \param nx
  \param nt
  \param x0
  \param dx
  \param dt
*/
void write_vtk_header(ofstream&  outfile, int nx, int nt, double x0,double dx, double dt)
{ 
  outfile << "# vtk DataFile Version 2.0"  << std::endl;
  outfile << "slip data"  << std::endl;
  outfile << "ASCII" << std::endl;
  outfile << "DATASET RECTILINEAR_GRID" << std::endl;
  outfile << "DIMENSIONS " << nx  << " " << nt << "  1 " << std::endl;
  outfile << "X_COORDINATES " << nx << " float"  << std::endl;
  for(int i=0;i<nx;i++){
    outfile << x0+double(i)*dx << " "; 
  }
  outfile << std::endl;
  outfile << "Y_COORDINATES " << nt << " float" << std::endl;
  for(int i=0;i<nt;i++){
    outfile << i*dt << " ";
  }
  outfile << std::endl;
  outfile << "Z_COORDINATES " << " 1  float "  << std::endl << "0" << std::endl;
  outfile << "POINT_DATA " <<  nt*nx << std::endl;
  outfile << "SCALARS displacement float" << std::endl;
  outfile << "LOOKUP_TABLE default" << std::endl;
}

/*!
  read a displacement file into a map

  \param filename
*/
map<double,double> read_file_to_map(const string& filename)
{
  map<double,double> res;

  // open input file
  ifstream infile(filename.c_str()); 
  // read data 
  while(!infile.eof()){
    Vec3 pos,disp;
    infile >> pos >> disp;  
    if(!infile.fail()){
      res.insert(make_pair(pos.X(),disp.X()));
     }
    }
  // close input file
  infile.close();

  return res;
}

/*!
  read line of a RAW_SERIES file into map

  \param infile the input file
  \pos_vec vector of particle x-positions 
*/
map<double,double> read_line_to_map(ifstream& infile, const vector<double>& pos_vec)
{
  map<double,double> res;

  // read data 
  int ndata=pos_vec.size();
  for(int i=0;i<ndata;i++){
    Vec3 disp;
    infile >> disp;
    res.insert(make_pair(pos_vec[i],disp.X()));
  }
  
  return res;
}

/*!
  get displacement at given point from map (linear interpolation)
  
  \param x 
  \param dispmap
*/
double get_disp_from_map_linear(const map<double,double>& dispmap,double x)
{
  double disp;
  map<double,double>::const_iterator it1=dispmap.upper_bound(x); // the one "above" (x>px)
  // check if there is one "below" (x<px)
  if(it1!=dispmap.begin()){
    map<double,double>::const_iterator it2=it1;
    it2--;
    disp=it2->second+((x-it2->first)/(it1->first-it2->first))*(it1->second-it2->second);
  } else {
    disp=it1->second;
  }

  return disp;
}

/*!
  convert a number of file pairs containing particle displacements above/below a 1D fault
  in a 2D medium to a vtk-file (rect. grid) of the relative displacement

  \param infilebase_up base file name for data above the fault 
  \param infilebase_down base file name for data below the fault 
  \param outfilename file name for the putput file
  \param t0
  \param nt
  \param dt
  \param tlen length of the output data in t-direction 
  \param x0 minimum x-value
  \param x1 maximum x-value
  \param nx nr. of grid points in x
*/
void slip2vtk2d(const string& infilebase_up,
		const string& infilebase_down,
		const string& outfilename,
		int t0,
		int nt,
		int dt,
		double tlen,
		double x0,
		double x1,
		int nx)
{
  std::cout << "slip2vtk2d"  << std::endl;
  map<double,double> upper_map,lower_map;
  ostringstream infilename_up,infilename_down;
  vector<double> disp_old(nx);
  // open output file
  ofstream outfile(outfilename.c_str());
  // write header
  write_vtk_header(outfile,nx,((nt-t0)/dt)-1 ,x0,(x1-x0)/double(nx),tlen/(double(nt)));
  // first ts -> fill up disp_old
  infilename_up << infilebase_up << "." << t0 << ".dat";  
  infilename_down << infilebase_down << "." << t0 << ".dat";  
  upper_map=read_file_to_map(infilename_up.str());
  lower_map=read_file_to_map(infilename_down.str());
  double dx=(x1-x0)/double(nx);
  for(int i=0;i<nx;i++){
    double px=x0+double(i)*dx; // current sample point
    double disp_up,disp_down;
    disp_up=get_disp_from_map_linear(upper_map,px);
    disp_down=get_disp_from_map_linear(lower_map,px);
    disp_old[i]=disp_up-disp_down;
  }
  // clean up
  lower_map.clear();
  upper_map.clear();
  infilename_up.seekp(0,ios_base::beg);
  infilename_down.seekp(0,ios_base::beg);
  
  // ---- further steps ----
  for(int t=t0+1;t<nt;t+=dt){
    
    // generate input filenames
    infilename_up << infilebase_up << "." << t0+t*dt << ".dat";  
    infilename_down << infilebase_down << "." << t0+t*dt << ".dat";  
    std::cout << infilename_up.str()  << std::endl;
    upper_map=read_file_to_map(infilename_up.str());
    lower_map=read_file_to_map(infilename_down.str());
    // --- calc 
    // get x-step
    double dx=(x1-x0)/double(nx);
    for(int i=0;i<nx;i++){
      double px=x0+double(i)*dx; // current sample point
      double disp_up,disp_down;
      disp_up=get_disp_from_map_linear(upper_map,px);
      disp_down=get_disp_from_map_linear(lower_map,px);
      double disp=disp_up-disp_down;
      outfile << disp-disp_old[i] << " ";
    }
    outfile << std::endl;
    // clean up
    lower_map.clear();
    upper_map.clear();
    infilename_up.seekp(0,ios_base::beg);
    infilename_down.seekp(0,ios_base::beg);
  }  
  // close output file
  outfile.close();
}

/*!
  convert a number of file pairs containing particle displacements above/below a 1D fault
  in a 2D medium to a vtk-file (rect. grid) of the relative displacement rate

  \param infilebase_up base file name for data above the fault 
  \param infilebase_down base file name for data below the fault 
  \param outfilename file name for the putput file
  \param t0
  \param nt
  \param dt
  \param tlen length of the output data in t-direction 
  \param x0 minimum x-value
  \param x1 maximum x-value
  \param nx nr. of grid points in x
*/
void slip2vtk2d_rate(const string& infilebase_up,
		const string& infilebase_down,
		const string& outfilename,
		int t0,
		int nt,
		int dt,
		double tlen,
		double x0,
		double x1,
		int nx)
{
  ostringstream infilename_up,infilename_down;
  map<double,double> upper_map,lower_map;
  vector<double> disp_old(nx);
 // open output file
  ofstream outfile(outfilename.c_str());
  // write header
  write_vtk_header(outfile,nx,((nt-t0)/dt)-1 ,x0,(x1-x0)/double(nx),tlen/(double(nt)));
  // first ts -> fill up disp_old
  infilename_up << infilebase_up << "." << t0 << ".dat";  
  infilename_down << infilebase_down << "." << t0 << ".dat";  
  upper_map=read_file_to_map(infilename_up.str());
  lower_map=read_file_to_map(infilename_down.str());
  double dx=(x1-x0)/double(nx);
  for(int i=0;i<nx;i++){
    double px=x0+double(i)*dx; // current sample point
    double disp_up,disp_down;
    disp_up=get_disp_from_map_linear(upper_map,px);
    disp_down=get_disp_from_map_linear(lower_map,px);
    disp_old[i]=disp_up-disp_down;
  }
  // clean up
  lower_map.clear();
  upper_map.clear();
  infilename_up.seekp(0,ios_base::beg);
  infilename_down.seekp(0,ios_base::beg);
  
  // ---- further steps ----
  for(int t=t0+1;t<nt;t+=dt){
    
    // generate input filenames
    infilename_up << infilebase_up << "." << t0+t*dt << ".dat";  
    infilename_down << infilebase_down << "." << t0+t*dt << ".dat";  
    upper_map=read_file_to_map(infilename_up.str());
    lower_map=read_file_to_map(infilename_down.str());
    // --- calc 
    // get x-step
    double dx=(x1-x0)/double(nx);
    for(int i=0;i<nx;i++){
      double px=x0+double(i)*dx; // current sample point
      double disp_up,disp_down;
      disp_up=get_disp_from_map_linear(upper_map,px);
      disp_down=get_disp_from_map_linear(lower_map,px);
      double disp=disp_up-disp_down;
      outfile << disp-disp_old[i] << " ";
      disp_old[i]=disp;
    }
    outfile << std::endl;
    // clean up
    lower_map.clear();
    upper_map.clear();
    infilename_up.seekp(0,ios_base::beg);
    infilename_down.seekp(0,ios_base::beg);
  }  
  // close output file
  outfile.close();
}

/*! 
  convert a pair of files containing particle displacements above/below a 1D fault
  and a pair of files containing the initial particle positions to a vtk-file (rect. grid) of 
  the relative displacement rate

  \param infile_up  file name for displacement data above the fault 
  \param infile_down base file name for displacement data below the fault 
  \param posfile_up  file name for position data above the fault 
  \param posfile_down base file name for position data below the fault 
  \param outfilename file name for the putput file
  \param t0
  \param nt
  \param dt
  \param tlen length of the output data in t-direction 
  \param x0 minimum x-value
  \param x1 maximum x-value
  \param nx nr. of grid points in x
*/
void slip2vtk2d_rate_rs(const string& infile_up, 
		      const string& infile_down, 
		      const string& posfile_up,
		      const string& posfile_down,  
		      const string& outfilename, 
		      int t0,int nt,int dt,double tlen,double x0,double x1,int nx)
{
  cout << "slip2vtk2d_rate_rs, displacement from " <<  infile_up << " , " <<  infile_down << " pos. from  "  << posfile_up << " , " << posfile_down << endl;

  // --- get position vector, nr of. particles
  vector<double> pos_up;
  vector<double> pos_down;
  int np_up,np_down;

  // open "up" pos file
  ifstream pos_up_file(posfile_up.c_str());
  // read data
  while(!pos_up_file.eof()){
    Vec3 v;
    pos_up_file >> v;
    if(!pos_up_file.eof()){
      pos_up.push_back(v.X());
    }
  }
  np_up=pos_up.size();
  // close file 
  pos_up_file.close();
  // open "down" pos file
  ifstream pos_down_file(posfile_down.c_str());
  // read data
  while(!pos_down_file.eof()){
    Vec3 v;
    pos_down_file >> v;
    if(!pos_down_file.eof()){
      pos_down.push_back(v.X());
    }
  }
  np_down=pos_down.size();
  // close file 
  pos_down_file.close();
  cout << "found " << np_up << " particles up and " << np_down << " particles down" << endl;
  // --- read displacement data , output slip --
  map<double,double> upper_map,lower_map;
  // open disp files
  ifstream disp_up_file(infile_up.c_str());
  ifstream disp_down_file(infile_down.c_str());
 // open output file
  ofstream outfile(outfilename.c_str());
  // write header
  write_vtk_header(outfile,nx,((nt-t0)/dt)-1 ,x0,(x1-x0)/double(nx),tlen/(double(nt)));
  // first ts -> fill up disp_old
  vector<double> disp_old(nx);
  for(int t=0;t<t0;t++){
    upper_map=read_line_to_map(disp_up_file,pos_up);
    lower_map=read_line_to_map(disp_down_file,pos_down);
    double dx=(x1-x0)/double(nx);
    for(int i=0;i<nx;i++){
      double px=x0+double(i)*dx; // current sample point
      double disp_up,disp_down;
      disp_up=get_disp_from_map_linear(upper_map,px);
      disp_down=get_disp_from_map_linear(lower_map,px);
      disp_old[i]=disp_up-disp_down;
    }
    // clean up
    lower_map.clear();
    upper_map.clear();
  }
  // ---- further steps ----
  for(int t=t0+1;t<nt;t++){
   upper_map=read_line_to_map(disp_up_file,pos_up);
   lower_map=read_line_to_map(disp_down_file,pos_down);
    // --- calc 
    // get x-step
   if(t%dt==0){
     double dx=(x1-x0)/double(nx);
     for(int i=0;i<nx;i++){
       double px=x0+double(i)*dx; // current sample point
       double disp_up,disp_down;
       disp_up=get_disp_from_map_linear(upper_map,px);
       disp_down=get_disp_from_map_linear(lower_map,px);
       double disp=disp_up-disp_down;
       outfile << disp-disp_old[i] << " ";
       disp_old[i]=disp;
     }
     outfile << std::endl;
   }
   // clean up
    lower_map.clear();
    upper_map.clear(); 
  }
  // close disp files
  disp_up_file.close();
  disp_down_file.close();
  outfile.close();
}
/*! 
  convert a pair of files containing particle displacements above/below a 1D fault
  and a pair of files containing the initial particle positions to a raw file (column per 
 data point) of  the relative displacement 

  \param infile_up  file name for displacement data above the fault 
  \param infile_down base file name for displacement data below the fault 
  \param posfile_up  file name for position data above the fault 
  \param posfile_down base file name for position data below the fault 
  \param outfilename file name for the putput file
  \param t0
  \param nt
  \param dt
  \param x0 minimum x-value
  \param x1 maximum x-value
  \param nx nr. of grid points in x
  \param shift shift output data by given abount (circular)
  \param mult interpolation steps between data points
*/
void slip2raw2d_rs(const string& infile_up, 
		      const string& infile_down, 
		      const string& posfile_up,
		      const string& posfile_down,  
		      const string& outfilename, 
		   int t0,int nt,int dt,double x0,double x1,int nx,int shift,int mult)
{
  cout << "slip2raw2d_rs, displacement from " <<  infile_up << " , " <<  infile_down << " pos. from  "  << posfile_up << " , " << posfile_down << " shift:" << shift << " mult: " << mult << endl;

  // --- get position vector, nr of. particles
  vector<double> pos_up;
  vector<double> pos_down;
  int np_up,np_down;

  // open "up" pos file
  ifstream pos_up_file(posfile_up.c_str());
  // read data
  while(!pos_up_file.eof()){
    Vec3 v;
    pos_up_file >> v;
    if(!pos_up_file.eof()){
      pos_up.push_back(v.X());
    }
  }
  np_up=pos_up.size();
  // close file 
  pos_up_file.close();
  // open "down" pos file
  ifstream pos_down_file(posfile_down.c_str());
  // read data
  while(!pos_down_file.eof()){
    Vec3 v;
    pos_down_file >> v;
    if(!pos_down_file.eof()){
      pos_down.push_back(v.X());
    }
  }
  np_down=pos_down.size();
  // close file 
  pos_down_file.close();
  cout << "found " << np_up << " particles up and " << np_down << " particles down" << endl;
  // --- read displacement data , output slip --
  map<double,double> upper_map,lower_map;
  // open disp files
  ifstream disp_up_file(infile_up.c_str());
  ifstream disp_down_file(infile_down.c_str());
 // open output file
  ofstream outfile(outfilename.c_str());
  // first ts -> fill up disp_init
  vector<double> disp_init(nx);
  vector<double> disp_old(nx);
  vector<double> disp_new(nx);
  for(int t=0;t<=t0;t++){
    upper_map=read_line_to_map(disp_up_file,pos_up);
    lower_map=read_line_to_map(disp_down_file,pos_down);
    double dx=(x1-x0)/double(nx);
    for(int i=0;i<nx;i++){
      int ix=(i+shift)%nx;
      double px=x0+double(i)*dx; // current sample point
      double disp_up,disp_down;
      disp_up=get_disp_from_map_linear(upper_map,px);
      disp_down=get_disp_from_map_linear(lower_map,px);
      disp_init[ix]=disp_up-disp_down;
      disp_new[ix]=disp_init[ix];
    }
  }
  // clean up
  lower_map.clear();
  upper_map.clear();
  // ---- further steps ----
  for(int t=t0+1;t<nt;t++){
   upper_map=read_line_to_map(disp_up_file,pos_up);
   lower_map=read_line_to_map(disp_down_file,pos_down);
    // --- calc 
    // get x-step
   if(t%dt==0){
     double dx=(x1-x0)/double(nx);
     for(int i=0;i<nx;i++){
       int ix=(i+shift)%nx;
       double px=x0+double(i)*dx; // current sample point
       double disp_up,disp_down;
       disp_up=get_disp_from_map_linear(upper_map,px);
       disp_down=get_disp_from_map_linear(lower_map,px);
       disp_old[ix]=disp_new[ix];
       disp_new[ix]=disp_up-disp_down;
    }
     // write to file
     for(int im=1;im<=mult;im++){
       for(int i=0;i<nx;i++){
	 double disp=disp_old[i]+(disp_new[i]-disp_old[i])*(double(im)/double(mult));
	 outfile << disp-disp_init[i] << " ";	 
       }
       outfile << std::endl;
     }
     // clean up
    lower_map.clear();
    upper_map.clear(); 
   }
  }
  // close disp files
  disp_up_file.close();
  disp_down_file.close();
  outfile.close();
}

/*! 
  convert a pair of files containing particle displacements above/below a 1D fault
  and a pair of files containing the initial particle positions to a raw file (column per 
   data point) of the relative velocity

  \param infile_up  file name for displacement data above the fault 
  \param infile_down base file name for displacement data below the fault 
  \param posfile_up  file name for position data above the fault 
  \param posfile_down base file name for position data below the fault 
  \param outfilename file name for the putput file
  \param t0
  \param nt
  \param dt
  \param x0 minimum x-value
  \param x1 maximum x-value
  \param nx nr. of grid points in x
*/
void slip2raw2d_rate_rs(const string& infile_up, 
			const string& infile_down, 
			const string& posfile_up,
			const string& posfile_down,  
			const string& outfilename, 
			int t0,int nt,int dt,double x0,double x1,int nx)
{
  cout << "slip2raw2d_rs, displacement from " <<  infile_up << " , " <<  infile_down << " pos. from  "  << posfile_up << " , " << posfile_down << endl;

  // --- get position vector, nr of. particles
  vector<double> pos_up;
  vector<double> pos_down;
  int np_up,np_down;

  // open "up" pos file
  ifstream pos_up_file(posfile_up.c_str());
  // read data
  while(!pos_up_file.eof()){
    Vec3 v;
    pos_up_file >> v;
    if(!pos_up_file.eof()){
      pos_up.push_back(v.X());
    }
  }
  np_up=pos_up.size();
  // close file 
  pos_up_file.close();
  // open "down" pos file
  ifstream pos_down_file(posfile_down.c_str());
  // read data
  while(!pos_down_file.eof()){
    Vec3 v;
    pos_down_file >> v;
    if(!pos_down_file.eof()){
      pos_down.push_back(v.X());
    }
  }
  np_down=pos_down.size();
  // close file 
  pos_down_file.close();
  cout << "found " << np_up << " particles up and " << np_down << " particles down" << endl;
  // --- read displacement data , output slip --
  map<double,double> upper_map,lower_map;
  // open disp files
  ifstream disp_up_file(infile_up.c_str());
  ifstream disp_down_file(infile_down.c_str());
 // open output file
  ofstream outfile(outfilename.c_str());
  // first ts -> fill up disp_old
  vector<double> disp_old(nx);
  for(int t=0;t<=t0;t++){
    upper_map=read_line_to_map(disp_up_file,pos_up);
    lower_map=read_line_to_map(disp_down_file,pos_down);
    double dx=(x1-x0)/double(nx);
    for(int i=0;i<nx;i++){
      double px=x0+double(i)*dx; // current sample point
      double disp_up,disp_down;
      disp_up=get_disp_from_map_linear(upper_map,px);
      disp_down=get_disp_from_map_linear(lower_map,px);
      disp_old[i]=disp_up-disp_down;
    }
    // clean up
    lower_map.clear();
    upper_map.clear();
  }
  // ---- further steps ----
  for(int t=t0+1;t<nt;t++){
   upper_map=read_line_to_map(disp_up_file,pos_up);
   lower_map=read_line_to_map(disp_down_file,pos_down);
    // --- calc 
    // get x-step
   if(t%dt==0){
     double dx=(x1-x0)/double(nx);
     for(int i=0;i<nx;i++){
       double px=x0+double(i)*dx; // current sample point
       double disp_up,disp_down;
       disp_up=get_disp_from_map_linear(upper_map,px);
       disp_down=get_disp_from_map_linear(lower_map,px);
       double disp=disp_up-disp_down;
       outfile << disp-disp_old[i] << " ";
       disp_old[i]=disp;
     }
     outfile << std::endl;
   }
   // clean up
    lower_map.clear();
    upper_map.clear(); 
  }
  // close disp files
  disp_up_file.close();
  disp_down_file.close();
  outfile.close();
}


/*!
  convert a number of file pairs containing particle displacements above/below a 1D fault
  in a 2D medium to a file which represents the time evolution of the displacement at a given
  point along the fault

  \param  infilebase_up base file name for data above the fault 
  \param infilebase_down base file name for data below the fault 
  \param outfilename file name for the putput file
  \param t0
  \param nt
  \param dt
  \param x0 location / sample point
*/
void slip_x2slip_t2d( const string& infilebase_up,
		      const string& infilebase_down,
		      const string& outfilename,
		      int t0,int nt ,int dt,double x0)
{
  ostringstream infilename_up,infilename_down;
  map<double,double> upper_map,lower_map;
 // open output file
  ofstream outfile(outfilename.c_str());
 // ---- further steps ----
  for(int t=t0+1;t<nt;t+=dt){
    
    // generate input filenames
    infilename_up << infilebase_up << "." << t0+t*dt << ".dat";  
    infilename_down << infilebase_down << "." << t0+t*dt << ".dat";  
    std::cout << infilename_up.str()  << std::endl;
    upper_map=read_file_to_map(infilename_up.str());
    lower_map=read_file_to_map(infilename_down.str());
    double disp_up=get_disp_from_map_linear(upper_map,x0);
    double disp_down=get_disp_from_map_linear(lower_map,x0);
    double disp=disp_up-disp_down;
    outfile << disp << std::endl;
    // clean up
    lower_map.clear();
    upper_map.clear();
    infilename_up.seekp(0,ios_base::beg);
    infilename_down.seekp(0,ios_base::beg);
  }  
  // close output file
  outfile.close();
}

/*!
  get the time for each trace when the displacement first exceeds a given threshold
*/
void slip2rf(const string& infilebase_up,
	     const string& infilebase_down,
	     const string& outfilename,
	     int t0,
	     int nt,
	     int dt,
	     double x0,
	     double x1,
	     int nx,
	     double thr)
{
  map<double,double> upper_map,lower_map;
  ostringstream infilename_up,infilename_down;
  vector<double> disp_old(nx);
  vector<double> disp_field(nx*(nt-t0)/dt);

  // first ts -> fill up disp_old
  infilename_up << infilebase_up << "." << t0 << ".dat";  
  infilename_down << infilebase_down << "." << t0 << ".dat";  
  upper_map=read_file_to_map(infilename_up.str());
  lower_map=read_file_to_map(infilename_down.str());
  double dx=(x1-x0)/double(nx);
  for(int i=0;i<nx;i++){
    double px=x0+double(i)*dx; // current sample point
    double disp_up,disp_down;
    disp_up=get_disp_from_map_linear(upper_map,px);
    disp_down=get_disp_from_map_linear(lower_map,px);
    disp_old[i]=disp_up-disp_down;
  }
  // clean up
  lower_map.clear();
  upper_map.clear();
  infilename_up.seekp(0,ios_base::beg);
  infilename_down.seekp(0,ios_base::beg);
  
  // ---- further steps ----
  for(int t=t0+1;t<nt;t+=dt){
    int it=(t-(t0+1))/dt;
    // generate input filenames
    infilename_up << infilebase_up << "." << t0+t<< ".dat";  
    infilename_down << infilebase_down << "." << t0+t << ".dat";  
    std::cout << infilename_up.str()  << std::endl;
    upper_map=read_file_to_map(infilename_up.str());
    lower_map=read_file_to_map(infilename_down.str());
    // --- calc 
    // get x-step
    double dx=(x1-x0)/double(nx);
    for(int i=0;i<nx;i++){
      double px=x0+double(i)*dx; // current sample point
      double disp_up,disp_down;
      disp_up=get_disp_from_map_linear(upper_map,px);
      disp_down=get_disp_from_map_linear(lower_map,px);
      double disp=disp_up-disp_down;
      disp_field[it*nx+i]=disp;
    }
    // clean up
    lower_map.clear();
    upper_map.clear();
    infilename_up.seekp(0,ios_base::beg);
    infilename_down.seekp(0,ios_base::beg);
  }
  // for each x, find  min y where disp > threshold
 // open output file
  ofstream outfile(outfilename.c_str());

  for(int i=0;i<nx;i++){
    int iy=1;
    bool found=false;
    while((iy<(nt-t0)/dt) && (!found)){
      //outfile << disp_field[iy*nx+i]-disp_field[(iy-1)*nx+i] << " ";
      found=(disp_field[iy*nx+i]-disp_field[(iy-1)*nx+i])>thr;
      iy++;
    }
  //outfile << endl;
    if(found){
      outfile << i << "  " << iy << endl; 
    }
  }

  // close output file
  outfile.close();
}
/*!
  get the time for each trace when the displacement first exceeds a given threshold 
  using velocity files as input
*/
void vel2rf(const string& infilename_up,
	     const string& infilename_down,
	     const string& outfilename,
	     int t0,
	     int nt,
	     int dt,
	    int ofs, 
	     double x0,
	     double x1,
	     int nx,
	    double px0,
	    double px1,
	    int npx,
	     double thr)
{
  map<double,double> upper_map,lower_map;
  vector<double> disp_old(nx);
  vector<double> disp_field(nx*(nt-t0)/dt);
  vector<double> posvec;

  // setup position vector
  double dpx=(px1-px0)/double(npx);
  for(int i=0;i<npx;i++){
    posvec.push_back(double(i)*dpx+px0);
  }
  // open files
  ifstream infile_up(infilename_up.c_str());
  ifstream infile_down(infilename_down.c_str());
  // first ts -> fill up disp_old
  // read stuff before start of delayed file
  if(ofs>0){
    //    std::cout << "ofs = " << ofs << std::endl;
    for(int i=0;i<ofs;i++){
      //      std::cout << "pre read " << i << std::endl; 
      upper_map=read_line_to_map(infile_up,posvec);  
      upper_map.clear();
    } 
  } else {
    //    std::cout << "ofs = " << ofs << std::endl;
    for(int i=0;i>ofs;i--){
      //      std::cout << "pre read " << i << std::endl; 
      lower_map=read_line_to_map(infile_down,posvec);
      lower_map.clear();
    }
  }  
  upper_map=read_line_to_map(infile_up,posvec);
  lower_map=read_line_to_map(infile_down,posvec);
  double dx=(x1-x0)/double(nx);
  for(int i=0;i<nx;i++){
    double px=x0+double(i)*dx; // current sample point
    double disp_up,disp_down;
    disp_up=get_disp_from_map_linear(upper_map,px);
    disp_down=get_disp_from_map_linear(lower_map,px);
    disp_old[i]=disp_up-disp_down;
    disp_field[i]=0.0;
  }
  // clean up
  lower_map.clear();
  upper_map.clear();
  
  // ---- further steps ----
  for(int t=t0+1;t<nt;t+=dt){
    int it=(t-(t0+1))/dt;
    //    std::cout << "read step " << t << "   [ " << it << " ]" << std::endl; 
    upper_map=read_line_to_map(infile_up,posvec);
    lower_map=read_line_to_map(infile_down,posvec);
    // --- calc 
    // get x-step
    double dx=(x1-x0)/double(nx);
    for(int i=0;i<nx;i++){
      double px=x0+double(i)*dx; // current sample point
      double disp_up,disp_down;
      disp_up=get_disp_from_map_linear(upper_map,px);
      disp_down=get_disp_from_map_linear(lower_map,px);
      double disp=disp_up-disp_down;
      disp_field[it*nx+i]=(disp_field[(it-1)*nx+i]+disp)-disp_old[i];
    }
    // clean up
    lower_map.clear();
    upper_map.clear();
 }
  infile_up.close();
  infile_down.close();
  // for each x, find  min y where disp > threshold
 // open output file
  ofstream outfile(outfilename.c_str());

  for(int i=0;i<nx;i++){
    int iy=1;
    bool found=false;
    //    std::cout << "calc rf for point " << i << std::endl;
    while((iy<(nt-t0)/dt) && (!found)){
      //      std::cout << "disp : " << disp_field[iy*nx+i] << std::endl;
      found=(disp_field[iy*nx+i])>thr;
      iy++;
    }
  //outfile << endl;
    if(found){
      outfile << i << "  " << iy << endl; 
    }
  }

  // close output file
  outfile.close();
}

/*!
  get total slip distribution at time t from posfiles & RAW_SERIES

  \param infile_up  file name for displacement data above the fault 
  \param infile_down base file name for displacement data below the fault 
  \param posfile_up  file name for position data above the fault 
  \param posfile_down base file name for position data below the fault 
  \param outfilename file name for the putput file
  \param t0 start time
  \param te end time 
  \param x0 minimum x-value
  \param x1 maximum x-value
  \param nx nr. of grid points in x
*/
void slip2d_total_rs(const string&infile_up, 
	                  const string& infile_down,
	                  const  string& posfile_up, 
	                  const string& posfile_down,
	                  const string& outfilename, 
		int t0, int te,double x0,double x1,int nx)
{
  cout << "slip2d_total_rs:  displacement from " <<  infile_up << " , " <<  infile_down << " pos. from  "  << posfile_up << " , " << posfile_down << endl;

  // --- get position vector, nr of. particles
  vector<double> pos_up;
  vector<double> pos_down;
  int np_up,np_down;

  // open "up" pos file
  ifstream pos_up_file(posfile_up.c_str());
  // read data
  while(!pos_up_file.eof()){
    Vec3 v;
    pos_up_file >> v;
    if(!pos_up_file.eof()){
      pos_up.push_back(v.X());
    }
  }
  np_up=pos_up.size();
  // close file 
  pos_up_file.close();
  // open "down" pos file
  ifstream pos_down_file(posfile_down.c_str());
  // read data
  while(!pos_down_file.eof()){
    Vec3 v;
    pos_down_file >> v;
    if(!pos_down_file.eof()){
      pos_down.push_back(v.X());
    }
  }
  np_down=pos_down.size();
  // close file 
  pos_down_file.close();
  cout << "found " << np_up << " particles up and " << np_down << " particles down" << endl;
  // --- read displacement data , output slip --
  map<double,double> upper_map,lower_map;
  // open disp files
  ifstream disp_up_file(infile_up.c_str());
  ifstream disp_down_file(infile_down.c_str());
 // open output file
  ofstream outfile(outfilename.c_str());
   // first ts -> fill up disp_old
  vector<double> disp_old(nx);
  upper_map=read_line_to_map(disp_up_file,pos_up);
  lower_map=read_line_to_map(disp_down_file,pos_down);
  double dx=(x1-x0)/double(nx);
  for(int i=0;i<nx;i++){
    double px=x0+double(i)*dx; // current sample point
    double disp_up,disp_down;
    disp_up=get_disp_from_map_linear(upper_map,px);
    disp_down=get_disp_from_map_linear(lower_map,px);
    disp_old[i]=disp_up-disp_down;
  }
  // clean up
  lower_map.clear();
  upper_map.clear();
  // ---- further steps ----
  for(int t=t0+1;t<te-1;t++){
    upper_map=read_line_to_map(disp_up_file,pos_up);
    lower_map=read_line_to_map(disp_down_file,pos_down);
   // clean up
    lower_map.clear();
    upper_map.clear(); 
  }
  upper_map=read_line_to_map(disp_up_file,pos_up);
  lower_map=read_line_to_map(disp_down_file,pos_down);
  double sum_disp=0.0;
  for(int i=0;i<nx;i++){
    double px=x0+double(i)*dx; // current sample point
    double disp_up,disp_down;
    disp_up=get_disp_from_map_linear(upper_map,px);
    disp_down=get_disp_from_map_linear(lower_map,px);
    double disp=disp_up-disp_down;
    sum_disp+=disp-disp_old[i];
    // cout << "disp at point: " << i << "  " << disp << " (up,down): " << disp_up << "  "  << disp_down << " old : " << disp_old[i] << endl;
    outfile << disp-disp_old[i] << std::endl;
  }
  cout << " integral slip : " << sum_disp/(double(nx)) << endl;
  // close disp files
  disp_up_file.close();
  disp_down_file.close();
  outfile.close();
}

/*! 
  calculate the moment rate function (actually the potency rate function) from
  a pair of files containing particle displacements above/below a 1D fault
  and a pair of files containing the initial particle positions

  \param infile_up  file name for displacement data above the fault 
  \param infile_down base file name for displacement data below the fault 
  \param posfile_up  file name for position data above the fault 
  \param posfile_down base file name for position data below the fault 
  \param outfilename file name for the putput file
  \param t0
  \param nt
  \param dt
  \param x0 minimum x-value
  \param x1 maximum x-value
  \param nx nr. of grid points in x
*/
void slip2momrate_rs(const string& infile_up, 
		     const string& infile_down, 
		     const string& posfile_up,
		     const string& posfile_down,  
		     const string& outfilename, 
		     int t0,int nt,int dt,double x0,double x1,int nx)
{
  cout << "slip2momrate_rs, displacement from " <<  infile_up << " , " <<  infile_down << " pos. from  "  << posfile_up << " , " << posfile_down << endl;

  // --- get position vector, nr of. particles
  vector<double> pos_up;
  vector<double> pos_down;
  int np_up,np_down;

  // open "up" pos file
  ifstream pos_up_file(posfile_up.c_str());
  // read data
  while(!pos_up_file.eof()){
    Vec3 v;
    pos_up_file >> v;
    if(!pos_up_file.eof()){
      pos_up.push_back(v.X());
    }
  }
  np_up=pos_up.size();
  // close file 
  pos_up_file.close();
  // open "down" pos file
  ifstream pos_down_file(posfile_down.c_str());
  // read data
  while(!pos_down_file.eof()){
    Vec3 v;
    pos_down_file >> v;
    if(!pos_down_file.eof()){
      pos_down.push_back(v.X());
    }
  }
  np_down=pos_down.size();
  // close file 
  pos_down_file.close();
  cout << "found " << np_up << " particles up and " << np_down << " particles down" << endl;


  // --- read displacement data , output slip --
  map<double,double> upper_map,lower_map;
  // open disp files
  ifstream disp_up_file(infile_up.c_str());
  ifstream disp_down_file(infile_down.c_str());
  // open output file
  ofstream outfile(outfilename.c_str());
  // first ts -> fill up disp_old
  vector<double> disp_old(nx);
  upper_map=read_line_to_map(disp_up_file,pos_up);
  lower_map=read_line_to_map(disp_down_file,pos_down);
  double dx=(x1-x0)/double(nx);
  for(int i=0;i<nx;i++){
    double px=x0+double(i)*dx; // current sample point
    double disp_up,disp_down;
    disp_up=get_disp_from_map_linear(upper_map,px);
    disp_down=get_disp_from_map_linear(lower_map,px);
    disp_old[i]=disp_up-disp_down;
  }
  // clean up
  lower_map.clear();
  upper_map.clear();
  // ---- further steps ---- 
  double sum_disp_old=0.0;
  for(int t=t0+dt;t<nt;t+=dt){
    for(int tt=t+1;tt<t+dt;tt++){
      //cout << "tt: " << tt << endl;
      upper_map=read_line_to_map(disp_up_file,pos_up);
      lower_map=read_line_to_map(disp_down_file,pos_down);
      // clean up
      lower_map.clear();
      upper_map.clear(); 
    }
    //cout << "t: " << t << endl;
    upper_map=read_line_to_map(disp_up_file,pos_up);
    lower_map=read_line_to_map(disp_down_file,pos_down);
    double sum_disp=0.0;
    for(int i=0;i<nx;i++){
      double px=x0+double(i)*dx; // current sample point
      double disp_up,disp_down;
      disp_up=get_disp_from_map_linear(upper_map,px);
      disp_down=get_disp_from_map_linear(lower_map,px);
      double disp=disp_up-disp_down;
      sum_disp+=disp-disp_old[i];
      // cout << "disp at point: " << i << "  " << disp << " (up,down): " << disp_up << "  "  << disp_down << " old : " << disp_old[i] << endl;
    }
    outfile << t << "  " << ((sum_disp-sum_disp_old)*((x1-x0)/double(nx)))/double(dt) << std::endl;
    sum_disp_old=sum_disp;
    // clean up
    lower_map.clear();
    upper_map.clear();
    

  }
  // disp files
  disp_up_file.close();
  disp_down_file.close();
  outfile.close();
 
}

/*!
  convert a number of file pairs containing particle displacements above/below a 1D fault
  in a 2D medium to a vtk-file (rect. grid) of the relative displacement

  \param infilebase_up base file name for data above the fault 
  \param infilebase_down base file name for data below the fault 
  \param outfilename file name for the putput file
  \param t0
  \param nt
  \param dt
  \param x0 minimum x-value
  \param x1 maximum x-value
  \param nx nr. of grid points in x
*/
void slip2raw2d(const string& infilebase_up,
		const string& infilebase_down,
		const string& outfilename,
		int t0,
		int nt,
		int dt,
		double x0,
		double x1,
		int nx)
{
  std::cout << "slip2raw2d t0,nt,dt: " << t0 << " " << nt << " " << dt <<  std::endl;  
  map<double,double> upper_map,lower_map;
  ostringstream infilename_up,infilename_down;
  vector<double> disp_old(nx);
  // open output file
  ofstream outfile(outfilename.c_str());
  // first ts -> fill up disp_old
  infilename_up << infilebase_up << "." << t0 << ".dat";  
  infilename_down << infilebase_down << "." << t0 << ".dat";  
  upper_map=read_file_to_map(infilename_up.str());
  lower_map=read_file_to_map(infilename_down.str());
  double dx=(x1-x0)/double(nx);
  for(int i=0;i<nx;i++){
    double px=x0+double(i)*dx; // current sample point
    double disp_up,disp_down;
    disp_up=get_disp_from_map_linear(upper_map,px);
    disp_down=get_disp_from_map_linear(lower_map,px);
    disp_old[i]=disp_up-disp_down;
  }
  // clean up
  lower_map.clear();
  upper_map.clear();
  infilename_up.seekp(0,ios_base::beg);
  infilename_down.seekp(0,ios_base::beg);
  
  // ---- further steps ----
  for(int t=t0+dt;t<nt;t+=dt){
    
    // generate input filenames
    infilename_up << infilebase_up << "." << t << ".dat";  
    infilename_down << infilebase_down << "." << t << ".dat";  
    std::cout << infilename_up.str()  << std::endl;
    upper_map=read_file_to_map(infilename_up.str());
    lower_map=read_file_to_map(infilename_down.str());
    // --- calc 
    // get x-step
    double dx=(x1-x0)/double(nx);
    for(int i=0;i<nx;i++){
      double px=x0+double(i)*dx; // current sample point
      double disp_up,disp_down;
      disp_up=get_disp_from_map_linear(upper_map,px);
      disp_down=get_disp_from_map_linear(lower_map,px);
      double disp=disp_up-disp_down;
      outfile << disp-disp_old[i] << " ";
    }
    outfile << std::endl;
    // clean up
    lower_map.clear();
    upper_map.clear();
    infilename_up.seekp(0,ios_base::beg);
    infilename_down.seekp(0,ios_base::beg);
  }  
  // close output file
  outfile.close();
}
