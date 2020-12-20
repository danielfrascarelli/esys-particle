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

/* extract connected aggregates (grains) from snapshots */
/* Options:
   -i FILENAME : input file name
   -oc FILENAME : output nr. or particles per grain
   -list FILENAME : output list of particle id/ grain id pairs
   -om FILENAME : output mass of grains
   -dm FILENAME BASE CUMULATIVE(0/1) : output grain mass distribution
   -dd FILENAME BASE CUMULATIVE(0/1) : output grain diameter distribution
   -ds FILENAME BASE SMIN: "sieve" distribution
   -p PERC: get sieve size (diameter) for x-th mass percentile
   -vtk FILENAME MINMASS : write out grains above a minimal mass as vtk-xml files  
   -tag MINTAG : ignore particles with tag<MINTAG
   -tagrange MINTAG MAXTAG : include only particles with MINTAG <= tag <= MAXTAG
   -diff FILENAME_1 FILENAME_2 : difference between to snapshots (req. -rmass or -fvtk) 
   -rmass FILENAME : in combination with -diff, write relative size of largest fragment
   -rmf FILENAME COUNT : in combination with -diff, write the sizes of all fragments, return nr. of fragments written
   -minmass MINMASS : with -rmf, writes only fragments of grains with m>MINMASS 
   -fvtk FILENAME EXP :
   -fv : write new grain id as fake vector field
   -profile FILENAME YMIN YMAX NBIN write y-profile of equiv. diameter
   -mprofile FILENAME YMIN YMAX NBIN SIZE_LIMIT write y-profile of matrix fraction
   -grid FILENAME XMIN XMAX YMIN YMAX ZMIN ZMAX CELLDIM
   -crossxy POS FILENAME W H FS : cross section in x-y plane, FS: 1 to filter singles, 0 not to filter
   -crossxz POS FILENAME W H FS : cross section in x-z plane
   -crossyz POS FILENAME W H FS : cross section in x-z plane
   -ppm : if cross section, save as PPM
   -abrasion FILENAME : output splitting/abrasion generated area. FILENAME: original (time step 0) data
   -btag TAG : only consider bonds with tag TAG
   -ocp FILENAME : write list of grain id / grain center position
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <limits>

// --- Project includes ---
#include "graph.h"
#include "readSnap.h"
#include "Frac.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ofstream;

int main(int argc, char** argv)
{
  string infilename;
  string infilename2;
  string countoutfilename;
  string massoutfilename;
  string countdistfilename;
  string massdistfilename;
  string sdistfilename;
  string rmassfilename;
  string vtkfilename;
  string listfilename;
  string rotlistfilename;
  string diadistfilename;
  string profilename;
  string crossfilename;
  string originalFileName;
  string posfilename;
  bool options_valid=true;
  double massdistbase=0.0;
  double minmass=0.0;
  double diadistbase=0.0;
  double sdistbase=0.0;
  double xmin=0.0,xmax=0.0,ymin=0.0,ymax=0.0,zmin=0.0,zmax=0.0;
  double celldim=0.0, slimit=0.0;
  double cross_pos=0.0;
  double matrix_cutoff=0.0;
  double split_ratio=0.0;
  double perc=0.0;
  float exp=0.0;
  int mintag=-1;
  int maxtag=std::numeric_limits<int>::max(); // defaults to largerst representable integer
  std::cerr << "maxtag=" << maxtag << std::endl;
  int rmf_count=0;
  int ret=0;
  int nbin=0;
  int cum=0;
  int cross_width=0, cross_height=0;
  int btag=-1;
  bool printcount=false;
  bool printmass=false;
  bool printcountdist=false;
  bool printmassdist=false;
  bool is_diff=false;
  bool writemassratio=false;
  bool writeallmass=false;
  bool writeallvtk=false;
  bool writegrains=false;
  bool writeidlist=false;
  bool writerotlist=false;
  bool writefracvtk=false;
  bool fakevectors=false;
  bool writemasstagged=false;
  bool printdiadist=false;
  bool writeprofile=false;
  bool writemprofile=false;
  bool writegrid=false;
  bool writematrixfraction=false;
  bool crossxy=false;
  bool crossxz=false;
  bool crossyz=false;
  bool write_as_ppm=false;
  bool cross_filter_singles=false;
  bool abrasion=false;
  bool writesdist=false;
  bool getperc=false;
  bool readgeoformat=false;
  bool writegraincenterpos=false;

  // process args
  int args_read=1;  
  while(args_read<argc){
    string option=string(argv[args_read]);
    if(option=="-i"){
      if(argc>args_read){
        infilename=string(argv[args_read+1]);
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-oc"){
      if(argc>args_read){
        countoutfilename=string(argv[args_read+1]);
        printcount=true;
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-list"){
      if(argc>args_read){
        listfilename=string(argv[args_read+1]);
        writeidlist=true;
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-rotlist"){
      if(argc>args_read){
        rotlistfilename=string(argv[args_read+1]);
        writerotlist=true;
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-om"){
      if(argc>args_read){
        massoutfilename=string(argv[args_read+1]);
        printmass=true;
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-dm"){
      if(argc>args_read+2){
        massdistfilename=string(argv[args_read+1]);
        massdistbase=atof(argv[args_read+2]);
        cum=atoi(argv[args_read+3]);
        printmassdist=true;
        args_read+=4;
      } else {
        options_valid=false;
      }
    }  else if(option=="-profile"){
      if(argc>args_read+3){
        profilename=string(argv[args_read+1]);
        ymin=atof(argv[args_read+2]);
        ymax=atof(argv[args_read+3]);
        nbin=atoi(argv[args_read+4]);
        writeprofile=true;
        args_read+=5;
      } else {
        options_valid=false;
      }
    } else if(option=="-mprofile"){
      if(argc>args_read+4){
        profilename=string(argv[args_read+1]);
        ymin=atof(argv[args_read+2]);
        ymax=atof(argv[args_read+3]);
        nbin=atoi(argv[args_read+4]);
        matrix_cutoff=atof(argv[args_read+5]);
        writemprofile=true;
        args_read+=6;
      } else {
        options_valid=false;
      }
    }  else if(option=="-grid"){
     if(argc>args_read+7){
        profilename=string(argv[args_read+1]);
        xmin=atof(argv[args_read+2]);
        xmax=atof(argv[args_read+3]);
        ymin=atof(argv[args_read+4]);
        ymax=atof(argv[args_read+5]);
        zmin=atof(argv[args_read+6]);
        zmax=atof(argv[args_read+7]);
        celldim=atof(argv[args_read+8]);
        writegrid=true;
        args_read+=9;
      } else {
        options_valid=false;
      }
    }  else if(option=="-matrix"){
      if(argc>args_read+8){
        profilename=string(argv[args_read+1]);
        xmin=atof(argv[args_read+2]);
        xmax=atof(argv[args_read+3]);
        ymin=atof(argv[args_read+4]);
        ymax=atof(argv[args_read+5]);
        zmin=atof(argv[args_read+6]);
        zmax=atof(argv[args_read+7]);
        celldim=atof(argv[args_read+8]);
        slimit=atof(argv[args_read+9]);
        writematrixfraction=true;
        args_read+=10;
      } else {
        options_valid=false;
      }
    } else if(option=="-dd"){
      if(argc>args_read+2){
        diadistfilename=string(argv[args_read+1]);
        diadistbase=atof(argv[args_read+2]);
        cum=atoi(argv[args_read+3]);
        printdiadist=true;
        args_read+=4;
      } else {
        options_valid=false;
      }
    } else if(option=="-ds"){
      if(argc>args_read+1){
        sdistfilename=string(argv[args_read+1]);
        sdistbase=atof(argv[args_read+2]);
        writesdist=true;
        args_read+=3;
      } else {
        options_valid=false;
      }
    } else if(option=="-p"){
      if(argc>args_read){
        perc=atof(argv[args_read+1]);
        getperc=true;
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-vtk"){
      if(argc>args_read+1){
        vtkfilename=string(argv[args_read+1]);
        minmass=atof(argv[args_read+2]);
        writegrains=true;
        args_read+=3;
      } else {
        options_valid=false;
      }
    } else if(option=="-avtk"){
      if(argc>args_read+1){
        vtkfilename=string(argv[args_read+1]);
        writeallvtk=true;
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-fvtk"){
      if(argc>args_read+1){
        vtkfilename=string(argv[args_read+1]);
        exp=atof(argv[args_read+2]);
        writefracvtk=true;
        args_read+=3;
      } else {
        options_valid=false;
      }
    } else if(option=="-dc"){
      if(argc>args_read){
        countdistfilename=string(argv[args_read+1]);
        printcountdist=true;
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-tag"){
      if(argc>args_read){
        mintag=atoi(argv[args_read+1]);
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-tagrange"){
      if(argc>args_read+1){
        mintag=atoi(argv[args_read+1]);
        maxtag=atoi(argv[args_read+2]);
	std::cerr << "maxtag=" << maxtag << std::endl;
        args_read+=3;
      } else {
        options_valid=false;
      }
    } else if(option=="-minmass"){
      if(argc>args_read){
        minmass=atof(argv[args_read+1]);
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-diff"){
      if(argc>args_read+1){
        infilename=string(argv[args_read+1]);
        infilename2=string(argv[args_read+2]);
        args_read+=3;
        is_diff=true;
      } else {
        options_valid=false;
      }
    } else if(option=="-rmass"){
      if(argc>args_read){
        rmassfilename=string(argv[args_read+1]);
        writemassratio=true;
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-rmf"){
      if(argc>args_read+1){
        rmassfilename=string(argv[args_read+1]);
        rmf_count=atoi(argv[args_read+2]);
        writeallmass=true;
        args_read+=3;
      } else {
        options_valid=false;
      }
    } else if(option=="-crossxy"){
      if(argc>args_read+8){
        cross_pos=atof(argv[args_read+1]);
        xmin=atof(argv[args_read+2]);
        xmax=atof(argv[args_read+3]);
        ymin=atof(argv[args_read+4]);
        ymax=atof(argv[args_read+5]);
        crossfilename=string(argv[args_read+6]);
        cross_width=atoi(argv[args_read+7]);
        cross_height=atoi(argv[args_read+8]);
        cross_filter_singles=(atoi(argv[args_read+9])==1);
        crossxy=true;
        args_read+=10;
      } else {
        options_valid=false;
      }
    }else if(option=="-crossxz"){
      if(argc>args_read+8){
        cross_pos=atof(argv[args_read+1]);
        xmin=atof(argv[args_read+2]);
        xmax=atof(argv[args_read+3]);
        ymin=atof(argv[args_read+4]);
        ymax=atof(argv[args_read+5]);
        crossfilename=string(argv[args_read+6]);
        cross_width=atoi(argv[args_read+7]);
        cross_height=atoi(argv[args_read+8]);
        cross_filter_singles=(atoi(argv[args_read+9])==1);
        crossxz=true;
        args_read+=10;
      } else {
        options_valid=false;
      }
    }else if(option=="-crossyz"){
      if(argc>args_read+8){
        cross_pos=atof(argv[args_read+1]);
        xmin=atof(argv[args_read+2]);
        xmax=atof(argv[args_read+3]);
        ymin=atof(argv[args_read+4]);
        ymax=atof(argv[args_read+5]);
        crossfilename=string(argv[args_read+6]);
        cross_width=atoi(argv[args_read+7]);
        cross_height=atoi(argv[args_read+8]);
        cross_filter_singles=(atoi(argv[args_read+9])==1);
        crossyz=true;
        args_read+=10;
      } else {
        options_valid=false;
      }
    }else if(option=="-rmft"){
      if(argc>args_read+1){
        rmassfilename=string(argv[args_read+1]);
        rmf_count=atoi(argv[args_read+2]);
        writeallmass=true;
        writemasstagged=true;
        args_read+=3;
      } else {
        options_valid=false;
      }
    } else if(option=="-fv"){
      fakevectors=true;
      args_read++;
    } else if(option=="-ppm"){
      write_as_ppm=true;
      args_read++;
    } else if(option=="-abrasion"){
      if(argc>args_read+2){
        originalFileName=string(argv[args_read+1]);
        split_ratio=atof(argv[args_read+2]);
        minmass=atof(argv[args_read+3]);
        abrasion=true;
        args_read+=4;
      } else {
        options_valid=false;
      }
    } else if(option=="-btag"){
      if(argc>args_read){
        btag=atoi(argv[args_read+1]);
	std::cerr << "bond tag: " << btag << std::endl;
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-ocp"){
      if(argc>args_read){
        posfilename=string(argv[args_read+1]);
        writegraincenterpos=true;
        args_read+=2;
      } else {
        options_valid=false;
      }
    } else if(option=="-geo"){
      readgeoformat=true;
      args_read++;
    } else {
      cerr << "Unknown option " << option << endl;
      options_valid = false;
      break;
    }
  }

  if(options_valid){
    // read data from snapshots into graph
    Graph g;
    if(readgeoformat){
	readGeo(infilename,mintag,maxtag,g,btag);
    } else {
	readSnap(infilename,mintag,maxtag,g,btag);
    }
    g.makeConnComp();
    if(is_diff){
      // differential between snapshots 
      Graph g2(readSnap(infilename2,mintag,maxtag,btag)); 
      g2.makeConnComp();
      Frac fr(g,g2);
      if(writemassratio){
        ofstream outfile(rmassfilename.c_str());
        fr.writeMassRatio(outfile,minmass);
        outfile.close();
      }
      if(writeallmass){
        ofstream outfile(rmassfilename.c_str());
        ret=fr.writeAllMass(outfile,minmass,rmf_count,writemasstagged);
        outfile.close();
      }
      if(writefracvtk){
        fr.writeAsVtk(vtkfilename,exp,g2,fakevectors);
      }
      if(abrasion){
        Graph orig(readSnap(originalFileName,mintag,maxtag));
        pair<double,double> abratio=fr.getSplitAbrasion(orig,split_ratio,minmass);
        std::cout << "Split: " << abratio.first << " Abrasion: " << abratio.second << std::endl;
      }
     } else {
      // single snapshot data
      if(printcount){
        ofstream outfile(countoutfilename.c_str());
        g.printGrainPCount(outfile);
        outfile.close();
      }
      if(printmass){
        ofstream outfile(massoutfilename.c_str());
        g.printGrainMass(outfile);
        outfile.close();
      }
      if(printcountdist){
        g.printGrainCountDist(countdistfilename);
      }
      if(printmassdist){
        g.printGrainMassDist(massdistfilename,massdistbase,cum);
      }
      if(printdiadist){
        g.printGrainDiamDist(diadistfilename,diadistbase,cum);
      }
      if(writesdist){
        g.printSieveDist(sdistfilename,sdistbase);
      }
      if(getperc){
        std::cout << "Diameter at " << perc << " percentile mass : " << g.getPercentile(perc) << std::endl;
      }
      if(writegrains){
        g.printGrainsAsVtk(vtkfilename,minmass);
      }
      if(writeallvtk){
	g.printAllAsVtk(vtkfilename);
      }
      if(writeidlist){
        g.printIdList(listfilename);
      }
      if(writerotlist){
        g.printRotList(rotlistfilename);
      }
      if(writeprofile){
        g.writeAvgGrainSizeProfile(profilename,ymin,ymax,nbin);
      }
      if(writemprofile){
        g.writeMatrixFractionProfile(profilename,ymin,ymax,nbin,matrix_cutoff);
      }
      if(writegrid){
        g.writeAvgGrainSizeGrid(profilename,xmin,xmax,ymin,ymax,zmin,zmax,celldim);
      }
      if(writematrixfraction){
        g.writeMatrixFractionGrid(profilename,xmin,xmax,ymin,ymax,zmin,zmax,celldim,slimit);
      }
      if(crossxy){
        g.printCrossSection(Vec3(0.0,0.0,cross_pos),Vec3(1.0,0.0,0.0),Vec3(0.0,1.0,0.0),crossfilename,cross_width,cross_height,xmin,xmax,ymin,ymax,write_as_ppm,cross_filter_singles);
      }
      if(crossxz){
        g.printCrossSection(Vec3(0.0,cross_pos,0.0),Vec3(1.0,0.0,0.0),Vec3(0.0,0.0,-1.0),crossfilename,cross_width,cross_height,xmin,xmax,ymin,ymax,write_as_ppm,cross_filter_singles);
      }
      if(crossyz){
        g.printCrossSection(Vec3(cross_pos,0.0,0.0),Vec3(0.0,0.0,-1.0),Vec3(0.0,1.0,0.0),crossfilename,cross_width,cross_height,xmin,xmax,ymin,ymax,write_as_ppm,cross_filter_singles);
      }
      if(writegraincenterpos){
	g.printGrainCenterPosition(posfilename);
      }
      ret=0;
    } 
  }else {
    ret=1;
  }

  return ret;
}
