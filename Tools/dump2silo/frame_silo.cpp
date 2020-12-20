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

#include "frame_silo.h"

// --- project includes ---
#include "vec3.h"

#include <silo.h>

// --- STL includes ---
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iterator>
#include <map>
#include <utility>
#include <sstream>

using std::vector;
using std::string;
using std::ifstream;
using std::istream_iterator;
using std::back_inserter;
using std::map;
using std::pair;
using std::make_pair;
using std::ostringstream;

typedef vector<float> FloatVec;
typedef vector<int> IntVec;

struct Arr3 {
    void reserve(int n) { x.reserve(n); y.reserve(n); z.reserve(n); }
    void push_back(const Vec3& v) {
        x.push_back(v.X());
        y.push_back(v.Y());
        z.push_back(v.Z());
    }
    FloatVec x;
    FloatVec y;
    FloatVec z;
};

struct nr_part{
  Vec3 pos;
  Vec3 vel;
  Vec3 force;
  double rad;
  double mass;
  int tag;
};

struct ParticleData {
    void append(int i, int t, float r, float m, const Vec3& p, const Vec3& v,
                const Vec3& f) {
        id.push_back(i);
        tag.push_back(t);
        rad.push_back(r);
        mass.push_back(m);
        pos.push_back(p);
        vel.push_back(v);
        force.push_back(f);
    }
    void append(int i, int t, float r, float m, const Vec3& p, const Vec3& v,
                const Vec3& f, const Vec3& av) {
        append(i, t, r, m, p, v, f);
        angvel.push_back(av);
    }

    void reserve(int n, bool rot=true) {
        id.reserve(n);
        tag.reserve(n);
        rad.reserve(n);
        mass.reserve(n);
        pos.reserve(n);
        vel.reserve(n);
        force.reserve(n);
        if (rot)
            angvel.reserve(n);
    }

    IntVec id;
    IntVec tag;
    FloatVec rad;
    FloatVec mass;
    Arr3 pos;
    Arr3 vel;
    Arr3 force;
    FloatVec q1,q2,q3,q4;
    Arr3 angvel;
};


DBfile* openSiloFile(const string& filename)
{
  DBfile* dbfile = NULL;

  // check if file exists and append data if so
  if (DBInqFile(filename.c_str()) > 0) {
    // open existing SILO file, try HDF5 format first
    dbfile = DBOpen(filename.c_str(), DB_HDF5, DB_APPEND);
    if (!dbfile)
      dbfile = DBOpen(filename.c_str(), DB_PDB, DB_APPEND);
  } else {
    // create new SILO file, prefer HDF5 format if available
    dbfile = DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_HDF5);
    if (!dbfile)
      dbfile = DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL, NULL, DB_PDB);
  }
  return dbfile;
}


void writeBonds(DBfile* dbfile, const ParticleData& data,
        const vector<pair<int,int> >& bonds, map<int,int>& id2idx)
{
    vector<int> nodes;
    vector<float> strain;

    for (vector<pair<int,int> >::const_iterator iter=bonds.begin(); iter!=bonds.end(); iter++) {
        int i1 = id2idx[iter->first];
        int i2 = id2idx[iter->second];
        nodes.push_back(i1);
        nodes.push_back(i2);
        Vec3 pos1(data.pos.x[i1], data.pos.y[i1], data.pos.z[i1]);
        Vec3 pos2(data.pos.x[i2], data.pos.y[i2], data.pos.z[i2]);
        double r1=data.rad[i1];
        double r2=data.rad[i2];
        double strn=((pos1-pos2).norm()-(r1+r2))/(r1+r2);
        strain.push_back((float)strn);
    }

    int nodesPerElement = 2;
    int numBonds = bonds.size();
    int elType = DB_ZONETYPE_BEAM;

    DBPutZonelist2(dbfile, "bonds", numBonds, 3, &nodes[0], nodesPerElement*numBonds, 0, 0, 0, &elType, &nodesPerElement, &numBonds, 1, NULL);

    DBPutUcdvar1(dbfile, "strain", "fullmesh", &strain[0], strain.size(), NULL, 0, DB_FLOAT, DB_ZONECENT, NULL);
}

int getVersion(const string& infilename)
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

vector<string> getFilenames(const string& infilename, int version)
{
  cout << "infilename : " << infilename << endl;
  ifstream headerfile(infilename.c_str()); 
  float dummy,xmax,ymax,zmax,xmin,ymin,zmin;
  vector<string> filenames;
  string dummystring;

  if(version==0){
    headerfile >> dummy >> dummy >> dummy;
    headerfile >> dummystring >> dummy;
  } else if ((version==1)||(version==2) || (version==3)) {
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

void readParticles(ParticleData& data, vector< pair<int,int> >& bonds,
        const string& infilename, bool rot)
{
    int version=getVersion(infilename);
    vector<string> filenames=getFilenames(infilename,version);

    for (vector<string>::iterator iter=filenames.begin(); iter!=filenames.end(); iter++) {
        cout << *iter << endl;
        ifstream datafile(iter->c_str());
        // get particles
        int npart;
        int id, tag;
        float rad, mass;
        Vec3 initpos, oldpos, pos, vel, force, angvel;
        datafile >> npart;
        data.reserve(npart, rot);
        if (version < 2) {
            for (int i=0; i<npart; i++) {
                datafile >> pos >> rad >> id >> tag >> mass >> initpos >> oldpos
                    >> vel >> force;
                if (rot) {
                    float q1, q2, q3, q4;
                    datafile >> q1 >> q2 >> q3 >> q4 >> angvel;
                    data.append(id, tag, rad, mass, pos, vel, force, angvel);
                } else {
                    data.append(id, tag, rad, mass, pos, vel, force);
                }
            } 
        } else {
            Vec3 circ_shift;
            for (int i=0; i<npart; i++) {
                datafile >> pos >> rad >> id >> tag >> mass >> initpos >> oldpos
                    >> vel >> force >> circ_shift;
                if (rot) {
                    float q1, q2, q3, q4;
                    datafile >> q1 >> q2 >> q3 >> q4 >> angvel;
                    data.append(id, tag, rad, mass, pos, vel, force, angvel);
                } else {
                    data.append(id, tag, rad, mass, pos, vel, force);
                }
            } 
        }  

        // get bonds
        int ngrp;
        int nbond;
        string type;

        datafile >> ngrp;
        if (version>0) {
            datafile >> type;
            if ((type=="Bonded") || (rot==true && type=="RotBonded")) {
                datafile >> nbond;
                for (int i=0; i<nbond; i++) {
                    int b1, b2;
                    float dummy;
              
                    datafile >> b1 >> b2 >> dummy;
                    bonds.push_back(make_pair(b1,b2));
                }
            }
        } else { // pre - V1 snapshot -> assume bondend pair IG
            datafile >> nbond;
            for (int i=0; i<nbond;i++) {
                int b1, b2;
                float dummy;
            
                datafile >> b1 >> b2 >> dummy;
                bonds.push_back(make_pair(b1,b2));
            }
        }

        datafile.close();
    }
}

void readParticlesSliced(ParticleData& data, vector< pair<int,int> >& bonds,
        const string& infilename, bool rot, double minZ, double maxZ)
{
    int version=getVersion(infilename);
    vector<string> filenames=getFilenames(infilename,version);

    for (vector<string>::iterator iter=filenames.begin(); iter!=filenames.end(); iter++) {
        cout << *iter << endl;
        ifstream datafile(iter->c_str());
        // get particles
        int npart;
        int id, tag;
        float rad, mass;
        Vec3 initpos, oldpos, pos, vel, force, angvel;
        datafile >> npart;
        data.reserve(npart, rot);
        if (version < 2) {
            for (int i=0; i<npart; i++) {
                datafile >> pos >> rad >> id >> tag >> mass >> initpos >> oldpos
                    >> vel >> force;
                // flatten in Z to middle of slice area
                pos.Z() = (minZ+maxZ)*0.5;
                if (rot) {
                    float q1, q2, q3, q4;
                    datafile >> q1 >> q2 >> q3 >> q4 >> angvel;
                    // only append if in slice
                    if ((pos.Z() > minZ) && (pos.Z() < maxZ)) {
                        data.append(id, tag, rad, mass, pos, vel, force, angvel);
                    }
                } else {
                    // only append if in slice
                    if ((pos.Z() > minZ) && (pos.Z() < maxZ)) {
                        data.append(id, tag, rad, mass, pos, vel, force);
                    }
                }
            } 
        } else {
            Vec3 circ_shift;
            for (int i=0; i<npart; i++) {
                datafile >> pos >> rad >> id >> tag >> mass >> initpos >> oldpos
                    >> vel >> force >> circ_shift;
                // flatten in Z to middle of slice area
                pos.Z() = (minZ+maxZ)*0.5;
                if (rot) {
                    float q1, q2, q3, q4;
                    datafile >> q1 >> q2 >> q3 >> q4 >> angvel;
                    // only append if in slice
                    if ((pos.Z() > minZ) && (pos.Z() < maxZ)) {
                        data.append(id, tag, rad, mass, pos, vel, force, angvel);
                    }
                } else {
                    // only append if in slice
                    if ((pos.Z() > minZ) && (pos.Z() < maxZ)) {
                        data.append(id, tag, rad, mass, pos, vel, force);
                    }
                }
            } 
        }  

        datafile.close();
    }
}

void saveSiloSnap(const string& infilename, const string& outfilename,
        int iframe, bool with_list, const string& listfilename, bool rot)
{
    ParticleData data;
    vector< pair<int,int> > bonds;
    map<int,int> id2idx;

    readParticles(data, bonds, infilename, rot);

    // generate reverse mapping between particle id and point idx
    int count=0;
    for (IntVec::const_iterator iter=data.id.begin(); iter!=data.id.end(); iter++)
    {
        id2idx.insert(make_pair(*iter, count));
        count++;
    }

    // open output file
    ostringstream silofilename;
    silofilename << outfilename << iframe << ".silo";
    DBfile* dbfile = openSiloFile(silofilename.str());

    if (!dbfile) {
        cerr << "Could not open output file " << silofilename.str() << endl;
        return;
    }

    float* coords[] = { &data.pos.x[0], &data.pos.y[0], &data.pos.z[0] };
    DBPutUcdmesh(dbfile, "fullmesh", 3, NULL, coords, data.pos.x.size(), bonds.size(), "bonds", NULL, DB_FLOAT, NULL);

    // write bonds 
    writeBonds(dbfile, data, bonds, id2idx);

    // FIXME: Have to decide whether to use pointmesh or not!
if (1) {
    DBPutUcdvar1(dbfile, "radius", "fullmesh", &data.rad[0], data.rad.size(), NULL, 0, DB_FLOAT, DB_NODECENT, NULL);
    DBPutUcdvar1(dbfile, "tag", "fullmesh", (float*)&data.tag[0], data.tag.size(), NULL, 0, DB_INT, DB_NODECENT, NULL);
    DBPutUcdvar1(dbfile, "id", "fullmesh", (float*)&data.id[0], data.id.size(), NULL, 0, DB_INT, DB_NODECENT, NULL);
    coords[0] = &data.vel.x[0];
    coords[1] = &data.vel.y[0];
    coords[2] = &data.vel.z[0];
    char* varnames[] = { "velx", "vely", "velz" };
    DBPutUcdvar(dbfile, "velocity", "fullmesh", 3, varnames, coords, data.vel.x.size(), NULL, 0, DB_FLOAT, DB_NODECENT, NULL);
    if (rot) {
        coords[0] = &data.angvel.x[0];
        coords[1] = &data.angvel.y[0];
        coords[2] = &data.angvel.z[0];
        char* varnames[] = { "angvelx", "angvely", "angvelz" };
        DBPutUcdvar(dbfile, "angular_velocity", "fullmesh", 3, varnames, coords, data.angvel.x.size(), NULL, 0, DB_FLOAT, DB_NODECENT, NULL);
    }

} else {

    DBPutPointmesh(dbfile, "mesh", 3, coords, data.pos.x.size(), DB_FLOAT, NULL);
    DBPutPointvar1(dbfile, "radius", "mesh", &data.rad[0], data.rad.size(), DB_FLOAT, NULL);
    DBPutPointvar1(dbfile, "tag", "mesh", (float*)&data.tag[0], data.tag.size(), DB_INT, NULL);
    DBPutPointvar1(dbfile, "id", "mesh", (float*)&data.id[0], data.id.size(), DB_INT, NULL);
    coords[0] = &data.vel.x[0];
    coords[1] = &data.vel.y[0];
    coords[2] = &data.vel.z[0];
    DBPutPointvar(dbfile, "velocity", "mesh", 3, coords, data.vel.x.size(), DB_FLOAT, NULL);

    if (rot) {
        coords[0] = &data.angvel.x[0];
        coords[1] = &data.angvel.y[0];
        coords[2] = &data.angvel.z[0];
        DBPutPointvar(dbfile, "angular_velocity", "mesh", 3, coords, data.angvel.x.size(), DB_FLOAT, NULL);
    }

    // if option set, read list file int ng_id_map;
    // current grain id
    if (with_list) {
        map<int,int> ng_id_map;
        vector<int> grId;
        ifstream listfile(listfilename.c_str());
        int id, nid;
        while (!listfile.eof()) {
            listfile >> id >> nid;
            ng_id_map[id]=nid;
        }
        listfile.close();

        for (IntVec::const_iterator iter=data.id.begin(); iter!=data.id.end(); iter++) {
            grId.push_back(ng_id_map[*iter]);
        }
        DBPutPointvar1(dbfile, "grainId", "mesh", (float*)&grId[0], grId.size(), DB_INT, NULL);
    }

}

    // close file
    DBClose(dbfile);
}


/*!
  writes a slice, doesn't write bonds
*/
void saveSiloSnapSliced(const string& infilename, const string& outfilename,
        int iframe, bool with_list, const string& listfilename,
        double slz_min, double slz_max)
{  
    cout << "saveSiloSnapSliced(" << slz_min << ", " << slz_max << ")" << endl;

    ParticleData data;
    vector< pair<int,int> > bonds;
    map<int,int> id2idx;

    readParticles(data, bonds, infilename, true);

    // generate reverse mapping between particle id and point idx
    int count=0;
    for (IntVec::const_iterator iter=data.id.begin(); iter!=data.id.end(); iter++)
    {
        id2idx.insert(make_pair(*iter, count));
        count++;
    }

    // open output file
    ostringstream silofilename;
    silofilename << outfilename << iframe << ".silo";
    DBfile* dbfile = openSiloFile(silofilename.str());

    if (!dbfile) {
        cerr << "Could not open output file " << silofilename.str() << endl;
        return;
    }

    float* coords[] = { &data.pos.x[0], &data.pos.y[0], &data.pos.z[0] };

    DBPutPointmesh(dbfile, "mesh", 3, coords, data.pos.x.size(), DB_FLOAT, NULL);
    DBPutPointvar1(dbfile, "radius", "mesh", &data.rad[0], data.rad.size(), DB_FLOAT, NULL);
    DBPutPointvar1(dbfile, "tag", "mesh", (float*)&data.tag[0], data.tag.size(), DB_INT, NULL);
    DBPutPointvar1(dbfile, "id", "mesh", (float*)&data.id[0], data.id.size(), DB_INT, NULL);
    coords[0] = &data.vel.x[0];
    coords[1] = &data.vel.y[0];
    coords[2] = &data.vel.z[0];
    DBPutPointvar(dbfile, "velocity", "mesh", 3, coords, data.vel.x.size(), DB_FLOAT, NULL);
    coords[0] = &data.angvel.x[0];
    coords[1] = &data.angvel.y[0];
    coords[2] = &data.angvel.z[0];
    DBPutPointvar(dbfile, "angular_velocity", "mesh", 3, coords, data.angvel.x.size(), DB_FLOAT, NULL);

    // if option set, read list file int ng_id_map;
    // current grain id
    if (with_list) {
        map<int,int> ng_id_map;
        vector<int> grId;
        ifstream listfile(listfilename.c_str());
        int id, nid;
        while (!listfile.eof()) {
            listfile >> id >> nid;
            ng_id_map[id]=nid;
        }
        listfile.close();

        for (IntVec::const_iterator iter=data.id.begin(); iter!=data.id.end(); iter++) {
            grId.push_back(ng_id_map[*iter]);
        }
        DBPutPointvar1(dbfile, "grainId", "mesh", (float*)&grId[0], grId.size(), DB_INT, NULL);
    }

    // close file
    DBClose(dbfile);
}

