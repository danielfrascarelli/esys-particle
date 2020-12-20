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

#ifndef __NTABLE_H
#define __NTABLE_H

#include "Foundation/console.h"

//--- project includes ---
#include "Foundation/vec3.h"
#include "ntable/src/dslice.h"
#include "ntable/src/handle.h"
#include "Geometry/Triangle.h"
#include "Geometry/AEdge.h"
#include "Model/FluidCell.h" //fluid contents

//--- STL includes ---
#include <utility>
#include <valarray>
#include <vector>
#include <list>

#if HAVE_CONFIG_H
#include <config.h>
#endif

/* #if HAVE_HASH_MAP */
/* #include <hash_map> */
/* #else */
#include <map>
/* #endif */

using std::valarray;
using std::vector;
using std::list;
using std::pair;
using std::make_pair;

//--- IO includes ---
#include <iostream>

//--- forward declarations ---
template<typename T>
class NTSlab;

template<typename T>
class NTBlock;

template<typename T>
class NTSlab_iter;

/*!
  \class NeighborTable
  \brief class for neighbor search
*/
template <typename T>
class NeighborTable
{
 public: // types
  typedef valarray<vector<typename list<T>::iterator> > arraytype;
  typedef vector<typename list<T>::iterator>            pointtype;
  typedef pair<int,int>                                 indextype;
  typedef list<pair<T*,T*> >                            pairlist;
  typedef list<T*>                                      particlelist;
/* #if HAVE_HASH_MAP */
/*   typedef std::hash_map<int, T*>                        IdParticleMap; */
/* #else */
  typedef std::map<int, T*>                             IdParticleMap;
/* #endif */
  typedef typename list<T>::iterator                    iterator;

  typedef list<pair<T*,CFluidCell*> >                   particlecelllist; //fluid contents

 private:
  // member variables
  list<T>       m_list;          //!< list of particles
  arraytype     m_array;         //!< search array
  IdParticleMap m_idParticleMap; //!< mapping between particle-id and particle-pointer
  Vec3          m_p0_global;     //!< minimum corner of global search space
  Vec3          m_pmax_global;   //!< maximum corner of global search space //fluid contents
  double        m_dim;           //!< grid spacing of search array
  double        m_alpha;         //!< padding factor (dim=2*rmax+alpha)
  int           m_global_idx;    //!< minimum corner index (x component)
  int           m_global_idy;    //!< minimum corner index (y component)
  int           m_global_idz;    //!< minimum corner index (z component)
  int           m_xsize;         //!< number of grid point, x direction
  int           m_ysize;         //!< number of grid point, y direction
  int           m_zsize;         //!< number of grid point, z direction
  bool          m_valid;
  Vec3          m_min_corner;    //!< minimum corner of the search array
  Vec3          m_max_corner;    //!< maximum corner of the search array

  /****fluid contents: begin****/
  bool          m_fluidexist;
  arraytype     m_cellarray;     //!< fluid cell search array
  double        m_xside;         //!< side length of fluid cells in x direction
  double        m_yside;         //!< side length of fluid cells in y direction
  double        m_zside;         //!< side length of fluid cells in z direction
  int           m_global_cellidx;//!< minimum corner index of fluid cell(x component)
  int           m_global_cellidy;//!< minimum corner index of fluid cell(y component)
  int           m_global_cellidz;//!< minimum corner index of fluid cell(z component)
  int           m_xcell;         //!< number of fluid cells, x direction
  int           m_ycell;         //!< number of fluid cells, y direction
  int           m_zcell;         //!< number of fluid cells, z direction

  double        m_Bw;            //!< bulk modulus of fluid
  double        m_Bp;            //!< bulk modulus of solid particle
  double        m_Mu;            //!< fluid viscosity
  double        m_adjust;        //!< adjusting factor between two time steps
  double        m_flowrate;      //!< inflow rate
  double        m_pressure;      //!< pressure gradient
  Vec3          m_inflow;        //!< inflow directions (0 - closed boundary, 1 - positive axis direction, -1 - negative axis direction)
  Vec3          m_outflow;       //!< outflow directions (0 - closed boundary, 1 - positive axis direction, -1 - negative axis direction)
  vector<CFluidCell> m_cells;
  /****fluid contents: end****/


  // private member functions
  void clear_search_array();
  void clear_cell_array(); //fluid contents
  int index(const Vec3&); // index from position
  int cellindex(const Vec3&); // fluid contents

  void addPairsToList(T_Handle<pairlist>,int,int);
  void addPairsToListLocal(T_Handle<pairlist>,int);
  void addPairsToListFlagged(T_Handle<pairlist>,int,int);
  void addPairsToListLocalFlagged(T_Handle<pairlist>,int);

  arraytype *array(){return &m_array;};

 public:
  //! Constructors
  NeighborTable();
  NeighborTable(int,int,int,double,double,const Vec3&,const Vec3&,int,int,int);
  ~NeighborTable();

  /****fluid contents: begin***/
  void initiateFluid(double,double,double,double,double,double,Vec3,Vec3,int,int,int,int,int,int);
  void calPorositySphere0();
  void calPorositySphere();
  void updateFluidcells();
  void calCoeffi(double,int);
  void setPressure(vector<pair<Vec3,double> >);
  void calVelocity();
  vector<CFluidCell> get_yz_cellslab(int);
  vector<CFluidCell> get_xz_cellslab(int);
  vector<CFluidCell> get_xy_cellslab(int);
  void set_yz_cellslab(int, vector<CFluidCell>);
  void set_xz_cellslab(int, vector<CFluidCell>);
  void set_xy_cellslab(int, vector<CFluidCell>);
  template<typename P> vector<pair<Vec3,P> > forAllInnerCellsGet(P (CFluidCell::*rdf)() const);
  template<typename P> vector<pair<Vec3,P> > forAllInnerCellsGetIndexed(P (CFluidCell::*rdf)() const);
  template<typename P> vector<P> forAllInnerCellsGetSum(P (CFluidCell::*rdf)() const);
  /****fluid contents: end***/

  void insert(const T&); //!< particle insertion
  void insert(iterator i,const T& data){insert(data);}; //!< STL compat. insert
  void build(); //!< build search array
  inline int index(int,int,int) const; // total index from x,y,z indices
  inline int cellindex(int,int,int) const; // fluid contents

  //!< iterators
  iterator begin(){return m_list.begin();};
  iterator end(){return m_list.end();};

  //!< number of particles at a given gridpoint
  unsigned int nparts_at_gridpoint(unsigned int idx) const
  {
    return (m_array[idx]).size();
  }

  // check in posn. is in inner part
  bool isInInner(const Vec3&);

  //!< dimensions
  int xsize(){return m_xsize;};
  int ysize(){return m_ysize;};
  int zsize(){return m_zsize;};
  int xcell(){return m_xcell;}; //fluid contents
  int ycell(){return m_ycell;}; //fluid contents
  int zcell(){return m_zcell;}; //fluid contents
  int size(){return m_list.size();};
  Vec3 base_point() const {return m_p0_global;};
  Vec3 max_point() const {return m_pmax_global;};
  int base_idx_x() const {return m_global_idx;};
  int base_idx_y() const {return m_global_idy;};
  int base_idx_z() const {return m_global_idz;};
  int base_cellidx_x() const {return m_global_cellidx;}; //fluid contents
  int base_cellidx_y() const {return m_global_cellidy;}; //fluid contents
  int base_cellidx_z() const {return m_global_cellidz;}; //fluid contents
  double dim(){return m_dim;}; //fluid contents

  //!< partial access functions
  NTSlab<T> xy_slab(int);
  NTSlab<T> xz_slab(int);
  NTSlab<T> yz_slab(int);
  NTBlock<T> block(int,int,int,int,int,int);
  NTBlock<T> block(const Vec3&,const Vec3&);
  NTBlock<T> inner();

  // access ops
  T* ptr(NeighborTable<T>::indextype);
  T& ref(NeighborTable<T>::indextype);
  T* ptr_by_id(int);
  T *getNearestPtr(const Vec3&);

  // erase particle (not meant to be used directly)
  void erase(NeighborTable<T>::indextype);

  // neigbor list creation
  T_Handle<pairlist> getFullList();
  T_Handle<pairlist> getNewList();
  T_Handle<particlelist> getParticlesAtPlane(const Vec3&,const Vec3&); // distance as parameter ??
  T_Handle<particlelist> getParticlesNearSphere(const Vec3&,const double&);
  T_Handle<particlelist> getParticlesNearTriangle(const Triangle&);
  T_Handle<particlelist> getParticlesNearEdge(const AEdge*);
  T_Handle<particlelist> getParticlesNearPoint(const Vec3&);
  T_Handle<particlelist> getAllParticles();
  T_Handle<particlecelllist> getParticleCellList();//fluid content
  T_Handle<particlecelllist> getNewParticleCellList();//fluid content

  //!output
  template <typename TT>
  friend std::ostream& operator<<(std::ostream &, const NeighborTable<TT> &);
};

#include "ntable/src/ntable.hpp"

#endif // __NTABLE_H
