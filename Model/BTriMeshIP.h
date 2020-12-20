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

#ifndef __BTRIMESHIP_H
#define __BTRIMESHIP_H

#include <string>

class BTriMeshIP
{
public:
  double k;
  double brk;

  BTriMeshIP();

  BTriMeshIP(
    const std::string &interactionName,
    const std::string &meshName,
    double normalK,
    double breakDistance
  );
  
  virtual ~BTriMeshIP();

  void setMeshName(const std::string &meshName);

  const std::string &getMeshName() const;

  void setName(const std::string &name);

  const std::string &getName() const;

  virtual std::string getTypeString() const;

private:
  std::string m_name;
  std::string m_meshName;
};

class MeshBuildPrms
{
public:
  MeshBuildPrms();

  virtual ~MeshBuildPrms();

  virtual std::string getTypeString() const = 0;

private:
};

class MeshTagBuildPrms : public MeshBuildPrms
{
public:
  int m_tag;
  int m_mask;

  MeshTagBuildPrms();

  MeshTagBuildPrms(int tag, int mask);

  virtual std::string getTypeString() const;
};

class MeshGapBuildPrms : public MeshBuildPrms
{
public:
  double m_maxGap;
  
  MeshGapBuildPrms();

  MeshGapBuildPrms(double maxGap);
  
  virtual std::string getTypeString() const;
};

#endif // __BTRIMESHIP_H
