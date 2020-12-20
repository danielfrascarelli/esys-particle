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


#ifndef CHECKPOINTCONTROLLER_H
#define CHECKPOINTCONTROLLER_H

// --- MPI includes ---
#include <mpi.h>

// --- Project includes ---
#include "Geometry/GeometryInfo.h"
#include "Foundation/StringUtil.h"
#include "Foundation/BoundingBox.h"

// --- STL includes ---
#include <string>

/**
 * Controls the issue of check-pointing commands to slave processes.
 */
class CheckPointController
{
public:
  /**
   * Default constructor.
   *
   * Default parameters cause isCheckPoint to return false.
   */
  CheckPointController();

  /**
   * Instantiates and sets the parameters which determine when
   * a check-point should occur.
   *
   * @see setCheckpointParams
   */
  CheckPointController(
    const std::string &fileNamePrefix,
    int beginTime,
    int endTime,
    int timeInterval,
    bool writeThroughMaster
  );

  /**
   *
   */
  virtual ~CheckPointController();

  /**
   * Determines whether specified argument is a check-point
   * time (see isCheckPoint). If currentTime is a check-point,
   * issues commands to slave processes to perform check-point.
   *
   * @param currentTime Parameter used to determine whether
   *                    a check point should occur.
   */
  virtual void performCheckPoint(int currentTime);
  virtual void performSnapShot(int currentTime);

  /**
   * Issues the check-point command to slave processes.
   *
   * @param currentTime The check-point time.
   */
  virtual void issueCheckPointCmd(int currentTime);
  virtual void issueCheckPointCmdWTM(int currentTime);
  virtual void issueSnapShotCmd(int currentTime);
  virtual void issueCheckPointLoadingCmd(const std::string&);

  /**
   * Returns whether specified time is a check-point time.
   *
   * @param time This value is checked against the checkpoint
   *             parameters.
   */
  bool isCheckPoint(int time);

  /**
   * Sets the parameters which determine when a check-point should occur.
   */
  void setCheckPointParams(
    const std::string &fileNamePrefix,
    int beginTime,
    int endTime,
    int timeInterval,
    bool writeThroughMaster,
    int precision=12
  );

  std::string getLatticeDataFileName(const std::string &fileNamePrefix, int timeStep, int rank, bool bin=false);

  esys::lsm::StringVector getLatticeDataFiles(int timeStep, int size);

  esys::lsm::StringVector getLatticeDataFiles(int timeStep);

  /**
   * Set 2-D information to true if the particle data are two-dimensional; 
   * otherwise set to false.
   */
  void set_is2d(bool do2d);

  /**
   * Set the LSMGeometry version for use in geometry files. 
   */
  void setLsmGeoVersion(float version);

  /**
   * Set the periodicity of the x, y and z dimensions. 
   */
  void setPeriodicDimensions(esys::lsm::BoolVector periodicDimensions);

  /**
   * Sets the spatial extent in which particles are tracked.
   */
  void setGeometryInfo(const esys::lsm::GeometryInfo &geoInfo);

  /**
   * Sets geometry info.
   */
  void setSpatialDomain(const esys::lsm::BoundingBox &bBox);

  /**
   * Returns geometry info.
   */
  esys::lsm::GeometryInfo getGeometryInfo() const;

  /**
   *
   */
  int getNumTimeSteps() const;
  
  /**
   * Sets the number of time steps.
   */
  void setNumTimeSteps(int numTimeSteps);

  /*!
    get time step size
  */
  double getTimeStepSize() const;
  
  /**
   * Sets the time step size.
   */
  void setTimeStepSize(double timeStepSize);

  void setPrecision(int precision){m_precision=precision;};

  bool spatialDomainHasBeenSet() const;
 

  MPI_Comm getMpiComm() const;

  void setMpiComm(MPI_Comm mpiComm);
  

protected:
  /**
   *
   */
  MPI_Comm m_mpiComm;


  /**
   * Prefix of check-point files.
   */
  std::string m_fileNamePrefix;

  /**
   *
   */
  int m_beginTime;

  /**
   *
   */
  int m_endTime;

  /**
   *
   */
  int m_timeInterval;

  /**
   *
   */
  esys::lsm::GeometryInfo m_geoInfo;

  /**
   *
   */
  int m_numTimeSteps;

  /**
   *
   */
  double m_timeStepSize;

  bool m_spatialDomainHasBeenSet;
  
  /*!
    If set, pipe all write operations through master process. Useful on installations
    where only proc. 0 can write to files
  */
  bool m_writeThroughMaster;  
  int m_precision; //! output precision (digits)
};

#endif
