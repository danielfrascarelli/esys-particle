#===========================================================
FILE(GLOB GEOPY_SRC "*.cpp")
ADD_LIBRARY(GeometryPy SHARED ${GEOPY_SRC})
SET_TARGET_PROPERTIES(GeometryPy PROPERTIES PREFIX "")

TARGET_LINK_LIBRARIES(GeometryPy FoundationPy Geometry BoostPythonUtil LsmPy ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

# Install targets hard-coded for python 3.8
INSTALL(TARGETS GeometryPy DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python3.8/dist-packages/esys/lsm/geometry)
INSTALL(FILES ${CMAKE_SOURCE_DIR}/Python/esys/lsm/geometry/__init__.py DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python3.8/dist-packages/esys/lsm/geometry)
INSTALL(FILES ${CMAKE_SOURCE_DIR}/Python/esys/lsm/geometry/SingleParticleGrainGen.py DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python3.8/dist-packages/esys/lsm/geometry)
INSTALL(FILES ${CMAKE_SOURCE_DIR}/Python/esys/lsm/geometry/SphericalGrainGen.py DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python3.8/dist-packages/esys/lsm/geometry)
