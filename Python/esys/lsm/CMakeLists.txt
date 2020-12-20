#===========================================================
ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/Python/esys/lsm/util")
ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/Python/esys/lsm/geometry")

FILE(GLOB LSMPY_SRC "*.cpp")
ADD_LIBRARY(LsmPy SHARED ${LSMPY_SRC})
TARGET_COMPILE_OPTIONS(LsmPy PUBLIC -O3 -march=native)
SET_TARGET_PROPERTIES(LsmPy PROPERTIES PREFIX "")

TARGET_LINK_LIBRARIES(LsmPy Parallel Foundation Geometry Model FoundationPy BoostPythonUtil ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})

# Install targets hard-coded for python 3.8
INSTALL(TARGETS LsmPy DESTINATION ${CMAKE_INSTALL_PREFIX}/lib//python3.8/dist-packages/esys/lsm/)
INSTALL(FILES ${CMAKE_SOURCE_DIR}/Python/esys/lsm/__init__.py DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python3.8/dist-packages/esys/lsm/)
INSTALL(FILES ${CMAKE_SOURCE_DIR}/Python/esys/lsm/Logging.py DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python3.8/dist-packages/esys/lsm/)

