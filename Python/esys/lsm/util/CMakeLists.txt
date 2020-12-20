#===========================================================

ADD_LIBRARY(libFoundationPy SHARED Vec3Py.cpp BoundingBoxPy.cpp BoundingSpherePy.cpp RngPy.cpp SetPy.cpp VectorPy.cpp QuaternionPy.cpp)
ADD_LIBRARY(FoundationPy SHARED ExportModuleUtil.cpp)
TARGET_COMPILE_OPTIONS(FoundationPy PUBLIC -O3 -march=native)
TARGET_COMPILE_OPTIONS(libFoundationPy PUBLIC -O3 -march=native)
SET_TARGET_PROPERTIES(libFoundationPy PROPERTIES PREFIX "")
SET_TARGET_PROPERTIES(FoundationPy PROPERTIES PREFIX "")

# configure the InstallInfo file
configure_file (
	"${CMAKE_SOURCE_DIR}/Python/esys/lsm/util/InstallInfo.cmake_tmpl"
	"${CMAKE_CURRENT_BINARY_DIR}/InstallInfo.py"
	)

TARGET_LINK_LIBRARIES(libFoundationPy ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} Foundation )
TARGET_LINK_LIBRARIES(FoundationPy ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} Foundation libFoundationPy)
INSTALL(TARGETS libFoundationPy DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/)
INSTALL(TARGETS FoundationPy DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python3.8/dist-packages/esys/lsm/util/)
INSTALL(FILES ${CMAKE_SOURCE_DIR}/Python/esys/lsm/util/__init__.py DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python3.8/dist-packages/esys/lsm/util)
INSTALL(FILES ${CMAKE_SOURCE_DIR}/Python/esys/lsm/util/FileUtil.py DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python3.8/dist-packages/esys/lsm/util)
#INSTALL(FILES ${CMAKE_SOURCE_DIR}/Python/esys/lsm/util/InstallInfo.py DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python3.8/dist-packages/esys/lsm/util)
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/InstallInfo.py DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python3.8/dist-packages/esys/lsm/util)
INSTALL(FILES ${CMAKE_SOURCE_DIR}/Python/esys/lsm/util/OptParse.py DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python3.8/dist-packages/esys/lsm/util)
INSTALL(FILES ${CMAKE_SOURCE_DIR}/Python/esys/lsm/util/pathSearcher.py DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python3.8/dist-packages/esys/lsm/util)
INSTALL(FILES ${CMAKE_SOURCE_DIR}/Python/esys/lsm/util/process.py DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python3.8/dist-packages/esys/lsm/util)

