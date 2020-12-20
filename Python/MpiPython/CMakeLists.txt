FILE(GLOB MODEL_SRC_LIB "*.cpp")
#===========================================================

ADD_EXECUTABLE(esysparticle ${MODEL_SRC_LIB})
TARGET_COMPILE_OPTIONS(esysparticle PUBLIC -O3 -march=native)
TARGET_LINK_LIBRARIES(esysparticle Foundation Parallel Geometry Model TmlMessage TmlType ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
INSTALL(TARGETS esysparticle DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
