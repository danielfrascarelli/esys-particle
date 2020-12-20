FILE(GLOB GEOMETRY_SRC_LIB "*.cpp")
#===========================================================
ADD_LIBRARY(Geometry SHARED ${GEOMETRY_SRC_LIB})
TARGET_COMPILE_OPTIONS(Geometry PUBLIC -O3 -march=native)
TARGET_LINK_LIBRARIES(Geometry Foundation ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
INSTALL(TARGETS Geometry DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
