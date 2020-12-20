FILE(GLOB SRC_LIB "*.cpp")
#===========================================================

ADD_LIBRARY(BoostPythonUtil SHARED ${SRC_LIB})
TARGET_COMPILE_OPTIONS(BoostPythonUtil PUBLIC -O3 -march=native)
TARGET_LINK_LIBRARIES(BoostPythonUtil ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
INSTALL(TARGETS BoostPythonUtil DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
