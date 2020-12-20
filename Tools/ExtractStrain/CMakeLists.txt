FILE(GLOB SRC "*.cpp")
#===========================================================
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/ntable/src)
ADD_EXECUTABLE(strainextract ${SRC})
TARGET_LINK_LIBRARIES(strainextract Foundation ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
INSTALL(TARGETS strainextract DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
