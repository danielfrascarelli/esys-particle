FILE(GLOB SRC "*.cpp")
#===========================================================
ADD_EXECUTABLE(fracextract ${SRC})
TARGET_LINK_LIBRARIES(fracextract Geometry)
INSTALL(TARGETS fracextract DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
