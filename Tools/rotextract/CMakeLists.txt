FILE(GLOB SRC "*.cpp")
#===========================================================
ADD_EXECUTABLE(rotextract ${SRC})
INSTALL(TARGETS rotextract DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
