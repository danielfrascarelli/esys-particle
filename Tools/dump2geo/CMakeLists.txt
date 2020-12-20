FILE(GLOB SRC "*.cpp")
#===========================================================
ADD_EXECUTABLE(dump2geo ${SRC})
INSTALL(TARGETS dump2geo DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
