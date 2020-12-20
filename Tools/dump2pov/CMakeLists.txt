FILE(GLOB SRC "*.cpp")
#===========================================================
ADD_EXECUTABLE(dump2pov ${SRC})
TARGET_LINK_LIBRARIES(dump2pov Parallel Model ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
INSTALL(TARGETS dump2pov DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
