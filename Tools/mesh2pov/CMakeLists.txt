FILE(GLOB SRC "*.cpp")
#===========================================================
ADD_EXECUTABLE(mesh2pov ${SRC})
TARGET_LINK_LIBRARIES(mesh2pov Parallel Model ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} )
INSTALL(TARGETS mesh2pov DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
