FILE(GLOB SRC "*.cpp")
if(${Boost_FOUND})
        message(${Boost_LIBRARIES})
endif(${Boost_FOUND})

#===========================================================
ADD_EXECUTABLE(dump2vtk ${SRC})
target_link_libraries(dump2vtk ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
INSTALL(TARGETS dump2vtk DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
