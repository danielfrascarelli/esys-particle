FILE(GLOB SRC "*.cpp")
#===========================================================

message("PYTHONLIBS:" ${PYTHON_LIBRARIES})

ADD_EXECUTABLE(raw2tostress ${SRC})
TARGET_LINK_LIBRARIES(raw2tostress Foundation ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
INSTALL(TARGETS raw2tostress DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
