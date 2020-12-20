FILE(GLOB SRC "*.cpp")
#===========================================================
ADD_EXECUTABLE(fcconv ${SRC})
TARGET_LINK_LIBRARIES(fcconv Foundation ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
INSTALL(TARGETS fcconv DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
