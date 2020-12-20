FILE(GLOB MODEL_SRC_LIB "*.cpp")
#===========================================================

ADD_LIBRARY(TmlComm SHARED ${MODEL_SRC_LIB})
TARGET_LINK_LIBRARIES(TmlComm Foundation ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
INSTALL(TARGETS TmlComm DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
