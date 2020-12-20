FILE(GLOB MODEL_SRC_LIB "*.cpp")
#===========================================================

ADD_LIBRARY(TmlMessage SHARED ${MODEL_SRC_LIB})
TARGET_LINK_LIBRARIES(TmlMessage ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
INSTALL(TARGETS TmlMessage DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
