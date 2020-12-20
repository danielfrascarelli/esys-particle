FILE(GLOB MODEL_SRC_LIB "src/*.cpp")
#===========================================================

ADD_LIBRARY(ppa SHARED ${MODEL_SRC_LIB})
TARGET_LINK_LIBRARIES(ppa ${Boost_LIBRARIES} ${PYTHON_LIBRARIES})
INSTALL(TARGETS ppa DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
