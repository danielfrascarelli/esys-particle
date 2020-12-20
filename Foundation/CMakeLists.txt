FILE(GLOB FOUNDATION_SRC_LIB "*.cpp")
#===========================================================
ADD_LIBRARY(Foundation SHARED ${FOUNDATION_SRC_LIB})
TARGET_LINK_LIBRARIES(Foundation ${Boost_LIBRARIES})
TARGET_COMPILE_OPTIONS(Foundation PUBLIC -O3 -march=native)
#${PYTHON_LIBRARIES})
INSTALL(TARGETS Foundation DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
