ADD_SUBDIRECTORY("${CMAKE_SOURCE_DIR}/Python/esys/lsm")
INSTALL(FILES ${CMAKE_SOURCE_DIR}/Python/esys/__init__.py DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python3.8/dist-packages/esys)
