ADD_EXECUTABLE(pop main.cc)
TARGET_LINK_LIBRARIES(pop pop_modeler)

ADD_EXECUTABLE(pooling pooling.cc)
TARGET_LINK_LIBRARIES(pooling pop_modeler pooling_modeler)
 
