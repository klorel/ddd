
FILE(
	GLOB
	POOLING_FILES
	./*.cc
	./*.h  
	./*.hxx
	)

ADD_LIBRARY(pooling_modeler STATIC ${POOLING_FILES} )
