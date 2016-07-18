
FILE(
	GLOB
	POP_FILES
	./*.cc
	./*.h  
	./*.hxx
	)

ADD_LIBRARY(pop_modeler STATIC ${POP_FILES} )
