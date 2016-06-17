

FILE(
	GLOB
	OPF_FILES
	./*.cc
	./*.h  
	./*.hxx
	)


ADD_LIBRARY(opf_modeler STATIC ${OPF_FILES} )
