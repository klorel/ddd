

FILE(
	GLOB
	QP_FILES
	./*.cc
	./*.h  
	./*.hxx
	)

ADD_LIBRARY(qp_modeler STATIC ${QP_FILES} )
