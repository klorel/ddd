FILE(
	GLOB
	EXE_FILES
	./*.cc
	./*.h  
	./*.hxx
	)


ADD_EXECUTABLE(qp ${EXE_FILES})
TARGET_LINK_LIBRARIES(qp qp_modeler opf_modeler)
