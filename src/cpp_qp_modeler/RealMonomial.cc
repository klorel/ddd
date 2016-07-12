#include "RealMonomial.h"

ComplexMonomialPtr  ComplexMonomial::ZeroPtr = ComplexMonomialPtr(new ComplexMonomial);


ComplexMonomialPtr ComplexMonomial::Build(PosInt id) {
	ComplexMonomialPtr result(new ComplexMonomial);
	result->z().alpha()[id] = 1;
	result->z().degree() = 1;
	return result;
}
ComplexMonomialPtr ComplexMonomial::BuildH(PosInt id) {
	ComplexMonomialPtr result(new ComplexMonomial);
	result->zH().alpha()[id] = 1;
	result->zH().degree() = 1;
	return result;
}

