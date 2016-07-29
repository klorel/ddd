#include "common.h"

#include "Problem.h"
#include "Constraint.h"


Index2 get_index(int i, int j){
	return std::make_pair(std::min(i, j), std::max(i, j));
}
bool isZero(Number value) {
	return std::fabs(value)< 1e-10;
}
bool isZero(NumberPtr const & value) {
	return isZero(*value);
}

bool isZero(ComplexNumber value) {
	return std::abs(value)< 1e-10;
}

std::string format(Number value){
	std::stringstream buffer;
	if (value > 0){
		buffer << "+";
		if (value != 1)
			buffer << value;
	}
	else if (isZero(value + 1)){
		buffer << "-";
	}
	else{
		buffer << value;
	}
	return buffer.str();
}
std::string format(ComplexNumber const & value) {
	std::stringstream buffer;
	if (!isZero(value)) {
		if(!isZero(value.real()) && !isZero(value.imag()))
			buffer << "("<<format(value.real())<< format(value.imag())<<"i)";
		else if (!isZero(value.real()))
			buffer << format(value.real());
		else if (!isZero(value.imag()))
			buffer << format(value.imag())<<"i";
	}
	else {
		buffer << "0";
	}
	return buffer.str();
}
//
//std::ostream & operator<<(std::ostream & stream, FunctionReal const & rhs) {
//	rhs.print(stream);
//	return stream;
//}
//std::ostream & operator<<(std::ostream & stream, FunctionComplex const & rhs) {
//	rhs.print(stream);
//	return stream;
//}

//
//std::ostream & operator<<(std::ostream & stream, Problem const & rhs) {
//	rhs.print(stream);
//	return stream;
//}
//std::ostream & operator<<(std::ostream & stream, Constraint const & rhs) {
//	rhs.print(stream);
//	return stream;
//}


std::ostream & printAlpha(std::ostream & stream, IntVector const & rhs){
	for (size_t i(0); i < rhs.size(); ++i)
		stream << rhs[i];
	return stream;
}

Number posInfinity() {
	return 1e20;
}
Number negInfinity() {
	return -1e20;
}
