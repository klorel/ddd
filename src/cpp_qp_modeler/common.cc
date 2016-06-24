#include "common.h"
#include "FunctionComplex.h"
#include "FunctionReal.h"
#include "Problem.h"
#include "Constraint.h"


Index2 get_index(int i, int j){
	return std::make_pair(std::min(i, j), std::max(i, j));
}
bool isZero(Number value){
	return std::fabs(value) < 1e-10;
}
bool isZero(NumberPtr const & value){
	return isZero(*value);
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
Number PosInfinity(){
	return 1e20;
}

Number NegInfinity(){
	return -1e20;
}


std::ostream & operator<<(std::ostream & stream, FunctionReal const & rhs) {
	rhs.print(stream);
	return stream;
}
std::ostream & operator<<(std::ostream & stream, FunctionComplex const & rhs) {
	rhs.print(stream);
	return stream;
}


std::ostream & operator<<(std::ostream & stream, Problem const & rhs) {
	rhs.print(stream);
	return stream;
}
std::ostream & operator<<(std::ostream & stream, Constraint const & rhs) {
	rhs.print(stream);
	return stream;
}


std::ostream & printAlpha(std::ostream & stream, IntVector const & rhs){
	for (size_t i(0); i < rhs.size(); ++i)
		stream << rhs[i];
	return stream;
}