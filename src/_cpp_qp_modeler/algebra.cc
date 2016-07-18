#include "algebra.h"
#include "FunctionReal.h"
#include "FunctionComplex.h"
#include "Constraint.h"

#define __DEBUG_ALGEBRA__ 0
#if __DEBUG_ALGEBRA__
#define __DEBUG_ALGEBRA_CALL__ std::cout <<"CALL "<<__FUNCTION__ <<"("<<lhs<<", "<<rhs<<")"<< std::endl
#define __DEBUG_ALGEBRA_END__ std::cout <<"END  "<<__FUNCTION__ <<"("<<lhs<<", "<<rhs<<") : "<<result<< std::endl

#else
#define __DEBUG_ALGEBRA_CALL__
#define __DEBUG_ALGEBRA_END__
#endif


FunctionReal operator+(FunctionReal const & rhs) {
	return 0+rhs;
}

FunctionReal operator-(FunctionReal const & rhs) {
	return 0 - rhs;
}


FunctionReal operator+(FunctionReal const & lhs, FunctionReal const & rhs) {
	__DEBUG_ALGEBRA_CALL__;
	FunctionReal result;
	result += lhs;
	result += rhs;
	__DEBUG_ALGEBRA_END__;
	return result;
}
FunctionReal operator-(FunctionReal const & lhs, FunctionReal const & rhs) {
	__DEBUG_ALGEBRA_CALL__;
	FunctionReal result;
	result += lhs;
	result -= rhs;
	__DEBUG_ALGEBRA_END__;
	return result;
}
FunctionReal operator*(FunctionReal const & lhs, FunctionReal const & rhs) {
	__DEBUG_ALGEBRA_CALL__;
	FunctionReal result;
	result += lhs;
	result *= rhs;
	__DEBUG_ALGEBRA_END__;
	return result;
}
FunctionReal operator/(FunctionReal const & lhs, FunctionReal const & rhs) {
	__DEBUG_ALGEBRA_CALL__;
	FunctionReal result;
	result += lhs;
	result /= rhs;
	__DEBUG_ALGEBRA_END__;
	return result;
}

Constraint operator<=(FunctionReal const & lhs, Number rhs) {
	Constraint result;
	result.f() = lhs;
	result.ub() = rhs;
	return result;

}
Constraint operator==(FunctionReal const & lhs, Number rhs) {
	Constraint result;
	result.f() = lhs;
	result.ub() = rhs;
	result.lb() = rhs;
	return result;

}
Constraint operator>=(FunctionReal const & lhs, Number rhs) {
	Constraint result;
	result.f() = lhs;
	result.lb() = rhs;
	return result;

}


Constraint operator<=(Number lhs, FunctionReal const & rhs) {
	Constraint result;
	result.f() = rhs;
	result.lb() = lhs;
	return result;

}
Constraint operator==(Number lhs, FunctionReal const & rhs) {
	Constraint result;
	result.f() = rhs;
	result.ub() = lhs;
	result.lb() = lhs;
	return result;

}
Constraint operator>=(Number lhs, FunctionReal const & rhs) {
	Constraint result;
	result.f() = rhs;
	result.ub() = lhs;
	return result;

}

FunctionComplex operator+(FunctionComplex const & rhs) {
	return 0 + rhs;
}
FunctionComplex operator-(FunctionComplex const & rhs) {
	return 0 - rhs;
}

FunctionComplex operator+(FunctionComplex const & lhs, FunctionComplex const & rhs) {
	__DEBUG_ALGEBRA_CALL__;
	FunctionComplex result;
	result += lhs;
	result += rhs;
	__DEBUG_ALGEBRA_END__;
	return result;
}
FunctionComplex operator-(FunctionComplex const & lhs, FunctionComplex const & rhs) {
	__DEBUG_ALGEBRA_CALL__;
	FunctionComplex result;
	result += lhs;
	result -= rhs;
	__DEBUG_ALGEBRA_END__;
	return result;
}
FunctionComplex operator*(FunctionComplex const & lhs, FunctionComplex const & rhs) {
	__DEBUG_ALGEBRA_CALL__;
	FunctionComplex result;
	result += lhs;
	result *= rhs;
	__DEBUG_ALGEBRA_END__;
	return result;
}
FunctionComplex operator/(FunctionComplex const & lhs, FunctionComplex const & rhs) {
	__DEBUG_ALGEBRA_CALL__;
	FunctionComplex result;	
	result += lhs;
	result /= rhs;
	__DEBUG_ALGEBRA_END__;
	return result;
}

#undef __DEBUG_ALGEBRA__