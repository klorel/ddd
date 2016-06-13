#pragma once

#include "common.h"

FunctionReal operator+(FunctionReal const &);
FunctionReal operator-(FunctionReal const &);

FunctionReal operator+(FunctionReal const &, FunctionReal const &);
FunctionReal operator-(FunctionReal const &, FunctionReal const &);
FunctionReal operator*(FunctionReal const &, FunctionReal const &);
FunctionReal operator/(FunctionReal const &, FunctionReal const &);




Constraint operator<=(FunctionReal const & lhs, Number rhs);
Constraint operator==(FunctionReal const & lhs, Number rhs);
Constraint operator>=(FunctionReal const & lhs, Number rhs);


Constraint operator<=(Number lhs, FunctionReal const & rhs);
Constraint operator==(Number lhs, FunctionReal const & rhs);
Constraint operator>=(Number lhs, FunctionReal const & rhs);



FunctionComplex operator+(FunctionComplex const &);
FunctionComplex operator-(FunctionComplex const &);

FunctionComplex operator+(FunctionComplex const &, FunctionComplex const &);
FunctionComplex operator-(FunctionComplex const &, FunctionComplex const &);
FunctionComplex operator*(FunctionComplex const &, FunctionComplex const &);
FunctionComplex operator/(FunctionComplex const &, FunctionComplex const &);