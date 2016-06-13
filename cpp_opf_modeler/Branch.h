#pragma once

#include "FunctionComplex.h"

class Branch
{
public:
public:
	Branch();
	~Branch();
public:
	ComplexNumbers _attributes;

	template<BranchAttributes attribute> ComplexNumber & get();
	template<BranchAttributes attribute> ComplexNumber const & get()const;
};

template<BranchAttributes attribute> inline ComplexNumber & Branch::get() {
	return _attributes[attribute];
}
template<BranchAttributes attribute> inline ComplexNumber const & Branch::get()const {
	return _attributes[attribute];
}
