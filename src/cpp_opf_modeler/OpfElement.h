#pragma once

#include "common_opf.h"
#include "FunctionComplex.h"

template<class _NAME_T_>
class ElementT
{
public:
public:
	ElementT();
	~ElementT();
public:
	ComplexNumbers _attributes;

	template<enum _NAME_T_::Attributes attribute> ComplexNumber & get();
	template<enum _NAME_T_::Attributes attribute> ComplexNumber const & get()const;
};

template<typename _NAME_T_>
inline ElementT<_NAME_T_>::ElementT():_attributes(_NAME_T_::NUM) {

}
template<typename _NAME_T_>
inline ElementT<_NAME_T_>::~ElementT() {

}
template<typename _NAME_T_>
template<enum _NAME_T_::Attributes attribute> inline ComplexNumber & ElementT<_NAME_T_>::get() {
	return _attributes[attribute];
}
template<typename _NAME_T_>
template<enum _NAME_T_::Attributes attribute> inline ComplexNumber const & ElementT<_NAME_T_>::get() const {
	return _attributes[attribute];
}