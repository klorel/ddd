#include "Constraint.h"


Constraint::Constraint() :_lb(new Number(negInfinity())), _ub(new Number(posInfinity())), _f()
{
}


Constraint::~Constraint()
{
}

Number const & Constraint::lb()const{
	return *_lb;
}
Number const & Constraint::ub()const{
	return *_ub;
}
FunctionReal const & Constraint::f()const{
	return _f;
}

Number & Constraint::lb(){
	return *_lb;
}
Number & Constraint::ub(){
	return *_ub;
}
FunctionReal & Constraint::f(){
	return _f;
}
void Constraint::print(std::ostream & stream)const {
	stream << lb() << " <= " << f() << " <= " << ub();
}
void Constraint::print(std::ostream & stream, Problem const & problem)const {
	if (lb() == ub())
		stream << lb() << " = ";
	else if (lb() > negInfinity())
		stream << lb() << " <= ";
	f().print(stream, problem);
	if (lb() != ub() && ub() < posInfinity())
		stream << " <= " << ub();

}
Sense Constraint::sense()const {
	if (lb() == ub())
		return EQ;
	else if (lb() > negInfinity() && ub() >= posInfinity())
		return GEQ;
	else if (ub() < posInfinity() && lb() <= negInfinity())
		return LEQ;
	else
		return RNG;
}