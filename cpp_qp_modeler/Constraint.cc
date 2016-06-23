#include "Constraint.h"


Constraint::Constraint() :_lb(new Number(NegInfinity())), _ub(new Number(PosInfinity())), _f()
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
	else if (lb() > -1e20)
		stream << lb() << " <= ";
	f().print(stream, problem);
	if (lb() != ub() && ub() < 1e20)
		stream << " <= " << ub();

}
void Constraint::addSparsityPattern(SparsityPattern & sparsityPattern)const {
	f().addSparsityPattern(sparsityPattern);
}