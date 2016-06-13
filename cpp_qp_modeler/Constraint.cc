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
	stream << lb() << " <= ";
	f().print(stream, problem);
	stream << " <= " << ub();

}