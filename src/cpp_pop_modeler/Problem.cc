#include "Problem.h"
#include "Constraint.h"


void Problem::get_all_monomial(ComplexMonomialPtr2Int & output)const {
	output.clear();
	minimize().get_all_monomial(output);
	for (auto const & ctr : _constraints) {
		ctr.f().get_all_monomial(output);
	}
	int i(-1);
	for (auto & kvp : output) {
		kvp.second = ++i;
		//std::cout << kvp.first << " | " << kvp.first->z().degree() << ", " << kvp.first->zH().degree() << std::endl;
	}
}

Problem::Problem()
{
}


Problem::~Problem()
{
}

std::string  & Problem::name(int key) {
	return *_varnames[key];
}

std::string const & Problem::name(int key)const {
	return *_varnames[key];
}
std::string & Problem::ctrname(int key) {
	return *_ctrnames[key];
}
std::string const &Problem::ctrname(int key)const {
	return *_ctrnames[key];
}

Constraint& Problem::ctr(int key) {
	return _constraints[key];
}

Constraint const & Problem::ctr(int key) const {
	return _constraints[key];
}

ComplexPolynomial Problem::variable(int i)const {
	return ComplexPolynomial::Build(i);
}
//
//ComplexPolynomial Problem::variable(std::string const &name, int i)const {
//	return ComplexPolynomial::Build(idvar(name, i));
//}
//ComplexPolynomial Problem::variable(std::string const &name, int i1, int i2)const {
//	return ComplexPolynomial::Build(idvar(name, i1, i2));
//}
//ComplexPolynomial Problem::variable(std::string const &name, IntPair const & i)const {
//	return variable(name, i.first, i.second);
//}
Constraint & Problem::ctr(std::string const &name, int i1) {
	return _constraints[idctr(name, i1)];
}
Constraint & Problem::ctr(std::string const &name, int i1, int i2) {
	return _constraints[idctr(name, i1, i2)];
}
Constraint & Problem::ctr(std::string const &name, IntPair const & i) {
	return ctr(name, i.first, i.second);
}

//FunctionComplex Problem::ivariable(std::string const &name, int i)const {
//	return variable(name + "_real", name + "_imag", i);
//}

//FunctionComplex Problem::variable(std::string const &name1, std::string const &name2, int i1)const{
//	FunctionComplex result;
//	result.real().add(id(name1, i1), 1, 1);
//	result.imag().add(id(name2, i1), 1, 1);
//	return result;
//}
//
//FunctionComplex Problem::variable(std::string const &name1, int i1, std::string const &name2, int i2)const{
//	FunctionComplex result;
//	result.real().add(id(name1, i1), 1, 1);
//	result.imag().add(id(name2, i2), 1, 1);
//	return result;
//}
//
//FunctionComplex Problem::variable(std::string const &name1, int i1, int i2)const{
//	FunctionComplex result;
//	result.real().add(id(name1, i1), 1, 1);
//	result.imag().add(id(name1, i2), 1, 1);
//	return result;
//}


IndexedPool const & Problem::varpool(std::string const & name)const {
	Str2Pool::const_iterator it(_varpools.find(name));
	if (it == _varpools.end()) {
		throw std::invalid_argument("var pool name was not found");
	}
	else {
		return *it->second;
	}
}

IndexedPool const & Problem::ctrpool(std::string const &name)const {
	Str2Pool::const_iterator it(_ctrpools.find(name));
	if (it == _ctrpools.end()) {
		throw std::invalid_argument("ctr pool name was not found");
	}
	else {
		return *it->second;
	}
}

IndexedPool const & Problem::newvarpool(std::string const & poolname, size_t ids) {
	return *newpool(_varpools, _varnames, poolname, ids);
}

IndexedPool const & Problem::newvarpool(std::string const & poolname, IntPairSet const & ids) {
	return *newpool(_varpools, _varnames, poolname, ids);
}

IndexedPool const & Problem::newvarpool(std::string const & poolname, IntSet const & ids) {
	return *newpool(_varpools, _varnames, poolname, ids);
}

IndexedPool const & Problem::newctrpool(std::string const & poolname, size_t ids) {
	IndexedPool const & result = *newpool(_ctrpools, _ctrnames, poolname, ids);
	_constraints.resize(_ctrnames.size());
	return result;
}

IndexedPool const & Problem::newctrpool(std::string const & poolname, IntSet const & ids) {
	IndexedPool const & result = *newpool(_ctrpools, _ctrnames, poolname, ids);
	_constraints.resize(_ctrnames.size());
	return result;
}

IndexedPool const & Problem::newctrpool(std::string const & poolname, IntPairSet const & ids) {
	IndexedPool const & result = *newpool(_ctrpools, _ctrnames, poolname, ids);
	_constraints.resize(_ctrnames.size());
	return result;
}


IndexedPoolPtr newpool(Str2Pool & pool, StrPtrVector & names, std::string const & poolname, size_t ids) {
	auto result = pool.insert({ poolname, IndexedPoolPtr(new IndexedPool1Dense(poolname, names, ids)) });
	if (!result.second)
		throw std::invalid_argument("var pool name was not created");
	return result.first->second;
}

IndexedPoolPtr newpool(Str2Pool & pool, StrPtrVector & names, std::string const & poolname, IntSet const & ids) {
	auto result = pool.insert({ poolname, IndexedPoolPtr(new IndexedPool1Sparse(poolname, names, ids)) });
	if (!result.second)
		throw std::invalid_argument("var pool name was not created");
	return result.first->second;
}

IndexedPoolPtr newpool(Str2Pool & pool, StrPtrVector & names, std::string const & poolname, IntPairSet const & ids) {
	auto result = pool.insert({ poolname, IndexedPoolPtr(new IndexedPool2Sparse(poolname, names, ids)) });
	if (!result.second)
		throw std::invalid_argument("var pool name was not created");
	return result.first->second;
}

//void Problem::ivariablePool(std::string const & poolname, size_t size) {
//	variablePool(poolname + "_real", size);
//	variablePool(poolname + "_imag", size);
//}

int Problem::idvar(std::string const & name, int i1)const {
	return id(_varpools, name, i1);
}

int Problem::idvar(std::string const & name, int i1, int i2)const {
	return id(_varpools, name, i1, i2);
}

int Problem::idctr(std::string const & name, int i1)const {
	return id(_ctrpools, name, i1);
}

int Problem::idctr(std::string const & name, int i1, int i2)const {
	return id(_ctrpools, name, i1, i2);
}

int id(Str2Pool const & pool, std::string const & name, int i) {
	Str2Pool::const_iterator it(pool.find(name));
	if (it == pool.end()) {
		throw std::invalid_argument("pool name was not created");
		return -1;
	}
	else {
		return it->second->id(i);
	}
}
int id(Str2Pool const & pool, std::string const & name, int i1, int i2) {
	Str2Pool::const_iterator it(pool.find(name));
	if (it == pool.end()) {
		throw std::invalid_argument("pool name was not created");
		return -1;
	}
	else {
		return it->second->id(i1, i2);
	}
}

void Problem::add(Constraint const &rhs) {
	_constraints.push_back(rhs);
	_ctrnames.push_back(StrPtr(new std::string(Str("ctr_", nctrs()))));
}


ComplexPolynomial & Problem::minimize() {
	return _minimize;
}

ComplexPolynomial const & Problem::minimize()const {
	return _minimize;
}

int Problem::nvars()const {
	return static_cast<int>(_varnames.size());
}

int Problem::nctrs()const {
	return static_cast<int>(_ctrnames.size());
}
void Problem::clear() {
	_ctrpools.clear();
	_ctrnames.clear();
	_constraints.clear();

	_varpools.clear();
	_varnames.clear();
	_minimize = 0;
}

std::ostream & operator<<(std::ostream & lhs, Problem const & rhs) {
	rhs.print(lhs);
	return lhs;
}

void operator<<(Problem & problem, ComplexPolynomial const & rhs) {
	problem.minimize() = rhs;
}

void Problem::print(std::ostream &stream)const {
	stream << "minimize ";
	minimize().print(stream, *this);
	//minimize().print(stream);
	stream << std::endl;
	for (size_t i(0); i < _constraints.size(); ++i) {
		//_constraints[i].print(stream << *_ctrnames[i] << " : ");
		_constraints[i].print(stream << *_ctrnames[i] << " : ", *this);
		stream << std::endl;
	}
}

//void Problem::addSparsityPattern(SparsityPattern & sparsityPattern)const {
//	sparsityPattern.resize(nvars());
//	//minimize().addSparsityPattern(sparsityPattern);
//	for (auto const & ctr : _constraints)
//		ctr.f().addSparsityPattern(sparsityPattern);
//}

//void Problem::addSupport(SparsityPattern & sparsityPattern)const {
//	sparsityPattern.assign(nctrs() + 1, IntSet());
//	minimize().addSupport(sparsityPattern[nctrs()]);
//	for (int i(0); i < _constraints.size(); ++i) {
//		_constraints[i].f().addSupport(sparsityPattern[i]);
//		++i;
//	}
//}
//

//void Problem::removeInequality() {
//	IntSet slacks;
//
//	for (int i(0); i < _constraints.size(); ++i) {
//		if (_constraints[i].sense() == RNG)
//			throw std::invalid_argument("in Problem::removeInequality, constraint is RNG");
//		if (_constraints[i].sense() != EQ) {
//			slacks.insert(i);
//		}
//	}
//	newvarpool("slack", slacks);
//	for (auto const ctr : slacks) {
//		FunctionReal s(variable("slack", ctr));
//		if (_constraints[ctr].sense() == GEQ) {
//			_constraints[ctr].f() -= s*s;
//			_constraints[ctr].ub() = _constraints[ctr].lb();
//		}
//		else {
//			_constraints[ctr].f() += s*s;
//			_constraints[ctr].lb() = _constraints[ctr].ub();
//		}
//	}
//
//}

//void Problem::amplExport(std::string const & rootname)const {
//	{
//		std::ofstream file((rootname + ".mod").c_str());
//		std::ostream & stream(file);
//		stream << "param N;" << std::endl;
//		stream << "param M;" << std::endl;
//
//		stream << "set IJ dimen 3;" << std::endl;
//		stream << "set I dimen 2;" << std::endl;
//
//		stream << "param Q{IJ};" << std::endl;
//		stream << "param L{I};" << std::endl;
//		stream << "param LB{1..M} default -1e20;" << std::endl;
//		stream << "param UB{1..M} default +1e20;" << std::endl;
//		stream << "set EQ  := {m in 1..M:abs(LB[m]-UB[m])<1e-10};" << std::endl;
//		stream << "set LEQ := {m in 1..M:LB[m]<=-1e20 && UB[m]< +1e20};" << std::endl;
//		stream << "set GEQ := {m in 1..M:LB[m]> -1e20 && UB[m]>=+1e20};" << std::endl;
//		stream << "set RNG := {m in 1..M:abs(LB[m]-UB[m])>1e-10 && LB[m]> -1e20 && UB[m]<+1e20};" << std::endl;
//
//		stream << "var x{0..N-1};" << std::endl;
//
//		stream << "minimize obj:sum{(0,i,j) in IJ}(Q[0,i,j]*x[i]*x[j])+sum{(0,i) in I}(L[0,i]*x[i]);" << std::endl;
//		stream << "subject to ctr_eq {m in  EQ}:sum{(m,i,j) in IJ}(Q[m,i,j]*x[i]*x[j])+sum{(m,i) in I}(L[m,i]*x[i])=LB[m];" << std::endl;
//		stream << "subject to ctr_geq{m in GEQ}:sum{(m,i,j) in IJ}(Q[m,i,j]*x[i]*x[j])+sum{(m,i) in I}(L[m,i]*x[i])>=LB[m];" << std::endl;
//		stream << "subject to ctr_leq{m in LEQ}:sum{(m,i,j) in IJ}(Q[m,i,j]*x[i]*x[j])+sum{(m,i) in I}(L[m,i]*x[i])<=UB[m];" << std::endl;
//		stream << "subject to ctr_rng{m in RNG}:LB[m]<=sum{(m,i,j) in IJ}(Q[m,i,j]*x[i]*x[j])+sum{(m,i) in I}(L[m,i]*x[i])<=UB[m];" << std::endl;
//		file.close();
//	}
//	{
//		std::ofstream file((rootname + ".dat").c_str());
//		file << "param N := " << nvars() << ";" << std::endl;
//		file << "param M := " << nctrs() << ";" << std::endl;
//		file << "param: IJ:Q := include " << rootname << "_q.txt;" << std::endl;
//		file << "param: I:L := include " << rootname << "_l.txt;" << std::endl;
//		file << "param LB := include " << rootname << "_lb.txt;" << std::endl;
//		file << "param UB := include " << rootname << "_ub.txt;" << std::endl;
//		file.close();
//	}
//	{
//		std::ofstream q((rootname + "_q.txt").c_str());
//		std::ofstream l((rootname + "_l.txt").c_str());
//		std::ofstream lb((rootname + "_lb.txt").c_str());
//		std::ofstream ub((rootname + "_ub.txt").c_str());
//
//		for (auto const & term : minimize().quadratic()) {
//			q << std::setw(8) << 0;
//			q << std::setw(8) << term.first.first;
//			q << std::setw(8) << term.first.second;
//			q << std::setw(30) << std::setprecision(20) << term.second;
//			q << std::endl;
//		}
//		for (auto const & term : minimize().linear()) {
//			l << std::setw(8) << 0;
//			l << std::setw(8) << term.first;
//			l << std::setw(30) << std::setprecision(20) << term.second;
//			l << std::endl;
//		}
//		for (size_t i(0); i < nctrs(); ++i) {
//			Constraint const & ctr(_constraints[i]);
//			FunctionReal const & f(ctr.f());
//			for (auto const & term : f.quadratic()) {
//				q << std::setw(8) << i + 1;
//				q << std::setw(8) << term.first.first;
//				q << std::setw(8) << term.first.second;
//				q << std::setw(30) << std::setprecision(20) << term.second;
//				q << std::endl;
//			}
//			for (auto const & term : f.linear()) {
//				l << std::setw(8) << i + 1;
//				l << std::setw(8) << term.first;
//				l << std::setw(30) << std::setprecision(20) << term.second;
//				l << std::endl;
//			}
//			if (ctr.lb() > -1e20) {
//				lb << std::setw(8) << i + 1;
//				lb << std::setw(30) << std::setprecision(20) << ctr.lb() - f.constant();
//				lb << std::endl;
//			}
//			if (ctr.ub() < +1e20) {
//				ub << std::setw(8) << i + 1;
//				ub << std::setw(30) << std::setprecision(20) << ctr.ub() - f.constant();
//				ub << std::endl;
//			}
//		}
//		q.close();
//		l.close();
//		lb.close();
//		ub.close();
//	}
//	{
//		std::ofstream file((rootname + ".run").c_str());
//		file << "reset;" << std::endl;
//		file << "model " << rootname << ".mod;" << std::endl;
//		file << "data " << rootname << ".dat;" << std::endl;
//		file << "option solver couenne;" << std::endl;
//		file << "solve;" << std::endl;
//	}
//}

