#include "Problem.h"
#include "Constraint.h"

Problem::Problem()
{
}


Problem::~Problem()
{
}

std::string  & Problem::name(int key){
	return *_varnames[key];
}

std::string const & Problem::name(int key)const{
	return *_varnames[key];
}
std::string & Problem::ctrname(int key){
	return *_ctrnames[key];
}
std::string const &Problem::ctrname(int key)const{
	return *_ctrnames[key];
}

Constraint& Problem::ctr(int key) {
	return _constraints[key];
}

Constraint const & Problem::ctr(int key) const{
	return _constraints[key];
}

FunctionReal Problem::variable(int i)const{
	FunctionReal result;
	result.add(i, 1);
	return result;
}

FunctionReal Problem::variable(std::string const &name, int i)const{
	FunctionReal result;
	result.add(idvar(name, i), 1, 1);
	return result;
}
FunctionReal Problem::variable(std::string const &name, int i1, int i2)const{
	FunctionReal result;
	result.add(idvar(name, i1, i2), 1, 1);
	return result;
}
FunctionReal Problem::variable(std::string const &name, IntPair const & i)const{
	return variable(name, i.first, i.second);
}
Constraint & Problem::ctr(std::string const &name, int i1){
	return _constraints[idctr(name, i1)];
}
Constraint & Problem::ctr(std::string const &name, int i1, int i2){
	return _constraints[idctr(name, i1, i2)];
}
Constraint & Problem::ctr(std::string const &name, IntPair const & i){
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

void Problem::newvarpool(std::string const & poolname, size_t ids) {
	newpool(_varpools, _varnames, poolname, ids);
}
void Problem::newvarpool(std::string const & poolname, IntPairSet const & ids) {
	newpool(_varpools, _varnames, poolname, ids);
}
void Problem::newvarpool(std::string const & poolname, IntSet const & ids) {
	newpool(_varpools, _varnames, poolname, ids);
}

void Problem::newctrpool(std::string const & poolname, size_t ids){
	newpool(_ctrpools, _ctrnames, poolname, ids);
	_constraints.resize(_ctrnames.size());
}
void Problem::newctrpool(std::string const & poolname, IntSet const & ids){
	newpool(_ctrpools, _ctrnames, poolname, ids);
	_constraints.resize(_ctrnames.size());
}
void Problem::newctrpool(std::string const & poolname, IntPairSet const & ids){
	newpool(_ctrpools, _ctrnames, poolname, ids);
	_constraints.resize(_ctrnames.size());
}


void newpool(Str2Pool & pool, StrPtrVector & names, std::string const & poolname, size_t ids){
	pool[poolname] = IndexedPoolPtr(new IndexedPool1Dense(poolname, names, ids));

}
void newpool(Str2Pool & pool, StrPtrVector & names, std::string const & poolname, IntSet const & ids){
	pool[poolname] = IndexedPoolPtr(new IndexedPool1Sparse(poolname, names, ids));

}
void newpool(Str2Pool & pool, StrPtrVector & names, std::string const & poolname, IntPairSet const & ids){
	pool[poolname] = IndexedPoolPtr(new IndexedPool2Sparse(poolname, names, ids));
}

//void Problem::ivariablePool(std::string const & poolname, size_t size) {
//	variablePool(poolname + "_real", size);
//	variablePool(poolname + "_imag", size);
//}
int Problem::idvar(std::string const & name, int i1)const{
	return id(_varpools, name, i1);
}

int Problem::idvar(std::string const & name, int i1, int i2)const{
	return id(_varpools, name, i1, i2);
}
int Problem::idctr(std::string const & name, int i1)const{
	return id(_ctrpools, name, i1);
}

int Problem::idctr(std::string const & name, int i1, int i2)const{
	return id(_ctrpools, name, i1, i2);
}

int id(Str2Pool const & pool, std::string const & name, int i){
	Str2Pool::const_iterator it(pool.find(name));
	if (it == pool.end()){
		throw std::invalid_argument("pool name was not created");
		return -1;
	}
	else{
		return it->second->id(i);
	}
}
int id(Str2Pool const & pool, std::string const & name, int i1, int i2){
	Str2Pool::const_iterator it(pool.find(name));
	if (it == pool.end()){
		throw std::invalid_argument("pool name was not created");
		return -1;
	}
	else{
		return it->second->id(i1, i2);
	}
}

void Problem::add(Constraint const &rhs){
	_constraints.push_back(rhs);
}


FunctionReal & Problem::minimize() {
	return _minimize;
}

FunctionReal const & Problem::minimize()const {
	return _minimize;
}

void operator<<(Problem & problem, FunctionReal const & rhs) {
	problem.minimize() = rhs.clone();
}

void Problem::print(std::ostream &stream)const {
	stream << "minimize ";
	minimize().print(stream, *this);
	stream << std::endl;
	for (size_t i(0); i < _constraints.size(); ++i) {
		_constraints[i].print(stream<<*_ctrnames[i]<<" : ", *this);
		stream << std::endl;
	}
}

void Problem::addSparsityPattern(SparsityPattern & sparsityPattern)const {
	sparsityPattern.resize(nvars());
	//minimize().addSparsityPattern(sparsityPattern);
	for (auto const & ctr : _constraints)
		ctr.f().addSparsityPattern(sparsityPattern);
}

void Problem::addSupport(SparsityPattern & sparsityPattern)const {
	sparsityPattern.assign(nctrs()+1, IntSet());
	minimize().addSupport(sparsityPattern[nctrs()]);	
	for (int i(0); i < _constraints.size(); ++i){
		_constraints[i].f().addSupport(sparsityPattern[i]);
		++i;
	}
}


int Problem::nvars()const{
	return static_cast<int>(_varnames.size());
}
int Problem::nctrs()const{
	return static_cast<int>(_ctrnames.size());
}

void Problem::removeInequality() {
	IntSet slacks;
	
	for (int i(0); i < _constraints.size(); ++i) {
		if(_constraints[i].sense()==RNG)
			throw std::invalid_argument("in Problem::removeInequality, constraint is RNG");
		if (_constraints[i].sense() != EQ) {
			slacks.insert(i);
		}
	}
	newvarpool("slack", slacks);
	for (auto const ctr : slacks) {
		FunctionReal s(variable("slack", ctr));
		if (_constraints[ctr].sense() == GEQ) {
			_constraints[ctr].f() -= s*s;
			_constraints[ctr].ub() = _constraints[ctr].lb();
		}
		else {
			_constraints[ctr].f() += s*s;
			_constraints[ctr].lb() = _constraints[ctr].ub();
		}
	}

}