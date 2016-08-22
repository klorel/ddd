/*
* Pooling.cpp
*
*  Created on: 4 sept. 2012
*      Author: manuel
*/

#include "Pooling.h"
#include "PolynomialOptimizationProblem.h"

Pooling::Pooling(size_t n, size_t k) {
	allocate(n, k);
}
void Pooling::allocate(size_t n, size_t k) {
	_n = n;
	_k = k;
	_nPlus.assign(n, IntSet());
	_nMoins.assign(n, IntSet());
	_b.assign(n, NumberVector(2, 0));
	_q.assign(n, NumberVector(k, 0));
	_c.assign(n, NumberVector(n, 0));
}
Pooling::~Pooling() {
}

Number Pooling::c(IntPair const & ij) const{
	return c(ij.first,ij.second);
}
Number Pooling::q(IntPair const &ij) const{
	return q(ij.first, ij.second);
}
Number Pooling::c(int i, int j) const {
	return _c[i][j];
}
Number Pooling::q(int i, int j) const {
	return _q[i][j];
}

Number Pooling::bL(int i) const {
	return _b[i][0];
}
Number Pooling::bU(int i) const {
	return _b[i][1];

}
void Pooling::newArc(int i, int j) {
	_nPlus[i].insert(j);
	_nMoins[j].insert(i);
}

void Pooling::Structure::build(Pooling const & pooling){
	S.clear();
	I.clear();
	T.clear();
	SI.clear();
	ST.clear();
	IT.clear();
	SIT.clear();
	IK.clear();
	TK.clear();
	bL.clear();
	bU.clear();
	for (int n(0); n < pooling._n; ++n) {
		if (pooling._nPlus[n].empty())
			T.insert(n);
		else if (pooling._nMoins[n].empty())
			S.insert(n);
		else
			I.insert(n);
	}

	std::vector<bool> isI(pooling._n, false);
	for (auto const& i : I){
		isI[i].flip();
		for (auto const & s : pooling._nMoins[i])
			SI.insert({ s, i });
		for (auto const &t : pooling._nPlus[i])
			IT.insert({ i, t });
		for (auto const & s : pooling._nMoins[i]){
			for (auto const &t : pooling._nPlus[i]){
				SIT.insert({ { s, i, }, t });
			}
		}
	}
	for (auto const & t : T){
		for (auto const & s : pooling._nMoins[t])
			if (!isI[s])
				ST.insert({ s, t });
	}
	for (int k(0); k < pooling._k; ++k) {
		for (auto const & i : I){
			if (pooling.q(i,k) != 0)
				IK.insert({ i, k });
		}
		for (auto const & t : T){
			if (pooling.q(t, k) != 0)
				TK.insert({ t, k });
		}
	}
	for (int n(0); n < pooling._n; ++n){
		if (pooling.bL(n) != 0){
			bL.insert(n);
		}
		if (pooling.bU(n) != 0){
			bU.insert(n);
		}
	}
}
void Pooling::get(Structure & result) const{
	result.build(*this);
}

void Pooling::read(std::string const & fileName) {
	std::ifstream file(fileName.c_str());
	if (file.good()) {
		int i(-1), j(-1);
		file >> i;
		file >> j;
		std::string header;
		allocate(i, j);
		while (file >> header) {
			if (header == "c") {
				file >> i;
				file >> j;
				file >> _c[i][j];
				newArc(i, j);
			}
			else if (header == "q") {
				file >> i;
				file >> j;
				file >> _q[i][j];
			}
			else if (header == "bu") {
				file >> i;
				file >> _b[i][1];
			}
			else if (header == "bl") {
				file >> i;
				file >> _b[i][0];
			}
		}
	}
	else {
		std::cout << "unable to open " << fileName << "\n";
	}
}
void Pooling::out(std::ostream & stream) const {
	stream << _n << " " << _k << "\n";
	for (size_t i(0); i < _n; ++i) {
		for (auto const & j : _nPlus[i])
			stream << "c " << i << " " << j << " " << _c[i][j] << "\n";
	}
	for (size_t i(0); i < _n; ++i)
		for (size_t k(0); k < _k; ++k)
			//			if (!isZero(_q.get(i, k)))
			stream << "q " << i << " " << k << " " << _q[i][k] << "\n";

	for (size_t i(0); i < _n; ++i)
		if (!isZero(_b[i][0]))
			stream << "bl " << i << " " << _b[i][0] << "\n";
	for (size_t i(0); i < _n; ++i)
		if (!isZero(_b[i][1]))
			stream << "bu " << i << " " << _b[i][1] << "\n";
}
void Pooling::outAmpl(std::ostream & stream) const {
	stream << "data;" << std::endl;
	stream << "param nNodes      := " << _n << ";" << std::endl;
	stream << "param nAttributes := " << _k << ";" << std::endl;
	stream << "param: Edges: c := " << std::endl;
	for (size_t i(0); i < _n; ++i) {
		for (auto const & j : _nPlus[i]) {
			stream << std::setw(4) << i + 1;
			stream << std::setw(4) << j + 1;
			stream << std::setw(25) << _c[i][j];
			stream << std::endl;
		}
	}
	stream << ";" << std::endl;
	stream << "param q := " << std::endl;
	for (size_t i(0); i < _n; ++i) {
		for (size_t k(0); k < _k; ++k) {
			if (!isZero(_q[i][k])) {
				stream << std::setw(4) << i + 1;
				stream << std::setw(4) << k + 1;
				stream << std::setw(25) << _q[i][k];
				stream << std::endl;
			}
		}
	}
	stream << ";" << std::endl;

	stream << "param bL := " << std::endl;
	for (size_t i(0); i < _n; ++i) {
		if (!isZero(_b[i][0])) {
			stream << std::setw(4) << i + 1;
			stream << std::setw(25) << _b[i][0];
			stream << std::endl;
		}
	}
	stream << ";" << std::endl;
	stream << "param bU := " << std::endl;
	for (size_t i(0); i < _n; ++i) {
		if (!isZero(_b[i][1])) {
			stream << std::setw(4) << i + 1;
			stream << std::setw(25) << _b[i][1];
			stream << std::endl;
		}
	}
	stream << ";" << std::endl;

}

size_t Pooling::n() const {
	return _n;
}

size_t Pooling::k() const {
	return _k;
}

IntSet const & Pooling::nPlus(int i) const {
	return _nPlus[i];
}
IntSet const & Pooling::nMoins(int i) const {
	return _nMoins[i];
}


void Pooling::pqFormulation(PolynomialOptimizationProblem & result)const{
	Structure structure;
	get(structure);

	IndexedPool const & x =  result.newvarpool("x", structure.ST);
	IndexedPool const & y = result.newvarpool("y", structure.SI);
	IndexedPool const & z = result.newvarpool("z", structure.IT);


	// balance
	result.newctrpool("balance", structure.I);
	for (auto const & i : structure.I){
		result.ctr("balance", i).lb() = 1;
		result.ctr("balance", i).ub() = 1;
	}
	for (auto const & si : structure.SI){
		result.ctr("balance", si.second).f() += y(si);
	}
	// attributes
	result.newctrpool("attribute", structure.TK);
	for (int k(0); k < _k; ++k){
		// <=0
		for (auto const & t : structure.T)
			result.ctr("attribute", t, k).ub() = 0;
		// (q[s,k]-q[t,k])*x[s,t]
		for (auto const & st : structure.ST){
			int const s(st.first);
			int const t(st.second);
			result.ctr("attribute", t, k).f() += (q(s, k) - q(t, k))* x(st);
		}
		// (q[s,k]-q[t,k])*y[s,i]*z[i,t]
		for (auto const & sit : structure.SIT){
			int const s(sit.first.first);
			int const i(sit.first.second);
			int const t(sit.second);
			result.ctr("attribute", t, k).f() += (q(s, k) - q(t, k)) * y(s,i) *z( i, t);
		}
	}
	// bL
	result.newctrpool("bL", structure.bL);
	for (auto const & n : structure.bL){
		result.ctr("bL", n).lb() = bL(n);
	}
	// bU
	result.newctrpool("bU", structure.bU);
	for (auto const & n : structure.bU){
		result.ctr("bU", n).ub() = bU(n);
	}
	// 
	for (auto const & st : structure.ST){
		int const s(st.first);
		int const t(st.second);
		if (bL(s) != 0){
			result.ctr("bL", s).f() += x(s,t);
		}
		if (bU(s) != 0){
			result.ctr("bU", s).f() += x(s,t);
		}
		if (bL(t) != 0){
			result.ctr("bL", t).f() += x(s, t);
		}
		if (bU(t) != 0){
			result.ctr("bU", t).f() += x(s, t);
		}
	}
	for (auto const & it : structure.IT){
		int const i(it.first);
		int const t(it.second);
		if (bL(i) != 0){
			result.ctr("bL", i).f() += z(i, t);
		}
		if (bU(i) != 0){
			result.ctr("bU", i).f() += z(i, t);
		}
		if (bL(t) != 0){
			result.ctr("bL", t).f() += z(i, t);
		}
		if (bU(t) != 0){
			result.ctr("bU", t).f() += z(i, t);
		}
	}
	for (auto const & sit : structure.SIT){
		int const s(sit.first.first);
		int const i(sit.first.second);
		int const t(sit.second);
		if (bL(s) != 0){
			result.ctr("bL", s).f() += y(s,i)*z(i, t);
		}
		if (bU(s) != 0){
			result.ctr("bU", s).f() += y(s, i)*z(i, t);
		}
	}
	// positivity
	result.newctrpool("x_pos", structure.ST);
	for (auto const & st : structure.ST){
		result.ctr("x_pos", st).lb() = 0;
		result.ctr("x_pos", st).f() += x(st);
	}
	result.newctrpool("y_pos", structure.SI);
	for (auto const & si : structure.SI){
		result.ctr("y_pos", si).lb() = 0;
		result.ctr("y_pos", si).f() += y(si);
	}
	result.newctrpool("z_pos", structure.IT);
	for (auto const & it : structure.IT){
		result.ctr("z_pos", it).lb() = 0;
		result.ctr("z_pos", it).f() += z(it);
	}
	// cost

	for (auto const & st : structure.ST){
		result.minimize() += c(st)*x(st);
	}
	for (auto const & sit : structure.SIT){
		int const s(sit.first.first);
		int const i(sit.first.second);
		int const t(sit.second);
		result.minimize() += (c(s, i) + c(i, t))*y(s,i)*z(i,t);
	}
}



void Pooling::pFormulation(PolynomialOptimizationProblem & result)const{
	Structure structure;
	get(structure);
	
	IntPairSet IJ;
	IJ.insert(structure.ST.begin(), structure.ST.end());
	IJ.insert(structure.SI.begin(), structure.SI.end());
	IJ.insert(structure.IT.begin(), structure.IT.end());
	IntSet N;
	N.insert(structure.S.begin(), structure.S.end());
	N.insert(structure.I.begin(), structure.I.end());
	N.insert(structure.T.begin(), structure.T.end());
	// flow
	IndexedPool const & f =  result.newvarpool("f", IJ);
	// positivity of flow
	result.newctrpool("f_pos", IJ);
	for (auto const & ij : IJ){
		result.ctr("f_pos", ij).f() += f(ij);
		result.ctr("f_pos", ij).lb() = 0;
	}
	// balance
	result.newctrpool("balance", structure.I);
	for (auto const & i : structure.I){
		result.ctr("balance", i).lb() = 0;
		result.ctr("balance", i).ub() = 0;
	}
	for (auto const & si : structure.SI){
		int const s(si.first);
		int const i(si.second);
		result.ctr("balance", i).f() += f(si);
	}
	for (auto const & it : structure.IT){
		int const i(it.first);
		int const t(it.second);
		result.ctr("balance", i).f() -= f(it);
	}
	// bU
	result.newctrpool("bU", structure.bU);
	for (auto const & n : structure.bU){
		result.ctr("bU", n).ub() = bU(n);
	}
	// 
	result.newctrpool("s_cost", structure.S);
	result.newctrpool("i_cost", structure.I);

	IndexedPool const & s_cost = result.newvarpool("s_cost", structure.S);
	IndexedPool const & i_cost = result.newvarpool("i_cost", structure.I);

	for (auto const s : structure.S){
		result.ctr("s_cost", s).lb() = 0;
		result.ctr("s_cost", s).ub() = 0;
		result.ctr("s_cost", s).f() -= s_cost(s);
		result.minimize() += s_cost(s);
	}
	for (auto const i : structure.I){
		result.ctr("i_cost", i).lb() = 0;
		result.ctr("i_cost", i).ub() = 0;
		result.ctr("i_cost", i).f() -= i_cost(i);
		result.minimize() += i_cost(i);
	}
	for (auto const & st : structure.ST)
		result.ctr("s_cost", st.first).f() += c(st)*f(st);
	for (auto const & it : structure.IT)
		result.ctr("i_cost", it.first).f() += c(it)*f(it);

	// flow 
	for (auto const & ij : IJ){
		int const i(ij.first);
		int const j(ij.second);
		// upper bound on i
		if (bU(i) != 0)
			result.ctr("bU", i).f() += f(ij);
		// upper bound on j
		if (bU(j) != 0)
			result.ctr("bU", j).f() += f(ij);
	}
	// attributes
	std::vector<bool> isI(_n, false);
	for (auto const i : structure.I)
		isI[i].flip();
	IntPairSet IK;
	for (int k(0); k < _k; ++k){
		for (auto const & i : structure.I){
			IK.insert({ i, k });
		}
	}
	IntPairSet NK;
	NK.insert(IK.begin(), IK.end());
	NK.insert(structure.TK.begin(), structure.TK.end());

	result.newctrpool("attribute", NK);
	IndexedPool const & t = result.newvarpool("t", IK);
	for (int k(0); k < _k; ++k){
		// <=0
		for (auto const & i : structure.I){
			result.ctr("attribute", i, k).lb() = 0;
			result.ctr("attribute", i, k).ub() = 0;
		}
		// sum q[i,k]*f[s,i] = t[i,k] sum f[i,t] = t[i,k] sum f[s,i]
		for (auto const & si : structure.SI){
			int const s(si.first);
			int const i(si.second);
			result.ctr("attribute", i, k).f() += q(s,k)*f(si);
			result.ctr("attribute", i, k).f() += f(si)*t(i,k);
		}
	}
	// attributes
	for (int k(0); k < _k; ++k){
		// <=0
		for (auto const & t : structure.T)
			result.ctr("attribute", t, k).ub() = 0;
		// (q[s,k]-q[t,k])*x[s,t]
		for (auto const & st : structure.ST){
			int const s(st.first);
			int const t(st.second);
			result.ctr("attribute", t, k).f() += (q(s, k) - q(t, k))* f(st);
		}
		// (q[s,k]-q[t,k])*y[s,i]*z[i,t]
		for (auto const & it : structure.IT){
			int const i(it.first);
			result.ctr("attribute", it.second, k).f() += t(i, k)*f(it);
			result.ctr("attribute", it.second, k).f() -= q(it.second,k)*f(it);
		}
	}

}

void Pooling::qFormulation(PolynomialOptimizationProblem & result)const{

}