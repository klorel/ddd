#pragma once

#include "common.h"
#include "ComplexPolynomial.h"
#include "PolynomialOptimizationProblem.h"

class Gencost;
typedef std::shared_ptr<Gencost> GencostPtr;
typedef std::vector<GencostPtr> GencostPtrVector;

class Gencost {
public:
	Number kind;
	Number startup;
	Number shutdown;
	Number n;
	NumberVector c;

	Gencost(std::string const &line);
};

class MatPowerData {
public:
public:
	Number baseMVA;
	enum Bus {
#define __MAT_POWER_DATA__(x) bus_ ## x,
#include "MatPowerBus.hxx"
#undef __MAT_POWER_DATA__
		BUS_MAX_SIZE
	};
	enum Gen {
#define __MAT_POWER_DATA__(x) gen_ ## x,
#include "MatPowerGen.hxx"
#undef __MAT_POWER_DATA__
		GEN_MAX_SIZE
	};
	enum Branch {
#define __MAT_POWER_DATA__(x) branch_ ## x,
#include "MatPowerBranch.hxx"
#undef __MAT_POWER_DATA__
		BRANCH_MAX_SIZE
	};
	template<class __DATA__>int MAX_SIZE() {
		return -1;
	}
	template<>int MAX_SIZE<Bus>() {
		return Bus::BUS_MAX_SIZE;
	}
	template<>int MAX_SIZE<Branch>() {
		return Branch::BRANCH_MAX_SIZE;
	}
	template<>int MAX_SIZE<Gen>() {
		return Gen::GEN_MAX_SIZE;
	}

	template<class __DATA__>
	NumberVectorPtr newData(std::string const & line) {
		std::stringstream buffer(line);
		NumberVectorPtr ptr(new NumberVector(MAX_SIZE<__DATA__>()));
		NumberVector & v(*ptr);
		for (int i(0); i < v.size(); ++i) {
			if (!(buffer >> v[i])) {
				std::cout << "fails to read " << i << " of " << line << std::endl;
			}
		}
		return ptr;
	}
	ComplexPolynomial get_shift(NumberVector const &)const;
	ComplexPolynomial get_ib(NumberVector const &)const;
	ComplexPolynomial get_z(NumberVector const &)const;
	ComplexPolynomial get_y(NumberVector const &)const;
	ComplexPolynomial get_11(NumberVector const &)const;
	ComplexPolynomial get_21(NumberVector const &)const;
	ComplexPolynomial get_12(NumberVector const &)const;
	ComplexPolynomial get_22(NumberVector const &)const;
public:
	void read_file(std::string const & file_name);
	void generate_opf(PolynomialOptimizationProblem & output)const;
public:
	NumberVectorPtrVector bus;
	NumberVectorPtrVector gen;
	NumberVectorPtrVector branch;
	GencostPtrVector gencost;
	IntVector bus2index;


};