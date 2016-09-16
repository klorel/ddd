#pragma once

#include "common.h"
#include "ComplexPolynomial.h"


class IndexedPool {
public:
	IndexedPool(std::string const & poolname, StrPtrVector & names, size_t size) :_first(static_cast<int>(names.size())), _size(size), _name(poolname) {
		names.resize(names.size() + _size);
	}
	virtual ~IndexedPool() {
	}
public:
	int _first;
	size_t _size;
	std::string _name;

	virtual void fill(StrPtrVector & names)const { throw std::logic_error("NOT IMPLEMENTED"); }
	virtual int size()const { return (int)_size; }
	virtual int first()const { return  _first; }

	virtual int id(int)const { throw std::logic_error("NOT IMPLEMENTED"); }
	virtual int id(int, int)const { throw std::logic_error("NOT IMPLEMENTED"); }
	virtual int id(int, int, int)const { throw std::logic_error("NOT IMPLEMENTED"); }
	virtual int id(int, int, int, int)const { throw std::logic_error("NOT IMPLEMENTED"); }

	ComplexPolynomial operator()(int key1) const { return ComplexPolynomial::Build(id(key1)); }
	ComplexPolynomial operator()(int key1, int key2)  const { return ComplexPolynomial::Build(id(key1, key2)); }
	ComplexPolynomial operator()(int key1, int key2, int key3)  const { return ComplexPolynomial::Build(id(key1, key2, key3)); }
	ComplexPolynomial operator()(int key1, int key2, int key3, int key4) const { return ComplexPolynomial::Build(id(key1, key2, key3, key4)); }

	ComplexPolynomial operator()(IntPair const & key)  const { return operator()(key.first, key.second); }
};

class IndexedPool1Dense : public IndexedPool {
public:
	IndexedPool1Dense(std::string const & poolname, StrPtrVector & names, size_t size) :IndexedPool(poolname, names, size) {
		if (_size > 1) {
			for (size_t i(0); i < _size; ++i)
			{
				std::stringstream buffer;
				buffer << poolname << "[" << i << "]";
				names[names.size() - _size + i] = std::make_shared<std::string>(buffer.str());
			}
		}
		else {
			std::stringstream buffer;
			buffer << poolname;
			names[names.size() - _size] = std::make_shared<std::string>(buffer.str());
		}
	}
	virtual ~IndexedPool1Dense() {

	}

	virtual int id(int i)const { return  _first + i; }
};

class IndexedPool1Sparse : public IndexedPool {
public:
	IndexedPool1Sparse(std::string const & poolname, StrPtrVector & names, IntSet const & ids) :IndexedPool(poolname, names, ids.size()) {
		int i(0);
		_id.clear();
		for (auto const & id : ids)
		{
			_id[id] = i;
			std::stringstream buffer;
			buffer << poolname << "[" << id << "]";
			names[names.size() - _size + i] = std::make_shared<std::string>(buffer.str());
			++i;
		}
	}
	virtual ~IndexedPool1Sparse() {

	}

	virtual int id(int i)const { return _first + _id.find(i)->second; }
protected:
	Int2Int _id;
};
class IndexedPool2Sparse : public IndexedPool {
public:
	IndexedPool2Sparse(std::string const & poolname, StrPtrVector & names, IntPairSet const & ids) :IndexedPool(poolname, names, ids.size()) {
		int i(0);
		_id.clear();
		for (auto const & id : ids)
		{
			_id[id] = i;
			std::stringstream buffer;
			buffer << poolname << "[" << id.first << ", " << id.second << "]";
			names[names.size() - _size + i] = std::make_shared<std::string>(buffer.str());
			++i;
		}
	}
	virtual ~IndexedPool2Sparse() {

	}

	virtual int id(int i, int j)const { return _first + _id.find({ i, j })->second; }
protected:
	IntPair2Int _id;
};

class IndexedPool2Square : public IndexedPool {
public:
	IndexedPool2Square(std::string const & poolname, StrPtrVector & names, int n, int id) :IndexedPool(poolname, names, n*(n + 1) / 2), _n(n), _id(id) {
		int i(0);
		for (int i(0); i < _n; ++i) {
			for (int j(i); j < _n; ++j)
			{
				std::stringstream buffer;
				buffer << poolname << "[" << i << ", " << j << "]";
				names[names.size() - _size + i] = std::make_shared<std::string>(buffer.str());
			}
		}
	}
	virtual ~IndexedPool2Square() {

	}

	virtual int id(int i, int j)const { return _first+i*(2 * _n - i + 1) / 2 + j; }

	int len()const { return _n; }
protected:
	int _n;
	int _id;
};
typedef std::shared_ptr<IndexedPool2Square> IndexedPool2SquarePtr;