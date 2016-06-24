#pragma once

#include "common.h"

class IndexedPool;

typedef std::shared_ptr<IndexedPool> IndexedPoolPtr;
typedef std::map<std::string, IndexedPoolPtr> Str2Pool;

class IndexedPool{
public:
	IndexedPool(std::string const & poolname, StrPtrVector & names, size_t size) :_first(static_cast<int>(names.size())), _size(size), _name(poolname){
		names.resize(names.size() + _size);
	}
	virtual ~IndexedPool(){
	}
public:
	int _first;
	size_t _size;
	std::string _name;

	virtual void fill(StrPtrVector & names)const { throw std::logic_error("NOT IMPLEMENTED"); }
	virtual int size()const { throw std::logic_error("NOT IMPLEMENTED"); }
	virtual int first()const { throw std::logic_error("NOT IMPLEMENTED"); }

	virtual int id(int)const{ throw std::logic_error("NOT IMPLEMENTED"); }
	virtual int id(int, int)const{ throw std::logic_error("NOT IMPLEMENTED"); }
	virtual int id(int, int, int)const{ throw std::logic_error("NOT IMPLEMENTED"); }
	virtual int id(int, int, int, int)const{ throw std::logic_error("NOT IMPLEMENTED"); }
};

class IndexedPool1Dense : public IndexedPool{
public:
	IndexedPool1Dense(std::string const & poolname, StrPtrVector & names, size_t size) :IndexedPool(poolname, names, size){
		for (size_t i(0); i < _size; ++i)
		{
			std::stringstream buffer;
			buffer << poolname << "[" << i << "]";
			names[names.size() - _size + i] = std::make_shared<std::string>(buffer.str());
		}
	}
	virtual ~IndexedPool1Dense(){

	}

	virtual int id(int i)const{ return  _first + i; }
};

class IndexedPool1Sparse : public IndexedPool{
public:
	IndexedPool1Sparse(std::string const & poolname, StrPtrVector & names, IntSet const & ids) :IndexedPool(poolname, names, ids.size()){
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
	virtual ~IndexedPool1Sparse(){

	}

	virtual int id(int i)const{ return _first+_id.find(i)->second; }
protected:
	Int2Int _id;
};
class IndexedPool2Sparse : public IndexedPool{
public:
	IndexedPool2Sparse(std::string const & poolname, StrPtrVector & names, IntPairSet const & ids) :IndexedPool(poolname, names, ids.size()){
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
	virtual ~IndexedPool2Sparse(){

	}

	virtual int id(int i, int j)const{ return _first+ _id.find({ i, j })->second; }
protected:
	IntPair2Int _id;
};
