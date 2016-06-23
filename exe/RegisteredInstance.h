#pragma once

#include "common.h"
#include "Pooling.h"

enum AvailableInstances {
#define REGISTER_INSTANCE(x) x,
#include "RegisteredInstance.hxx"
#undef REGISTER_INSTANCE
	SIZE
};

class RegisteredInstance: public Pooling {
public:
public:
	void setName(AvailableInstances id);
public:
	void out() const;
	std::string const & name() const {
		return _name;
	}
private:
	std::string _name;
public:
	RegisteredInstance(AvailableInstances id);
	RegisteredInstance(size_t i);
	virtual ~RegisteredInstance() {
	}
public:
	static void PrintAvailable(std::ostream & stream = std::cout);

};

inline void RegisteredInstance::PrintAvailable(std::ostream & stream) {
	stream << "Available instances : \n";
	for (size_t i(0); i < AvailableInstances::SIZE; ++i) {
		AvailableInstances id(static_cast<AvailableInstances>(i));
		RegisteredInstance instance(id);
		stream << std::setw(3) << i;
		stream << " : ";
		stream << std::setw(30) << std::left << instance.name() << std::right;
		stream << "\n";
	}
	stream << "<exe> <id of selected instance> \n";

}
class Info {
public:
	static std::string const InstancesPath;
};

inline void RegisteredInstance::out() const {
	std::cout <<"Instance name is "<<_name<<"\n";
	std::cout << "Data were read from " << Info::InstancesPath + _name + ".data" << "\n";
	Pooling::out();
}
inline RegisteredInstance::RegisteredInstance(AvailableInstances id) {
	setName(id);
	read(Info::InstancesPath + _name + ".data");
}

inline RegisteredInstance::RegisteredInstance(size_t i) {
	AvailableInstances const id(static_cast<AvailableInstances>(i));
	setName(id);
	read(Info::InstancesPath + _name + ".data");
}

inline void RegisteredInstance::setName(AvailableInstances id) {
	switch (id) {
#define REGISTER_INSTANCE(x) case x:_name =#x;break;
#include "RegisteredInstance.hxx"
#undef REGISTER_INSTANCE
	default:
		std::cout << id << std::endl;
		throw std::invalid_argument("UNKNOWN INSTANCE");
		break;
	}
}
