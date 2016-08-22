#include "MatPowerData.h"

Gencost::Gencost(std::string const &line) {
	std::stringstream buffer(line);
	buffer >> kind;
	buffer >> startup;
	buffer >> shutdown;
	buffer >> n;
	c.assign((int)n, 0);
	for (int i(0); i < n; ++i) {
		if (!(buffer >> c[i])) {
			std::cout << "could not read gencost c[" << i << "] in : " << line << std::endl;
		}
	}
}


void MatPowerData::read_file(std::string const & file_name) {
	std::ifstream file(file_name.c_str());
	std::string line;
	int max_bus(0);
	while (std::getline(file, line)) {
		if (line.find("mpc.bus_name") != std::string::npos) {
		}
		else if (line.find("mpc.bus") != std::string::npos) {
			std::getline(file, line);
			while (line[0] != ']') {
				NumberVectorPtr ptr = newData<Bus>(line);
				bus.push_back(ptr);
				max_bus = std::max(max_bus, (int)(*ptr)[0]);
				std::getline(file, line);
			};
		}
		else if (line.find("mpc.branch") != std::string::npos) {
			std::getline(file, line);
			while (line[0] != ']') {
				NumberVectorPtr ptr = newData<Branch>(line);
				branch.push_back(ptr);
				std::getline(file, line);
			};
		}
		else if (line.find("mpc.gencost") != std::string::npos) {
			std::getline(file, line);
			while (line[0] != ']') {
				GencostPtr ptr(new Gencost(line));
				gencost.push_back(ptr);
				std::getline(file, line);
			}
		}
		else if (line.find("mpc.gen") != std::string::npos) {
			std::getline(file, line);
			while (line[0] != ']') {
				NumberVectorPtr ptr = newData<Gen>(line);
				gen.push_back(ptr);
				std::getline(file, line);
			};
		}
		else if (line.find("mpc.baseMVA") != std::string::npos) {
			std::stringstream buffer(line);
			buffer >> line;
			buffer >> line;
			if (!(buffer >> baseMVA))std::cout << "unable to read baseMVA" << std::endl;
		}
	}
	std::cout << "max_bus is " << max_bus << std::endl;
	bus2index.assign(max_bus + 1, -1);
	for (int i(0); i < bus.size(); ++i) {
		bus2index[(int)(*bus[i])[0]] = i;
	}
}

void MatPowerData::generate_opf(PolynomialOptimizationProblem & output)const {
	output.clear();
	int const nBus((int)bus.size());
	int const nBranch((int)branch.size());
	int const nGen((int)gen.size());
	int const nGenCost((int)gencost.size());

	IndexedPool	const & v = output.newvarpool("V", nBus);

	IndexedPool	const & balance = output.newctrpool("balance", nBus);
	IndexedPool	const & vmin = output.newctrpool("vmin", nBus);
	IndexedPool	const & vmax = output.newctrpool("vmax", nBus);

	for (int i(0); i < nBus; ++i) {
		NumberVector const & data(*bus[i]);
		// vmin
		output.ctr(vmin.id(i)).f() += v(i).conjugate()*v(i);
		output.ctr(vmin.id(i)).lb() = data[Bus::bus_Vmin] * data[Bus::bus_Vmin];
		// vmax
		output.ctr(vmax.id(i)).f() += v(i).conjugate()*v(i);
		output.ctr(vmax.id(i)).ub() = data[Bus::bus_Vmax] * data[Bus::bus_Vmax];
		// balance
		output.ctr(balance.id(i)).lb() = -ComplexNumber(data[Bus::bus_Pd], data[Bus::bus_Qd]);
		output.ctr(balance.id(i)).ub() = output.ctr(balance.id(i)).lb();
	}
	for (int i(0); i < nBranch; ++i) {
		NumberVector const & data(*branch[i]);
		int const or (bus2index[(int)data[Branch::branch_fbus]]);
		int const ex(bus2index[(int)data[Branch::branch_tbus]]);

		ComplexPolynomial current_or = baseMVA*(get_11(data)*v(or ) + get_12(data)*v(ex));
		ComplexPolynomial current_ex = baseMVA*(get_21(data)*v(or ) + get_22(data)*v(ex));
		output.ctr(balance.id(or )).f() += current_or.conjugate() * v(or );
		output.ctr(balance.id(ex)).f() += current_ex.conjugate() * v(ex);
		//if(or==3)
		//std::exit(0);
	}
	for (int i(0); i < nGen; ++i) {
		NumberVector const & data(*gen[i]);
		int const bus(bus2index[(int)data[Gen::gen_bus]]);
		output.ctr(balance.id(bus)).lb() = ComplexNumber(data[Gen::gen_Pmin], data[Gen::gen_Qmin]);
		output.ctr(balance.id(bus)).ub() = ComplexNumber(data[Gen::gen_Pmax], data[Gen::gen_Qmax]);
	}
	for (int i(0); i < nGenCost; ++i) {
		Gencost const & data(*gencost[i]);
		if (data.kind != 2) {
			throw std::invalid_argument("gencost type is not 2");
		}
		int const bus(bus2index[(int)(*gen[i])[Gen::gen_bus]]);
		for (int j(0); j < data.n; ++j) {
			if (data.c[j] != 0 && j > 0) {
				ComplexPolynomial term(1);
				for (int k(0); k < data.n - 1 - j; ++k) {
					term *= output.ctr(balance.id(bus)).f();
				}
				output.minimize() += data.c[j] * term;
				//output.minimize().print(std::cout, output);
				//std::cout << std::endl;
			}
		}
	}
	//std::cout << output << std::endl;
}
ComplexPolynomial MatPowerData::get_shift(NumberVector const & data)const {
	ComplexPolynomial result;
	if (data[Branch::branch_ratio] != 0) {
		result += ComplexPolynomial(data[Branch::branch_ratio] * std::cos(data[Branch::branch_angle]), data[Branch::branch_ratio] * std::sin(data[Branch::branch_angle]));
	}
	else {
		result += 1;
	}
	return result;
}
ComplexPolynomial MatPowerData::get_z(NumberVector const & data)const {
	return data[Branch::branch_r] + ComplexPolynomial::i()* data[Branch::branch_x];
}

ComplexPolynomial MatPowerData::get_ib(NumberVector const & data)const {
	return data[Branch::branch_b] * ComplexPolynomial::i();
}

ComplexPolynomial MatPowerData::get_y(NumberVector const & data)const {
	return 1 / get_z(data);
}

ComplexPolynomial MatPowerData::get_11(NumberVector const & data)const {
	ComplexPolynomial tau_2 = get_shift(data).conjugate()*get_shift(data);
	return get_22(data) / tau_2;
}

ComplexPolynomial MatPowerData::get_22(NumberVector const & data)const {
	return get_y(data) + 0.5*get_ib(data);
}

ComplexPolynomial MatPowerData::get_21(NumberVector const & data)const {
	return -get_y(data) / (get_shift(data));
}

ComplexPolynomial MatPowerData::get_12(NumberVector const & data)const {
	return -get_y(data) / (get_shift(data).conjugate());
}