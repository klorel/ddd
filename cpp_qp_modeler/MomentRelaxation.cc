#include "MomentRelaxation.h"


MomentRelaxation::MomentRelaxation(size_t nvariables, size_t order):_nvariables(nvariables), _order(order) {

}
MomentRelaxation::~MomentRelaxation() {

}


double MomentRelaxation::size()const {
	return 1 + std::pow(_nvariables, _order);
}