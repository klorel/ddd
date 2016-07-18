
#include "common.h"
#include "Problem.h"

int main(int argc, char** argv) {
	Problem pop;
	IndexedPool const & z = pop.newvarpool("x", 10);
	pop.add(z(2) * 3 <= 1);
	std::cout << pop << std::endl;
	return 0;
}