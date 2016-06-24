#pragma once

#include "common.h"

class Branch;
class Bus;


class branch {
public:
	enum Attributes {
		Z, YOR, YOEX, PSTOR, PSTEX, NUM
	};
};

class bus {
public:
	enum Attributes {
		SD, VC, BOUNDS, NUM
	};
};