#pragma once

#include "GENERIC_DEFINITIONS.h"
#include "CONVENTIONAL_MACROS.h"
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <ctime>

class RandomNumberGenerator
{
public:
	boost::random::mt19937 &generator_;
	boost::random::uniform_int_distribution<> &dist_;

public:
	RandomNumberGenerator(const int seed)
		: generator_(*(new boost::random::mt19937(40+seed))), dist_(*(new boost::random::uniform_int_distribution<>(0,(int)1e6)))
	{}

	~RandomNumberGenerator(void)
	{
		delete &generator_;
		delete &dist_;
	}

	int getProbabilityInteger(const T& probability)	const // getProbabilityInteger(3.5) returns 3 + (0 or 1)
	{
		return (int)probability + ((probability - (T)((int)probability)) > getNumber() ? 1 : 0);
	}

	T getNumber() const		// 0 ~ 1
	{
		return (T)(dist_(generator_)*1e-6);
	}

	TV getVector() const	// (0~1, 0~1, 0~1)
	{
		return TV(getNumber(), getNumber(), getNumber());
	}

	TV getUnitVector() const
	{
		// for spherical coordinates, see http://en.wikipedia.org/wiki/Spherical_coordinate_system
		const static T TWO_PI = (T)2 * (T)PI;

		const T theta = getNumber()*(T)PI;
		const T phi = getNumber()*TWO_PI;

		return TV(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
	}
};

