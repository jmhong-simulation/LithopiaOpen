// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "GENERIC_DEFINITIONS.h"
#include "CUBIC_POLYNOMIAL.h"
#include "DataStructure/Array1D.h"

class ParametricCurveSegment
{
public:
	CubicPolynomial x_, y_, z_;

public:
	TV getPosition(const T t) const
	{
		return TV(x_.getValue(t), y_.getValue(t), z_.getValue(t));
	}

	T getLength(const T t0, const T t1) const
	{
		const TV x0 = getPosition(t0), x1 = getPosition(t1);

		return (T)(x1 - x0).getMagnitudeDouble();
	}

	T getLength(const T t0, const T t1, const int num_segments) const
	{
		const T dt = (t1 - t0) / (T)num_segments;

		T length = (T)0;
		for (T t = t0; t < t1; t += dt)
		{
			length += getLength(t0, t0 + dt);
		}

		return length;
	}
};