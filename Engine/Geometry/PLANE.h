// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "GENERIC_DEFINITIONS.h"

class PLANE
{
public:
	TV n_, p_;	// use d_ instead of p_
	T d_;

	PLANE()
	{}

	PLANE(const TV& _n, const TV& _p)
		: n_(_n), p_(_p)
	{
		d_ = dotProduct(n_, p_);
	}

	const T GetSignedDistance(const TV x) const
	{
		return dotProduct(x - p_, n_);
	}

	const T GetDistance(const TV x) const
	{
		return ABS(dotProduct(x - p_, n_));
	}

	const TV intersectRay(const TV& p0, const TV& p1) const {

		assert(false);

		TV v = (p0 - p1).normalized();
		T t = -(dotProduct(p0, n_) + d_) / dotProduct(v, n_);

		return p0 + v*t;
	}

	const bool intersectRay(const TV& p0, const TV& p1, TV& out) const {

		assert(false);

		TV v = (p0 - p1).normalized();
		T t = -(dotProduct(p0, n_) + d_) / dotProduct(v, n_);

		out =  p0 + v*t;

		TV cen = n_ * -d_;
		T dot0 = dotProduct(n_, p0 - cen);
		T dot1 = dotProduct(n_, p1 - cen);

		if (dot0 * dot1 < 0.0f) return true;

		return false;
	}
};