// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "../GENERIC_DEFINITIONS.h"
#include "../DataStructure/Vector4D.h"

class CubicPolynomial
{
public:
	T a_, b_, c_, d_;	// follow notations in http://en.wikipedia.org/wiki/Cubic_function

public:
	CubicPolynomial()
	{}

	~CubicPolynomial()
	{}

	void InitializeZeroOneDomain(const Vector4D<T> curve)
	{
		InitializeZeroOneDomain(curve.x_, curve.y_, curve.z_, curve.w_);
	}

	void InitializeZeroOneDomain(const T f0, const T fdot0, const T f1, const T fdot1)
	{
		// special case when x0 = 0, x1 = 1
		d_ = f0;
		c_ = fdot0;
		a_ = (T)2 * f0 - (T)2 * f1 + fdot0 + fdot1;
		b_ = (T)-3 * f0 + (T)3 * f1 - (T)2 * fdot0 - fdot1;
	}

	const T getValue(const T x) const
	{
		return a_*POW3(x) + b_*POW2(x) + c_*x + d_;
	}

	const T Get1stOrderDerivative(const T x) const
	{
		return (T)3 * a_*POW2(x) + (T)2 * b_*x + c_;
	}

	const T Get2ndOrderDerivative(const T x) const
	{
		return (T)6 * a_ * x + (T)2 * b_;
	}
};