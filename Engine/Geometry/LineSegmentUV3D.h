// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "GENERIC_DEFINITIONS.h"

class LineSegmentUV3D {

public:

	TV3 p0_, p1_;
	TV2 u0_, u1_;

	float v0_, v1_;

	int ix_;

	LineSegmentUV3D()
	{}

	LineSegmentUV3D(const TV3& _p0, const TV3& _p1, const int _ix)
		: p0_(_p0), p1_(_p1), ix_(_ix)
	{}
	LineSegmentUV3D(const TV3& _p0, const TV3& _p1, const TV2& _u0, const TV2& _u1, const int _ix)
		: p0_(_p0), p1_(_p1), u0_(_u0), u1_(_u1), ix_(_ix)
	{}
	LineSegmentUV3D(const TV3& _p0, const TV3& _p1, const TV2& _u0, const TV2& _u1, const float _v0, const float _v1, const int _ix)
		: p0_(_p0), p1_(_p1), u0_(_u0), u1_(_u1), v0_(_v0), v1_(_v1), ix_(_ix)
	{}

	void Initialize(const TV3& _p0, const TV3& _p1, const int _ix) {

		p0_ = _p0;
		p1_ = _p1;

		ix_ = _ix;
	}

	void Initialize(const TV3& _p0, const TV3& _p1, const TV2& _u0, const TV2& _u1, const int _ix) {

		p0_ = _p0;
		p1_ = _p1;

		u0_ = _u0;
		u1_ = _u1;

		ix_ = _ix;
	}

	void Initialize(const TV3& _p0, const TV3& _p1, const TV2& _u0, const TV2& _u1, const float _v0, const float _v1, const int _ix) {

		p0_ = _p0;
		p1_ = _p1;

		u0_ = _u0;
		u1_ = _u1;

		v0_ = _v0;
		v1_ = _v1;

		ix_ = _ix;
	}

	T getDistance(const TV& v) {

		const T mag = (T)(p1_ - p0_).getMagnitudeDouble();
		const T t = dotProduct(v - p0_, (p1_ - p0_).normalizedDouble());

		if (t <= 0.0f) return (v - p0_).getMagnitude();

		if (t >= mag) return (v - p1_).getMagnitude();

		return (T)(v - ((p1_ - p0_).normalizedDouble()*t + p0_)).getMagnitudeDouble();
	}

	T getDistance(TV2& v) {

		const T mag = (T)(u1_ - u0_).getMagnitude();
		const T t = dotProduct(v - u0_, (u1_ - u0_).getNormalized());

		if (t <= 0.0f) return (v - u0_).getMagnitude();

		if (t >= mag) return (v - u1_).getMagnitude();

		return (T)(v - ((u1_ - u0_).getNormalized()*t + u0_)).getMagnitude();

	}

	T getDistance(TV2 v, const int wid, const int hei) {

		v.x_ *= (float)(wid - 1);
		v.y_ *= (float)(hei - 1);

		return getDistance(v);
	}
};