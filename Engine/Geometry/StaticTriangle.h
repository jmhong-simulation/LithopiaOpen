// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "GENERIC_DEFINITIONS.h"
#include "DataStructure/LinkedArray.h"
#include "DataStructure/DynamicArray.h"
#include "DataStructure/Vector3D.h"

class StaticTriangle
{
public:
	TV v0_, v1_, v2_;

public:
	StaticTriangle()
	{}

	StaticTriangle(const TV& _v0, const TV& _v1, const TV& _v2)
		: v0_(_v0), v1_(_v1), v2_(_v2)
	{}

	TV getClosestPointFromLine(const TV& location, const TV& x1, const TV& x2) const;
	TV getBarycentricCoordinates(const TV& location, const TV& x1, const TV& x2, const TV& x3) const;
	TV getClosestPosition(const TV& location) const;
	T  getDistance(const TV& location) const;
	T  getArea() const;
	T  getAspectRatio() const;

	TV getNormal() const {

		return crossProduct(v1_ - v0_, v2_ - v0_).normalized();

	}

	bool checkOnTriangle(const TV& p, const float rad) {

		TV3 cen = (v0_ + v1_ + v2_) / 3.0f;

		TV3 n = crossProduct(v1_ - v0_, v2_ - v0_).normalized();
		float h = dotProduct(n, p - v0_);

		TV3 n0 = crossProduct(n, (v1_ - v0_).normalized()).normalized();
		if (dotProduct(v0_ - cen, n0) < 0.0f) {
			n0 *= -1.0f;
		}

		TV3 n1 = crossProduct(n, (v2_ - v1_).normalized()).normalized();
		if (dotProduct(v1_ - cen, n1) < 0.0f) {
			n1 *= -1.0f;
		}

		TV3 n2 = crossProduct(n, (v0_ - v2_).normalized()).normalized();
		if (dotProduct(v2_ - cen, n2) < 0.0f) {
			n2 *= -1.0f;
		}

		float d0 = dotProduct(n0, p - v0_);
		float d1 = dotProduct(n1, p - v1_);
		float d2 = dotProduct(n2, p - v2_);

		if (d0 <= 1e-04 && d1 <= 1e-04 && d2 <= 1e-04) {

			if (h >= -rad && h <= rad) {

				return true;

			}

		}

		return false;
	}

	void getXMinMax(T& x_min, T& x_max) const;
	void getYMinMax(T& y_min, T& y_max) const;
	void getZMinMax(T& z_min, T& z_max) const;


	template<class EE> const EE Interpolate(const TV& q, const TV& p0, const TV& p1, const EE& v0, const EE& v1) const {

		float dot = dotProduct(q - p0, p1 - p0) / dotProduct(p1 - p0, p1 - p0);
		dot = CLAMP(dot, 0.0f, 1.0f);
		return v0*(1.0f-dot) + v1*dot;
	}

	template<class EE> const EE Interpolate(const TV& p, const EE arr[3]) const	{

		float dot = dotProduct(p - v0_, v1_ - v0_) / (dotProduct(v1_ - v0_, v1_ - v0_));
		dot = CLAMP(dot, 0.0f, 1.0f);

		TV v00 = v0_ + (v1_ - v0_) * dot;
		EE e00 = arr[0] * (1.0f - dot) + arr[1] * dot;


		dot = dotProduct(p - v0_, v2_ - v0_) / (dotProduct(v2_ - v0_, v2_ - v0_));
		dot = CLAMP(dot, 0.0f, 1.0f);

		TV v11 = v0_ + (v2_ - v0_) * dot;
		EE e11 = arr[0] * (1.0f - dot) + arr[2] * dot;

		
		dot = dotProduct(p - v00, v11 - v00) / dotProduct(v11 - v00, v11 - v00);
		dot = CLAMP(dot, 0.0f, 1.0f);


		EE result = e00*(1.0f - dot) + e11*dot;

		return result;
	}
};