#pragma once

#include "GENERIC_DEFINITIONS.h"
#include "Geometry/BOX_2D.h"

class Triangle2D
{
public:
	TV2 v0_, v1_, v2_;

	Triangle2D(const TV2& _v0, const TV2& _v1, const TV2& _v2)
		: v0_(_v0), v1_(_v1), v2_(_v2)
	{}

	T getArea()
	{
		// http://www.mathopenref.com/coordtrianglearea.html

		return ABS(getSignedArea());
	}

	T getSignedArea() const
	{
		return (v0_.x_ * (v1_.y_ - v2_.y_) + v1_.x_ * (v2_.y_ - v0_.y_) + v2_.x_ * (v0_.y_ - v1_.y_)) * (T)0.5;		// negative means CW 
	}

	T getSignedDoubleArea() const	// positive means CCW, negative means CW
	{
		return (v0_.x_ * (v1_.y_ - v2_.y_) + v1_.x_ * (v2_.y_ - v0_.y_) + v2_.x_ * (v0_.y_ - v1_.y_));
	}

	BOX_2D<T> getAABB() const;

	TV3 getBarycentricCoordinates(const TV2& location) const;

	bool isInside(const TV2& p) const;

	int getInsideFlag(const TV2& p) const;	// return 0 if outside, 1 if CCW and inside, -1 if CW and inside
};