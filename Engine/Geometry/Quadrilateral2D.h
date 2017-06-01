// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "GENERIC_DEFINITIONS.h"
#include "Triangle2D.h"
#include "BOX_2D.h"

typedef class Quadrilateral2D
{
public:
	TV2 v0_, v1_, v2_, v3_;		// CCW order

	Quadrilateral2D(const TV2& _v0, const TV2& _v1, const TV2& _v2, const TV2& _v3)
		: v0_(_v0), v1_(_v1), v2_(_v2), v3_(_v3)
	{}

	T getArea()
	{
		return Triangle2D(v0_, v1_, v3_).getArea() + Triangle2D(v3_, v1_, v2_).getArea();		// sum of two triangles
	}

	bool isInside(const TV2& p)
	{
		//http://stackoverflow.com/questions/5922027/how-to-determine-if-a-point-is-within-a-quadrilateral

		if (getAreaOfSubtriangles(p) > getArea()) return false;
		else return true;
	}

	T getAreaOfSubtriangles(const TV2& p)
	{
		//http://stackoverflow.com/questions/5922027/how-to-determine-if-a-point-is-within-a-quadrilateral for inside test

		return Triangle2D(v0_, p, v3_).getArea() + Triangle2D(v0_, v1_, p).getArea() + Triangle2D(p, v1_, v2_).getArea() + Triangle2D(v3_, p, v2_).getArea();
	}

	BOX_2D<T> getAABB() const;

	int getInsideFlag(const TV2& p) const
	{
		return Triangle2D(v0_, v1_, v3_).getInsideFlag(p) + Triangle2D(v3_, v1_, v2_).getInsideFlag(p);
	}

} Quadrangle2D, Tetragon2D;

