#pragma once

#include "GENERIC_DEFINITIONS.h"
#include "RAY.h"
#include "Sphere3D.h"
#include "PLANE.h"

namespace Intersection
{
	// ray r = p + td, |d| = 1
	// returns intersecting t and intersection point q
	int checkRaySphere(const TV3& p, const TV& d, const Sphere3D& sphere, T& t, TV& q);

	int checkRaySphere(const RAY& ray, const Sphere3D& sphere, T& t, TV& intersection_point);

	int checkLinePlane(const TV& a, const TV& b, const PLANE& p, T& t, TV& q);	// q = intersection point
}