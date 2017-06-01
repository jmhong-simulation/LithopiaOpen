#include "RAY.h"

const TV RAY::GetPosition(const T& t) const
{
	return origin_ + direction_*t_end_;
}

// http://www.scratchapixel.com/lessons/3d-basic-lessons/lesson-9-ray-triangle-intersection/m-ller-trumbore-algorithm/
const bool RAY::CheckTriangleIntersection(const TV& v0, const TV& v1, const TV& v2, const RAY& segment) const
{
	const TV edge1 = v1 - v0;
	const TV edge2 = v2 - v0;
	const TV pvec  = crossProduct(segment.direction_, edge2);
	const T	 det   = dotProduct(edge1, pvec);

	if (det == 0) return false;			// the direction is parallel to the triangle

	const T  invDet = 1.0f / det;

	const TV tvec = segment.origin_ - v0;

	const T u = dotProduct(tvec, pvec) * invDet;

	if (u < 0 || u > 1) return false;

	const TV qvec = crossProduct(tvec, edge1);

	const T v = dotProduct(segment.direction_, qvec) * invDet;

	if (v < 0 || u + v > 1) return false;

	const T t = dotProduct(edge2, qvec) * invDet;

	if (t < 0 || t > segment.t_end_) return false;

	return true;
}

//Code for Ray Sphere Intersection :
//rayStart is camera position
//rayDirection is direction of the ray i.e. (s - e)
const T RAY::GetSphereIntersection(const TV& center, const T& radius) const
{
	T t1 = -1;
	T t2 = -1;
	T t = -1;

	//temporary == e-c
	TV temporary = origin_ - center;

	T b = 2 * (dotProduct(direction_, temporary));
	T a = dotProduct(direction_, direction_);
	T c = dotProduct(temporary, temporary) - (radius * radius);
	T disc;
	disc = b*b - 4 * a*c;
	if (disc < 0){
		return -1;	// no intersection
	}
	else if (disc >= 0)	// two intersections
	{
		T discriminent = sqrt(disc);
		t1 = (-b + discriminent) / (2 * a);
		t2 = (-b - discriminent) / (2 * a);

		if (t1 >= 0 && t2 >= 0) return MIN2(t1, t2);
		else
		{
			if (t1 >= 0) return t1;
			if (t2 >= 0) return t2;
		}
	}
	else{		// one intersection
		T discriminent = sqrt(disc);

		t1 = (-b + discriminent) / (2 * a);

		if (t1 >= 0) return t1;
	}

	return -1;
}
