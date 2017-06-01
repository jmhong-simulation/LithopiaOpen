#include "Intersection.h"

namespace Intersection
{
	// ray r = p + td, |d| = 1
	// returns intersecting t and intersection point q
	int checkRaySphere(const TV3& p, const TV& d, const Sphere3D& sphere, T& t, TV& q)
	{
		const TV m = p - sphere.center_;
		const T b = dotProduct(m, d);
		const T c = dotProduct(m, m) - sphere.radius_ * sphere.radius_;

		// exit if r's origin outside s (c > 0) and r pointing away from s (b > 0)
		if (c > (T)0 && b > (T)0) return 0;	// no collision

		T discr = b*b - c;

		// a negative discriminant corresponds to ray missing sphere
		if (discr < (T)0) return 0;

		// ray now found to intersect sphere, compute smallest t value of intersection
		t = -b - sqrt(discr);

		// if t is negative, ray started inside sphere so clamp t to zero
		if (t < 0.0f) t = 0.0f;
		q = q + t * d;
		return 1;
	}

	int checkRaySphere(const RAY& ray, const Sphere3D& sphere, T& t, TV& intersection_point)
	{
		return checkRaySphere(ray.origin_, ray.direction_, sphere, t, intersection_point);
	}

	int checkLinePlane(const TV& a, const TV& b, const PLANE& p, T& t, TV& q)	// q = intersection point
	{
		// compute the t value for the directed line ab intersecting the plane
		const TV ab = b - a;

		t = (p.d_ - dotProduct(p.n_, a)) / dotProduct(p.n_, ab);

		// if t in [0..1] compute and return intersection point
		if (t >= (T)0 && t <= (T)1)
		{
			q = a + t * ab;
			return 1;
		}

		// else no intersection
		return 0;
	}
}