#pragma once

#include "../GENERIC_DEFINITIONS.h"

namespace UV_TO_XYZ		// 0 <= u, v <= 1
{
	const TV GetCylinderXYZ(const Vector2D<T>& uv, const T& upper_radius_, const T& bottom_radius_, const T& width_);

	const TV GetPrismXYZ(const Vector2D<T>& uv, const int& num_sides_, const T& upper_radius_, const T& bottom_radius_, const T& width_);

	const TV GetPrismXYZLinearCorner(const Vector2D<T>& uv, const int& num_sides_, const float& corner_angle, const T& upper_radius_, const T& bottom_radius_, const T& width_);	// corner_angle is degree from user input.
};
