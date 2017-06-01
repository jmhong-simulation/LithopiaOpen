#include "UV_TO_XYZ.h"

namespace UV_TO_XYZ		// 0 <= u, v <= 1
{
	const TV GetCylinderXYZ(const Vector2D<T>& uv, const T& upper_radius_, const T& bottom_radius_, const T& width_)
	{
		const T theta = 2.0f*PI*uv.v_;
		const T radius = bottom_radius_ + (upper_radius_ - bottom_radius_)*uv.u_;

		return TV(uv.u_*width_, radius*cos(theta), -radius*sin(theta));
	}

	const TV GetPrismXYZ(const Vector2D<T>& uv, const int& num_sides_, const T& upper_radius_, const T& bottom_radius_, const T& width_)
	{
		const T side_length = (T)1 / (T)num_sides_;
		const int side_index = (int)(uv.v_ / side_length);

		if (side_length < 0){
			std::cout << "GetPrismXYZ Error" << std::endl; exit(-1);
		}

		if (side_length >= num_sides_){
			std::cout << "GetPrismXYZ Error" << std::endl; exit(-1);
		}

		const T v_0 = side_length*(T)side_index;
		const T v_1 = side_length*(T)(side_index + 1);

		const T alpha = (uv.v_ - v_0) / side_length;

		const Vector2D<T> uv0(uv.u_, v_0);
		const Vector2D<T> uv1(uv.u_, v_1);

		const TV XYZ0 = GetCylinderXYZ(uv0, upper_radius_, bottom_radius_, width_);
		const TV XYZ1 = GetCylinderXYZ(uv1, upper_radius_, bottom_radius_, width_);

		return XYZ0*((T)1 - alpha) + XYZ1*alpha;
	}

	const TV GetPrismXYZLinearCorner(const Vector2D<T>& uv, const int& num_sides_, const float& corner_angle, const T& upper_radius_, const T& bottom_radius_, const T& width_)	// corner_angle is degree from user input.
	{
		const T side_length = (T)1 / (T)num_sides_;
		const T corner_length = corner_angle / (T)360;
		const int side_index = (int)(uv.v_ / side_length);

		if (side_length < 0)
		{
			std::cout << "GetPrismXYZ Error" << std::endl; exit(-1);
		}

		if (side_length >= num_sides_)
		{
			std::cout << "GetPrismXYZ Error" << std::endl; exit(-1);
		}

		const T v_0 = side_length*(T)side_index;

		const T relative_v = (uv.v_ - v_0);

		if (relative_v <= side_length - corner_length)		// non-corner side
		{
			const T v_1 = side_length*(T)(side_index + 1) - corner_length;

			const T alpha = (uv.v_ - v_0) / (side_length - corner_length);

			const Vector2D<T> uv0(uv.u_, v_0);
			const Vector2D<T> uv1(uv.u_, v_1);

			const TV XYZ0 = GetCylinderXYZ(uv0, upper_radius_, bottom_radius_, width_);
			const TV XYZ1 = GetCylinderXYZ(uv1, upper_radius_, bottom_radius_, width_);

			return XYZ0*((T)1 - alpha) + XYZ1*alpha;
		}
		else
		{
			const T v_0 = side_length*(T)(side_index + 1) - corner_length;	// new v_0
			const T v_1 = side_length*(T)(side_index + 1);

			const T alpha = (uv.v_ - v_0) / corner_length;

			const Vector2D<T> uv0(uv.u_, v_0);
			const Vector2D<T> uv1(uv.u_, v_1);

			const TV XYZ0 = GetCylinderXYZ(uv0, upper_radius_, bottom_radius_, width_);
			const TV XYZ1 = GetCylinderXYZ(uv1, upper_radius_, bottom_radius_, width_);

			return XYZ0*((T)1 - alpha) + XYZ1*alpha;
		}		
	}
};
