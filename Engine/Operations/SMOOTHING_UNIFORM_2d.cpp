// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "SMOOTHING_UNIFORM_2D.h"
#include "../CONVENTIONAL_MACROS.h"

void SMOOTHING_UNIFORM_2D::SmoothIClampJRepeat(const GridUniform2D& grid, Array2D<T>& height_field, const int smoothing_repeat)
{
	// height smoothing
	Array2D<T> height_field_temp;
	Array2D<bool> fixed;
	grid.InitializeNodeArray(height_field_temp);
	grid.InitializeNodeArray(fixed);

	for (int j = 0; j < fixed.j_res_; ++j)
	for (int i = 0; i < fixed.i_res_; ++i)
	{
		// 		if (height(i, j) == 1.0f) fixed(i, j) = true;
		//		else fixed(i, j) = false;
		fixed(i, j) = false;
	}

	const T coeff_fixed = 0.00f;
	const T coeff = 0.5f;

	for (int s = 0; s < smoothing_repeat; ++s)
	{
		for (int j = 0; j < height_field.j_res_ - 1; ++j)
		for (int i = 0; i < height_field.i_res_; ++i)
		{
			T sum = 0.0f;

			sum += height_field(MIN2(height_field.i_end_, i + 1), j);

			if (j + 1 >= height_field.j_end_) sum += height_field(i, 0);
			else sum += height_field(i, j + 1);

			sum += height_field(MAX2(0, i - 1), j);

			if (j - 1 < 0) sum += height_field(i, height_field.j_end_ - 1);
			else sum += height_field(i, j - 1);

			sum *= 0.25f;

			if (fixed(i, j) == true)
				height_field_temp(i, j) = height_field(i, j)*(1.0f - coeff_fixed) + sum*coeff_fixed;
			else
				height_field_temp(i, j) = height_field(i, j)*(1.0f - coeff) + sum*coeff;
		}

		for (int i = 0; i < height_field.i_res_; ++i)
		{
			height_field_temp(i, height_field.j_end_) = height_field_temp(i, 0);
		}

		height_field.copyFrom(height_field_temp);
	}
}


void SMOOTHING_UNIFORM_2D::SmoothIRepeatJClamp(const GridUniform2D& grid, Array2D<T>& height_field, const int smoothing_repeat)
{
	// height smoothing
	Array2D<T> height_field_temp;
	Array2D<bool> fixed;
	grid.InitializeNodeArray(height_field_temp);
	grid.InitializeNodeArray(fixed);

	for (int j = 0; j < fixed.j_res_; ++j)
	for (int i = 0; i < fixed.i_res_; ++i)
	{
// 		if (height(i, j) == 1.0f) fixed(i, j) = true;
//		else fixed(i, j) = false;
		fixed(i, j) = false;
	}

	const T coeff_fixed = 0.00f;
	const T coeff = 0.5f;

	for (int s = 0; s < smoothing_repeat; ++s)
	{
		for (int j = 0; j < height_field.j_res_; ++j)
		for (int i = 0; i < height_field.i_res_; ++i)
		{
			T sum = 0.0f;

			sum += height_field.getIRepeatedJClamped(i + 1, j);
			sum += height_field.getIRepeatedJClamped(i - 1, j);
			sum += height_field.getIRepeatedJClamped(i, j + 1);
			sum += height_field.getIRepeatedJClamped(i, j - 1);

			sum *= 0.25f;

			if (fixed(i, j) == true)
				height_field_temp(i, j) = height_field(i, j)*(1.0f - coeff_fixed) + sum*coeff_fixed;
			else
				height_field_temp(i, j) = height_field(i, j)*(1.0f - coeff) + sum*coeff;
		}

// 		for (int i = 0; i < height_field.i_res_; ++i)
// 		{
// 			height_field_temp(i, height_field.j_end_) = height_field_temp(i, 0);
// 		}

		height_field.copyFrom(height_field_temp);
	}
}



void SMOOTHING_UNIFORM_2D::SmoothIRepeatJClampCell(const GridUniform2D& grid, Array2D<T>& height_field, const int smoothing_repeat)
{
	// height smoothing
	Array2D<T> height_field_temp;
	Array2D<bool> fixed;
	grid.InitializeCellArray(height_field_temp);
	grid.InitializeCellArray(fixed);

	for (int j = 0; j < fixed.j_res_; ++j)
		for (int i = 0; i < fixed.i_res_; ++i)
		{
			// 		if (height(i, j) == 1.0f) fixed(i, j) = true;
			//		else fixed(i, j) = false;
			fixed(i, j) = false;
		}

	const T coeff_fixed = 0.00f;
	const T coeff = 0.5f;

	for (int s = 0; s < smoothing_repeat; ++s)
	{
		for (int j = 0; j < height_field.j_res_; ++j)
			for (int i = 0; i < height_field.i_res_; ++i)
			{
				T sum = 0.0f;

				sum += height_field.getIRepeatedJClamped(i + 1, j);
				sum += height_field.getIRepeatedJClamped(i - 1, j);
				sum += height_field.getIRepeatedJClamped(i, j + 1);
				sum += height_field.getIRepeatedJClamped(i, j - 1);

				sum *= 0.25f;

				if (fixed(i, j) == true)
					height_field_temp(i, j) = height_field(i, j)*(1.0f - coeff_fixed) + sum*coeff_fixed;
				else
					height_field_temp(i, j) = height_field(i, j)*(1.0f - coeff) + sum*coeff;
			}

		// 		for (int i = 0; i < height_field.i_res_; ++i)
		// 		{
		// 			height_field_temp(i, height_field.j_end_) = height_field_temp(i, 0);
		// 		}

		height_field.copyFrom(height_field_temp);
	}
}

void SMOOTHING_UNIFORM_2D::SmoothDirichletIJ(const GridUniform2D& grid, Array2D<T>& height_field, const int smoothing_repeat)
{
	// height smoothing
	Array2D<T> height_field_temp;
	Array2D<bool> fixed;
	grid.InitializeNodeArray(height_field_temp);
	grid.InitializeNodeArray(fixed);

	for (int j = 0; j < grid.j_res_; ++j)
	for (int i = 0; i < grid.i_res_; ++i)
	{
		// 		if (height(i, j) == 1.0f) fixed(i, j) = true;
		// 		else fixed(i, j) = false;
		fixed(i, j) = false;
	}

	const T coeff_fixed = 0.01f;
	const T coeff = 0.5f;
	for (int s = 0; s < smoothing_repeat; ++s)
	{
		for (int j = 0; j < grid.j_res_; ++j)
		for (int i = 0; i < grid.i_res_; ++i)
		{
			T sum = 0.0f;

			if (height_field.isValid(i + 1, j)) sum += height_field(i + 1, j);
			if (height_field.isValid(i, j + 1)) sum += height_field(i, j + 1);
			if (height_field.isValid(i - 1, j)) sum += height_field(i - 1, j);
			if (height_field.isValid(i, j - 1)) sum += height_field(i, j - 1);

			sum *= 0.25f;

			if (fixed(i, j) == true)
			{
				height_field_temp(i, j) = height_field(i, j)*(1.0f - coeff_fixed) + sum*coeff_fixed;
			}
			else
			{
				height_field_temp(i, j) = height_field(i, j)*(1.0f - coeff) + sum*coeff;
			}
		}

		height_field.copyFrom(height_field_temp);

		//		std::cout << "itr " << std::endl;
	}
}

void SMOOTHING_UNIFORM_2D::SmoothSignedDistance(const GridUniform2D& grid, Array2D<T>& height_field, const int smoothing_repeat, const T th)
{
	const T dx = 1.0f;
	const T h((T)1), hh(h*h), hh2(hh*(T)2), hh3(hh*(T)3);

	// height field to level set
	for (int j = 0; j < grid.j_res_; ++j)
	for (int i = 0; i < grid.i_res_; ++i)
	{
		height_field(i, j) -= th;
	}

	// height smoothing
	Array2D<T> height_field_temp;
	Array2D<bool> fixed;
	grid.InitializeNodeArray(height_field_temp);
	grid.InitializeNodeArray(fixed);	

	for (int s = 0; s < smoothing_repeat; ++s)
	{
		// re-distancing interfacial cells
		for (int j = 0; j < grid.j_res_; ++j)
		for (int i = 0; i < grid.i_res_; ++i)
		{
			if (fixed(i, j) == true) continue;

			T phi = height_field(i, j);
			T up, down, left, right;

			if (height_field.isValid(i + 1, j)) right = height_field(i + 1, j); else right = phi;
			if (height_field.isValid(i - 1, j)) left = height_field(i - 1, j); else left = phi;
			if (height_field.isValid(i, j + 1)) up = height_field(i, j + 1); else up = phi;
			if (height_field.isValid(i, j - 1)) down = height_field(i, j - 1); else down = phi;

			T a1, a2;
			T new_phi;

			if (phi >(T)0)
			{
				// See definitions of u^h_{x min}, u^h_{y min} in page 605.
				T a = MIN2(left, right);
				T b = MIN2(up, down);

				INCREASING_SORT2(a, b, a1, a2);

				new_phi = a1 + h;

				if (new_phi > a2) new_phi = 0.5f*(a1 + a2 + sqrt(hh2 - POW2(a1 - a2)));
				
				phi = MIN2(new_phi, phi);
			}
			else // phi <= (T)0
			{
				// See Equation (2.8) in page 607 for the implementation of negative parts.
				T a = MAX2(left, right);
				T b = MAX2(up, down);

				INCREASING_SORT2(a, b, a1, a2);

				new_phi = a2 - h;

				if (new_phi < a2) new_phi = 0.5f*(a2 + a1 - sqrt(hh2 - POW2(a2 - a1)));

				phi = MAX2(new_phi, phi);
			}

			height_field_temp(i, j) = phi;

//			if (height_field(i, j) != phi) std::cout << height_field(i, j) << " to " << phi << std::endl;
		}

		height_field.copyFrom(height_field_temp);
	}

	// height field to level set
	for (int j = 0; j < grid.j_res_; ++j)
	for (int i = 0; i < grid.i_res_; ++i)
	{
		height_field(i, j) += th;

// 		if (height_field(i, j) > th) height_field(i, j) = 1.0f;
// 		else height_field(i, j) = 0.0f;

//		height_field(i, j) = CLAMP(height_field(i, j), (T)0, (T)1);
	}
}