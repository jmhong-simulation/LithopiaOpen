// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "GridUniform1D.h"
#include "Array1D.h"
#include "Geometry/PARAMETRIC_CURVE_SEGMENT.h"

class NonlinearMapping1D
{
public:
	GridUniform1D grid_;
	Array1D<T> u2t_node_array_;	// maps texture coordinate to value coordinate
	Array1D<T> vel_node_array_;	// node velocity for iteration
	Array1D<T> length_cell_array_;

	void initialize(const int num_samples, const T u_min, const T u_max)
	{
		grid_.initialize(0, num_samples, u_min, u_max);

		grid_.initializeNodeArray(u2t_node_array_);

		grid_.initializeNodeArray(vel_node_array_);

		grid_.initializeCellArray(length_cell_array_);

		for (int i = 0; i < grid_.getNumAllNodes(); i++)
		{
			u2t_node_array_[i] = grid_.getNodePosition(i);		// initial condition

			vel_node_array_[i] = (T)0;
		}
	}

	T getT(const T& u)
	{
		return grid_.getLinearInterpolationNode(u2t_node_array_, u);
	}

	void normalizeByLength(const ParametricCurveSegment& curve)
	{
		const int num_repeat = 1000000;
		const T cfl = (T) 0.025;

		for (int repeat = 0; repeat < num_repeat; repeat++)
		{
			T length_sum = (T)0;
			for (int i = 0; i < grid_.getNumAllCells(); i++)
			{
				const T length = curve.getLength(u2t_node_array_[i], u2t_node_array_[i + 1]);

				length_cell_array_[i] = length;

				length_sum += length;
			}

			//			std::cout << length_cell_array_ << std::endl;

			for (int i = 0; i < grid_.getNumAllCells(); i++)
			{
				length_cell_array_[i] /= length_sum;
			}

			// 			length_cell_array_.accumulateFromLeft();
			// 
			// 			for (int i = 0; i < grid_.getNumAllCells(); i++)
			// 			{
			// 				u2t_node_array_[i+1] = length_cell_array_[i];
			// 			}

			//			std::cout << u2t_node_array_ << std::endl;

			for (int i = 1; i < grid_.getNumAllNodes() - 1; i++)
			{
				const T vel = length_cell_array_[i] - length_cell_array_[i - 1];

				vel_node_array_[i] = vel;
			}

			// update nodes
			for (int i = 0; i < grid_.getNumAllNodes(); i++)
			{
				u2t_node_array_[i] += vel_node_array_[i] * cfl;
			}
		}

		std::cout << length_cell_array_ << std::endl;
	}

	void normalizeByLengthOverX(const ParametricCurveSegment& curve)
	{
		const int num_repeat = 1000000;
		const T cfl = (T) 0.025;

		for (int repeat = 0; repeat < num_repeat; repeat++)
		{
			T length_sum = (T)0;
			for (int i = 0; i < grid_.getNumAllCells(); i++)
			{
//				const T length = curve.getLength(u2t_node_array_[i], u2t_node_array_[i + 1]);
				const T length = curve.getLength(u2t_node_array_[i], u2t_node_array_[i + 1]) / curve.x_.getValue((u2t_node_array_[i + 1] + u2t_node_array_[i]) * (T)0.5);

				length_cell_array_[i] = length;

				length_sum += length;
			}
			
//			std::cout << length_cell_array_ << std::endl;

			for (int i = 0; i < grid_.getNumAllCells(); i++)
			{
				length_cell_array_[i] /= length_sum;
			}

// 			length_cell_array_.accumulateFromLeft();
// 
// 			for (int i = 0; i < grid_.getNumAllCells(); i++)
// 			{
// 				u2t_node_array_[i+1] = length_cell_array_[i];
// 			}

//			std::cout << u2t_node_array_ << std::endl;

			for (int i = 1; i < grid_.getNumAllNodes() - 1; i++)
			{
				const T vel = length_cell_array_[i] - length_cell_array_[i - 1];

				vel_node_array_[i] = vel;
			}

			// update nodes
			for (int i = 0; i < grid_.getNumAllNodes(); i++)
			{
				u2t_node_array_[i] += vel_node_array_[i] * cfl;
			}
		}		

		std::cout << length_cell_array_ << std::endl;
	}

	void normalizeBySOXRArea(const ParametricCurveSegment& curve)
	{
		//		std::cout << u2t_node_array_ << std::endl;

		T area_sum = (T)0;

		for (int i = 0; i < grid_.getNumAllCells(); i++)
		{
			const T length = curve.getLength(u2t_node_array_[i], u2t_node_array_[i + 1]);
			const T middle_radius = curve.getLength(u2t_node_array_[i], u2t_node_array_[i + 1]) * curve.x_.getValue((u2t_node_array_[i] + u2t_node_array_[i + 1]) * (T)0.5);
			const T area = length / middle_radius;

			area_sum += area;

			length_cell_array_[i] = area;
		}

//		std::cout << length_cell_array_ << std::endl;

		for (int i = 0; i < grid_.getNumAllCells(); i++)
		{
			length_cell_array_[i] = length_cell_array_[i] / area_sum * (grid_.x_max_ - grid_.x_min_);
		}

//		std::cout << length_cell_array_ << std::endl;

		length_cell_array_.accumulateFromLeft();

		//		std::cout << length_cell_array_ << std::endl;

		for (int i = 0; i < grid_.getNumAllCells(); i++)
		{ 
			u2t_node_array_[i + 1] = (u2t_node_array_[i + 1] + length_cell_array_[i]) * (T)0.5;
		}

		std::cout << u2t_node_array_ << std::endl;
	}
};