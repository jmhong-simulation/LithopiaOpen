// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "Parallelism/Multithreading.h"
#include "GridUniform2D.h"
#include "Array2D.h"

class ArrayTools
{
public:
	template<class TT>
	static void fillGhostCells(MT* mt, const int thread_id, const GridUniform2D& grid, const GridUniform2D& grid_ghost, Array2D<TT>& arr)
	{
		BEGIN_GRID_ITERATION_2D(grid_ghost)
		{
			arr(i, j) = arr(grid.ClampedIndex(i, j));		//TODO: don't need to copy non-ghost cells
		}
		END_GRID_ITERATION_2D;
	}

	template<class TT>
	static void fillGhostCells(const GridUniform2D& grid, const GridUniform2D& grid_ghost_, Array2D<TT>& arr)
	{
		for (int j = grid_ghost_.j_start_; j <= grid_ghost_.j_end_; j++)
			for (int i = grid_ghost_.i_start_; i <= grid_ghost_.i_end_; i++)
			{
				arr(i, j) = arr(grid.ClampedIndex(i, j));		//TODO: don't need to copy non-ghost cells
			}
	}

	template<class TT>
	static void copyArray(MT* mt, const int thread_id, const Array2D<TT>& from, Array2D<TT>& to)
	{
		BEGIN_ONE_THREAD_WORK
		{
			if (from.getNumAllValues() != to.getNumAllValues())
				to.initialize(from.i_start_, from.j_start_, from.i_res_, from.j_res_, false);
		}
		END_ONE_THREAD_WORK;

		BEGIN_1D_ITERATION(from.getNumAllValues())
		{
			to.values_[p] = from.values_[p];
		}
		END_1D_ITERATION;
	}

	template<class TT>
	static void smoothLaplacian(MT* mt, const int thread_id, const T& alpha, const int& repeat, Array2D<TT>& arr)
	{
		for (int r = 0; r < repeat; r++)
		{
			BEGIN_ARRAY_ITERATION_2D(arr)
			{
				const T average = (arr.getClamped(i - 1, j) + arr.getClamped(i + 1, j) + arr.getClamped(i, j - 1) + arr.getClamped(i, j + 1)) * (T)0.25;

				arr(i, j) = arr(i, j) * ((T)1 - alpha) + average * alpha;
			}
			END_ARRAY_ITERATION_2D;
		}
	}

	template<class TT>
	static void clampMinMaxZeroOne(MT* mt, const int thread_id, const TT& min, const TT& max, Array2D<TT>& arr)
	{
		BEGIN_ARRAY_ITERATION_2D(arr)
		{
			T value = arr(i, j);

			if (value < min) value = 0.0f;
			if (value > max) value = 1.0f;

			arr(i, j) = value;
		}
		END_ARRAY_ITERATION_2D;
	}


	template<class TT>
	static void addAllValues(MT* mt, const int thread_id, const T v_add, Array2D<TT>& arr)
	{
		BEGIN_1D_ITERATION(arr.getNumAllValues())
		{
			arr.values_[p] += v_add;
		}
		END_1D_ITERATION;
	}
};
