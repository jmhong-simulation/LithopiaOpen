// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "GENERIC_DEFINITIONS.h"
#include "DataStructure/Array1D.h"

class GridUniform1D
{
public:
	// grid resolution
	int i_res_;

	// start indices
	int i_start_;

	// end indices
	int i_end_;

	// grid domain (excluding ghost cells)
	T x_min_, x_max_;

	// grid spacing
	T dx_;

	// Inverse of grid spacing
	T one_over_dx_;

	// Inverse of grid spacing
	T one_over_2dx_;

public:// constructors and a destructor
	GridUniform1D()
	{}

	GridUniform1D(const int& i_start_input, const int& i_res_input, const T& x_min_input, const T& x_max_input)
	{
		initialize(i_start_input, i_res_input, x_min_input, x_max_input);
	}

	GridUniform1D(const GridUniform1D& grid_input)
	{
		initialize(grid_input.i_start_, grid_input.i_res_, grid_input.x_min_, grid_input.x_max_);
	}

	~GridUniform1D(void)
	{}

	void initialize(const int& i_start_input, const int& i_res_input, const T& x_min_input, const T& x_max_input);

	T getCellCenter(const int& i) const;
	T getNodePosition(const int& i) const;

	int getNumAllCells() const;
	int getNumAllNodes() const;


	template<class TT>
	T getLinearInterpolationNode(const Array1D<TT>& arr, const T pos)
	{
		const int i0 = (int)floor((pos - x_min_)*one_over_dx_) + i_start_;

		const T q0 = (pos - getCellCenter(i0)) * one_over_dx_;
	
		const T v0 = arr[i0];
		const T v1 = arr[i0 + 1];

		return v0 * ((T)1 - q0) + v1 * q0;
	}

	template<class TT>
	void initializeCellArray(Array1D<TT>& cell_array) const
	{
		cell_array.initialize(i_res_, false);
	}

	template<class TT>
	void initializeNodeArray(Array1D<TT>& node_array) const
	{
		node_array.initialize(i_res_ + 1, false);
	}
};