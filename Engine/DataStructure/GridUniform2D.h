// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include <iostream>
#include "Array1D.h"
#include "Array2D.h"
#include "CONVENTIONAL_MACROS.h"
#include "GENERIC_DEFINITIONS.h"
#include "Parallelism/MultiThreading.h"
#include "Geometry/BOX_2D.h"

class GridUniform2D
{
public:
	// grid resolution
	int i_res_, j_res_;

	// start indices
	int i_start_, j_start_;

	// end indices
	int i_end_, j_end_;

	// grid domain (excluding ghost cells)
	T x_min_, y_min_, x_max_, y_max_;

	// grid spacing
	T dx_, dy_;

	// Inverse of grid spacing
	T one_over_dx_, one_over_dy_;

	// Inverse of grid spacing
	T one_over_2dx_, one_over_2dy_;

public:// constructors and a destructor
	GridUniform2D()
	{}

	GridUniform2D(const int& i_start_input, const int& j_start_input, const int& i_res_input,const int& j_res_input, const T& x_min_input,const T& y_min_input, const T& x_max_input,const T& y_max_input)
	{
		Initialize(i_start_input, j_start_input, i_res_input, j_res_input, x_min_input, y_min_input, x_max_input, y_max_input);
	}

	GridUniform2D(const TV2_INT& ij_start_input, const TV2_INT& ij_res_input, const TV2& xy_min_input, const TV2& xy_max_input)
	{
		Initialize(ij_start_input.i_, ij_start_input.j_, ij_res_input.i_, ij_res_input.j_, xy_min_input.x_, xy_min_input.y_, xy_max_input.x_, xy_max_input.y_);
	}

	GridUniform2D(const GridUniform2D& grid_input)
	{
		Initialize(grid_input.i_start_, grid_input.j_start_, grid_input.i_res_, grid_input.j_res_, grid_input.x_min_, grid_input.y_min_, grid_input.x_max_, grid_input.y_max_);
	}

	~GridUniform2D(void)
	{}

public:
	void Initialize(const int& i_start_input, const int& j_start_input, const int& i_res_input,const int& j_res_input, const T& x_min_input,const T& y_min_input, const T& x_max_input,const T& y_max_input)
	{
		i_start_ = i_start_input;
		j_start_ = j_start_input;

		i_res_ = i_res_input;
		j_res_ = j_res_input;

		x_min_ = x_min_input;
		y_min_ = y_min_input;

		x_max_ = x_max_input;
		y_max_ = y_max_input;
		
		i_end_ = i_start_ + i_res_ - 1;
		j_end_ = j_start_ + j_res_ - 1;

		dx_ = (x_max_-x_min_)/(T)i_res_;
		dy_ = (y_max_-y_min_)/(T)j_res_;

		one_over_dx_ = (T)1/dx_;
		one_over_dy_ = (T)1/dy_;

		one_over_2dx_ = (T)0.5/dx_;
		one_over_2dy_ = (T)0.5/dy_;
	}

	void Initialize(const GridUniform2D& grid_input)
	{
		Initialize(grid_input.i_start_, grid_input.j_start_, grid_input.i_res_, grid_input.j_res_, grid_input.x_min_, grid_input.y_min_, grid_input.x_max_, grid_input.y_max_);
	}

	void InitializeDualGrid(GridUniform2D& dual_grid_input)
	{
		dual_grid_input.Initialize(i_start_, j_start_, i_res_-1, j_res_-1, x_min_+(T)0.5*dx_, y_min_+(T)0.5*dy_, x_max_-(T)0.5*dx_, y_max_-(T)0.5*dy_);	
	}

	template<class TT>
	void InitializeCellArray(Array2D<TT>& cell_array) const
	{
		cell_array.initialize(i_start_, j_start_, i_res_, j_res_, false);
	}

	template<class TT>
	void InitializeNodeArray(Array2D<TT>& node_array) const
	{
		node_array.initialize(i_start_, j_start_, i_res_ + 1, j_res_ + 1, false);
	}

	template<class TT>
	void InitializeUFaceArray(Array2D<TT>& x_edge_array) const		// U face 
	{
		x_edge_array.initialize(i_start_, j_start_, i_res_ + 1, j_res_, false);
	}

	template<class TT>
	void InitializeVFaceArray(Array2D<TT>& y_edge_array) const		// U face 
	{
		y_edge_array.initialize(i_start_, j_start_, i_res_, j_res_ + 1, false);
	}

public:
	TV2_INT ClampedIndex(const int& i, const int& j) const
	{
		return TV2_INT(CLAMP(i,i_start_,i_end_), CLAMP(j,j_start_,j_end_));
	}

	TV2_INT ClampedIndex(const TV2_INT& ix) const
	{
		return TV2_INT(CLAMP(ix.i_,i_start_,i_end_), CLAMP(ix.j_,j_start_,j_end_));
	}

	int Get1DIndex(const int& i, const int& j) const
	{
		return (i-i_start_) + (j-j_start_)*i_res_;
	}

	int ClampedIndexI(const int& i) const
	{
		return CLAMP(i,i_start_,i_end_);
	}

	int ClampedIndexJ(const int& j) const
	{
		return CLAMP(j,j_start_,j_end_);
	}

public:
	TV2 GetCellCenter(const int& i, const int& j) const
	{
		return TV2(x_min_ + ((T)0.5+(T)(i-i_start_))*dx_, y_min_ + ((T)0.5+(T)(j-j_start_))*dy_);
	}

	TV2 GetCellCenter(const TV2_INT& ix) const
	{
		return TV2(x_min_ + ((T)0.5+(T)(ix.i_-i_start_))*dx_, y_min_ + ((T)0.5+(T)(ix.j_-j_start_))*dy_);
	}

	void Center(const int& i, const int& j, TV2& position) const
	{
		position.x_ = x_min_ + ((T)0.5+(T)(i-i_start_))*dx_;
		position.y_ = y_min_ + ((T)0.5+(T)(j-j_start_))*dy_;
	}

	void Center(const TV2_INT& ix, TV2& position) const
	{
		position.x_ = x_min_ + ((T)0.5+(T)(ix.i_-i_start_))*dx_;
		position.y_ = y_min_ + ((T)0.5+(T)(ix.j_-j_start_))*dy_;
	}

	TV2_INT Cell(const TV2& position) const // return a cell index contains the position
	{
		//TODO: check!!!
		return TV2_INT((int)((position.x_-x_min_)*one_over_dx_)+i_start_, (int)((position.y_-y_min_)*one_over_dy_)+j_start_);
	}

	TV2_INT getLeftBottomCell(const TV2& position) const //TODO check performance compared to no return value functions.
	{
		return TV2_INT(i_start_ + (int)floor((position.x_ - x_min_)*one_over_dx_ - (T)0.5),
			j_start_ + (int)floor((position.y_ - y_min_)*one_over_dy_ - (T)0.5));
	}

	TV2_INT getRightTopCell(const TV2& position) const
	{
		return TV2_INT(i_start_ + (int)floor((position.x_ - x_min_)*one_over_dx_ + (T)0.5),
			j_start_ + (int)floor((position.y_ - y_min_)*one_over_dy_ + (T)0.5));
	}

	void Cell(const TV2& position, int& i, int& j) const // return a cell index contains the position
	{
		i = (int)((position.x_-x_min_)*one_over_dx_);
		j = (int)((position.y_-y_min_)*one_over_dy_);

		assert(i >= i_start_ && i <= i_end_);
		assert(j >= j_start_ && j <= j_end_);
	}

	TV2_INT ClampedCell(const TV2& position) const
	{
		TV2_INT index((int)((position.x_-x_min_)*one_over_dx_), (int)((position.y_-y_min_)*one_over_dy_));//TODO: check!!!

		index.i_ = CLAMP(index.i_, i_start_, i_end_);
		index.j_ = CLAMP(index.j_, j_start_, j_end_);

		return index;
	}

	TV2_INT ClampedCell(const TV2& position, TV2_INT& index) const
	{
		index.assign((int)((position.x_-x_min_)*one_over_dx_), (int)((position.y_-y_min_)*one_over_dy_));//TODO: check!!!

		index.i_ = CLAMP(index.i_, i_start_, i_end_);
		index.j_ = CLAMP(index.j_, j_start_, j_end_);

		return index;
	}

	void ClampedCell(const TV2& position, int& i, int& j) const
	{
		i = (int)((position.x_-x_min_)*one_over_dx_);
		j = (int)((position.y_-y_min_)*one_over_dy_);

		i = CLAMP(i, i_start_, i_end_);
		j = CLAMP(j, j_start_, j_end_);
	}

	BOX_2D<int> getIXBox(const BOX_2D<T>& aabb)
	{
		return BOX_2D<int>(getLeftBottomCell(aabb.getMin()), getRightTopCell(aabb.getMax()));
	}	

	void LeftBottomCell(const TV2& position, int& i, int& j) const
	{
		i = (int)((position.x_-x_min_)*one_over_dx_ - (T)0.5);//TODO: check if this works with 64 bit compiler. maybe we need to use floor instead of int casting.
		j = (int)((position.y_-y_min_)*one_over_dy_ - (T)0.5);
	}

	void LeftBottomCell(const TV2& position, TV2_INT& ix) const
	{
		return LeftBottomCell(position, ix.i_, ix.j_);
	}

	bool Inside(const TV2& position) const
	{
		if(position.x_ <= x_min_) return false;
		else if(position.x_ >= x_max_) return false;
		else if(position.y_ <= y_min_) return false;
		else if(position.y_ >= y_max_) return false;
		return true;
	}

	bool Inside(const TV2& position, const T& width) const
	{
		if(position.x_ <= x_min_+width) return false;
		else if(position.x_ >= x_max_-width) return false;
		else if(position.y_ <= y_min_+width) return false;
		else if(position.y_ >= y_max_-width) return false;
		return true;
	}

	bool Inside(const int& i, const int& j) const
	{
		if(i < i_start_) return false;
		else if(i > i_end_) return false;
		else if(j < j_start_) return false;
		else if(j > j_end_) return false;
		else return true;
	}

	bool Inside(const TV2_INT& ix) const
	{
		if(ix.i_ < i_start_) return false;
		else if(ix.i_ > i_end_) return false;
		else if(ix.j_ < j_start_) return false;
		else if(ix.j_ > j_end_) return false;
		else return true;
	}

	bool Inside(const TV2_INT& ix, const int& inner_width) const
	{
		if(ix.i_ < i_start_ + inner_width) return false;
		else if(ix.i_ > i_end_ - inner_width) return false;
		else if(ix.j_ < j_start_ + inner_width) return false;
		else if(ix.j_ > j_end_ - inner_width) return false;
		else return true;
	}

	bool InsideI(const int& i) const
	{
		if(i < i_start_) return false;
		else if(i > i_end_) return false;
		return true;
	}

	bool InsideJ(const int& j) const
	{
		if(j < j_start_) return false;
		else if(j > j_end_) return false;
		return true;
	}

	bool isBoundary(const int& i, const int& j) const
	{
		if (i == i_start_) return true;
		else if (i == i_end_) return true;
		else if (j == j_start_) return true;
		else if (j == j_end_) return true;

		return false;
	}

	GridUniform2D getEnlarged(const int& width) const
	{
		return GridUniform2D(i_start_-width, j_start_-width, i_res_+2*width, j_res_+2*width, x_min_-(T)width*dx_, y_min_-(T)width*dy_, x_max_+(T)width*dx_, y_max_+(T)width*dy_);
	}

	GridUniform2D getNodeGrid() const
	{
		return GridUniform2D(i_start_, j_start_, i_res_ + 1, j_res_ + 1, x_min_ - (T)0.5*dx_, y_min_ - (T)0.5*dy_, x_max_ + (T)0.5*dx_, y_max_ + (T)0.5*dy_);
	}

	GridUniform2D getCellGrid() const
	{
		return GridUniform2D(i_start_, j_start_, i_res_ - 1, j_res_ - 1, x_min_ + (T)0.5*dx_, y_min_ + (T)0.5*dy_, x_max_ - (T)0.5*dx_, y_max_ - (T)0.5*dy_);
	}

	void Enlarge(const int& width)
	{
		Initialize(i_start_-width, j_start_-width, i_res_+2*width, j_res_+2*width, x_min_-(T)width*dx_, y_min_-(T)width*dy_, x_max_+(T)width*dx_, y_max_+(T)width*dy_);
	}

	void operator = (const GridUniform2D& grid_input)
	{
		Initialize(grid_input);
	}

	void SplitInHeight(const int& num_threads, Array1D<GridUniform2D>& partial_grids)
	{
		partial_grids.initialize(num_threads);

		int quotient = j_res_ / num_threads;
		int remainder = j_res_ % num_threads;
		int j_start_p = j_start_;
		T y_min_p = y_min_;

		for(int i=0; i<num_threads; i++)
		{
			int y_height = i < remainder ? quotient+1 : quotient;

			T y_max_p = dy_ * y_height + y_min_p;

			partial_grids[i].Initialize(i_start_, j_start_p, i_res_, y_height, x_min_, y_min_p, x_max_, y_max_p);

			j_start_p += y_height;
			y_min_p = y_max_p;
		}
	}

	int GetNumAllNodes() const
	{
		return (i_res_ + 1)*(j_res_ + 1);
	}

	int getNumAllCells() const
	{
		return i_res_ * j_res_;
	}

	const Vector2D<T> GetNodeUV(const int& i, const int& j) const
	{
		return Vector2D<T>((T)i / (T)i_res_, (T)j / (T)j_res_);	// For node indices (i, j), x_res_ is correct. (Not x_res_ - 1)
	}

	const int GetCell1DIndex(const int& i, const int& j) const // 2D index to 1D Array1D index
	{
		assert(i >= i_start_ && i <= i_end_);
		assert(j >= j_start_ && j <= j_end_);

		return (i - i_start_) + (j - j_start_)*i_res_;//TODO: try pointer operation for optimization (check performance)
	}

	const TV2 GetNode(const int& i, const int& j) const
	{
		return TV2(x_min_ + (T)(i - i_start_)*dx_, y_min_ + (T)(j - j_start_)*dy_);
	}

	template<class TT> const TT GetClampedLinearInterpolationCell(const Array2D<TT>& arr_cell, const TV2& position) const
	{
		const T& posx = position.x_;
		const T& posy = position.y_;

		const int i0 = (int)floor((posx - x_min_)*one_over_dx_ - (T)0.5) + i_start_;
		const int j0 = (int)floor((posy - y_min_)*one_over_dy_ - (T)0.5) + j_start_;

		const T a = (posx - (x_min_ + ((T)0.5 + (T)(i0 - i_start_))*dx_))*one_over_dx_;
		const T b = (posy - (y_min_ + ((T)0.5 + (T)(j0 - j_start_))*dy_))*one_over_dy_;

		const T	w00 = ((T)1 - a)*((T)1 - b), w10 = a*((T)1 - b), w01 = ((T)1 - a)*b, w11 = a*b;

		const TT v00 = arr_cell.getClamped(i0, j0), v10 = arr_cell.getClamped(i0 + 1, j0), v01 = arr_cell.getClamped(i0, j0 + 1), v11 = arr_cell.getClamped(i0 + 1, j0 + 1);

		return v00 * w00 + v10 * w10 + v01 * w01 + v11 * w11;
	}

	template<class TT> const TT getLinearInterpolationCell(const Array2D<TT>& arr_cell, const TV2& position) const
	{
		//TODO: fix!!!!

		const T& posx = position.x_;
		const T& posy = position.y_;

		const int i0 = (int)floor((posx - x_min_)*one_over_dx_ - (T)0.5) + i_start_;
		const int j0 = (int)floor((posy - y_min_)*one_over_dy_ - (T)0.5) + j_start_;

		TV2 q0 = GetCellCenter(i0, j0);
		TV2 q1 = GetCellCenter(i0 + 1, j0);
		TV2 q2 = GetCellCenter(i0, j0 + 1);
		TV2 q3 = GetCellCenter(i0 + 1, j0 + 1);

		TT& v0 = arr_cell(i0, j0);
		TT& v1 = arr_cell(i0 + 1, j0);
		TT& v2 = arr_cell(i0, j0 + 1);
		TT& v3 = arr_cell(i0 + 1, j0 + 1);

		T nx = posx - q0.x_;
		T ny = posy - q0.y_;

		TT vv0 = v0 + (v1 - v0) * nx / dx_;
		TT vv1 = v3 + (v2 - v3) * nx / dx_;

		return vv0 + (vv1 - vv0) * ny / dy_;
	}

	template<class TT> const TT GetRepeatedLinearInterpolationCell(const Array2D<TT>& arr_cell, const TV2& position) const
	{
		T posx = position.x_;
		T posy = position.y_;

		T cxMin = x_min_ + dx_*0.55f;
		T cyMin = y_min_ + dy_*0.55f;
		T cxMax = x_max_ - dx_*0.55f;
		T cyMax = y_max_ - dy_*0.55f;

		// clamp and repeat
		if (posx < cxMin) {
			posx = cxMax - ABS(cxMin - posx);
		}
		else if (posx > cxMax) {
			posx = cxMin + ABS(posx - cxMax);
		}
		
		if (posy < cyMin) {
			posy = cyMax - ABS(cyMin - posy);
		}
		else if (posy > cyMax) {
			posy = cyMin + ABS(posy - cyMax);
		}

		//
		const int i0 = (int)floor((posx - x_min_)*one_over_dx_ - (T)0.5) + i_start_;
		const int j0 = (int)floor((posy - y_min_)*one_over_dy_ - (T)0.5) + j_start_;

		TV2 q0 = GetCellCenter(i0, j0);
		TV2 q1 = GetCellCenter(i0+1, j0);
		TV2 q2 = GetCellCenter(i0, j0+1);
		TV2 q3 = GetCellCenter(i0+1, j0+1);

		TT& v0 = arr_cell(i0, j0);
		TT& v1 = arr_cell(i0+1, j0);
		TT& v2 = arr_cell(i0, j0+1);
		TT& v3 = arr_cell(i0+1, j0+1);

		T nx = posx - q0.x_;
		T ny = posy - q0.y_;
		
		TT vv0 = v0 + (v1 - v0) * nx / dx_;
		TT vv1 = v3 + (v2 - v3) * nx / dx_;
		
		return vv0 + (vv1 - vv0) * ny / dy_;
	}
};

inline std::ostream&
operator << (std::ostream& output,const GridUniform2D& grid)
{
	return output << "GRID_UNIFORM_3D ["
				  << "Resolution = "<<grid.i_res_<<" "<<grid.j_res_
		          << " Index range = ("<<grid.i_start_<<" "<<grid.j_start_<<") to ("<< grid.i_end_<<" "<<grid.j_end_<<") "
				  << " Range = ("<< grid.x_min_<<" "<<grid.y_min_<<") to ("<<grid.x_max_<<" "<<grid.y_max_<<") "
				  << " DX = " <<grid.dx_<<" "<<grid.dy_<<"]";
}

