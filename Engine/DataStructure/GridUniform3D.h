// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "../GENERIC_DEFINITIONS.h"
#include "../CONVENTIONAL_MACROS.h"
//#include "SCRIPT_READER.h"
#include "Array1D.h"
#include "../Geometry/BOX_3D.h"
#include "Array3D.h"
#include <iostream>

class GridUniform3D
{
public:
	typedef BOX_3D<T> BOX;

	union{// grid resolution
		struct{int i_res_, j_res_, k_res_;};
		int res_[3];};

	union{// start indices
		struct{int i_start_, j_start_, k_start_;};
		int ix_start_[3];};

	union{// end indices
		struct{int i_end_, j_end_, k_end_;};
		int ix_end_[3];};

	union{// grid domain (excluding ghost cells)
		struct{T x_min_, y_min_, z_min_, x_max_, y_max_, z_max_;};
		struct{T min_[3], max_[3];};};

	union{// grid spacing
		struct{T dx_, dy_, dz_;};
		T dw_[3];};

	union{// Inverse of grid spacing
		struct{T one_over_dx_, one_over_dy_, one_over_dz_;};
		T one_over_dw_[3];};

	union{// Inverse of grid spacing
		struct{T one_over_2dx_, one_over_2dy_, one_over_2dz_;};
		T one_over_2dw_[3];};

	// speed up constants
	int ij_res_, ijk_res_;

public:// constructors and a destructor
	GridUniform3D();
	GridUniform3D(const GridUniform3D& grid_input);
	GridUniform3D(const BOX_3D<T>& bb, const TV_INT& res);
	GridUniform3D(const TV_INT& ijk_start_input, const TV_INT& ijk_res_input, const TV& xyz_min_input, const TV& xyz_max_input);
	GridUniform3D(const int& i_start_input, const int& j_start_input, const int& k_start_input, const int& i_res_input, const int& j_res_input, const int& k_res_input, const T& x_min_input, const T& y_min_input, const T& z_min_input, const T& x_max_input, const T& y_max_input, const T& z_max_input);

	~GridUniform3D(void);

	void initialize(const int& i_start_input, const int& j_start_input, const int& k_start_input, const int& i_res_input, const int& j_res_input, const int& k_res_input, const T& x_min_input, const T& y_min_input, const T& z_min_input, const T& x_max_input, const T& y_max_input, const T& z_max_input);
	void initialize(const TV_INT& start_input, const TV_INT& res_input, const TV& min_input, const TV& max_input);
	void initialize(const GridUniform3D& grid_input);
	void initialize(const GridUniform3D& grid_input, const T& resolution_scale);
	void initialize(const GridUniform3D& grid_input, const TV_INT& resolution_scale);
	void initialize(const BOX& box_input, const int& max_res);
	void initialize(const BOX& box_input, const T& dx_input);
	void initialize(const BOX& box_input, const int& i_res, const int& j_res, const int& k_res);
	void InitializeDualGrid(GridUniform3D& dual_grid_input);

	template<class TT>
	void initializeCenterArray(Array3D<TT>& _array)
	{
		_array.Initialize(i_start_, j_start_, k_start_, i_res_, j_res_, k_res_);
	}

	template<class TT>
	void initializeUArray(Array3D<TT>& _array)
	{
		_array.Initialize(i_start_, j_start_, k_start_, i_res_+1, j_res_, k_res_);
	}

	template<class TT>
	void initializeVArray(Array3D<TT>& _array)
	{
		_array.Initialize(i_start_, j_start_, k_start_, i_res_, j_res_ + 1, k_res_);
	}

	template<class TT>
	void initializeWArray(Array3D<TT>& _array)
	{
		_array.Initialize(i_start_, j_start_, k_start_, i_res_, j_res_, k_res_ + 1);
	}

// 	void InitializeFromBlock(const SCRIPT_BLOCK& block)
// 	{
// 		const T resolution_scale = block.GetFloat("resolution_scale", (T)1);
// 		const TV_INT start = block.GetInt3("start_indices");
// 		const TV_INT res = block.GetInt3("base_grid_resolution");
// 		const TV_INT res_scaled((int)((T)res.x_*resolution_scale), (int)((T)res.y_*resolution_scale), (int)((T)res.z_*resolution_scale));
// 		const TV min = block.GetVector3("base_grid_min");
// 		const TV max = block.GetVector3("base_grid_max");
// 
// 		Initialize(start.i_, start.j_, start.k_, res_scaled.i_, res_scaled.j_, res_scaled.k_, min.i_, min.j_, min.k_, max.i_, max.j_, max.k_);
// 	}

public:
	TV_INT getClampedIndex(const int& i, const int& j, const int& k) const;
	TV_INT getClampedIndex(const TV_INT& ix) const;

	void getClampStartEndIndices(int& l_start, int& m_start, int& n_start, int& l_end, int& m_end, int& n_end);

	int get1DIndex(const int& i, const int& j, const int& k) const;

	TV_INT get3DIndex(const int& index_1d) const;

	int getClampedIndexI(const int& i) const;
	int getClampedIndexJ(const int& j) const;
	int getClampedIndexK(const int& k) const;

	TV getCellCenter(const int& i, const int& j, const int& k) const;
	TV getCellCenter(const TV_INT& ix) const;

	void getCenter(const int& i, const int& j, const int& k, TV& position);
	void getCenter(const TV_INT& ix, TV& position) const;
	void getCell(const TV& position, TV_INT& index) const; // return a cell index contains the position
	void getCell(const TV& position, int& i, int& j, int& k) const; // return a cell index contains the position
	void getCell(const T& x, const T& y, const T& z, int& i, int& j, int& k) const; // return a cell index contains the position
	TV_INT getCell(const TV& position) const; // return a cell index contains the position

	void getClampedCell(const TV& position, TV_INT& index) const;
	void getClampedCell(const TV& position, int& i, int& j, int& k) const;
	TV_INT getClampedCell(const TV& position) const;

	void getLeftBottomCell(const TV& position, int& i, int& j, int& k) const;
	void getLeftBottomCell(const T& x, const T& y, const T& z, int& i, int& j, int& k) const;
	void getLeftBottomCell(const TV& position, TV_INT& ix) const;
	TV_INT getLeftBottomCell(const TV& position) const; //TODO check performance compared to no return value functions.
	
	TV_INT getRightTopCell(const TV& position) const;
	
	int getNumAllCells() const;

	bool isInside(const TV& position) const;
	bool isInside(const TV& position, const T& width) const;
	bool isInside(const int& i, const int& j, const int& k) const;
	bool isInside(const TV_INT& ix) const;
	bool isInside(const TV_INT& ix, const int& inner_width) const;
	bool isInsideI(const int & i);	// check if i_start_ <= i <= i_end_
	bool isInsideJ(const int & j);	// check if j_start_ <= j <= j_end_
	bool isInsideK(const int & k);	// check if k_start_ <= k <= k_end_
	
	GridUniform3D getEnlarged(const int& width) const;
	GridUniform3D getUFaceGrid() const;
	GridUniform3D getVFaceGrid() const;
	GridUniform3D getWFaceGrid() const;
	GridUniform3D getXEdgeGrid() const;	// getCellCenter returns the centers of x-direction edges
	GridUniform3D getYEdgeGrid() const;
	GridUniform3D getZEdgeGrid() const;
	GridUniform3D getCellGrid() const;	// returns cell grid when this is a node grid
	GridUniform3D getPartialGrid(const BOX_3D<int>& partial_ix_box) const;

	void enlarge(const int& width);
	void extend(const int width);
	void translate(const TV& deviation);
	void scaleToBox(const BOX& box);
	void scaleRes(const T scale);

	void operator = (const GridUniform3D& grid_input);
	bool operator == (const GridUniform3D& grid);	// for debugging only (not optimized)
	
	TV getMin() const;
	TV getMax() const;
	TV getCellMin(const int& i, const int& j, const int& k) const;
	TV getCellMax(const int& i, const int& j, const int& k) const;

	BOX_3D<T> getMimMax() const;
	BOX_3D<T> getCellMinMax(const int& i, const int& j, const int& k) const;
	BOX_3D<int> getIXBox() const;

	void splitInXDirection(const int& num_threads, Array1D<GridUniform3D>& partial_grids);// make partial grids by splitting this grid in X direction
	void splitInYDirection(const int& num_threads, Array1D<GridUniform3D>& partial_grids);// make partial grids by splitting this grid in Y direction
	void splitInZDirection(const int& num_threads, Array1D<GridUniform3D>& partial_grids);// make partial grids by splitting this grid in Z direction
	void splitInMaxDirection(const int& num_threads, Array1D<GridUniform3D>& partial_grids);// make partial grids by splitting direction automatically
	void splitInZDirectionPartially(const int& num_threads, Array1D<GridUniform3D>& partial_grids);

	int recommendMaxMultigridLevel(const int& lowest_level_res);

	T	getTrilinearWeight(const TV& deviation) const;
};// end of class GRID_UNIFORM_3D

template<class TT>
TT interpolateTrilinear(const GridUniform3D& grid, const Array3D<TT>& arr, const TV &position)
{
	T m = (T)1e-8;
	int i0, j0, k0;

	T x = CLAMP(position.x_, grid.x_min_+m, grid.x_max_-m);
	T y = CLAMP(position.y_, grid.y_min_+m, grid.y_max_-m);
	T z = CLAMP(position.z_, grid.z_min_+m, grid.z_max_-m);

	grid.getLeftBottomCell(x, y, z, i0, j0, k0);
	const TV p = grid.getCellCenter(i0,j0,k0);
	const int ix = arr.Get1DIndex(i0, j0, k0);

	const TT &v000 = arr.values_[ix];
	const TT &v100 = arr.values_[ix+1];
	const TT &v010 = arr.values_[ix+arr.i_res_];
	const TT &v001 = arr.values_[ix+arr.ij_res_];
	const TT &v101 = arr.values_[ix+1+arr.ij_res_];
	const TT &v011 = arr.values_[ix+arr.i_res_+arr.ij_res_];
	const TT &v110 = arr.values_[ix+1+arr.i_res_];
	const TT &v111 = arr.values_[ix+1+arr.i_res_+arr.ij_res_];

	x = (x - p.x_)/grid.dx_;
	y = (y - p.y_)/grid.dy_;
	z = (z - p.z_)/grid.dz_;

	const T w000 = ((T)1-x)*((T)1-y)*((T)1-z);
	const T w100 = x*((T)1-y)*((T)1-z);
	const T w010 = ((T)1-x)*y*((T)1-z);
	const T w001 = ((T)1-x)*((T)1-y)*z;
	const T w101 = x*((T)1-y)*z;
	const T w011 = ((T)1-x)*y*z;
	const T w110 = x*y*((T)1-z);
	const T w111 = x*y*z;

	assert((w000+w100+w010+w001+w101+w011+w110+w111) < (T)1.1);

	return v000*w000 + v100*w100 + v010*w010 + v001*w001 + v101*w101 + v011*w011 + v110*w110 + v111*w111;
}

inline std::ostream&
operator << (std::ostream& output,const GridUniform3D& grid)
{
// 	return output << "GRID_UNIFORM_3D ["
// 				  << "Resolution = "<<grid.i_res_<<" "<<grid.j_res_<<" "<<grid.k_res_
// 		          << " Index range = ("<<grid.i_start_<<" "<<grid.j_start_<<" "<<grid.k_start_<<") to ("<< grid.i_end_<<" "<<grid.j_end_<<" "<<grid.k_end_<<") "
// 				  << " Range = ("<< grid.x_min_<<" "<<grid.y_min_<<" "<<grid.z_min_<<") to ("<<grid.x_max_<<" "<<grid.y_max_<<" "<<grid.z_max_ <<") "
// 				  << " DX = " <<grid.dx_<<" "<<grid.dy_<<" "<<grid.dz_<<"]";
	return output <<" Index range = ("<<grid.i_start_<<" "<<grid.j_start_<<" "<<grid.k_start_<<") to ("<< grid.i_end_<<" "<<grid.j_end_<<" "<<grid.k_end_<<") ";
}