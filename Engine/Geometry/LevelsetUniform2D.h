// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "CONVENTIONAL_MACROS.h"
#include "GENERIC_DEFINITIONS.h"
#include "Parallelism/MultiThreading.h"
#include "DataStructure/GridUniform2D.h"
#include "DataStructure/Array2D.h"
#include "DataStructure/ArrayTools.h"

class LevelsetUniform2D
{
public:
	GridUniform2D grid_;
	GridUniform2D grid_ghost_;

	Array2D<T> phi_;
	Array2D<T> phi_ghost_;
	Array2D<T> phi_backup_;	// for convergence check
	Array2D<T> curvature_;
	Array2D<TV2> normal_;
	Array2D<bool> fixed_;	// "is fixed?" flag for narrow banded fast sweeping

	int sweep_direction_;

public:
	LevelsetUniform2D()
		: sweep_direction_(0)
	{}

	~LevelsetUniform2D()
	{}

	void initialize(const GridUniform2D& _grid, const int& ghost_width = 3)
	{
		grid_ = _grid;
		grid_ghost_ = grid_.getEnlarged(ghost_width);

		grid_ghost_.InitializeCellArray(phi_);
		grid_ghost_.InitializeCellArray(phi_ghost_);
		grid_ghost_.InitializeCellArray(normal_);
		grid_ghost_.InitializeCellArray(curvature_);
		grid_ghost_.InitializeCellArray(fixed_);

		phi_.assignAllValues(grid_.dx_*(T)3);

		sweep_direction_ = 0;
	}

	int getGhostWidth() const
	{
		return (grid_ghost_.i_res_ - grid_.i_res_) / 2;
	}

	void initialize(const LevelsetUniform2D& _levelset)
	{
		initialize(_levelset.grid_, _levelset.getGhostWidth());

		for (int p = 0; p < _levelset.grid_ghost_.getNumAllCells(); p++)
		{
			phi_.values_[p] = _levelset.phi_.values_[p];
		}
	}

	void initialize(MT* mt, const int& thread_id, const LevelsetUniform2D& _levelset)
	{
		BEGIN_ONE_THREAD_WORK
		{
			initialize(_levelset.grid_, getGhostWidth());
		}
		END_ONE_THREAD_WORK;

		BEGIN_1D_ITERATION(_levelset.grid_ghost_.getNumAllCells())
		{
			phi_.values_[p] = _levelset.phi_.values_[p];
		}
		END_1D_ITERATION;
	}

// 	void initializeFromScalarField(const Array2D<T>& height_map, const T th, const bool& high_to_inside)
// 	{
// 		if (high_to_inside == false)
// 		{
// 			for (int i = grid_.i_start_; i <= grid_.i_end_; i++)
// 				for (int j = grid_.j_start_; j <= grid_.j_end_; j++)
// 				{
// 					phi_(i, j) = (height_map(i, j) - th) * grid_.dx_;
// 				}
// 		}
// 		else
// 		{
// 			for (int i = grid_.i_start_; i <= grid_.i_end_; i++)
// 				for (int j = grid_.j_start_; j <= grid_.j_end_; j++)
// 				{
// 					phi_(i, j) = -(height_map(i, j) - th) * grid_.dx_;
// 				}
// 		}
// 
// 		ArrayTools::fillGhostCells(grid_, grid_ghost_, phi_);
// 	}

	void initializeFromScalarField(MultiThreading* mt, const int thread_id, const Array2D<T>& height_map, const T th, const bool& high_to_inside)
	{
		if (high_to_inside == false)
		{
			BEGIN_GRID_ITERATION_2D(grid_)
			{
				phi_(i, j) = (height_map(i, j) - th) * grid_.dx_;
			}
			END_GRID_ITERATION_2D;
		}
		else
		{
			BEGIN_GRID_ITERATION_2D(grid_)
			{
				phi_(i, j) = -(height_map(i, j) - th) * grid_.dx_;
			}
			END_GRID_ITERATION_2D;
		}

		fillGhostCells(mt, thread_id, phi_);

		reinitializeFSM(mt, thread_id, 8);
	}

	TV2 getGradient(const int i, const int j)
	{
		TV2 gradient;

		gradient.x_ = (phi_(i + 1, j) - phi_(i - 1, j))*grid_.one_over_2dx_;
		gradient.y_ = (phi_(i, j + 1) - phi_(i, j - 1))*grid_.one_over_2dy_;

		return gradient;
	}

	void reinitializeInterfacialCells(MT* mt, const int& thread_id);

	void fixInterfacialCells(MT* mt, const int& thread_id);	

	void reinitializeFSM(MT* mt, const int& thread_id, const int& itr = 2);
	void reinitializeFSMThreaded(const int& itr = 2);

//	void CuthillFastSweepMethod(const int& thread_id, const int& direction);

	void sweep(const int& i_start, const int& j_start, const int& i_end, const int& j_end);
//	void CuthillSweep(const int& thread_id, const int& i_start, const int& j_start);

// 	template<class TT> void FastSweepingExtrapolation(const int& thread_id, FIELD_UNIFORM_2D<TT>& value, FIELD_UNIFORM_2D<TT>& value_ghost, LevelsetUniform2D& object_levelset);
// 	template<class TT> void ExtrapolatingSweep(FIELD_UNIFORM_2D<TT>& value_ghost, LevelsetUniform2D& object_levelset, const int& i_start, const int& j_start, const int& i_end, const int& j_end);

	void computeNormals(MT* mt, const int& thread_id);
	void computeNormalThreaded(MT* mt);
	void computeCurvatures(MT* mt, const int& thread_id);
	T	 computeCurvature(const int& i, const int& j);

	T   getSignedDistance(const TV2& position) const
	{
		return grid_ghost_.getLinearInterpolationCell(phi_, position);
	}

	TV2 getNormal(const TV2& position) const
	{
		return grid_ghost_.getLinearInterpolationCell(normal_, position);
	}

	void getUnitNormal(const TV2& position, TV2& normal_output) const
	{
		normal_output = grid_ghost_.getLinearInterpolationCell(normal_, position).getNormalized();
	}

	TV2 getUnitNormal(const TV2& position) const
	{
		return grid_ghost_.getLinearInterpolationCell(normal_, position).getNormalized();
	}

	void smoothLaplacian(MT* mt, const int thread_id)
	{
		ArrayTools::copyArray(mt, thread_id, phi_, phi_ghost_);

		// smooth Laplacian phi, this helps smooth interfacial values
		BEGIN_GRID_ITERATION_2D(grid_)
		{
			const T average = (phi_ghost_(i - 1, j) + phi_ghost_(i + 1, j) + phi_ghost_(i, j - 1) + phi_ghost_(i, j + 1)) * (T)0.25;

			phi_(i, j) = average;
		}
		END_GRID_ITERATION_2D;
	}

	void advanceNormalFlowMain(MT* mt, const int thread_id, const T distance, const T CFL, const T beta)
	{
		T   dt = grid_.dx_ * CFL;
		const int num_steps = (int)(ABS(distance) / dt);
		T   dt_rest = ABS(distance) - (T)num_steps * dt;

		if(distance < (T)0)
		{
			dt = -dt;
			dt_rest = -dt_rest;
		}

		for (int n = 0; n < num_steps; n++)
		{
			advanceNormalFlow(mt, thread_id, dt, beta);

//			smoothLaplacian(mt, thread_id);

			reinitializeFSM(mt, thread_id, 8);
		}

		advanceNormalFlow(mt, thread_id, dt_rest, beta);

//		smoothLaplacian(mt, thread_id);

		reinitializeFSM(mt, thread_id, 8);
	}


	void advanceNormalFlow(MT* mt, const int thread_id, const T alpha, const T beta)
	{
		fillGhostCells(mt, thread_id, phi_);

		if (beta != (T)0)
		{
			computeNormals(mt, thread_id);
			computeCurvatures(mt, thread_id);
		}

		BEGIN_GRID_ITERATION_2D(grid_)
		{
			phi_(i, j) -= (alpha - curvature_(i,j) * beta);
		}
		END_GRID_ITERATION_2D;
	}

	T advanceShrinkWrapping(MT* mt, const int thread_id, const T CFL, const T beta, const Array2D<TV2>& velocity_field)
	{
		ArrayTools::copyArray(mt, thread_id, phi_, phi_backup_);

		fillGhostCells(mt, thread_id, phi_);

		computeNormals(mt, thread_id);

		computeCurvatures(mt, thread_id);

		const T dt = CFL * grid_.dx_;

		BEGIN_GRID_ITERATION_2D(grid_)
		{
			phi_(i, j) -= (dt * dotProduct(normal_(i, j), velocity_field(i, j)) - curvature_(i, j) *beta);
		}
		END_GRID_ITERATION_2D;

		smoothLaplacian(mt, thread_id);

		reinitializeFSM(mt, thread_id, 8);

		//TODO: converge condition
		T residual = (T)0;
		int count = 0;
		BEGIN_GRID_ITERATION_2D(grid_)
		{
			if (isInterfacial(i, j) == true)
			{
				residual += ABS(phi_backup_(i, j) - phi_(i, j));
				count++;
			}
		}
		END_GRID_ITERATION_2D;

		residual = mt->syncSum(thread_id, residual);
		count = mt->syncSum(thread_id, count);

		BEGIN_ONE_THREAD_WORK
		{
			phi_backup_.freeMemory();
		}
		END_ONE_THREAD_WORK;

		if (count == 0) return (T)0;
		else return residual * grid_.one_over_dx_ / (T)count;
	}

	template<class TT>
	void fillGhostCells(MT* mt, const int thread_id, Array2D<TT>& arr)
	{
		BEGIN_GRID_ITERATION_2D(grid_ghost_)
		{
			arr(i, j) = arr(grid_.ClampedIndex(i, j));		//TODO: don't need to copy non-ghost cells
		}
		END_GRID_ITERATION_2D;
	}

	void updateInterfacialPhi(const int i, const int j);

	bool isInterfacial(const int i, const int j) const
	{
		const int ix = phi_.get1DIndex(i, j);
		const int i_res = phi_.i_res_;

		const T* phi_values = phi_.values_;

		if (phi_.values_[ix] > (T)0)
		{
			//TODO: use pointer operations to access neighbor phis.
			if (phi_values[ix + 1] <= (T)0) return true;
			else if (phi_values[ix - 1] <= (T)0) return true;
			else if (phi_values[ix + i_res] <= (T)0) return true;
			else if (phi_values[ix - i_res] <= (T)0) return true;
		}
		else // if(phi(i,j) <= (T)0)
		{
			if (phi_values[ix + 1] > (T)0) return true;
			else if (phi_values[ix - 1] > (T)0) return true;
			else if (phi_values[ix + i_res] > (T)0) return true;
			else if (phi_values[ix - i_res] > (T)0) return true;
		}

		return false;
	}

	TV getPhiPos(const int i, const int j) const
	{
		const TV2 pos2 = grid_ghost_.GetCellCenter(i, j);
		return TV3(pos2.x_, pos2.y_, phi_(i, j));
	}
};
