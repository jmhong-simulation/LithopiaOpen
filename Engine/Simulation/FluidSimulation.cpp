// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "FluidSimulation.h"
#include "Utilities/RandomNumberGenerator.h"
#include "Utilities/ScriptBlock.h"
#include "Utilities/ScriptReader.h"
#include "Utilities/ScriptUtility.h"
#include "GL/GL_TOOLS.h"

void FluidSimulation::initialize()
{	
	ScriptReader script;
	script.readFile("script.sim");

	ScriptBlock &sb(script.getBlock("GeneralOptions"));

	RandomNumberGenerator random(123);

	gravity_ = sb.getValue("gravity", Vector3D<float>());

	// initialize grid

	ScriptUtility::initialize(sb, grid_);
	grid_.scaleRes(sb.getValue("res_scale", float()));

	grid_ghost_ = grid_.getEnlarged(3);		// num_ghost_width_ = 3
	u_grid_ghost_ = grid_ghost_.getUFaceGrid();
	v_grid_ghost_ = grid_ghost_.getVFaceGrid();
	w_grid_ghost_ = grid_ghost_.getWFaceGrid();

	u_grid_ghost_.initializeCenterArray(u_velocity_array_);
	v_grid_ghost_.initializeCenterArray(v_velocity_array_);
	w_grid_ghost_.initializeCenterArray(w_velocity_array_);

	u_grid_ghost_.initializeCenterArray(u_velocity_array_temp_);
	v_grid_ghost_.initializeCenterArray(v_velocity_array_temp_);
	w_grid_ghost_.initializeCenterArray(w_velocity_array_temp_);

	initializeSourceListFromScript(script.getBlock("SourceList"));

	sorted_particles_.initialize(grid_);
	sorted_particles_.sortThread();

	grid_ghost_.initializeCenterArray(water_levelset_);
	water_levelset_.AssignAllValues(grid_ghost_.dx_);

	initializeObjectListFromScript(*script.head_block_);

	// tag inside level set (marching cubes test only)
	for (int p = 0; p < sorted_particles_.num_all_pts_; p++)
	{
		water_levelset_(grid_ghost_.getCell(sorted_particles_.particle_data_arrays_.pts_pos_array_[p])) = -grid_ghost_.dx_;
	}

// 	for (int k = grid_ghost_.k_start_; k <= grid_ghost_.k_end_; k++)
// 	for (int j = grid_ghost_.j_start_; j <= grid_ghost_.j_end_; j++)
// 	for (int i = grid_ghost_.i_start_; i <= grid_ghost_.i_end_; i++)
// 	{
// 		if (k > 10) continue;
// 
// 		water_levelset_(i, j, k) = -grid_ghost_.dx_;
// 	}

	// Poisson solver parameters
	num_iteration_ = 0;
	max_iteration_ = 1000;
	tolerance_ = (T)1e-8;
	sqr_tolerance_ = tolerance_*tolerance_;
}

void FluidSimulation::initializeSourceListFromScript(ScriptBlock& super_block)
{
	num_pts_per_cell_per_second_ = super_block.getValue("num_pts_per_cell_per_second", int());

	source_velocity_ = super_block.getValue("velocity", TV());

	source_box_ = ScriptUtility::initializeBox3D(super_block);
}

void FluidSimulation::initializeObjectListFromScript(ScriptBlock& super_block)
{
	ScriptBlock &object_list_block = super_block.getBlock("ObjectList");

	fixed_solid_list_.initialize((int)object_list_block.children_.size());

	int ix = 0;
	for (std::list<ScriptBlock*>::iterator itr = object_list_block.children_.begin(); itr != object_list_block.children_.end(); ++itr)
	{
		ScriptBlock &object_block = *(*itr);

		fixed_solid_list_[ix] = ScriptUtility::initializeBox3D(object_block);

		ix++;
	}	
}

void FluidSimulation::initializeLinearSystem(MT* mt, const int& thread_id, const GridUniform3D& grid_ghost, Array3D<int>& boundary_array, CSRMatrix<D>& A_matrix, VectorND<D>& x_vector, VectorND<D>& b_vector)
{
	int full_ix(0);
	BEGIN_BOX_ITERATION_ARR_3D_SYNC(grid_ghost_)
	{
		if (boundary_array.values_[arr_ix] >= 0) full_ix++;
	}
	END_GRID_ITERATION_Z_SYNC;

	// correct full cell index for multi-threading and count number of non-zero elements of A matrix
	int nnz(0);		// number of non-zero values
	int i_res = grid_ghost.i_res_;
	int ij_res = grid_ghost.i_res_ * grid_ghost.j_res_;
	
	const TV2_INT ix_range = mt->getIncrementalRange(thread_id, full_ix);
	int start_full_ix = ix_range.t_min_, end_full_ix = ix_range.t_max_;
	BEGIN_BOX_ITERATION_ARR_3D_SYNC(grid_ghost_)
	{
		if (boundary_array.values_[arr_ix] >= 0)
		{
			boundary_array.values_[arr_ix] = start_full_ix++;

			nnz++;
			if (boundary_array.values_[arr_ix + 1] >= 0) nnz++;
			if (boundary_array.values_[arr_ix - 1] >= 0) nnz++;
			if (boundary_array.values_[arr_ix + i_res] >= 0) nnz++;
			if (boundary_array.values_[arr_ix - i_res] >= 0) nnz++;
			if (boundary_array.values_[arr_ix + ij_res] >= 0) nnz++;
			if (boundary_array.values_[arr_ix - ij_res] >= 0) nnz++;
		}
	}
	END_GRID_ITERATION_Z_SYNC;

	assert(start_full_ix - 1 == end_full_ix);

	full_ix = mt->syncSum(thread_id, full_ix);
	nnz = mt->syncSum(thread_id, nnz);

	ONE_THREAD_WORK(A_matrix.Initialize(mt->num_threads_, full_ix, nnz));
	ONE_THREAD_WORK(x_vector.Initialize(full_ix, true));
	ONE_THREAD_WORK(b_vector.Initialize(full_ix));
}

TV FluidSimulation::getCollisionVelocity(const TV& vel, const TV& normal, const T& res, const T& fric_dt)
{
//	std::cout << normal << std::endl;

	const T n_com = dotProduct(vel, normal);
	const TV vel_n = n_com*normal;
	const TV vel_t = vel - vel_n;

	assert(fric_dt <= (T)1);

	return vel_t * ((T)1 - fric_dt) + vel_n * (n_com < 0 ? -res : (T)1);
}

void FluidSimulation::setupBoundaryCondition(MT* mt, const int& thread_id)
{
	// smoke
// 	BEGIN_BOX_ITERATION_ARR_3D_SYNC(grid_ghost_)
// 	{
// 		if (grid_.isInside(i, j, k) == false) bc_array_.values_[arr_ix] = BC_DIRICHLET;
// 		else bc_array_.values_[arr_ix] = BC_FULL;
// 	}
// 	END_GRID_ITERATION_Z_SYNC;

	// water
	const T fric_dt = (T)0.01;
	const T res = (T)0;
//	const BOX_3D<T> obstacle((T)0.5, (T)0, (T)0, (T)0.8, (T)1, (T)1);
	int* num_pts_array_ = (int*)sorted_particles_.num_pts_array_;
	BEGIN_BOX_ITERATION_ARR_3D_SYNC(grid_ghost_)
	{
		if (grid_.isInside(i, j, k) == false)
		{
//			bc_array_.values_[arr_ix] = BC_NEUMANN;

			if (j < grid_.j_start_)
			{
				bc_array_.values_[arr_ix] = BC_NEUMANN;

				v_velocity_array_(i, j, k) = 0;
				v_velocity_array_(i, j + 1, k) = 0;

				v_velocity_array_(i, j, k) = 0;
				v_velocity_array_(i, j + 1, k) = 0;
				
				w_velocity_array_(i, j, k) = 0;
				w_velocity_array_(i, j, k + 1) = 0;
			}
			else
				bc_array_.values_[arr_ix] = BC_DIRICHLET;

			// wall velocity (TODO: move to elsewhere)	
// 
// 			u_velocity_array_(i, j, k) = 0;
// 			u_velocity_array_(i + 1, j, k) = 0;
// 
// 			v_velocity_array_(i, j, k) = 0;
// 			v_velocity_array_(i, j + 1, k) = 0;
// 
// 			w_velocity_array_(i, j, k) = 0;
// 			w_velocity_array_(i, j, k + 1) = 0;
		}
		else
		{
			const int num_pt_this_cell = num_pts_array_[grid_.get1DIndex(i, j, k)];	// num_pts_array uses grid_, bc_array uses grid_ghost.

			if (num_pt_this_cell > 5) bc_array_.values_[arr_ix] = BC_FULL;
			else bc_array_.values_[arr_ix] = BC_DIRICHLET;
		}
	}
	END_GRID_ITERATION_Z_SYNC;

	for (int oix = 0; oix < fixed_solid_list_.num_elements_; oix++)
	{
		const BOX_3D<T>& object = fixed_solid_list_[oix];
		const BOX_3D<int> ix_box = BOX_3D<int>(grid_ghost_.getLeftBottomCell(object.GetMin()), grid_ghost_.getRightTopCell(object.GetMax()));
		const GridUniform3D partial_grid = grid_ghost_.getPartialGrid(ix_box);

		{
			const GridUniform3D &_grid = partial_grid;

			const T half_dx = (T)0.5*_grid.dx_;
			const TV2_INT k_range = mt->getParallelRange(thread_id, _grid.k_start_, _grid.k_end_);
			const int i_start(_grid.i_start_), i_end(_grid.i_end_), j_start(_grid.j_start_), j_end(_grid.j_end_), k_start(k_range.t_min_), k_end(k_range.t_max_);
			int i, j, k, arr_ix, u_ix, v_ix, w_ix;
			TV cell_center;
			for (k = k_start; k <= k_end; k ++)
				for (j = j_start; j <= j_end; ++j)
					for (i = i_start, 
						arr_ix = _grid.get1DIndex(i, j, k), 
						u_ix = u_grid_ghost_.get1DIndex(i, j, k),
						v_ix = v_grid_ghost_.get1DIndex(i, j, k),
						w_ix = w_grid_ghost_.get1DIndex(i, j, k),
						cell_center = grid_ghost_.getCellCenter(i,j,k);
						i <= i_end; i++, arr_ix++, u_ix++, v_ix++, w_ix++, cell_center.x_ += grid_ghost_.dx_)
					{
						if (object.getSignedDistance(grid_.getCellCenter(i, j, k)) <= (T)0)
						{
							bc_array_.values_[arr_ix] = BC_NEUMANN;
							// object velocity (TODO: move to elsewhere)

							//TODO friction, move to elsewhere for many object, use tangent friction
							u_velocity_array_.values_[u_ix] = getCollisionVelocity(TV(u_velocity_array_.values_[u_ix], 0, 0), object.getNormal(cell_center - TV(half_dx, (T)0, (T)0), grid_.dx_), res, fric_dt).x_;
							u_velocity_array_.values_[u_ix + 1] = getCollisionVelocity(TV(u_velocity_array_.values_[u_ix + 1], 0, 0), object.getNormal(cell_center + TV(half_dx, (T)0, (T)0), grid_.dx_), res, fric_dt).x_;

							v_velocity_array_.values_[v_ix] = getCollisionVelocity(TV(0, v_velocity_array_.values_[v_ix], 0), object.getNormal(cell_center - TV((T)0, half_dx, (T)0), grid_.dy_), res, fric_dt).y_;
							v_velocity_array_.values_[v_ix + v_grid_ghost_.i_res_] = getCollisionVelocity(TV(0, v_velocity_array_.values_[v_ix + v_grid_ghost_.i_res_], 0), object.getNormal(cell_center + TV((T)0, half_dx, (T)0), grid_.dy_), res, fric_dt).y_;

							w_velocity_array_.values_[w_ix] = getCollisionVelocity(TV(0, 0, w_velocity_array_(i, j, k)), object.getNormal(cell_center - TV((T)0, (T)0, half_dx), grid_.dz_), res, fric_dt).z_;
							w_velocity_array_.values_[w_ix + w_grid_ghost_.ij_res_] = getCollisionVelocity(TV(0, 0, w_velocity_array_.values_[w_ix + w_grid_ghost_.ij_res_]), object.getNormal(cell_center + TV((T)0, (T)0, half_dx), grid_.dz_), res, fric_dt).z_;
						}
					}
			mt->sync();
		}
	}
}

void FluidSimulation::project(MT* mt, const int thread_id, const T dt)
{
	BEGIN_ONE_THREAD_WORK
	{
		grid_ghost_.initializeCenterArray(bc_array_);
		grid_ghost_.initializeCenterArray(div_array_);
		grid_ghost_.initializeCenterArray(p_array_);
	}
	END_ONE_THREAD_WORK;

	setupBoundaryCondition(mt, thread_id);

	BEGIN_BOX_ITERATION_ARR_3D_SYNC(grid_ghost_)
	{
		T div = u_velocity_array_(i + 1, j, k) - u_velocity_array_(i, j, k);
		div += v_velocity_array_(i, j + 1, k) - v_velocity_array_(i, j, k);
		div += w_velocity_array_(i, j, k + 1) - w_velocity_array_(i, j, k);

//		div_array_(i, j, k) = - div * grid_ghost_.one_over_dx_;	// dx = dy = dz
		div_array_(i, j, k) = div;
	}
	END_GRID_ITERATION_Z_SYNC;

	BEGIN_PARALLEL_FOR(p, 0, grid_ghost_.getNumAllCells() - 1)
	{
		p_array_.values_[p] = (T)0;
	}
	END_PARALLEL_FOR;

	initializeLinearSystem(mt, thread_id, grid_ghost_, bc_array_, A_matrix_, x_vector_, b_vector_);
	setupLinearSystem(mt, thread_id, grid_ghost_, bc_array_, div_array_, p_array_);
	solveLinearSystemCPUCG(mt, thread_id, A_matrix_, x_vector_, b_vector_);

	BEGIN_BOX_ITERATION_ARR_3D_SYNC(grid_ghost_)
	{
		if (bc_array_.values_[arr_ix] >= 0) p_array_.values_[arr_ix] = (T)x_vector_[bc_array_.values_[arr_ix]];

		//		if (pressure_array.values_[arr_ix] * dt_dt > safe_dist) pressure_array_[ix] = (T)0;
	}
	END_GRID_ITERATION_Z_SYNC;

	updateVelocityByPressureGradient(mt, thread_id);
}

void FluidSimulation::solveLinearSystemCPUCG(MT* mt, const int& thread_id, CSRMatrix<D>& A_matrix, VectorND<D>& x_vector, VectorND<D>& b_vector)
{
	BEGIN_ONE_THREAD_WORK
	{
		mt->splitRange(0, A_matrix.N_ - 1);
	}
	END_ONE_THREAD_WORK;

	const int N(x_vector.num_dimension_), start_ix(mt->start_ix_1D_[thread_id]), end_ix(mt->end_ix_1D_[thread_id]);

	ONE_THREAD_WORK(res_.Initialize(N));
	ONE_THREAD_WORK(p_.Initialize(N));
	ONE_THREAD_WORK(Ap_.Initialize(N));
	ONE_THREAD_WORK(num_iteration_ = 0;);

	D *rval(res_.values_), *pval(p_.values_), *Apval(Ap_.values_), *xval(x_vector.values_);

	D alpha, res_old, res_new;

	A_matrix.ComputeResidual(mt, thread_id, x_vector, b_vector, res_);

	for (int i = start_ix; i <= end_ix; i++)
	{
		p_.values_[i] = res_.values_[i];
	}
	mt->sync();

	res_old = dotProduct(mt, thread_id, res_, p_);

	if (res_old < sqr_tolerance_)
	{
		num_iteration_ = 0;
		return;
	}

	num_iteration_ = 0;
	while (num_iteration_ < max_iteration_)
	{
		A_matrix.Multiply(mt, thread_id, p_, Ap_);

		alpha = res_old / dotProduct(mt, thread_id, p_, Ap_);

		for (int i = start_ix; i <= end_ix; i++)
		{
			xval[i] += alpha * pval[i];
			rval[i] -= alpha * Apval[i];
		}
		mt->sync();

		res_new = dotProduct(mt, thread_id, res_, res_);

		if (res_new < sqr_tolerance_) break;			// L2 norm

		for (int i = start_ix; i <= end_ix; i++)
		{
			const D k = res_new / res_old;

			pval[i] = res_.values_[i] + k*p_.values_[i];
		}
		mt->sync();

		res_old = res_new;

		BEGIN_ONE_THREAD_WORK
		{
			num_iteration_++;
		}
		END_ONE_THREAD_WORK;
	}

	mt->sync();

	BEGIN_ONE_THREAD_WORK
	{
		residual_ = sqrt(res_new);
		//std::cout<<"Iteration = "<<num_iteration_<<" residual = "<<residual_<<std::endl;
	}
	END_ONE_THREAD_WORK;
}

void FluidSimulation::updateVelocityByPressureGradient(MT* mt, const int& thread_id)
{
	T *pressure_array_ = p_array_.values_;		//TODO: use double
	int *boundary_array_ = bc_array_.values_;
//	TV *velocity_array_ = vedata_.velocity_array_;
	const int i_start = grid_.i_start_, i_end = grid_.i_end_;
	const int j_start = grid_.j_start_, j_end = grid_.j_end_;
	const int cell_num_i = grid_.i_res_;

	const int i_res = grid_ghost_.i_res_;
	const int ij_res = grid_ghost_.i_res_ * grid_ghost_.j_res_;

	// update velocity
//	T p_ijk, p_ip, p_im, p_jp, p_jm, p_kp, p_km;
	BEGIN_BOX_ITERATION_ARR_3D_SYNC(grid_ghost_)
	{
		if (boundary_array_[arr_ix] < 0) continue;

//		const T p_over_dx = pressure_array_[arr_ix] * grid_ghost_.one_over_dx_;	// assumes dx = dy = dz
		const T p_over_dx = pressure_array_[arr_ix];

		u_velocity_array_(i, j, k) -= p_over_dx;
		u_velocity_array_(i + 1, j, k) += p_over_dx;

		v_velocity_array_(i, j, k) -= p_over_dx;
		v_velocity_array_(i, j + 1, k) += p_over_dx;

		w_velocity_array_(i, j, k) -= p_over_dx;
		w_velocity_array_(i, j, k + 1) += p_over_dx;
/*
		p_ijk = pressure_array_[ix];

		if (boundary_array_[ix + 1] >= 0) p_ip = pressure_array_[ix + 1];
		else if (boundary_array_[ix + 1] <= BC_NEUMANN)	 p_ip = p_ijk;
		else if (boundary_array_[ix + 1] == BC_DIRICHLET) p_ip = -p_ijk;
		else { assert(false); }

		if (boundary_array_[ix - 1] >= 0) p_im = pressure_array_[ix - 1];
		else if (boundary_array_[ix - 1] <= BC_NEUMANN)   p_im = p_ijk;
		else if (boundary_array_[ix - 1] == BC_DIRICHLET) p_im = -p_ijk;
		else { assert(false); }


		if (boundary_array_[ix + i_res] >= 0) p_jp = pressure_array_[ix + i_res];
		else if (boundary_array_[ix + i_res] <= BC_NEUMANN)	 p_jp = p_ijk;
		else if (boundary_array_[ix + i_res] == BC_DIRICHLET) p_jp = -p_ijk;
		else { assert(false); }

		if (boundary_array_[ix - i_res] >= 0) p_jm = pressure_array_[ix - i_res];
		else if (boundary_array_[ix - i_res] <= BC_NEUMANN)	 p_jm = p_ijk;
		else if (boundary_array_[ix - i_res] == BC_DIRICHLET) p_jm = -p_ijk;
		else { assert(false); }


		if (boundary_array_[ix + ij_res] >= 0) p_kp = pressure_array_[ix + ij_res];
		else if (boundary_array_[ix + ij_res] <= BC_NEUMANN)	  p_kp = p_ijk;
		else if (boundary_array_[ix + ij_res] == BC_DIRICHLET) p_kp = -p_ijk;
		else { assert(false); }

		if (boundary_array_[ix - ij_res] >= 0) p_km = pressure_array_[ix - ij_res];
		else if (boundary_array_[ix - ij_res] <= BC_NEUMANN)	  p_km = p_ijk;
		else if (boundary_array_[ix - ij_res] == BC_DIRICHLET) p_km = -p_ijk;
		else { assert(false); }

		TV3 acc = TV3((T)0.5*(p_ip - p_im), (T)0.5*(p_jp - p_jm), (T)0.5*(p_kp - p_km));

		velocity_array_[ix] -= acc;
		//data_.mpm_force_array_[ix] = TV2((T)-0.5*(p_ip-p_im), (T)-0.5*(p_jp-p_jm));			
*/
	}
	END_GRID_ITERATION_Z_SYNC;
}

void FluidSimulation::setupLinearSystem(MT* mt, const int& thread_id, const GridUniform3D& grid_ghost, Array3D<int>& boundary_array, Array3D<T>& divergence_array, Array3D<T>& pressure_array)
{
	BEGIN_ONE_THREAD_WORK
	{
		A_matrix_.PrepareForMultithreading(mt);
	}
	END_ONE_THREAD_WORK;

	int i_res = grid_ghost.i_res_;
	int ij_res = grid_ghost.i_res_ * grid_ghost.j_res_;

	BEGIN_BOX_ITERATION_ARR_3D_SYNC(grid_ghost_)
	{
		if (boundary_array.values_[arr_ix] < 0) continue;

		const int bc_ijk = boundary_array.values_[arr_ix];

		int coeff_ijk = 0; // for optimization, inv_dxdx is multiplied at the end.

		// if neighbor is full cell
		if (boundary_array.values_[arr_ix + 1] >= 0)
		{
			coeff_ijk += 1;
			A_matrix_.AssignValue(thread_id, bc_ijk, boundary_array.values_[arr_ix + 1], -1);
		}

		if (boundary_array.values_[arr_ix - 1] >= 0)
		{
			coeff_ijk += 1;
			A_matrix_.AssignValue(thread_id, bc_ijk, boundary_array.values_[arr_ix - 1], -1);
		}

		if (boundary_array.values_[arr_ix + i_res] >= 0)
		{
			coeff_ijk += 1;
			A_matrix_.AssignValue(thread_id, bc_ijk, boundary_array.values_[arr_ix + i_res], -1);
		}

		if (boundary_array.values_[arr_ix - i_res] >= 0)
		{
			coeff_ijk += 1;
			A_matrix_.AssignValue(thread_id, bc_ijk, boundary_array.values_[arr_ix - i_res], -1);
		}

		if (boundary_array.values_[arr_ix + ij_res] >= 0)
		{
			coeff_ijk += 1;
			A_matrix_.AssignValue(thread_id, bc_ijk, boundary_array.values_[arr_ix + ij_res], -1);
		}

		if (boundary_array.values_[arr_ix - ij_res] >= 0)
		{
			coeff_ijk += 1;
			A_matrix_.AssignValue(thread_id, bc_ijk, boundary_array.values_[arr_ix - ij_res], -1);
		}

		if (boundary_array.values_[arr_ix + 1] == BC_DIRICHLET) coeff_ijk += 1;
		if (boundary_array.values_[arr_ix - 1] == BC_DIRICHLET) coeff_ijk += 1;
		if (boundary_array.values_[arr_ix + i_res] == BC_DIRICHLET) coeff_ijk += 1;
		if (boundary_array.values_[arr_ix - i_res] == BC_DIRICHLET) coeff_ijk += 1;
		if (boundary_array.values_[arr_ix + ij_res] == BC_DIRICHLET) coeff_ijk += 1;
		if (boundary_array.values_[arr_ix - ij_res] == BC_DIRICHLET) coeff_ijk += 1;

		if (coeff_ijk == 0){ assert(false); }// no connected FULL and DIR cells, all Neumann condition (null space)

		A_matrix_.AssignValue(thread_id, bc_ijk, bc_ijk, (T)coeff_ijk);

		b_vector_[bc_ijk] = -divergence_array.values_[arr_ix];// Note: We solve -Ax = -b
		//x_vector_[bc_ijk] = pressure_array_[arr_ix];
		x_vector_[bc_ijk] = 0;
	}
	END_GRID_ITERATION_Z_SYNC;

/*
	if (use_cuda_)
	{
		BEGIN_HEAD_THREAD_WORK
		{
			GPUCGSolve(A_matrix_.row_ptr_, A_matrix_.column_index_, A_matrix_.values_, x_vector_.values_, b_vector_.values_, x_vector_.num_dimension_, x_vector_.num_dimension_, A_matrix_.nz_, 0);
		}
		END_HEAD_THREAD_WORK;
	}
	else
	{
		CPUCGSolve(thread_id, A_matrix_, x_vector_, b_vector_);
	}
*/
}

void FluidSimulation::updateFaceVelocityHelper(const TV& particle_pos, const T& vel_compoment, const GridUniform3D& face_grid, const int &i, const int& j, const int& k, const Array3D<T>& face_vel_array, const Array3D<T>& face_vel_array_temp)
{
	const TV deviation = particle_pos - face_grid.getCellCenter(i, j, k);

	if (deviation.x_ > grid_ghost_.dx_ || deviation.y_ > grid_ghost_.dy_ || deviation.z_ > grid_ghost_.dz_) return;

	const T weight = face_grid.getTrilinearWeight(deviation);

	face_vel_array(i, j, k) += vel_compoment * weight;
	face_vel_array_temp(i, j, k) += weight;
}

void FluidSimulation::addNewParticlesFromSourceList(MT* mt, const int& thread_id)
{
//	for (int oix = 0; oix < fixed_solid_list_.num_elements_; oix++)
	{
		RandomNumberGenerator &random(mt->getRandom(thread_id));

		const TV source_velocity = source_velocity_;

//		const BOX_3D<T>& object = fixed_solid_list_[oix];
		const BOX_3D<T>& object = source_box_;
		const BOX_3D<int> ix_box = BOX_3D<int>(grid_ghost_.getLeftBottomCell(object.GetMin()), grid_ghost_.getRightTopCell(object.GetMax()));
		const GridUniform3D partial_grid = grid_ghost_.getPartialGrid(ix_box);

		{
			const GridUniform3D &_grid = partial_grid;

			const T half_dx = (T)0.5*_grid.dx_;
			const TV2_INT k_range = mt->getParallelRange(thread_id, _grid.k_start_, _grid.k_end_);
			const int i_start(_grid.i_start_), i_end(_grid.i_end_), j_start(_grid.j_start_), j_end(_grid.j_end_), k_start(k_range.t_min_), k_end(k_range.t_max_);
			int i, j, k, arr_ix, u_ix, v_ix, w_ix;
			TV cell_center;
			for (k = k_start; k <= k_end; k++)
				for (j = j_start; j <= j_end; ++j)
					for (i = i_start,
						arr_ix = _grid.get1DIndex(i, j, k),
						u_ix = u_grid_ghost_.get1DIndex(i, j, k),
						v_ix = v_grid_ghost_.get1DIndex(i, j, k),
						w_ix = w_grid_ghost_.get1DIndex(i, j, k),
						cell_center = grid_ghost_.getCellCenter(i, j, k);
						i <= i_end; i++, arr_ix++, u_ix++, v_ix++, w_ix++, cell_center.x_ += grid_ghost_.dx_)
					{
						const int num_pts_to_seed = random.getProbabilityInteger((float)num_pts_per_cell_per_second_ * dt_);

						for (int p = 0; p < num_pts_to_seed; p++)
						{
							const TV pos = grid_ghost_.getCellMin(i, j, k) + random.getVector() * grid_ghost_.dx_;	//TODO: dx, dy, dz

							if (object.isInside(pos) == true)
							{
								sorted_particles_.particle_data_arrays_.addNewParticle(pos, source_velocity);
							}
						}
					}
			mt->sync();
		}
	}
}

void FluidSimulation::getGridVelocityFromParticles(MT* mt, const int& thread_id)
{
	// particle velocity to grid
	mt->setArray(thread_id, u_grid_ghost_.getNumAllCells(), (T)0, u_velocity_array_.values_);
	mt->setArray(thread_id, v_grid_ghost_.getNumAllCells(), (T)0, v_velocity_array_.values_);
	mt->setArray(thread_id, w_grid_ghost_.getNumAllCells(), (T)0, w_velocity_array_.values_);

	mt->setArray(thread_id, u_grid_ghost_.getNumAllCells(), (T)0, u_velocity_array_temp_.values_);
	mt->setArray(thread_id, v_grid_ghost_.getNumAllCells(), (T)0, v_velocity_array_temp_.values_);
	mt->setArray(thread_id, w_grid_ghost_.getNumAllCells(), (T)0, w_velocity_array_temp_.values_);

	TV* pts_pos_array_ = sorted_particles_.particle_data_arrays_.pts_pos_array_;
	TV* pts_vel_array_ = sorted_particles_.particle_data_arrays_.pts_vel_array_;
	int* pts_start_ix_array_ = (int*)sorted_particles_.pts_start_ix_array_;
	int* num_pts_array_ = (int*)sorted_particles_.num_pts_array_;

// 	const TV2_INT k_range = mt->getParallelRange(thread_id, grid_.k_start_, grid_.k_end_);
// 	const BOX_3D<int> partial_box = grid_.getIXBox().getZResized(k_range.t_min_, k_range.t_max_);
// 	const GridUniform3D partial_grid = grid_.getPartialGrid(partial_box);

	// try red green

	BEGIN_BOX_ITERATION_ARR_3D_SYNC(grid_)
	{
		int pt_count = pts_start_ix_array_[arr_ix];
		const int num_pt_this_cell = num_pts_array_[arr_ix];

		if (num_pt_this_cell == 0) continue;

		for (int pt = 0; pt < num_pt_this_cell; pt++)
		{
			TV &position = pts_pos_array_[pt_count];
			TV &velocity = pts_vel_array_[pt_count];

			// x velocity of this particle influences u velocity of 10 faces
			updateFaceVelocityHelper(position, velocity.x_, u_grid_ghost_, i, j, k, u_velocity_array_, u_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.x_, u_grid_ghost_, i, j - 1, k, u_velocity_array_, u_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.x_, u_grid_ghost_, i, j + 1, k, u_velocity_array_, u_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.x_, u_grid_ghost_, i, j, k - 1, u_velocity_array_, u_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.x_, u_grid_ghost_, i, j, k + 1, u_velocity_array_, u_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.x_, u_grid_ghost_, i + 1, j, k, u_velocity_array_, u_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.x_, u_grid_ghost_, i + 1, j - 1, k, u_velocity_array_, u_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.x_, u_grid_ghost_, i + 1, j + 1, k, u_velocity_array_, u_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.x_, u_grid_ghost_, i + 1, j, k - 1, u_velocity_array_, u_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.x_, u_grid_ghost_, i + 1, j, k + 1, u_velocity_array_, u_velocity_array_temp_);

			updateFaceVelocityHelper(position, velocity.y_, v_grid_ghost_, i, j, k, v_velocity_array_, v_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.y_, v_grid_ghost_, i - 1, j, k, v_velocity_array_, v_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.y_, v_grid_ghost_, i + 1, j, k, v_velocity_array_, v_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.y_, v_grid_ghost_, i, j, k - 1, v_velocity_array_, v_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.y_, v_grid_ghost_, i, j, k + 1, v_velocity_array_, v_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.y_, v_grid_ghost_, i, j + 1, k, v_velocity_array_, v_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.y_, v_grid_ghost_, i - 1, j + 1, k, v_velocity_array_, v_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.y_, v_grid_ghost_, i + 1, j + 1, k, v_velocity_array_, v_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.y_, v_grid_ghost_, i, j + 1, k - 1, v_velocity_array_, v_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.y_, v_grid_ghost_, i, j + 1, k + 1, v_velocity_array_, v_velocity_array_temp_);

			updateFaceVelocityHelper(position, velocity.z_, w_grid_ghost_, i, j, k, w_velocity_array_, w_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.z_, w_grid_ghost_, i - 1, j, k, w_velocity_array_, w_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.z_, w_grid_ghost_, i + 1, j, k, w_velocity_array_, w_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.z_, w_grid_ghost_, i, j - 1, k, w_velocity_array_, w_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.z_, w_grid_ghost_, i, j - 1, k, w_velocity_array_, w_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.z_, w_grid_ghost_, i, j, k + 1, w_velocity_array_, w_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.z_, w_grid_ghost_, i - 1, j, k + 1, w_velocity_array_, w_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.z_, w_grid_ghost_, i + 1, j, k + 1, w_velocity_array_, w_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.z_, w_grid_ghost_, i, j - 1, k + 1, w_velocity_array_, w_velocity_array_temp_);
			updateFaceVelocityHelper(position, velocity.z_, w_grid_ghost_, i, j - 1, k + 1, w_velocity_array_, w_velocity_array_temp_);

			// TODO multithreading error

			pt_count++;
		}
	}
	END_GRID_ITERATION_Z_SYNC;

	BEGIN_1D_ITERATION(u_grid_ghost_.getNumAllCells())
	{
		const T weight = u_velocity_array_temp_.values_[p];

		if (weight >= 1e-8) u_velocity_array_.values_[p] /= weight;
		else u_velocity_array_.values_[p] = (T)0;
	}
	END_1D_ITERATION;

	BEGIN_1D_ITERATION(v_grid_ghost_.getNumAllCells())
	{
		const T weight = v_velocity_array_temp_.values_[p];

		if (weight >= 1e-8) v_velocity_array_.values_[p] /= weight;
		else v_velocity_array_.values_[p] = (T)0;
	}
	END_1D_ITERATION;

	BEGIN_1D_ITERATION(w_grid_ghost_.getNumAllCells())
	{
		const T weight = w_velocity_array_temp_.values_[p];

		if (weight >= 1e-8) w_velocity_array_.values_[p] /= weight;
		else w_velocity_array_.values_[p] = (T)0;
	}
	END_1D_ITERATION;
}

void FluidSimulation::advectParticles(MT* mt, const int& thread_id)
{
	const T flip_coef = (T)0.99;

	// using sorted particles
	TV* pts_pos_array_ = sorted_particles_.particle_data_arrays_.pts_pos_array_;
	TV* pts_vel_array_ = sorted_particles_.particle_data_arrays_.pts_vel_array_;
	int* pts_start_ix_array_ = (int*)sorted_particles_.pts_start_ix_array_;
	int* num_pts_array_ = (int*)sorted_particles_.num_pts_array_;

	BEGIN_BOX_ITERATION_ARR_3D_SYNC(grid_)
	{
		int pt_count = pts_start_ix_array_[arr_ix];
		const int num_pt_this_cell = num_pts_array_[arr_ix];

		if (num_pt_this_cell == 0) continue;

		for (int pt = 0; pt < num_pt_this_cell; pt++)
		{
			TV &position = pts_pos_array_[pt_count];
			TV &velocity = pts_vel_array_[pt_count];

			const TV interpolated_velocity_difference(interpolateTrilinear(u_grid_ghost_, u_velocity_array_temp_, position), interpolateTrilinear(v_grid_ghost_, v_velocity_array_temp_, position), interpolateTrilinear(w_grid_ghost_, w_velocity_array_temp_, position));
			const TV interpolated_velocity(interpolateTrilinear(u_grid_ghost_, u_velocity_array_, position), interpolateTrilinear(v_grid_ghost_, v_velocity_array_, position), interpolateTrilinear(w_grid_ghost_, w_velocity_array_, position));

			velocity += interpolated_velocity_difference;

			velocity = velocity * flip_coef + interpolated_velocity*((T)1 - flip_coef);
			position += interpolated_velocity * dt_;

			velocity += gravity_*dt_;

			// wall collision (move to elsewhere)
			const T push_back_thickness = (T)1e-6;
			const T restitution = (T)0.5;

// 			if (position.x_ < grid_.x_min_)
// 			{
// 				if(velocity.x_ < (T)0) velocity.x_ *= -restitution;
// 				position.x_ = grid_.x_min_ + push_back_thickness;
// 			}
// 			else if (position.x_ > grid_.x_max_)
// 			{
// 				if (velocity.x_ > (T)0) velocity.x_ *= -restitution;
// 				position.x_ = grid_.x_max_ - push_back_thickness;
// 			}
			
			if (position.y_ < grid_.y_min_)
			{
				if (velocity.y_ < (T)0) velocity.y_ *= -restitution;
				position.y_ = grid_.y_min_ + push_back_thickness;
			}
// 			else if (position.y_ > grid_.y_max_)
// 			{
// 				if (velocity.y_ > (T)0) velocity.y_ *= -restitution;
// 				position.y_ = grid_.y_max_ - push_back_thickness;
// 			}			
// 
// 			if (position.z_ < grid_.z_min_)
// 			{
// 				if (velocity.z_ < (T)0) velocity.z_ *= -restitution;
// 				position.z_ = grid_.z_min_ + push_back_thickness;
// 			}
// 			else if (position.z_ > grid_.z_max_)
// 			{
// 				if (velocity.z_ > (T)0) velocity.z_ *= -restitution;
// 				position.z_ = grid_.z_max_ - push_back_thickness;
// 			}

			pt_count++;
		}
	}
	END_GRID_ITERATION_Z_SYNC;
}

void FluidSimulation::advectEulerianVelocity(MT* mt, const int& thread_id)
{
	// Eulerian advection
	BEGIN_BOX_ITERATION_ARR_3D_SYNC(grid_ghost_)
	{
		u_velocity_array_temp_(i, j, k) = u_velocity_array_(i, j, k);
		v_velocity_array_temp_(i, j, k) = v_velocity_array_(i, j, k);
		w_velocity_array_temp_(i, j, k) = w_velocity_array_(i, j, k);
	}
	END_GRID_ITERATION_Z_SYNC;

	BEGIN_BOX_ITERATION_ARR_3D_SYNC(u_grid_ghost_)
	{
		if (grid_.isInside(i, j, k) == false) continue;	//TODO: remove

		const TV pos = u_grid_ghost_.getCellCenter(i, j, k);
		const T u_vel = interpolateTrilinear(u_grid_ghost_, u_velocity_array_temp_, pos);
		const T v_vel = interpolateTrilinear(v_grid_ghost_, v_velocity_array_temp_, pos);
		const T w_vel = interpolateTrilinear(w_grid_ghost_, w_velocity_array_temp_, pos);
		const TV vel(u_vel, v_vel, w_vel);
		const TV pos_back = pos - vel*dt_;

		u_velocity_array_(i, j, k) = interpolateTrilinear(u_grid_ghost_, u_velocity_array_temp_, pos_back);
	}
	END_GRID_ITERATION_Z_SYNC;

	BEGIN_BOX_ITERATION_ARR_3D_SYNC(v_grid_ghost_)
	{
		if (grid_.isInside(i, j, k) == false) continue;	//TODO: remove

		const TV pos = v_grid_ghost_.getCellCenter(i, j, k);
		const T u_vel = interpolateTrilinear(u_grid_ghost_, u_velocity_array_temp_, pos);
		const T v_vel = interpolateTrilinear(v_grid_ghost_, v_velocity_array_temp_, pos);
		const T w_vel = interpolateTrilinear(w_grid_ghost_, w_velocity_array_temp_, pos);
		const TV vel(u_vel, v_vel, w_vel);
		const TV pos_back = pos - vel*dt_;

		v_velocity_array_(i, j, k) = interpolateTrilinear(v_grid_ghost_, v_velocity_array_temp_, pos_back);
	}
	END_GRID_ITERATION_Z_SYNC;

	BEGIN_BOX_ITERATION_ARR_3D_SYNC(w_grid_ghost_)
	{
		if (grid_.isInside(i, j, k) == false) continue;	//TODO: remove

		const TV pos = w_grid_ghost_.getCellCenter(i, j, k);
		const T u_vel = interpolateTrilinear(u_grid_ghost_, u_velocity_array_temp_, pos);
		const T v_vel = interpolateTrilinear(v_grid_ghost_, v_velocity_array_temp_, pos);
		const T w_vel = interpolateTrilinear(w_grid_ghost_, w_velocity_array_temp_, pos);
		const TV vel(u_vel, v_vel, w_vel);
		const TV pos_back = pos - vel*dt_;

		w_velocity_array_(i, j, k) = interpolateTrilinear(w_grid_ghost_, w_velocity_array_temp_, pos_back);
	}
	END_GRID_ITERATION_Z_SYNC;
}

void FluidSimulation::saveVelocity(MT* mt, const int& thread_id)
{
	// backup before velocity
	BEGIN_1D_ITERATION(u_velocity_array_temp_.ijk_res_)
	{
		u_velocity_array_temp_.values_[p] = u_velocity_array_.values_[p];
	}
	END_1D_ITERATION;

	BEGIN_1D_ITERATION(v_velocity_array_temp_.ijk_res_)
	{
		v_velocity_array_temp_.values_[p] = v_velocity_array_.values_[p];
	}
	END_1D_ITERATION;

	BEGIN_1D_ITERATION(w_velocity_array_temp_.ijk_res_)
	{
		w_velocity_array_temp_.values_[p] = w_velocity_array_.values_[p];
	}
	END_1D_ITERATION;
}

void FluidSimulation::saveVelocityDifference(MT* mt, const int& thread_id)
{
	BEGIN_1D_ITERATION(u_velocity_array_temp_.ijk_res_)
	{
		u_velocity_array_temp_.values_[p] = u_velocity_array_.values_[p] - u_velocity_array_temp_.values_[p];
	}
	END_1D_ITERATION;

	BEGIN_1D_ITERATION(v_velocity_array_temp_.ijk_res_)
	{
		v_velocity_array_temp_.values_[p] = v_velocity_array_.values_[p] - v_velocity_array_temp_.values_[p];
	}
	END_1D_ITERATION;

	BEGIN_1D_ITERATION(w_velocity_array_temp_.ijk_res_)
	{
		w_velocity_array_temp_.values_[p] = w_velocity_array_.values_[p] - w_velocity_array_temp_.values_[p];
	}
	END_1D_ITERATION;
}

void FluidSimulation::updateOneStep(MT* mt, const int thread_id)
{
	// advance one frame
	const T dt = 0.001f;
	const int num_substeps = 10;
	
	for (int substep = 0; substep < num_substeps; substep++) advanceOneTimeStep(mt, thread_id, dt);
}

void FluidSimulation::advanceOneTimeStep(MT* mt, const int thread_id, const T dt)
{
	ONE_THREAD_WORK(dt_ = dt);

	advectParticles(mt, thread_id);

	addNewParticlesFromSourceList(mt, thread_id);

	sorted_particles_.sort(mt, thread_id);

	ONE_THREAD_WORK(std::cout << "Num pts " << sorted_particles_.num_all_pts_ << std::endl;);

	getGridVelocityFromParticles(mt, thread_id);

	saveVelocity(mt, thread_id);

	project(mt, thread_id, dt);

	saveVelocityDifference(mt, thread_id);
}

void FluidSimulation::updateDataDepot(MT* mt, const int thread_id, DataDepot* data_depot)
{
	BEGIN_ONE_THREAD_WORK
	{
		data_depot->reset();
	}
	END_ONE_THREAD_WORK;

	static ColoredParticlesData *flip_particles(nullptr);

	BEGIN_ONE_THREAD_WORK
	{
		flip_particles = new ColoredParticlesData;
		flip_particles->name_ = std::string("flip_particles");
		flip_particles->point_size_ = 2.0f;
		flip_particles->position_.initialize(sorted_particles_.num_all_pts_);
		flip_particles->color_.initialize(sorted_particles_.num_all_pts_);

		data_depot->colored_particles_list_.pushBack(flip_particles);
	}
	END_ONE_THREAD_WORK;

	const glm::vec4 slow_color(0.0f, 0.0f, 1.0f, 0.0f), fast_color(0.9f, 0.9f, 0.9f, 0.0f);
	const float max_vel = 4.0f;

	BEGIN_1D_ITERATION(sorted_particles_.num_all_pts_)
	{
		flip_particles->position_[p] = sorted_particles_.particle_data_arrays_.pts_pos_array_[p];

		const float alpha = MIN2(1.0f, sorted_particles_.particle_data_arrays_.pts_vel_array_[p].getMagnitude() / max_vel);

		flip_particles->color_[p] = slow_color * (1.0f - alpha) + fast_color * alpha;
	}
	END_1D_ITERATION;

	static LinesData *lines_data_temp = nullptr;

	BEGIN_ONE_THREAD_WORK
	{
		lines_data_temp = new LinesData;

		data_depot->lines_list_.pushBack(lines_data_temp);

		lines_data_temp->color_ = glm::vec4(0, 0, 0, 0);

		GL_TOOLS::AddCubeEdges(grid_.getMimMax(), lines_data_temp->vertices_);

		lines_data_temp = new LinesData;

		data_depot->lines_list_.pushBack(lines_data_temp);

//		std::cout << "adding line data" << std::endl;

		lines_data_temp->color_ = glm::vec4(1, 0, 0, 0);

		for (int i = 0; i < fixed_solid_list_.num_elements_; i++)
			GL_TOOLS::AddCubeEdges(fixed_solid_list_[i], lines_data_temp->vertices_);
	}
	END_ONE_THREAD_WORK;

	static PhongTrianglesData *phong_triangles_temp = nullptr;

	BEGIN_ONE_THREAD_WORK
	{
		phong_triangles_temp = new PhongTrianglesData;

		data_depot->phong_triangles_list_.pushBack(phong_triangles_temp);
	}
	END_ONE_THREAD_WORK;

	// generate water levelset
	BEGIN_1D_ITERATION(grid_ghost_.getNumAllCells())
	{
		water_levelset_.values_[p] = grid_ghost_.dx_;
	}
	END_1D_ITERATION;

	BEGIN_1D_ITERATION(sorted_particles_.num_all_pts_)
	{
		water_levelset_(grid_ghost_.getCell(sorted_particles_.particle_data_arrays_.pts_pos_array_[p])) = -grid_ghost_.dx_;
	}
	END_1D_ITERATION;

	static MarchingCubesAlgorithm mc_;
	static StaticTriangularSurface water_surface_;

	// polygonize water surface
	mc_.polygonize(mt, thread_id, grid_, water_levelset_, water_surface_);

	// copy water surface data
	BEGIN_ONE_THREAD_WORK
	{
//		water_surface_.

		water_surface_.copyRenderingData(phong_triangles_temp->positions_, phong_triangles_temp->normals_);	// TODO: multithreading
	}
	END_ONE_THREAD_WORK;
}