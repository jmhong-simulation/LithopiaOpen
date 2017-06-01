// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "SortedParticlesUniform3D.h"

template<class T_PARTICLE_DATA_ARRAYS>
void SortedParticlesUniform3D<T_PARTICLE_DATA_ARRAYS>::initialize(GridUniform3D& grid_input)
{
	grid_.initialize(grid_input);
//	domain_.InitializeMultithreading();

	num_all_pts_ = 0;

	num_pts_array_ = new std::atomic<int>[grid_.getNumAllCells()];
	pts_start_ix_array_ = new int [grid_.getNumAllCells()];

	for (int i = 0; i < grid_.getNumAllCells(); ++i) num_pts_array_[i] = 0;
	for (int i = 0; i < grid_.getNumAllCells(); ++i) pts_start_ix_array_[i] = -1;

	const int reserve_pts_num = 10000000;
	particle_data_arrays_.initialize(reserve_pts_num);
	particle_data_arrays_temp_.initialize(reserve_pts_num);

	pts_id_list_ = new std::atomic<int>[reserve_pts_num];
	pts_index_list_ = new std::atomic<int>[reserve_pts_num];
}

template<class T_PARTICLE_DATA_ARRAYS>
void SortedParticlesUniform3D<T_PARTICLE_DATA_ARRAYS>::prepareForSorting(MT* mt, const int& thread_id)
{
	ONE_THREAD_WORK(count_pts_num_ = 0;);
	ONE_THREAD_WORK(count_pts_start_ix_ = 0;);
	ONE_THREAD_WORK(num_all_pts_ = particle_data_arrays_.num_added_pts);

	BEGIN_BOX_ITERATION_ARR_3D_SYNC(grid_)
	{
		num_pts_array_[arr_ix] = 0;
		pts_start_ix_array_[arr_ix] = 0;
	}
	END_GRID_ITERATION_Z_SYNC;
}

template<class T_PARTICLE_DATA_ARRAYS>
void SortedParticlesUniform3D<T_PARTICLE_DATA_ARRAYS>::buildNumPtsArray(MT* mt, const int& thread_id)
{
	BEGIN_1D_ITERATION(num_all_pts_)
	{
		const TV& pts_pos = particle_data_arrays_.pts_pos_array_[p];
//		const TV_INT cell_index((int)pts_pos.x_, (int)pts_pos.y_, (int)pts_pos.z_);	
		const TV_INT cell_index(grid_.getCell(pts_pos));

		if(grid_.isInside(cell_index) == false) continue;

		const int arr_ix = grid_.get1DIndex(cell_index.i_, cell_index.j_, cell_index.k_);

		const int id = atomic_fetch_add(&num_pts_array_[arr_ix], 1);		

		pts_id_list_[p] = id;

		atomic_fetch_add(&count_pts_num_, 1);
	}
	END_1D_ITERATION;
}

template<class T_PARTICLE_DATA_ARRAYS>
void SortedParticlesUniform3D<T_PARTICLE_DATA_ARRAYS>::buildPtsStartIxArray(MT* mt, const int& thread_id)
{
	ONE_THREAD_WORK(count_pts_start_ix_ = 0;);

	int count_num_pts = 0;
	BEGIN_BOX_ITERATION_ARR_3D_SYNC(grid_)
	{
		pts_start_ix_array_[arr_ix] = count_num_pts;

		count_num_pts += num_pts_array_[arr_ix];		
	}
	END_GRID_ITERATION_Z_SYNC;

	const Vector2D<int> range = mt->getIncrementalRange(thread_id, count_num_pts);

	BEGIN_BOX_ITERATION_ARR_3D_SYNC(grid_)
	{
//		const int num = num_pts_array_[arr_ix];
//		pts_start_ix_array_[arr_ix] = atomic_fetch_add(&count_pts_start_ix_, num);

		pts_start_ix_array_[arr_ix] += range.t_min_;
	}
	END_GRID_ITERATION_Z_SYNC;
}

template<class T_PARTICLE_DATA_ARRAYS>
void SortedParticlesUniform3D<T_PARTICLE_DATA_ARRAYS>::copyParticles(MT* mt, const int& thread_id)
{
//	BEGIN_1D_ITERATION(num_all_pts_)
	{const TV2_INT k_range = mt->getParallelRange(thread_id, 0, num_all_pts_ - 1);
	const int _p_start(k_range.t_min_), _p_end(k_range.t_max_);
	for (int p = _p_start; p <= _p_end; p++)
	{
		const TV& pts_pos = particle_data_arrays_.pts_pos_array_[p];
		//const TV_INT cell_index((int)pts_pos.x_, (int)pts_pos.y_, (int)pts_pos.z_);	
		const TV_INT cell_index(grid_.getCell(pts_pos));

		if(grid_.isInside(cell_index) == false) continue;

		const int arr_ix = grid_.get1DIndex(cell_index.i_, cell_index.j_, cell_index.k_);

		const int start_ix = pts_start_ix_array_[arr_ix];
		const int id = pts_id_list_[p];

		pts_index_list_[start_ix+id] = p;
	}
//	END_1D_ITERATION;
	mt->sync(); }

	BEGIN_1D_ITERATION(count_pts_num_)
	{
		particle_data_arrays_temp_.copyFrom(p, pts_index_list_[p], particle_data_arrays_);
	}
	END_1D_ITERATION;

	BEGIN_ONE_THREAD_WORK
	{
		particle_data_arrays_.swapArray(particle_data_arrays_temp_);

		num_all_pts_ = (int)count_pts_num_;
		particle_data_arrays_.num_added_pts = (int)count_pts_num_;		
	}
	END_ONE_THREAD_WORK;
}

template<class T_PARTICLE_DATA_ARRAYS>
void SortedParticlesUniform3D<T_PARTICLE_DATA_ARRAYS>::sortThread()
{
	MultiThreading mt;
	mt.runWithID(&SortedParticlesUniform3D<T_PARTICLE_DATA_ARRAYS>::sort, this);
	mt.joinAll();
}

template<class T_PARTICLE_DATA_ARRAYS>
void SortedParticlesUniform3D<T_PARTICLE_DATA_ARRAYS>::sort(MT* mt, const int& thread_id)
{
	prepareForSorting(mt, thread_id);
	buildNumPtsArray(mt, thread_id);
	buildPtsStartIxArray(mt, thread_id);
	copyParticles(mt, thread_id);
}

template<class T_PARTICLE_DATA_ARRAYS>
void SortedParticlesUniform3D<T_PARTICLE_DATA_ARRAYS>::finalizeSorting(MT* mt, const int& thread_id)
{

}

template class SortedParticlesUniform3D<ParticleDataPV3D>;
