// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include <atomic>
#include "DataStructure/Array1D.h"
#include "DataStructure/GridUniform3D.h"
#include "Parallelism/MultiThreading.h"
//#include "SIM_DOMAIN_UNIFORM_3D.h"
#include "ParticleDataPV3D.h"

template<class T_PARTICLE_DATA_ARRAYS>
class SortedParticlesUniform3D
{
public:
//	SIM_DOMAIN_UNIFORM_3D domain_;
	GridUniform3D grid_;

	T_PARTICLE_DATA_ARRAYS particle_data_arrays_;		// reallocated at the head of rebuilding step for exact number of particles to save memory.
	T_PARTICLE_DATA_ARRAYS particle_data_arrays_temp_;	// deleted at the end of rebuilding step to save memory.

	std::atomic<int>* pts_id_list_;
	std::atomic<int>* pts_index_list_;

	std::atomic<int>* num_pts_array_;		// the number of particles contained in each cell	//TODO: define as integer, use casted as atomic<int>
	int* pts_start_ix_array_;	// the start index of the particle data arrays of the first particle in each cell

	int num_all_pts_;			// number of all particles. newly created particles are not included until data structure is rebuilt.

	std::atomic<int>  count_pts_num_;
	std::atomic<int>  count_pts_start_ix_;

public:
	SortedParticlesUniform3D()
		: num_pts_array_(0), pts_start_ix_array_(0)
	{}

	~SortedParticlesUniform3D()
	{}

public:
	void initialize(GridUniform3D& grid_input);
	void sortThread();

	void sort(MT* mt, const int& thread_id);
	void prepareForSorting(MT* mt, const int& thread_id);
	void buildNumPtsArray(MT* mt, const int& thread_id);
	void buildPtsStartIxArray(MT* mt, const int& thread_id);
	void copyParticles(MT* mt, const int& thread_id);
	void finalizeSorting(MT* mt, const int& thread_id);
};