// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once 

#include <atomic>
#include "DataStructure/Array1D.h"
#include "GENERIC_DEFINITIONS.h"

class ParticleDataPV3D	// position and velocity
{
public:
	std::atomic<int> num_added_pts;
	std::atomic<int> max_num_pts;

	TV *pts_pos_array_;
	TV *pts_vel_array_;

public:
	ParticleDataPV3D();

	void initialize(const int& max_num_particles_input);

	void addNewParticle(const TV& position, const TV& velocity);

	void copyFrom(const int& this_index, const int& source_index, const ParticleDataPV3D& source_arrays);

	void swapArray(ParticleDataPV3D& source_array);
};
