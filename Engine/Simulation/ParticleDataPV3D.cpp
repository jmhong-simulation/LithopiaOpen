// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "ParticleDataPV3D.h"

ParticleDataPV3D::ParticleDataPV3D()
	: pts_pos_array_(0), pts_vel_array_(0)
{}

void ParticleDataPV3D::initialize(const int& max_num_particles_input)
{
	num_added_pts = 0;
	max_num_pts = max_num_particles_input;

	SAFE_DELETE_ARRAY(pts_pos_array_);
	SAFE_DELETE_ARRAY(pts_vel_array_);

	pts_pos_array_ = new TV[max_num_pts];
	pts_vel_array_ = new TV[max_num_pts];
}

void ParticleDataPV3D::addNewParticle(const TV& position, const TV& velocity)
{
	assert(num_added_pts < max_num_pts);

	int ix = atomic_fetch_add(&num_added_pts, 1);

	pts_pos_array_[ix] = position;
	pts_vel_array_[ix] = velocity;
}

void ParticleDataPV3D::copyFrom(const int& this_index, const int& source_index, const ParticleDataPV3D& source_arrays)
{
	assert(this_index < max_num_pts);
	assert(source_index < max_num_pts);

	pts_pos_array_[this_index] = source_arrays.pts_pos_array_[source_index];
	pts_vel_array_[this_index] = source_arrays.pts_vel_array_[source_index];
}

void ParticleDataPV3D::swapArray(ParticleDataPV3D& source_array)
{
	SWAP(pts_pos_array_, source_array.pts_pos_array_, TV*);
	SWAP(pts_vel_array_, source_array.pts_vel_array_, TV*);
}