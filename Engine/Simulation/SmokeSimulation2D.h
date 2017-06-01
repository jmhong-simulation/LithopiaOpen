// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "SortedParticlesUniform3D.h"
#include "DataStructure/Array3D.h"
#include "DataStructure/CSRMatrix.h"
#include "DataStructure/GridUniform3D.h"
#include "DataStructure/VectorND.h"
#include "Geometry/MarchingCubesAlgorithm.h"
#include "Parallelism/MultiThreading.h"
#include "Utilities/ScriptReader.h"
#include "Framework/DataDepot.h"
#include "Framework/TaskManager.h"

#define BC_NULL			 1
#define BC_FULL			 0
#define BC_DIRICHLET	-1
#define BC_NEUMANN		-2		//NOTE: -2 means wall or mid cell Neumann conditions. -3 means rigid body -4 means rigid body 1, and so on.

class SmokeSimulation2D : public TaskManager
{
public:
	T dt_, time_;
	TV gravity_;

	int num_pts_per_cell_per_second_;
	TV source_velocity_;
	BOX_3D<T> source_box_;
	Array1D<BOX_3D<T>> fixed_solid_list_;

	Array1D<RandomNumberGenerator> random_;

	GridUniform3D grid_;
	GridUniform3D grid_ghost_;
	GridUniform3D u_grid_ghost_, v_grid_ghost_, w_grid_ghost_;

	SortedParticlesUniform3D<ParticleDataPV3D> sorted_particles_;

	Array3D<T> u_velocity_array_, u_velocity_array_temp_;		// both of them uses grid_ghost_
	Array3D<T> v_velocity_array_, v_velocity_array_temp_;		// both of them uses grid_ghost_
	Array3D<T> w_velocity_array_, w_velocity_array_temp_;		// both of them uses grid_ghost_

	Array3D<T> water_levelset_;

	Array3D<int> bc_array_;		// boundary condition
	Array3D<T> div_array_;		// divergence
	Array3D<T> p_array_;		// pressure

	CSRMatrix<D> A_matrix_;
	VectorND<D>  x_vector_;
	VectorND<D>  b_vector_;

	int num_iteration_;
	int max_iteration_;

	//TODO: initialize these
	D tolerance_;
	D sqr_tolerance_;
	D residual_;

	bool use_cuda_;
	int dev_id_;

	VectorND<D> res_, p_, Ap_;	// res_: residual vector;

public:
	SmokeSimulation2D(DataDepot* data)
		: TaskManager(data)
	{}

	void initialize();
	void updateOneStep(MT* mt, const int thread_id);
	void advanceOneTimeStep(MT* mt, const int thread_id, const T dt);

	void addNewParticlesFromSourceList(MT* mt, const int& thread_id);
	void advectParticles(MT* mt, const int& thread_id);
	void advectEulerianVelocity(MT* mt, const int& thread_id);
	void getGridVelocityFromParticles(MT* mt, const int& thread_id);
	void initializeLinearSystem(MT* mt, const int& thread_id, const GridUniform3D& grid_ghost, Array3D<int>& boundary_array, CSRMatrix<D>& A_matrix, VectorND<D>& x_vector, VectorND<D>& b_vector);
	void project(MT* mt, const int thread_id, const T dt);		// get divergence-free (or divergence constrained) velocity field
	void saveVelocity(MT* mt, const int& thread_id);
	void saveVelocityDifference(MT* mt, const int& thread_id);
	void setupBoundaryCondition(MT* mt, const int& thread_id);
	void setupLinearSystem(MT* mt, const int& thread_id, const GridUniform3D& grid_ghost, Array3D<int>& boundary_array, Array3D<T>& divergence_array, Array3D<T>& pressure_array);
	void solveLinearSystemCPUCG(MT* mt, const int& thread_id, CSRMatrix<D>& A_matrix, VectorND<D>& x_vector, VectorND<D>& b_vector);
	void updateFaceVelocityHelper(const TV& center_ijk, const T& vel_compoment, const GridUniform3D& face_grid, const int &i, const int& j, const int& k, const Array3D<T>& face_vel_array, const Array3D<T>& face_vel_array_temp);
	void updateVelocityByPressureGradient(MT* mt, const int& thread_id);
	TV	 getCollisionVelocity(const TV& vel, const TV& normal, const T& res, const T& fric_dt);
	void updateDataDepot(MT* mt, const int thread_id, DataDepot* data_depot);

	void initializeSourceListFromScript(ScriptBlock& sb);
	void initializeObjectListFromScript(ScriptBlock& sb);
};