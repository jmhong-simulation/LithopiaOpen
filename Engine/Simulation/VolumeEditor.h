// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include <iostream>
#include "Framework/TaskManager.h"
#include "DataStructure/GridUniform2D.h"
#include "DataStructure/GridUniform3D.h"
#include "Geometry/MarchingCubesAlgorithm.h"
#include "Geometry/StaticTriangularSurface.h"

class VolumeEditor : public TaskManager
{
public:
	GridUniform3D grid_;
	Array3D<T>    density_;
	Array3D<T>	  phi_;

	GridUniform2D grid_x_, grid_y_, grid_z_;	// x -> yz plane, etc.
	Array2D<T>	  phi_x_, phi_y_, phi_z_;

	VolumeEditor(DataDepot* data_)
		: TaskManager(data_)
	{}

	void initializeFromUltrasonography(const char* filename)
	{
		grid_.initialize(0, 0, 0, 300, 300, 300, 0, 0, 0, 1, 1, 1);
		grid_.initializeCenterArray(phi_);
		grid_.initializeCenterArray(density_);

		grid_x_.Initialize(0, 0, grid_.j_res_, grid_.k_res_, 0, 0, 1, 1);
		grid_y_.Initialize(0, 0, grid_.k_res_, grid_.i_res_, 0, 0, 1, 1);
		grid_z_.Initialize(0, 0, grid_.i_res_, grid_.j_res_, 0, 0, 1, 1);

		grid_x_.InitializeCellArray(phi_x_);
		grid_y_.InitializeCellArray(phi_y_);
		grid_z_.InitializeCellArray(phi_z_);

		phi_x_.assignAllValues(0);
		phi_y_.assignAllValues(0);
		phi_z_.assignAllValues(0);

		// read raw data
		{
			using namespace std;
			ifstream file(filename, ios::in | ios::binary);

			float min = 1e8, max = -1e8;
			for (int k = grid_.k_start_; k <= grid_.k_end_; k++)
				for (int i = grid_.i_start_; i <= grid_.i_end_; i++)
					for (int j = grid_.j_start_; j <= grid_.j_end_; j++)
					{
						char temp;

						file.read((char*)&temp, sizeof(temp));						

						density_(i, j, k) = (float)((int)temp) / (float)128.0f;		// 0 ~ 1 density

						min = MIN2(density_(i, j, k), min);
						max = MAX2(density_(i, j, k), max);
					}

			std::cout << "Raw data min max "<< min << " " << max << std::endl;

			file.close();
		}

		const T th = -0.0f;

		// converting raw data to level set (density field)
		{
			for (int k = grid_.k_start_; k <= grid_.k_end_; k++)
				for (int i = grid_.i_start_; i <= grid_.i_end_; i++)
					for (int j = grid_.j_start_; j <= grid_.j_end_; j++)
					{
						phi_(i, j, k) = density_(i, j, k) - th;
					}
		}
	}

	void updateOneStep(MT* mt, const int thread_id)
	{
		static int frame = 0;

		// converting raw data to level set (density field)

		const T one_over_weight_sum = (T)1 / ((grid_.one_over_dx_ + grid_.one_over_dy_ + grid_.one_over_dz_) * 2.0f);

		for (int r = 0; r < 5; r++)
		{
			const T alpha = 0.5;

			BEGIN_GRID_ITERATION_3D(grid_)
			{
				const float average = (phi_.getClamped(i - 1, j, k) * grid_.one_over_dx_ + phi_.getClamped(i + 1, j, k) * grid_.one_over_dx_
					+ phi_.getClamped(i, j - 1, k) * grid_.one_over_dy_ + phi_.getClamped(i, j + 1, k) * grid_.one_over_dy_
					+ phi_.getClamped(i, j, k - 1) * grid_.one_over_dz_ + phi_.getClamped(i, j, k + 1) * grid_.one_over_dz_) * one_over_weight_sum;

				phi_(i, j, k) = average * alpha + (1.0f - alpha)*phi_(i, j, k);
			}
			END_GRID_ITERATION_3D;
		}

// 		BEGIN_GRID_ITERATION_3D(grid_)
// 		{
// 			phi_(i, j, k) = (density_(i, j, k) - 0.5f);
// 		}
// 		END_GRID_ITERATION_3D;
	}

	void updateDataDepot(MT* mt, const int thread_id, DataDepot* data_depot)
	{
		BEGIN_ONE_THREAD_WORK
		{
			data_depot->reset();
		}
		END_ONE_THREAD_WORK;

		static MarchingCubesAlgorithm mc_;
		static StaticTriangularSurface water_surface_;

		static PhongTrianglesData *phong_triangles_temp = nullptr;

		BEGIN_ONE_THREAD_WORK
		{
			phong_triangles_temp = new PhongTrianglesData;

			data_depot->phong_triangles_list_.pushBack(phong_triangles_temp);
		}
		END_ONE_THREAD_WORK;

		// polygonize water surface
		mc_.polygonize(mt, thread_id, grid_, phi_, water_surface_);

		// copy water surface data
		BEGIN_ONE_THREAD_WORK
		{
			//		water_surface_.
			water_surface_.findAdjacentTrianglesOfVertices();
			water_surface_.determineFaceAveragedVertexNormals();

			water_surface_.copyRenderingData(phong_triangles_temp->positions_, phong_triangles_temp->normals_);	// TODO: multithreading
		}
		END_ONE_THREAD_WORK;

		return;

		static ColoredParticlesData *phi_center_pts(nullptr);

		{
			GridUniform2D &grid_plane = grid_x_;
			Array2D<T> &phi_plane = phi_x_;

			BEGIN_ONE_THREAD_WORK
			{
				phi_center_pts = new ColoredParticlesData;
				phi_center_pts->name_ = std::string("flip_particles");
				phi_center_pts->point_size_ = 4.0f;
				phi_center_pts->position_.initialize(grid_plane.getNumAllCells());
				phi_center_pts->color_.initialize(grid_plane.getNumAllCells());

				for (int j = grid_plane.j_start_; j <= grid_plane.j_end_; j++)
					for (int i = grid_plane.i_start_; i <= grid_plane.i_end_; i++)
					{
						const TV2 pos = grid_plane.GetCellCenter(i, j);
						const int ix1d = phi_plane.get1DIndex(i, j);
						phi_center_pts->position_.values_[ix1d] = TV(0.0f, pos.x_, pos.y_);
						const T scale = phi_plane(i, j);

						phi_center_pts->color_.values_[ix1d] = glm::vec4(scale, scale, scale, 0.0f);

						// 					if (levelset_.phi_(i, j) > 0.0f) 
						// 					else phi_center_pts->color_.values_[ix1d] = glm::vec4(1, 0, 0, 0);
					}

				data_depot->colored_particles_list_.pushBack(phi_center_pts);
			}
			END_ONE_THREAD_WORK;
		}

		{
			GridUniform2D &grid_plane = grid_y_;
			Array2D<T> &phi_plane = phi_y_;

			BEGIN_ONE_THREAD_WORK
			{
				phi_center_pts = new ColoredParticlesData;
				phi_center_pts->name_ = std::string("flip_particles");
				phi_center_pts->point_size_ = 4.0f;
				phi_center_pts->position_.initialize(grid_plane.getNumAllCells());
				phi_center_pts->color_.initialize(grid_plane.getNumAllCells());

				for (int j = grid_plane.j_start_; j <= grid_plane.j_end_; j++)
					for (int i = grid_plane.i_start_; i <= grid_plane.i_end_; i++)
					{
						const TV2 pos = grid_plane.GetCellCenter(i, j);
						const int ix1d = phi_plane.get1DIndex(i, j);
						phi_center_pts->position_.values_[ix1d] = TV(pos.y_, 0.0f, pos.x_);
						const T scale = phi_plane(i, j);

						phi_center_pts->color_.values_[ix1d] = glm::vec4(scale, scale, scale, 0.0f);

						// 					if (levelset_.phi_(i, j) > 0.0f) 
						// 					else phi_center_pts->color_.values_[ix1d] = glm::vec4(1, 0, 0, 0);
					}

				data_depot->colored_particles_list_.pushBack(phi_center_pts);
			}
			END_ONE_THREAD_WORK;
		}

		{
			GridUniform2D &grid_plane = grid_z_;
			Array2D<T> &phi_plane = phi_z_;

			BEGIN_ONE_THREAD_WORK
			{
				phi_center_pts = new ColoredParticlesData;
				phi_center_pts->name_ = std::string("flip_particles");
				phi_center_pts->point_size_ = 4.0f;
				phi_center_pts->position_.initialize(grid_plane.getNumAllCells());
				phi_center_pts->color_.initialize(grid_plane.getNumAllCells());

				for (int j = grid_plane.j_start_; j <= grid_plane.j_end_; j++)
					for (int i = grid_plane.i_start_; i <= grid_plane.i_end_; i++)
					{
						const TV2 pos = grid_plane.GetCellCenter(i, j);
						const int ix1d = phi_plane.get1DIndex(i, j);
						phi_center_pts->position_.values_[ix1d] = TV(pos.x_, pos.y_, 0.0f);
						const T scale = phi_plane(i, j);

						phi_center_pts->color_.values_[ix1d] = glm::vec4(scale, scale, scale, 0.0f);

						// 					if (levelset_.phi_(i, j) > 0.0f) 
						// 					else phi_center_pts->color_.values_[ix1d] = glm::vec4(1, 0, 0, 0);
					}

				data_depot->colored_particles_list_.pushBack(phi_center_pts);
			}
			END_ONE_THREAD_WORK;
		}
	}

};