// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "LevelsetEditor2D.h"
#include "GL/GL_TOOLS.h"
#include "Geometry/Sphere2D.h"

LevelsetEditor2D::LevelsetEditor2D(DataDepot* data_)
	: TaskManager(data_)
{}

void LevelsetEditor2D::initialize()
{
	grid_.Initialize(0, 0, 100, 100, 0, 0, 1, 1);
	levelset_.initialize(grid_);

//	Sphere2D object(TV2(0.5, 0.5), 0.3);
	BOX_3D<T> object(TV3((T)0.5, (T)0.5, (T)0.0), (T)0.3);

	for (int i = levelset_.grid_ghost_.i_start_; i <= levelset_.grid_ghost_.i_end_; i++)
		for (int j = levelset_.grid_ghost_.j_start_; j <= levelset_.grid_ghost_.j_end_; j++)
		{
			const TV2 p = grid_.GetCellCenter(i, j);
			if (object.getSignedDistance(TV(p.x_, p.y_, 0)) <= (T)0) levelset_.phi_(i, j) = -grid_.dx_;
			else levelset_.phi_(i, j) = grid_.dy_;
		}

	work_threads_.runWithID(&LevelsetUniform2D::reinitializeFSM, &levelset_, 8);
	work_threads_.joinAll();
}

void LevelsetEditor2D::updateOneStep(MT* mt, const int thread_id)
{
	for (int i = 0; i < 100; i++)
		levelset_.advanceNormalFlow(mt, thread_id, 0.01f, 0.000001f);		
}

void LevelsetEditor2D::updateDataDepot(MT* mt, const int thread_id, DataDepot* data_depot)
{
	BEGIN_ONE_THREAD_WORK
	{
		data_depot->reset();
	}
	END_ONE_THREAD_WORK;

	static ColoredParticlesData *phi_center_pts(nullptr);

	BEGIN_ONE_THREAD_WORK
	{
		phi_center_pts = new ColoredParticlesData;
		phi_center_pts->name_ = std::string("flip_particles");
		phi_center_pts->point_size_ = 4.0f;
		phi_center_pts->position_.initialize(levelset_.grid_ghost_.getNumAllCells());
		phi_center_pts->color_.initialize(levelset_.grid_ghost_.getNumAllCells());

		for (int j = levelset_.grid_ghost_.j_start_; j <= levelset_.grid_ghost_.j_end_; j++)
			for (int i = levelset_.grid_ghost_.i_start_; i <= levelset_.grid_ghost_.i_end_; i++)
			{
				const TV2 pos = levelset_.grid_ghost_.GetCellCenter(i, j);
				const int ix1d = levelset_.phi_.get1DIndex(i, j);
				phi_center_pts->position_.values_[ix1d] = TV(pos.x_, pos.y_, levelset_.phi_(i, j));
				if (levelset_.phi_(i, j) > 0.0f) phi_center_pts->color_.values_[ix1d] = glm::vec4(0, 0, 1, 0);
				else phi_center_pts->color_.values_[ix1d] = glm::vec4(1, 0, 0, 0);
			}

		data_depot->colored_particles_list_.pushBack(phi_center_pts);
	}
	END_ONE_THREAD_WORK;

	static LinesData *lines_data_temp = nullptr;

	BEGIN_ONE_THREAD_WORK
	{
		lines_data_temp = new LinesData;

		data_depot->lines_list_.pushBack(lines_data_temp);

		lines_data_temp->color_ = glm::vec4(0, 0, 0, 0);

		BOX_3D<T> box(grid_.x_min_, grid_.y_min_, 0.0f, grid_.x_max_, grid_.y_max_, 0.0f);

		GL_TOOLS::AddCubeEdges(box, lines_data_temp->vertices_);
	}
	END_ONE_THREAD_WORK;
}