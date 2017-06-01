// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "Photo2_5D.h"
#include "Utilities/RandomNumberGenerator.h"
#include "Utilities/ScriptBlock.h"
#include "Utilities/ScriptReader.h"
#include "Utilities/ScriptUtility.h"
#include "GL/GL_TOOLS.h"
#include "Geometry/Sphere3D.h"
#include "Geometry/Intersection.h"
#include "Geometry/LITHOPHANE.h"

void Photo2_5D::initializeFromScript(const char *script_filename)
{
	ScriptReader script;
	script.readFile(script_filename);

	ScriptBlock &sb(*script.head_block_);

	image_.ReadFileAndFitHeight((sb.getValue("input_path", std::string()) + sb.getValue("input_filename", std::string())).c_str(), 300);

	const float width = 3.0f;
	const float height = width / (float)image_.data_.i_res_ * (float)image_.data_.j_res_;

	grid_.Initialize(0, 0, image_.data_.i_res_, image_.data_.j_res_, 0, 0, width, height);

	grid_.InitializeCellArray(height_);
	grid_.InitializeCellArray(vertex_flag_);

	height_.assignAllValues((T)0);
	vertex_flag_.assignAllValues((int)0);

	vertex_flag_(grid_.i_start_, grid_.j_start_) = -1;
//	vertex_flag_(grid_.i_start_, grid_.j_end_) = -1;
	vertex_flag_(grid_.i_end_, grid_.j_start_) = -1;
//	vertex_flag_(grid_.i_end_, grid_.j_end_) = -1;
}

void Photo2_5D::resetMLS()
{
	//TODO: remove all MLS constraints

	for (int j = height_.j_start_; j <= height_.j_end_; j++)
		for (int i = height_.i_start_; i <= height_.i_end_; i++)
		{
			if(vertex_flag_(i, j) == -1)
				mls_.AddConstraint(grid_.GetCellCenter(i, j), height_(i, j));
		}

	const int i_center = grid_.i_res_ / 2;

	mls_.addLineConstraint(grid_.GetCellCenter(i_center + 20, 100), grid_.GetCellCenter(i_center + 20, 200), grid_.dx_*(T)10, grid_.dx_*(T)40, 0.5);

	mls_.addLineConstraint(grid_.GetCellCenter(i_center - 0, 100), grid_.GetCellCenter(i_center - 0, 200), grid_.dx_*(T)00, grid_.dx_*(T)00, 0.5);

	// top 
	mls_.addLineConstraint(grid_.GetCellCenter(grid_.i_start_, grid_.j_end_), grid_.GetCellCenter(grid_.i_end_, grid_.j_end_), grid_.dx_*(T)0, grid_.dx_*(T)0, 0.5);

	// bottom
	mls_.addLineConstraint(grid_.GetCellCenter(grid_.i_start_, grid_.j_start_ + 50), grid_.GetCellCenter(grid_.i_end_, grid_.j_start_ + 50), grid_.dx_*(T)0, grid_.dx_*(T)0, 0.5);

	// left
	mls_.addLineConstraint(grid_.GetCellCenter(grid_.i_start_, grid_.j_start_), grid_.GetCellCenter(grid_.i_start_, grid_.j_end_), grid_.dx_*(T)0, grid_.dx_*(T)0, 0.5);

	// right
	mls_.addLineConstraint(grid_.GetCellCenter(grid_.i_end_, grid_.j_start_), grid_.GetCellCenter(grid_.i_end_, grid_.j_end_), grid_.dx_*(T)0, grid_.dx_*(T)0, 0.5);

// 	mls_.AddConstraint(grid_.GetCellCenter(i_center, 200), grid_.dx_*(T)40, 0.5);
// 	mls_.AddConstraint(grid_.GetCellCenter(i_center, 100), grid_.dx_*(T)20, 0.5);
//	mls_.AddConstraint(grid_.GetCellCenter(grid_.i_res_/2, grid_.j_res_/2), grid_.dx_*(T)20);

	//TODO: multithreading
	for (int j = height_.j_start_; j <= height_.j_end_; j++)
		for (int i = height_.i_start_; i <= height_.i_end_; i++)
		{
			height_(i, j) = mls_.getScalar(grid_.GetCellCenter(i, j));
		}
}

BOX_3D<T> Photo2_5D::getAABB()
{
	return BOX_3D<T>(grid_.x_min_, grid_.y_min_, -0.5f, grid_.x_max_, grid_.y_max_, 0.5f);
}

void Photo2_5D::updateOneStep(MT* mt, const int thread_id)
{
}

void Photo2_5D::smoothHeight(const int& repeat)
{
	const T alpha = (T)0.5;

	for (int r = 0; r < repeat; r++)
	{
		for (int j = height_.j_start_; j <= height_.j_end_; j++)
			for (int i = height_.i_start_; i <= height_.i_end_; i++)
			{
				if (vertex_flag_(i,j) != 0) continue;

				const T average = (height_.getClamped(i - 1, j) + height_.getClamped(i + 1, j) + height_.getClamped(i, j - 1) + height_.getClamped(i, j + 1)) * (T)0.25;

				height_(i, j) = height_(i, j) * ((T)1 - alpha) + average * alpha;
			}
	}
}

void Photo2_5D::smoothHeightParallel(MT* mt, const int thread_id, const int& repeat)
{
	const T alpha = (T)0.5;

	for (int r = 0; r < repeat; r++)
	{
		BEGIN_GRID_ITERATION_2D(grid_)
		{
			if (vertex_flag_(i, j) != 0) continue;

			const T average = (height_.getClamped(i - 1, j) + height_.getClamped(i + 1, j) + height_.getClamped(i, j - 1) + height_.getClamped(i, j + 1)) * (T)0.25;

			height_(i, j) = height_(i, j) * ((T)1 - alpha) + average * alpha;
		}
		END_GRID_ITERATION_2D;
	}
}

/*
void Photo2_5D::receiveMessage(ScriptParameterList& message_parser)
{
	if (message_parser.isValid("pointer_x") == true)
	{
		pointer_.x_ = message_parser.getValue("pointer_x", float());
	}

	if (message_parser.isValid("pointer_y") == true)
	{
		pointer_.y_ = message_parser.getValue("pointer_y", float());
	}

	std::cout << pointer_.x_ << " "<<pointer_.y_<<" "<<pointer_.z_ << std::endl;
}
*/

void Photo2_5D::pickVertex(const RAY& ray)
{
	const T radius = grid_.dx_ * (T)2;
	
	T min_t = 1e8;

	for(int j = grid_.j_start_; j <= grid_.j_end_; j ++)
		for (int i = grid_.i_start_; i <= grid_.i_end_; i++)
		{
			const Sphere3D sphere(TV(grid_.GetCellCenter(i, j), height_(i, j)), radius);

			T t;
			TV intersection_point;

			if (Intersection::checkRaySphere(ray, sphere, t, intersection_point) > 0)
			{
				if (t < min_t)
				{
					selected_index_ = TV2_INT(i, j);
					min_t = t;

//					std::cout << "Picked " << i << " " << j << std::endl;
				}
			}
		}	

	vertex_flag_(selected_index_) = -1;	// fix vertex
}

void Photo2_5D::updateDataDepot(MT* mt, const int thread_id, DataDepot* data_depot)
{
	BEGIN_ONE_THREAD_WORK
	{
		updateDataDepotSingle(data_depot);
	}
	END_ONE_THREAD_WORK;
}


void Photo2_5D::updateDataDepotSingle(DataDepot* data_depot)
{
	data_depot->reset();

	static ColoredParticlesData *color_particles(nullptr);

	{
		color_particles = new ColoredParticlesData;
		color_particles->name_ = std::string("photo_2_5D");
		color_particles->point_size_ = 4.0f;
		color_particles->position_.initialize(grid_.getNumAllCells());
		color_particles->color_.initialize(grid_.getNumAllCells());

		for (int j = grid_.j_start_; j <= grid_.j_end_; j++)
			for (int i = grid_.i_start_; i <= grid_.i_end_; i++)
			{
				const TV2 pos = grid_.GetCellCenter(i, j);
				const int ix1d = grid_.Get1DIndex(i, j);

				color_particles->position_.values_[ix1d] = TV(pos.x_, pos.y_, height_(i,j));

				const T r = image_.data_(i, j).r_ / 255.0f;
				const T g = image_.data_(i, j).g_ / 255.0f;
				const T b = image_.data_(i, j).b_ / 255.0f;

				glm::vec4 color;

				if (vertex_flag_(i, j) == 0) color = glm::vec4(r, g, b, 1.0f);
				else if (vertex_flag_(i, j) == -1) color = glm::vec4(1, 0, 0, 1);
				else if (vertex_flag_(i, j) == -2) color = glm::vec4(0, 0, 1, 1);

				color_particles->color_.values_[ix1d] = color;
			}

		data_depot->colored_particles_list_.pushBack(color_particles);
	}

 	static PhongTrianglesData *phong_triangles_temp = nullptr;
 	phong_triangles_temp = new PhongTrianglesData;
 	data_depot->phong_triangles_list_.pushBack(phong_triangles_temp);
	
	LITHOPHANE litho_maker;
	StaticTriangularSurface surface;
	litho_maker.InitializePlane(height_, surface, grid_.dx_, 1.0f, grid_.x_max_ - grid_.x_min_);
	surface.translate(TV(grid_.x_max_, 0, 0));
	surface.copyRenderingData(phong_triangles_temp->positions_, phong_triangles_temp->normals_);

	static LinesData *lines_data_temp = nullptr;

	lines_data_temp = new LinesData;
	data_depot->lines_list_.pushBack(lines_data_temp);

	lines_data_temp->color_ = glm::vec4(1, 0, 0, 0);

	for (int j = grid_.j_start_; j <= grid_.j_end_; j++)
		for (int i = grid_.i_start_; i <= grid_.i_end_; i++)
		{
			const TV2 pos = grid_.GetCellCenter(i, j);

			if (vertex_flag_(i, j) < 0)
			{
				GL_TOOLS::addCircleLines(TV(pos.x_, pos.y_, height_(i, j)), (T)2.0f*grid_.dx_, 10, lines_data_temp->vertices_);
			}
			else if (vertex_flag_(i, j) == -2)
			{
				GL_TOOLS::addCircleLines(TV(pos.x_, pos.y_, height_(i, j)), (T)2.0f*grid_.dx_, 10, lines_data_temp->vertices_);
			}
		}
// 
// 	for (int j = grid_.j_start_; j < grid_.j_end_; j++)
// 		for (int i = grid_.i_start_; i < grid_.i_end_; i++)
// 		{
// 
// 		}
// 
// 		BEGIN_ONE_THREAD_WORK
// 		{
// 			phong_triangles_temp = new PhongTrianglesData;
// 
// 			data_depot->phong_triangles_list_.pushBack(phong_triangles_temp);
// 
// 			// 		phong_triangles_temp->positions_.Reset();
// 			// 		phong_triangles_temp->normals_.Reset();
// 
// 			//	LevelsetUniform2D &levelset_(levelset_image_);
// 
// 			cutter_surface_.use_face_normal_ = true;
// 			cutter_surface_.CopyRenderingData(phong_triangles_temp->positions_, phong_triangles_temp->normals_);
// 		}
// 		END_ONE_THREAD_WORK;
//	}

/*
	if (selected_index_.x_ != -1)
	{
		{
			color_particles = new ColoredParticlesData;
			color_particles->name_ = std::string("photo_2_5D_pointer");
			color_particles->point_size_ = 20.0f;
			color_particles->position_.initialize(1);
			color_particles->color_.initialize(1);

			const TV2 pointer = grid_.GetCellCenter(selected_index_);


			for (int j = grid_.j_start_; j <= grid_.j_end_; j++)
				for (int i = grid_.i_start_; i <= grid_.i_end_; i++)
				{
					const TV2 pos = grid_.GetCellCenter(i, j);
					const int ix1d = grid_.Get1DIndex(i, j);

					color_particles->position_.values_[ix1d] = TV(pos.x_, pos.y_, height_(i, j));

					const T r = image_.data_(i, j).r_ / 255.0f;
					const T g = image_.data_(i, j).g_ / 255.0f;
					const T b = image_.data_(i, j).b_ / 255.0f;

					color_particles->color_.values_[ix1d] = glm::vec4(r, g, b, 1.0f);
				}

			//		std::cout << pointer_.x_ << " " << pointer_.y_ << std::endl;

			data_depot->colored_particles_list_.pushBack(color_particles);
		}
	}*/
}