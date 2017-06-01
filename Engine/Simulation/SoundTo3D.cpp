// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "SoundTo3D.h"
#include "Utilities/RandomNumberGenerator.h"
#include "Utilities/ScriptBlock.h"
#include "Utilities/ScriptReader.h"
#include "Utilities/ScriptUtility.h"
#include "GL/GL_TOOLS.h"
#include "Geometry/Sphere3D.h"
#include "Geometry/Intersection.h"
#include "Geometry/LITHOPHANE.h"

void SoundTo3D::initializeFromScript(const char *script_filename)
{
	ScriptReader script;
	script.readFile(script_filename);

	ScriptBlock &sb(*script.head_block_);

	//TODO: script or automatically determine grid resolution
	grid_.Initialize(0, 0, 320/2 - 1, 320, 0, 0, 1, 2);
	grid_.InitializeCellArray(height_);
	current_j_ = 0;

	audio_spectrum_.initialize(320/2 - 1);

	height_.assignAllValues(0);

// 	InputTest input;
// 	input.show();
// 
// 	return a.exec();

// 	image_.ReadFileAndFitHeight((sb.getValue("input_path", std::string()) + sb.getValue("input_filename", std::string())).c_str(), 300);
// 
// 	const float width = 3.0f;
// 	const float height = width / (float)image_.data_.i_res_ * (float)image_.data_.j_res_;
// 
// 	grid_.Initialize(0, 0, image_.data_.i_res_, image_.data_.j_res_, 0, 0, width, height);
// 	grid_.InitializeCellArray(height_);
// 
// 	height_.assignAllValues((T)0);
}

void SoundTo3D::updateHeight(const Array1D<T>& audio_spectrum)
{
//	std::cout << audio_spectrum << std::endl;

// 	for (int i = 0; i < grid_.i_res_; i++)
// 		height_(i, 0) = audio_spectrum[i];

	for (int i = 0; i < grid_.i_res_; i++)
		height_(i, current_j_) = CLAMP(log(audio_spectrum[i]) / 5, 0, 1);

	current_j_++;

	if (current_j_ == grid_.j_res_)
	{
		current_j_ = 0;

		// write file
		IMAGE_2D image;
		image.Initialize(grid_.i_res_, grid_.j_res_);

		for (int j = 0; j < grid_.j_res_; j++)
			for (int i = 0; i < grid_.i_res_; i++)
			{
				image.data_(i, j).SetGray(1.0f - height_(i, j));
			}

		image.WriteBMP24("soundimage.bmp");
	}
}

BOX_3D<T> SoundTo3D::getAABB()
{
	return BOX_3D<T>(grid_.x_min_, grid_.y_min_, -2.0f, grid_.x_max_, grid_.y_max_, 2.0f);
}

SoundTo3D::~SoundTo3D()
{
	IMAGE_2D image;
	image.Initialize(grid_.i_res_, grid_.j_res_);

	for (int j = 0; j < grid_.j_res_; j++)
	for (int i = 0; i < grid_.i_res_; i++)
	{
		image.data_(i, j).SetGray(height_(i, j));
	}

	image.WriteBMP24("soundimage.bmp");
}

void SoundTo3D::updateOneStep(MT* mt, const int thread_id)
{
}

/*
void SoundTo3D::receiveMessage(ScriptParameterList& message_parser)
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

void SoundTo3D::updateDataDepot(MT* mt, const int thread_id, DataDepot* data_depot)
{
	BEGIN_ONE_THREAD_WORK
	{
		updateDataDepotSingle(data_depot);
	}
	END_ONE_THREAD_WORK;
}

void SoundTo3D::updateDataDepotSingle(DataDepot* data_depot)
{
	data_depot->reset();

	static LinesData *lines_data_temp = nullptr;

	{
		lines_data_temp = new LinesData;

		data_depot->lines_list_.pushBack(lines_data_temp);

		lines_data_temp->color_ = glm::vec4(1, 0, 0, 0);

//		int j = current_j_;
		for (int j = grid_.j_start_; j <= grid_.j_end_; j++)
			for (int i = grid_.i_start_; i <= grid_.i_end_; i++)
			{
				const TV2 pos = grid_.GetCellCenter(i, j);
				const int ix1d = grid_.Get1DIndex(i, j);

				lines_data_temp->vertices_.PushBack() = TV(pos.x_, pos.y_, 0);
				lines_data_temp->vertices_.PushBack() = TV(pos.x_, pos.y_, height_(i, j));
			}
	}
/*
	static ColoredParticlesData *color_particles(nullptr);

	{
		color_particles = new ColoredParticlesData;
		color_particles->name_ = std::string("sound_to_3d");
		color_particles->point_size_ = 4.0f;
		color_particles->position_.initialize(grid_.getNumAllCells());
		color_particles->color_.initialize(grid_.getNumAllCells());

		for (int j = grid_.j_start_; j <= grid_.j_end_; j++)
		for (int i = grid_.i_start_; i <= grid_.i_end_; i++)
		{
			const TV2 pos = grid_.GetCellCenter(i, j);
			const int ix1d = grid_.Get1DIndex(i, j);

			color_particles->position_.values_[ix1d] = TV(pos.x_, pos.y_, height_.values_[ix1d]);

			glm::vec4 color(1, 0, 0, 1);
			
			color_particles->color_.values_[ix1d] = color;
		}
		*/
// 		for (int j = grid_.j_start_; j <= grid_.j_end_; j++)
// 			for (int i = grid_.i_start_; i <= grid_.i_end_; i++)
// 			{
// 				const TV2 pos = grid_.GetCellCenter(i, j);
// 				const int ix1d = grid_.Get1DIndex(i, j);
// 
// 				color_particles->position_.values_[ix1d] = TV(pos.x_, pos.y_, height_(i,j));
// 
// 				const T r = image_.data_(i, j).r_ / 255.0f;
// 				const T g = image_.data_(i, j).g_ / 255.0f;
// 				const T b = image_.data_(i, j).b_ / 255.0f;
// 
// 				glm::vec4 color(1, 0, 0, 1);				
// 
// 				color_particles->color_.values_[ix1d] = color;
// 			}
// 
//		data_depot->colored_particles_list_.pushBack(color_particles);
//	}

//  	static PhongTrianglesData *phong_triangles_temp = nullptr;
//  	phong_triangles_temp = new PhongTrianglesData;
//  	data_depot->phong_triangles_list_.pushBack(phong_triangles_temp);
// 	
// 	LITHOPHANE litho_maker;
// 	StaticTriangularSurface surface;
// 	litho_maker.InitializePlane(height_, surface, grid_.dx_, 1.0f, grid_.x_max_ - grid_.x_min_);
// 	surface.translate(TV(grid_.x_max_, 0, 0));
// 	surface.copyRenderingData(phong_triangles_temp->positions_, phong_triangles_temp->normals_);

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