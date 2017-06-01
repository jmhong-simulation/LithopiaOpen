// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "CookieCutterMaker.h"
#include "GL/GL_TOOLS.h"
#include "Geometry/Sphere2D.h"
#include "Image/IMAGE_2D.h"
#include "Geometry/LITHOPHANE.h"
#include "Geometry/LineSegment2D.h"
#include "DataStructure/ArrayTools.h"
#include "Operations/SMOOTHING_UNIFORM_2D.h"
#include "Utilities/ScriptReader.h"

CookieCutterMaker::CookieCutterMaker(DataDepot* data_)
	: TaskManager(data_), step_(0), shrink_wrapping_steps_(0)
{}

void CookieCutterMaker::initializeFromScript(const char *script_filename)
{
	ScriptReader script;
	script.readFile(script_filename);

	ScriptBlock &sb(*script.head_block_);

	// parameters from script
	const T image_width = sb.getValue("width", float());
	cutter_padding_ = sb.getValue("cutter_padding", float());
	const int max_height = sb.getValue("max_height", int());
	
	output_path_ = sb.getValue("input_path", std::string());
	const std::string input_stamp_filename = sb.getValue("input_path", std::string()) + sb.getValue("input_stamp_filename", std::string());
	const std::string input_cutter_filename = sb.getValue("input_path", std::string()) + sb.getValue("input_cutter_filename", std::string());
	output_filename_prefix_ = sb.getValue("output_filename_prefix", std::string());

	extend_smoothing_ = sb.getValue("extend_smoothing", int());
	litho_smoothing_ = sb.getValue("litho_smoothing", int());
	stamp_extend_ = sb.getValue("stamp_extend", float());
	cutter_thickness_ = sb.getValue("cutter_thickness", float());
	stamp_cutter_spacing_ = sb.getValue("stamp_cutter_spacing", float());

	cutter_base_thickness_ = sb.getValue("cutter_base_thickness", float());
	cutter_print_thickness_ = sb.getValue("cutter_print_thickness", float());
	stamp_base_thickness_ = sb.getValue("stamp_base_thickness", float());
	stamp_print_thickness_ = sb.getValue("stamp_print_thickness", float());
	edge_smoothing_radius_ = sb.getValue("edge_smoothing_radius", float());

	input_stamp_image_.ReadFileAndFitHeight(input_stamp_filename, max_height);
	input_cutter_image_.ReadFileAndFitHeight(input_cutter_filename, max_height);

	const T dx = image_width / (T)input_stamp_image_.res_x_;
	extend_ = (int)ceil((cutter_padding_ + stamp_extend_ + cutter_thickness_ + stamp_cutter_spacing_) / dx * (T)1.5);

	grid_.Initialize(0, 0, input_stamp_image_.res_x_ + extend_ * 2, input_stamp_image_.res_y_ + extend_ * 2, 0, 0, dx*(T)(input_stamp_image_.res_x_ + extend_ * 2), dx*(T)(input_stamp_image_.res_y_ + extend_ * 2));

	levelset_.initialize(grid_);

	levelset_.grid_ghost_.InitializeCellArray(height_map_);
	height_map_.assignAllValues(0.0f);

	// copy input image to the center of height map. note the reflectLeftRight
	for (int j = 0; j < input_cutter_image_.res_y_; j++)
		for (int i = 0; i < input_cutter_image_.res_x_; i++)
		{
			const T height = input_cutter_image_.data_(input_cutter_image_.res_x_ - i - 1, j).GetReversedGray();	// black to 1, white to 0

			height_map_(i + extend_, j + extend_) = height;
		}

	ArrayTools::fillGhostCells(levelset_.grid_, levelset_.grid_ghost_, height_map_);

	// flood filling level set
	for (int j = grid_.j_start_; j <= grid_.j_end_; j ++)
		for (int i = grid_.i_start_; i <= grid_.i_end_; i++)
		{
			if (height_map_(i, j) >= 0.9f) levelset_.phi_(i, j) = -1.0f;
			else levelset_.phi_(i, j) = 0.0f;
		}

	levelset_.phi_.FloodFill(0, 0, -1.0f, 1.0f);

	for (int j = grid_.j_start_; j <= grid_.j_end_; j++)
		for (int i = grid_.i_start_; i <= grid_.i_end_; i++)
		{
			if (levelset_.phi_(i, j) == 1.0f) levelset_.phi_(i, j) = (T)0.5*grid_.dx_;
			else levelset_.phi_(i, j) = -(T)0.5*grid_.dx_;
		}

	work_threads_.runWithID(&LevelsetUniform2D::reinitializeFSM, &levelset_, 8);
	work_threads_.joinAll();

	levelset_.phi_.smoothLaplacian(0.5f, 3);

	work_threads_.runWithID(&LevelsetUniform2D::reinitializeFSM, &levelset_, 8);
	work_threads_.joinAll();

	work_threads_.runWithID(&CookieCutterMaker::updateOneStep, this);
	work_threads_.joinAll();
}

void CookieCutterMaker::makeCutterPart(MT* mt, const int thread_id)
{
	const float base_thickness = cutter_base_thickness_;
	const float front_print_thickness = cutter_print_thickness_;
	const float back_print_thickness = 0.0f;
	const int	smoothing = 0;

	static Array2D<T> shape_image, litho_image;

	BEGIN_ONE_THREAD_WORK
	{
		grid_.InitializeCellArray(shape_image);
		grid_.InitializeCellArray(litho_image);

		shape_image.assignAllValues(1);
		litho_image.assignAllValues(1);
	}
	END_ONE_THREAD_WORK;

	// space between stamp and cutter
	ArrayTools::addAllValues(mt, thread_id, -stamp_cutter_spacing_, levelset_.phi_);
	levelset_.reinitializeFSM(mt, thread_id, 8);

	BEGIN_GRID_ITERATION_2D(grid_)
	{
		if (levelset_.phi_(i, j) <= (T)0) shape_image(i, j) = 0.0f;
	}
	END_GRID_ITERATION_2D;

	// move 2 mm forward (cutter thickness)
	ArrayTools::addAllValues(mt, thread_id, -cutter_thickness_, levelset_.phi_);
	levelset_.reinitializeFSM(mt, thread_id, 8);

	BEGIN_GRID_ITERATION_2D(grid_)
	{
		if (levelset_.phi_(i, j) > (T)0) litho_image(i, j) = (T)0;
	}
	END_GRID_ITERATION_2D;

	ArrayTools::addAllValues(mt, thread_id, -cutter_padding_, levelset_.phi_);
	levelset_.reinitializeFSM(mt, thread_id, 8);

	BEGIN_GRID_ITERATION_2D(grid_)
	{
		if (levelset_.phi_(i, j) > (T)0) shape_image(i, j) = 0.0f;
	}
	END_GRID_ITERATION_2D;

	ArrayTools::smoothLaplacian(mt, thread_id, 0.5f, 10, shape_image);
	ArrayTools::smoothLaplacian(mt, thread_id, 0.5f, litho_smoothing_, litho_image);

//	cutter_surface_.reset(); //TODO implement

	static LITHOPHANE lithophane_maker;

	BEGIN_ONE_THREAD_WORK
	{
		lithophane_maker.makeFreeformStamp(shape_image, litho_image, cutter_surface_, base_thickness, front_print_thickness, back_print_thickness, grid_.x_max_ - grid_.x_min_, smoothing);
	}
	END_ONE_THREAD_WORK;
}

void CookieCutterMaker::makeStampPart(MT* mt, const int thread_id)
{
	// move 1 mm forward
	ArrayTools::addAllValues(mt, thread_id, -stamp_extend_, levelset_.phi_);
	levelset_.reinitializeFSM(mt, thread_id, 8);

	static Array2D<T> shape_image;

	BEGIN_ONE_THREAD_WORK
	{
		grid_.InitializeCellArray(shape_image);
	}
	END_ONE_THREAD_WORK;

	// make shape image (TODO: make it smoother)
	BEGIN_ONE_THREAD_WORK
	{
		for (int j = grid_.j_start_; j <= grid_.j_end_; ++j)
		for (int i = grid_.i_start_; i <= grid_.i_end_; ++i)
			shape_image(i, j) = CLAMP(-levelset_.phi_(i, j) * grid_.one_over_dx_ + (T)0.5, (T)0, (T)1);

		shape_image.smoothLaplacian((T)0.5, 10);

		const T alpha = 0.5f, min = 0.05f, max = 0.95f;

		height_map_.smoothMaxClamped(alpha, max, extend_smoothing_);
		height_map_.clampMinMaxZeroOne((T)0.5, (T)0.5);
		height_map_.smoothLaplacian(alpha, litho_smoothing_);

		ArrayTools::fillGhostCells(levelset_.grid_, levelset_.grid_ghost_, height_map_);
	}
	END_ONE_THREAD_WORK;

	static LITHOPHANE lithophane_maker;

	BEGIN_ONE_THREAD_WORK
	{
		lithophane_maker.makeFreeformStamp(shape_image, height_map_, stamp_surface_, stamp_base_thickness_, stamp_print_thickness_, 0.0f, grid_.x_max_ - grid_.x_min_, 0);
//		lithophane_maker.makeFreeformStampDoubleSided(shape_image, height_map_, stamp_surface_, stamp_base_thickness_, stamp_print_thickness_ - edge_smoothing_radius_, 0.0f, grid_.x_max_ - grid_.x_min_, 0);
	}
	END_ONE_THREAD_WORK;

	if (edge_smoothing_radius_ > (T)0)	// simplification
	{
		BEGIN_ONE_THREAD_WORK
		{
			const int num_v = stamp_surface_.vertex_positions_.num_elements_;
			const int num_tri = stamp_surface_.triangles_.num_elements_;

//			while (true) if (stamp_surface_.removeShortEdges(edge_smoothing_radius_ * (T)1.2, 1e8) == 0) break;
			
			stamp_surface_.findAdjacentTrianglesOfVertices();
			stamp_surface_.findEdgeTrianglesOfTriangles();

			for(int i = 1; i < 10; i ++)
				while (true) if (stamp_surface_.removeShortEdges(grid_.dx_ * (T)i, (T)1e-4) == 0) break;

			std::cout << "# vertices reduced to " << (T)stamp_surface_.vertex_positions_.num_elements_ / (T)num_v * (T)100 << "%, # triangles reduces to " <<
				(T)stamp_surface_.triangles_.num_elements_ / (T)num_tri * (T)100 << "%" << std::endl;
		}
		END_ONE_THREAD_WORK;
	}

	if (false)	// plane cut test
	{
		BEGIN_ONE_THREAD_WORK
		{
			const PLANE plane(TV(0, 0, 1), TV(0, 0, stamp_base_thickness_ + stamp_print_thickness_ - edge_smoothing_radius_));

			stamp_surface_.findAdjacentTrianglesOfVertices();
			stamp_surface_.findEdgeTrianglesOfTriangles();

			stamp_surface_.applyPlaneCutSuvdivision(plane);

			stamp_surface_.findAdjacentTrianglesOfVertices();
			stamp_surface_.determineFaceAveragedVertexNormals();			
		}
		END_ONE_THREAD_WORK;
	}

	if (edge_smoothing_radius_ > (T)0)
	{
		const T iso_value = (T)0.5; // hard coded in lithophane_maker.makeFreeformStamp

		BEGIN_ONE_THREAD_WORK
		{
			levelset_stamp_.initialize(grid_, 3);

//			stamp_surface_.ApplySubdivision(0);
		}
		END_ONE_THREAD_WORK;
/*
		levelset_stamp_.initializeFromScalarField(mt, thread_id, height_map_, iso_value, true);
 		levelset_stamp_.reinitializeInterfacialCells(mt, thread_id);
 		levelset_stamp_.reinitializeFSM(mt, thread_id, 4);
// 		for(int repeat = 0; repeat < 10; repeat ++) levelset_stamp_.smoothLaplacian(mt, thread_id);
// 		levelset_stamp_.reinitializeFSM(mt, thread_id, 4);
		levelset_stamp_.computeNormals(mt, thread_id);

		BEGIN_1D_ITERATION(stamp_surface_.vertex_positions_.num_elements_)
		{
			if (stamp_surface_.vertex_positions_.values_[p].z_ >= stamp_base_thickness_ + stamp_print_thickness_ - (T)1e-8)
			{
				const T phi = levelset_stamp_.getSignedDistance(TV2(stamp_surface_.vertex_positions_.values_[p].x_, stamp_surface_.vertex_positions_.values_[p].y_));

//				if (phi <= (T)0)
				{
//					const T edge_thickness = CLAMP(-phi, (T)0, edge_smoothing_radius_);
//					const T edge_thickness = ABS(phi);

//					stamp_surface_.vertex_positions_.values_[p].z_ = stamp_base_thickness_ + stamp_print_thickness_ + edge_thickness;

					// using level set normals
					TV2 pos = TV2(stamp_surface_.vertex_positions_.values_[p].x_, stamp_surface_.vertex_positions_.values_[p].y_);

					pos -= levelset_stamp_.getUnitNormal(pos) * edge_smoothing_radius_;

					stamp_surface_.vertex_positions_.values_[p].x_ = pos.x_;
					stamp_surface_.vertex_positions_.values_[p].y_ = pos.y_;
				}
			}
		}
		END_1D_ITERATION;
*/

		BEGIN_ONE_THREAD_WORK
		{
			stamp_surface_.findAdjacentTrianglesOfVertices();
			stamp_surface_.determineFaceAveragedVertexNormals();
		}
		END_ONE_THREAD_WORK;
	}
}

BOX_3D<T> CookieCutterMaker::getAABB()
{
	return BOX_3D<T>(grid_.x_min_, grid_.y_min_, -0.5f, grid_.x_max_, grid_.y_max_, 0.5f);
}

void CookieCutterMaker::updateOneStep(MT* mt, const int thread_id)
{
	makeStampPart(mt, thread_id);

//	makeCutterPart(mt, thread_id);

// 	BEGIN_ONE_THREAD_WORK
// 	{
// 		cutter_surface_.WriteSTL((output_path_ + output_filename_prefix_ + std::string("_cutter.stl")).c_str());
// 		std::cout << "End writing cutter part stl" << std::endl;
// 
// 		stamp_surface_.WriteSTL((output_path_ + output_filename_prefix_ + std::string("_stamp.stl")).c_str());
// 		std::cout << "End writing stamp part stl" << std::endl;
// 	}
// 	END_ONE_THREAD_WORK;

	pause();
}

void CookieCutterMaker::writeFiles()
{
	cutter_surface_.writeSTL((output_path_ + output_filename_prefix_ + std::string("_cutter.stl")).c_str());
	std::cout << "End writing cutter part stl" << std::endl;

	stamp_surface_.writeSTL((output_path_ + output_filename_prefix_ + std::string("_stamp.stl")).c_str());
	std::cout << "End writing stamp part stl" << std::endl;

 	stamp_surface_.writeOBJ((output_path_ + output_filename_prefix_ + std::string("_stamp.obj")).c_str());
 	std::cout << "End writing stamp part obj" << std::endl;
}

void CookieCutterMaker::updateDataDepot(MT* mt, const int thread_id, DataDepot* data_depot)
{
	BEGIN_ONE_THREAD_WORK
	{
		data_depot->reset();
	}
	END_ONE_THREAD_WORK;

// 	if (step_ <= 2)
// 	{
// 		static PhongTrianglesData *phong_triangles_temp = nullptr;
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
// 			LevelsetUniform2D &levelset_(levelset_);
// 
// 			for (int j = grid_.j_start_; j < grid_.j_end_; j++)
// 				for (int i = grid_.i_start_; i < grid_.i_end_; i++)
// 				{
// 					phong_triangles_temp->positions_.PushBack() = levelset_.getPhiPos(i, j);
// 					phong_triangles_temp->positions_.PushBack() = levelset_.getPhiPos(i + 1, j);
// 					phong_triangles_temp->positions_.PushBack() = levelset_.getPhiPos(i, j + 1);
// 
// 					phong_triangles_temp->positions_.PushBack() = levelset_.getPhiPos(i + 1, j);
// 					phong_triangles_temp->positions_.PushBack() = levelset_.getPhiPos(i + 1, j + 1);
// 					phong_triangles_temp->positions_.PushBack() = levelset_.getPhiPos(i, j + 1);
// 
// 					TV2 gradient;
// 
// 					gradient = -levelset_.getGradient(i, j);
// 					phong_triangles_temp->normals_.PushBack() = TV(gradient.x_, gradient.y_, 1).getSafeNormalized();
// 
// 					gradient = -levelset_.getGradient(i + 1, j);
// 					phong_triangles_temp->normals_.PushBack() = TV(gradient.x_, gradient.y_, 1).getSafeNormalized();
// 
// 					gradient = -levelset_.getGradient(i, j + 1);
// 					phong_triangles_temp->normals_.PushBack() = TV(gradient.x_, gradient.y_, 1).getSafeNormalized();
// 
// 					gradient = -levelset_.getGradient(i + 1, j);
// 					phong_triangles_temp->normals_.PushBack() = TV(gradient.x_, gradient.y_, 1).getSafeNormalized();
// 
// 					gradient = -levelset_.getGradient(i + 1, j + 1);
// 					phong_triangles_temp->normals_.PushBack() = TV(gradient.x_, gradient.y_, 1).getSafeNormalized();
// 
// 					gradient = -levelset_.getGradient(i, j + 1);
// 					phong_triangles_temp->normals_.PushBack() = TV(gradient.x_, gradient.y_, 1).getSafeNormalized();
// 				}
// 		}
// 		END_ONE_THREAD_WORK;
// 	}

	static PhongTrianglesData *phong_triangles_temp = nullptr;

//	if (step_ >= 3)
	{
		BEGIN_ONE_THREAD_WORK
		{
			phong_triangles_temp = new PhongTrianglesData;

			data_depot->phong_triangles_list_.pushBack(phong_triangles_temp);

			// 		phong_triangles_temp->positions_.Reset();
			// 		phong_triangles_temp->normals_.Reset();

			//	LevelsetUniform2D &levelset_(levelset_image_);
			stamp_surface_.use_face_normal_ = true;
			stamp_surface_.copyRenderingData(phong_triangles_temp->positions_, phong_triangles_temp->normals_);
		}
		END_ONE_THREAD_WORK;

		BEGIN_ONE_THREAD_WORK
		{
			phong_triangles_temp = new PhongTrianglesData;

			data_depot->phong_triangles_list_.pushBack(phong_triangles_temp);

			// 		phong_triangles_temp->positions_.Reset();
			// 		phong_triangles_temp->normals_.Reset();

			//	LevelsetUniform2D &levelset_(levelset_image_);

			cutter_surface_.use_face_normal_ = true;
			cutter_surface_.copyRenderingData(phong_triangles_temp->positions_, phong_triangles_temp->normals_);
		}
		END_ONE_THREAD_WORK;
	}
}