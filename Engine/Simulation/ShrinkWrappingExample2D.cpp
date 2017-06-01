// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "ShrinkWrappingExample2D.h"
#include "GL/GL_TOOLS.h"
#include "Geometry/Sphere2D.h"
#include "Image/IMAGE_2D.h"
#include "Geometry/LITHOPHANE.h"
#include "Geometry/LineSegment2D.h"
#include "DataStructure/ArrayTools.h"
#include "Operations/SMOOTHING_UNIFORM_2D.h"
#include "Utilities/ScriptReader.h"

ShrinkWrappingExample2D::ShrinkWrappingExample2D(DataDepot* data_)
	: TaskManager(data_), step_(0), shrink_wrapping_steps_(0)
{}

void ShrinkWrappingExample2D::initializeFromScript(const char *script_filename)
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

	input_stamp_image_.ReadFileAndFitHeight(input_stamp_filename, max_height);
	input_cutter_image_.ReadFileAndFitHeight(input_cutter_filename, max_height);

	const T dx = image_width / (T)input_stamp_image_.res_x_;
	extend_ = (int)ceil((cutter_padding_ + stamp_extend_ + cutter_thickness_ + stamp_cutter_spacing_) / dx * (T)1.5);

	grid_.Initialize(0, 0, input_stamp_image_.res_x_ + extend_ * 2, input_stamp_image_.res_y_ + extend_ * 2, 0, 0, dx*(T)(input_stamp_image_.res_x_ + extend_ * 2), dx*(T)(input_stamp_image_.res_y_ + extend_ * 2));

	levelset_.initialize(grid_);
	levelset_image_.initialize(grid_);

	levelset_.grid_ghost_.InitializeCellArray(height_map_);
	height_map_.assignAllValues(0.0f);

	// copy input image to the center of height map. note the reflectLeftRight
	for (int j = 0; j < input_cutter_image_.res_y_; j++)
		for (int i = 0; i < input_cutter_image_.res_x_; i++)
		{
			const T height = input_cutter_image_.data_(input_cutter_image_.res_x_ - i - 1, j).GetReversedGray();	// black to 1, white to 0

			height_map_(i + extend_, j + extend_) = height;
		}

	ArrayTools::fillGhostCells(levelset_image_.grid_, levelset_image_.grid_ghost_, height_map_);

	height_map_.clampMinMaxZeroOne((T)0.1, (T)0.6);

	// prepare for potential field (height field to levelset image
	work_threads_.runWithID(&LevelsetUniform2D::initializeFromScalarField, &levelset_image_, height_map_, (T)0.5, false);
	work_threads_.joinAll();

	work_threads_.runWithID(&LevelsetUniform2D::computeNormals, &levelset_image_);
	work_threads_.joinAll();

	prepareForShrinkWrapping();
}

void ShrinkWrappingExample2D::prepareForShrinkWrapping()
{
	// initialize shrink wrapping levelset
	for (int i = levelset_.grid_ghost_.i_start_; i <= levelset_.grid_ghost_.i_end_; i++)
		for (int j = levelset_.grid_ghost_.j_start_; j <= levelset_.grid_ghost_.j_end_; j++)
		{
			const TV2 p = grid_.GetCellCenter(i, j);

			if (grid_.Inside(TV2_INT(i, j), 3) == true) levelset_.phi_(i, j) = -grid_.dx_ * (T)0.5;
			else
				levelset_.phi_(i, j) = grid_.dx_  * (T)0.5;
		}

	LinkedArray<TV2> points_temp;
	points_temp.PushBack() = TV2(grid_.GetCellCenter(grid_.i_start_, grid_.j_start_));
	points_temp.PushBack() = TV2(grid_.GetCellCenter(grid_.i_end_, grid_.j_start_));
	points_temp.PushBack() = TV2(grid_.GetCellCenter(grid_.i_end_, grid_.j_end_));
	points_temp.PushBack() = TV2(grid_.GetCellCenter(grid_.i_start_, grid_.j_end_));
	points_temp.CopyToArray(contour_.vertex_positions_);

	LinkedArray<TV2_INT> lines_temp;
	lines_temp.PushBack() = TV2_INT(0, 1);
	lines_temp.PushBack() = TV2_INT(1, 2);
	lines_temp.PushBack() = TV2_INT(2, 3);
	lines_temp.PushBack() = TV2_INT(3, 0);
	lines_temp.CopyToArray(contour_.lines_);

	// split
	contour_.max_dl_ = grid_.dx_ * 1.0f;
	contour_.min_dl_ = grid_.dx_ * 1.0f * 0.25f;
	contour_.split();

	work_threads_.runWithID(&LevelsetUniform2D::reinitializeFSM, &levelset_, 8);
	work_threads_.joinAll();

	work_threads_.runWithID(&LevelsetUniform2D::computeNormals, &levelset_);
	work_threads_.joinAll();
}

void ShrinkWrappingExample2D::makeCutterPart(MT* mt, const int thread_id)
{
	const float base_thickness = cutter_base_thickness_;
	const float front_print_thickness = cutter_print_thickness_;
	const float back_print_thickness = 0.0f;
	//const float width = std::atof(argv[11]);
	const int smoothing = 0;

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
	levelset_.advanceNormalFlowMain(mt, thread_id, stamp_cutter_spacing_, 0.25f, 0.000001f);

	BEGIN_GRID_ITERATION_2D(grid_)
	{
		if (levelset_.phi_(i, j) <= (T)0) shape_image(i, j) = 0.0f;
	}
	END_GRID_ITERATION_2D;

	// move 2 mm forward (cutter thickness)
	levelset_.advanceNormalFlowMain(mt, thread_id, cutter_thickness_, 0.25f, 0.000001f);

	BEGIN_GRID_ITERATION_2D(grid_)
	{
		if (levelset_.phi_(i, j) > (T)0) litho_image(i, j) = (T)0;
	}
	END_GRID_ITERATION_2D;

	levelset_.advanceNormalFlowMain(mt, thread_id, cutter_padding_, 0.25f, 0.000001f);

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

void ShrinkWrappingExample2D::makeStampPart(MT* mt, const int thread_id)
{
	// move 1 mm forward
	levelset_.advanceNormalFlowMain(mt, thread_id, stamp_extend_, 0.25f, 0.000001f);

	const float base_thickness = stamp_base_thickness_;
	const float front_print_thickness = stamp_print_thickness_;
	const float back_print_thickness = 0.0f;
	const int smoothing = 0;

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
		lithophane_maker.makeFreeformStamp(shape_image, height_map_, stamp_surface_, base_thickness, front_print_thickness, back_print_thickness, grid_.x_max_ - grid_.x_min_, smoothing);
	}
	END_ONE_THREAD_WORK;
}

BOX_3D<T> ShrinkWrappingExample2D::getAABB()
{
	return BOX_3D<T>(grid_.x_min_, grid_.y_min_, -0.5f, grid_.x_max_, grid_.y_max_, 0.5f);
}

void ShrinkWrappingExample2D::findInnerline(MT* mt, const int thread_id)
{
	levelset_.advanceNormalFlowMain(mt, thread_id, 2.0f, 0.25f, 0.000001f);

	BEGIN_ONE_THREAD_WORK
	{
		inner_line_.initialize(levelset_);		//TODO: multithread
	}
	END_ONE_THREAD_WORK;
}

void ShrinkWrappingExample2D::updateOneStep(MT* mt, const int thread_id)
{
	findOutline(mt, thread_id);
}

void ShrinkWrappingExample2D::findOutlineSubstep(MT* mt, const int thread_id, const T max_dx)
{
	// update contour
	BEGIN_ONE_THREAD_WORK
	{
		contour_.split();
		contour_.removeShortLines();
	}
	END_ONE_THREAD_WORK;

	BEGIN_ONE_THREAD_WORK
	{
		contour_.vertex_velocities_.initialize(contour_.vertex_positions_.num_elements_);
	}
	END_ONE_THREAD_WORK;

	T max_vel_mag = 0;
	T phi_abs_sum = 0;

	// find particle velocity
	BEGIN_1D_ITERATION(contour_.vertex_positions_.num_elements_)
	{
		const TV2 point = contour_.vertex_positions_[p];

		const T phi = levelset_image_.grid_ghost_.GetClampedLinearInterpolationCell(levelset_image_.phi_, point);

		contour_.vertex_velocities_[p] = levelset_image_.grid_ghost_.GetClampedLinearInterpolationCell(levelset_image_.normal_, point).getSafeNormalized() * -phi;

		const T mag = contour_.vertex_velocities_[p].getMagnitude();

		max_vel_mag = MAX2(mag, max_vel_mag);

		phi_abs_sum += ABS(phi);
	}
	END_1D_ITERATION;

	phi_abs_sum = mt->syncSum(thread_id, phi_abs_sum);

	mt->syncMax(thread_id, max_vel_mag);

	T dt = max_dx / max_vel_mag;

	// removing intersecting lines (not necessary now)
	// 		BEGIN_ONE_THREAD_WORK
	// 		{
	// 			const int removed = contour_.removeIntersectingLines(dt);
	// 
	// 			if(removed > 0) std::cout << "# of intersecting lines " << removed << std::endl;
	// 
	// 			contour_.vertex_velocities_.initialize(contour_.vertex_positions_.num_elements_);
	// 		}
	// 		END_ONE_THREAD_WORK;
	// 
	// 		// find particle velocity
	// 		BEGIN_1D_ITERATION(contour_.vertex_positions_.num_elements_)
	// 		{
	// 			const TV2 point = contour_.vertex_positions_[p];
	// 
	// 			const T phi = levelset_image_.grid_ghost_.GetClampedLinearInterpolationCell(levelset_image_.phi_, point);
	// 
	// 			contour_.vertex_velocities_[p] = levelset_image_.grid_ghost_.GetClampedLinearInterpolationCell(levelset_image_.normal_, point).getSafeNormalized() * -phi;
	// 
	// 			const T mag = contour_.vertex_velocities_[p].getMagnitude();
	// 
	// 			max_vel_mag = MAX2(mag, max_vel_mag);
	// 
	// 			phi_abs_sum += ABS(phi);
	// 		}
	// 		END_1D_ITERATION;
	// 
	// 		phi_abs_sum = mt->syncSum(thread_id, phi_abs_sum);
	// 
	// 		mt->syncMax(thread_id, max_vel_mag);
	// 
	// 		dt = max_dx / max_vel_mag;

	if (max_vel_mag > (T)0)
	{
		// update topology
		BEGIN_1D_ITERATION(contour_.lines_.num_elements_)
		{
			Quadrilateral2D advance_quad = contour_.getAdvanceQuad(p, dt);

			BOX_2D<int> itr_box = levelset_.grid_ghost_.getIXBox(advance_quad.getAABB());

			for (int j = itr_box.j_start_; j <= itr_box.j_end_; j++)
				for (int i = itr_box.i_start_; i <= itr_box.i_end_; i++)
				{
					int inside_flag = advance_quad.getInsideFlag(levelset_.grid_ghost_.GetCellCenter(i, j));

					if (inside_flag > 0) levelset_.phi_(i, j) = -grid_.dx_ * (T)0.5;
					else if (inside_flag < 0) levelset_.phi_(i, j) = grid_.dx_ * (T)0.5;

					// 						if (inside_flag > 0) levelset_.phi_(i, j) = 1e8;		// TODO: use large positive float
					// 						else if (inside_flag < 0) levelset_.phi_(i, j) = -1e8;	// TODO: use large negative float
				}
		}
		END_1D_ITERATION;

		// advect particles
		BEGIN_1D_ITERATION(contour_.vertex_positions_.num_elements_)
		{
			contour_.vertex_positions_[p] += contour_.vertex_velocities_[p] * dt;
		}
		END_1D_ITERATION;
	}
	else return;

	//	break;	//TODO: stop condition
}

void ShrinkWrappingExample2D::findOutline(MT* mt, const int thread_id)
{
	if (step_ == 0)	// shrink wrapping step
	{
		const T CFL = 0.5f;
		const T max_dx = CFL*grid_.dx_;	// dt should be small enough
		const int repeat_max = (int)((T)MAX2(grid_.i_res_, grid_.j_res_) / CFL * (T)1.5);

		for (int repeat = 0; repeat < 30; repeat++)
		{
			findOutlineSubstep(mt, thread_id, max_dx);

			BEGIN_ONE_THREAD_WORK
			{
				shrink_wrapping_steps_++;
			}
			END_ONE_THREAD_WORK;
		}

		BEGIN_ONE_THREAD_WORK
		{
			if (shrink_wrapping_steps_ > repeat_max) step_++;		// go to next step
		}
		END_ONE_THREAD_WORK;
	}
	else if (step_ == 1)
	{
		BEGIN_1D_ITERATION(contour_.lines_.num_elements_)
		{
			const LineSegment2D line_seg(contour_.getP0(p), contour_.getP1(p));

			const BOX_2D<int> itr_box = levelset_.grid_ghost_.getIXBox(line_seg.getAABB()).getExtended(1);

			//TODO: atomic
			for (int j = itr_box.j_start_; j <= itr_box.j_end_; j++)
				for (int i = itr_box.i_start_; i <= itr_box.i_end_; i++)
				{
					const T dist = line_seg.getDistance(levelset_.grid_ghost_.GetCellCenter(i, j));
					const T abs_phi = MIN2(ABS(levelset_.phi_(i, j)), dist);

					if (levelset_.phi_(i, j) <= (T)0) levelset_.phi_(i, j) = -abs_phi;
					else levelset_.phi_(i, j) = abs_phi;
				}
		}
		END_1D_ITERATION;

		levelset_.reinitializeFSM(mt, thread_id, 8);
		//	levelset_.fixInterfacialCells(mt, thread_id);

		BEGIN_ONE_THREAD_WORK	//TODO: parallelize
		{
			height_map_.assignAllValues(0.0f);

			for (int j = 0; j < input_stamp_image_.res_y_; j++)
				for (int i = 0; i < input_stamp_image_.res_x_; i++)
				{
					T height = input_stamp_image_.data_(input_stamp_image_.res_x_ - i - 1, j).GetReversedGray();	// black to 1, white to 0

					height_map_(i + extend_, j + extend_) = height;
				}

			ArrayTools::fillGhostCells(levelset_image_.grid_, levelset_image_.grid_ghost_, height_map_);

			height_map_.clampMinMaxZeroOne((T)0.1, (T)0.6);
		}
		END_ONE_THREAD_WORK;

		ArrayTools::clampMinMaxZeroOne(mt, thread_id, 0.1f, 0.6f, height_map_);

		BEGIN_ONE_THREAD_WORK
		{
			step_++;
		}
		END_ONE_THREAD_WORK;
	}
	else if (step_ == 2)
	{
		makeStampPart(mt, thread_id);

		BEGIN_ONE_THREAD_WORK
		{
			step_++;
		}
		END_ONE_THREAD_WORK;
	}
	else if (step_ == 3)
	{
		makeCutterPart(mt, thread_id);

		BEGIN_ONE_THREAD_WORK
		{
			step_++;
		}
		END_ONE_THREAD_WORK;
	}
	else if (step_ == 4)
	{
		pause();

		// write files and end
		BEGIN_ONE_THREAD_WORK
		{
			cutter_surface_.writeSTL((output_path_ + output_filename_prefix_ + std::string("_cutter.stl")).c_str());
			std::cout << "End writing cutter part stl" << std::endl;

			stamp_surface_.writeSTL((output_path_ + output_filename_prefix_ + std::string("_stamp.stl")).c_str());
			std::cout << "End writing stamp part stl" << std::endl;
		}
		END_ONE_THREAD_WORK;
	}
}

void ShrinkWrappingExample2D::updateDataDepot(MT* mt, const int thread_id, DataDepot* data_depot)
{
	BEGIN_ONE_THREAD_WORK
	{
		data_depot->reset();
	}
	END_ONE_THREAD_WORK;

	static ColoredParticlesData *contour_pts(nullptr);

	// dynamic contour vertices
	if (step_ <= 1)
	{
		// 		BEGIN_ONE_THREAD_WORK
		// 		{
		// 			contour_pts = new ColoredParticlesData;
		// 			contour_pts->name_ = std::string("contour particles");
		// 			contour_pts->point_size_ = 3.0f;
		// 			contour_pts->position_.initialize(contour_.vertex_positions_.num_elements_);
		// 			contour_pts->color_.initialize(contour_.vertex_positions_.num_elements_);
		// 
		// 			int count = 0;
		// 			for (int p = 0; p < contour_.vertex_positions_.num_elements_; p++)
		// 			{
		// 				const TV2 &point = contour_.vertex_positions_[p];
		// 				contour_pts->position_.values_[count] = TV(point.x_, point.y_, 0.0f);
		// 				contour_pts->color_.values_[count] = glm::vec4(1, 0, 0, 0);
		// 				count++;
		// 			}
		// 
		// 			data_depot->colored_particles_list_.pushBack(contour_pts);
		// 		}
		// 		END_ONE_THREAD_WORK;

		static LinesData *lines_data_temp = nullptr;

		BEGIN_ONE_THREAD_WORK
		{
			lines_data_temp = new LinesData;

			data_depot->lines_list_.pushBack(lines_data_temp);

			lines_data_temp->color_ = glm::vec4(0, 0, 0, 0);

			BOX_3D<T> box(grid_.x_min_, grid_.y_min_, 0.0f, grid_.x_max_, grid_.y_max_, 0.0f);

			GL_TOOLS::AddCubeEdges(box, lines_data_temp->vertices_);

			// drawing level set normals
			// 		for(int j = grid_.j_start_; j <= grid_.j_end_; j ++)
			// 			for (int i = grid_.i_start_; i <= grid_.i_end_; i++)
			// 			{
			// 				const TV2 p0 = grid_.GetCellCenter(i, j);
			// 				const TV2 p1 = p0 + levelset_.normal_(i, j) * (T)0.5 * grid_.dx_;
			// 
			// 				lines_data_temp->vertices_.PushBack() = TV3(p0.x_, p0.y_, -0.5f);
			// 				lines_data_temp->vertices_.PushBack() = TV3(p1.x_, p1.y_, -0.5f);
			// 			}

			for (int l = 0; l < contour_.lines_.num_elements_; l++)
			{
				const TV2 &v0 = contour_.vertex_positions_[contour_.lines_[l].v0_];
				const TV2 &v1 = contour_.vertex_positions_[contour_.lines_[l].v1_];

				lines_data_temp->vertices_.PushBack() = TV3(v0.x_, v0.y_, 0.0f);
				lines_data_temp->vertices_.PushBack() = TV3(v1.x_, v1.y_, 0.0f);
			}
		}
		END_ONE_THREAD_WORK;
	}

	// 	static ColoredParticlesData *phi_center_pts(nullptr);
	// 
	// 	BEGIN_ONE_THREAD_WORK
	// 	{
	// 		phi_center_pts = new ColoredParticlesData;
	// 		phi_center_pts->name_ = std::string("flip_particles");
	// 		phi_center_pts->point_size_ = 4.0f;
	// 		phi_center_pts->position_.initialize(levelset_.grid_ghost_.getNumAllCells());
	// 		phi_center_pts->color_.initialize(levelset_.grid_ghost_.getNumAllCells());
	// 
	// 		for (int j = levelset_.grid_ghost_.j_start_; j <= levelset_.grid_ghost_.j_end_; j++)
	// 			for (int i = levelset_.grid_ghost_.i_start_; i <= levelset_.grid_ghost_.i_end_; i++)
	// 			{
	// 				const TV2 pos = levelset_.grid_ghost_.GetCellCenter(i, j);
	// 				const int ix1d = levelset_.phi_.get1DIndex(i, j);
	// 
	// 				phi_center_pts->position_.values_[ix1d] = TV(pos.x_, pos.y_, levelset_.phi_(i, j));
	// 
	// 				if (levelset_.phi_(i, j) > 0.0f) phi_center_pts->color_.values_[ix1d] = glm::vec4(0, 0, 1, 0);
	// 				else phi_center_pts->color_.values_[ix1d] = glm::vec4(1, 0, 0, 0);
	// 			}
	// 
	// 		data_depot->colored_particles_list_.pushBack(phi_center_pts);
	// 	}
	// 	END_ONE_THREAD_WORK;


	if (step_ <= 2 && step_ >= 1)
	{
		static PhongTrianglesData *phong_triangles_temp = nullptr;

		BEGIN_ONE_THREAD_WORK
		{
			phong_triangles_temp = new PhongTrianglesData;

			data_depot->phong_triangles_list_.pushBack(phong_triangles_temp);

			// 		phong_triangles_temp->positions_.Reset();
			// 		phong_triangles_temp->normals_.Reset();

			LevelsetUniform2D &levelset_(levelset_image_);

			for (int j = grid_.j_start_; j < grid_.j_end_; j++)
				for (int i = grid_.i_start_; i < grid_.i_end_; i++)
				{
					phong_triangles_temp->positions_.PushBack() = levelset_.getPhiPos(i, j);
					phong_triangles_temp->positions_.PushBack() = levelset_.getPhiPos(i + 1, j);
					phong_triangles_temp->positions_.PushBack() = levelset_.getPhiPos(i, j + 1);

					phong_triangles_temp->positions_.PushBack() = levelset_.getPhiPos(i + 1, j);
					phong_triangles_temp->positions_.PushBack() = levelset_.getPhiPos(i + 1, j + 1);
					phong_triangles_temp->positions_.PushBack() = levelset_.getPhiPos(i, j + 1);

					TV2 gradient;

					gradient = -levelset_.getGradient(i, j);
					phong_triangles_temp->normals_.PushBack() = TV(gradient.x_, gradient.y_, 1).getSafeNormalized();

					gradient = -levelset_.getGradient(i + 1, j);
					phong_triangles_temp->normals_.PushBack() = TV(gradient.x_, gradient.y_, 1).getSafeNormalized();

					gradient = -levelset_.getGradient(i, j + 1);
					phong_triangles_temp->normals_.PushBack() = TV(gradient.x_, gradient.y_, 1).getSafeNormalized();

					gradient = -levelset_.getGradient(i + 1, j);
					phong_triangles_temp->normals_.PushBack() = TV(gradient.x_, gradient.y_, 1).getSafeNormalized();

					gradient = -levelset_.getGradient(i + 1, j + 1);
					phong_triangles_temp->normals_.PushBack() = TV(gradient.x_, gradient.y_, 1).getSafeNormalized();

					gradient = -levelset_.getGradient(i, j + 1);
					phong_triangles_temp->normals_.PushBack() = TV(gradient.x_, gradient.y_, 1).getSafeNormalized();
				}
		}
		END_ONE_THREAD_WORK;
	}

	static PhongTrianglesData *phong_triangles_temp = nullptr;

	if (step_ >= 3)
	{
		BEGIN_ONE_THREAD_WORK
		{
			phong_triangles_temp = new PhongTrianglesData;

			data_depot->phong_triangles_list_.pushBack(phong_triangles_temp);

			// 		phong_triangles_temp->positions_.Reset();
			// 		phong_triangles_temp->normals_.Reset();

			//	LevelsetUniform2D &levelset_(levelset_image_);

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

			cutter_surface_.copyRenderingData(phong_triangles_temp->positions_, phong_triangles_temp->normals_);
		}
		END_ONE_THREAD_WORK;
	}
}