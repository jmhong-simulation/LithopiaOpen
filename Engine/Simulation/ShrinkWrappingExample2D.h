// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "Framework/TaskManager.h"
#include "Geometry/LevelsetUniform2D.h"
#include "Geometry/StaticTriangularSurface.h"
#include "Image/IMAGE_2D.h"
#include "Geometry/DynamicContour2D.h"

class ShrinkWrappingExample2D : public TaskManager
{
public:
	int step_;
	int shrink_wrapping_steps_;

	ShrinkWrappingExample2D(DataDepot* data_);

	Array2D<T>	height_map_;
	DynamicContour2D contour_;

	T cutter_padding_;

	int extend_smoothing_;
	int litho_smoothing_;

	float stamp_extend_;
	float cutter_thickness_;
	float stamp_cutter_spacing_;

	float cutter_base_thickness_;
	float cutter_print_thickness_;

	float stamp_base_thickness_;
	float stamp_print_thickness_;

	std::string output_path_;
	std::string output_filename_prefix_;

	StaticTriangularSurface cutter_surface_, stamp_surface_;

	GridUniform2D grid_;
	LevelsetUniform2D levelset_;
	LevelsetUniform2D levelset_image_;	// for shrink wrapping potential field

	LevelsetUniform2D inner_line_;

	IMAGE_2D input_stamp_image_;
	IMAGE_2D input_cutter_image_;
	int extend_;

	void initializeFromScript(const char *script_filename);

	void updateOneStep(MT* mt, const int thread_id);

	void updateDataDepot(MT* mt, const int thread_id, DataDepot* data_depot);

	void findOutline(MT* mt, const int thread_id);
	void findOutlineSubstep(MT* mt, const int thread_id, const T max_dx);

	void findInnerline(MT* mt, const int thread_id);

	void makeStampPart(MT* mt, const int thread_id);

	void makeCutterPart(MT* mt, const int thread_id);

	void prepareForShrinkWrapping();

	BOX_3D<T> getAABB();
};