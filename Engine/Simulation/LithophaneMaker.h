// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include <iostream>
#include "Framework/TaskManager.h"
#include "Geometry/StaticTriangularSurface.h"
#include "Geometry/SHELL_OF_REVOLUTION.h"
#include "Utilities/ScriptParameterList.h"

class LithophaneMaker : public TaskManager	 // currently levelset normal flow test
{
public:
	StaticTriangularSurface surface_, outer_base_surface_;

	LithophaneMaker(DataDepot* data_);

	std::string output_path_, output_filename_prefix_;

//	void initialize();
	void initializeFromScript(const std::string script_filename);
	void updateOneStep(MT* mt, const int thread_id){};
	void updateDataDepot(MT* mt, const int thread_id, DataDepot* data_depot);

	BOX_3D<T> getAABB();

	void writeFiles();
};