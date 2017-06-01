// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "Framework/TaskManager.h"
#include "Geometry/LevelsetUniform2D.h"

class LevelsetEditor2D : public TaskManager	 // currently levelset normal flow test
{
public:
	GridUniform2D grid_;
	LevelsetUniform2D levelset_;

	LevelsetEditor2D(DataDepot* data_);

	void initialize();

	void updateOneStep(MT* mt, const int thread_id);

	void updateDataDepot(MT* mt, const int thread_id, DataDepot* data_depot);
};