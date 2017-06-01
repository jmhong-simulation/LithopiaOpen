// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "DataStructure/GridUniform3D.h"
#include "DataStructure/VectorND.h"
#include "Parallelism/MultiThreading.h"
#include "Utilities/ScriptReader.h"
#include "Framework/DataDepot.h"
#include "Framework/TaskManager.h"
#include "Image/IMAGE_2D.h"
#include "DataStructure/GridUniform2D.h"

class SoundTo3D : public TaskManager
{
public:
	IMAGE_2D image_;
	GridUniform2D grid_;
	Array2D<T> height_;

	Array1D<T> audio_spectrum_;

public:
	SoundTo3D(DataDepot* data)
		: TaskManager(data)
	{}

	~SoundTo3D();

	void initializeFromScript(const char *script_filename);
	void updateOneStep(MT* mt, const int thread_id);
	void updateDataDepot(MT* mt, const int thread_id, DataDepot* data_depot);
	void updateDataDepotSingle(DataDepot* data_depot);
	BOX_3D<T> getAABB();

	int current_j_;

	void updateHeight(const Array1D<T>& audio_spectrum);

//	void receiveMessage(ScriptParameterList& message_parser);
};