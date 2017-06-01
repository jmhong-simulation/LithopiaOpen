// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "SortedParticlesUniform3D.h"
#include "DataStructure/Array3D.h"
#include "DataStructure/CSRMatrix.h"
#include "DataStructure/GridUniform3D.h"
#include "DataStructure/VectorND.h"
#include "Geometry/MarchingCubesAlgorithm.h"
#include "Parallelism/MultiThreading.h"
#include "Utilities/ScriptReader.h"
#include "Framework/DataDepot.h"
#include "Framework/TaskManager.h"
#include "Image/IMAGE_2D.h"
#include "DataStructure/GridUniform2D.h"
#include "Geometry/RAY.h"
#include "Math/MovingLeastSquares2D.h"

class Photo2_5D : public TaskManager
{
public:
	IMAGE_2D image_;
	GridUniform2D grid_;
	Array2D<T> height_;
	Array2D<int> vertex_flag_;

//	TV3 pointer_;

	TV2_INT selected_index_;

	MovingLeastSquares2D mls_;

public:
	Photo2_5D(DataDepot* data)
		: TaskManager(data), selected_index_(-1, -1), mls_(0, 1e-8)
	{}

	void initializeFromScript(const char *script_filename);
	void updateOneStep(MT* mt, const int thread_id);
	void updateDataDepot(MT* mt, const int thread_id, DataDepot* data_depot);
	void updateDataDepotSingle(DataDepot* data_depot);
	BOX_3D<T> getAABB();

//	void receiveMessage(ScriptParameterList& message_parser);

	void pickVertex(const RAY& ray);

	void smoothHeight(const int& number);
	void smoothHeightParallel(MT* mt, const int thread_id, const int& number);

	void resetMLS();
};