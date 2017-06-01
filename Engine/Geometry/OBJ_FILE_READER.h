// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include <glm/glm.hpp>
#include "DataStructure/Array1D.h"
#include "DataStructure/LinkedArray.h"
#include "DataStructure/Vector3D.h"
#include "GENERIC_DEFINITIONS.h"

class OBJ_FILE_READER
{
public:
	LinkedArray<TV> pos_stack_;
	LinkedArray<TV2> uv_stack_;
	LinkedArray<TV3> normal_stack_;

	LinkedArray<TV_INT> ix_stack_;
	LinkedArray<TV_INT> uv_ix_stack_;
	LinkedArray<TV_INT> nor_ix_stack_;

	float x_min_, x_max_, y_min_, y_max_, z_min_, z_max_;	// OOBB

	bool use_cout;
	bool cout_time_;

	OBJ_FILE_READER()
		: use_cout(false), cout_time_(false)
	{}

	void ReadOBJ(const char *filename);
	void WriteOBJ(const char* filename);

	const glm::vec3 GetCenterAABB() const;
	const float GetScaleAABB() const;
	const glm::vec3 GetScaleVecAABB() const;

	void GetPositionArray(Array1D<TV3>& pos_arr);
	void GetUVCoordArray(Array1D<TV2>& uv_arr, Array1D<TV3_INT>& tri_arr);
	void GetNormalArray(Array1D<TV3>& normal_arr);
	void GetTriangleArray(Array1D<TV3_INT>& tri_arr);
};
