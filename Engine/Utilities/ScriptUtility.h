// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "ScriptReader.h"
#include "Geometry/BOX_3D.h"
#include "DataStructure/GridUniform3D.h"

namespace ScriptUtility
{
	void initialize(ScriptBlock& sb, GridUniform3D& grid);

	BOX_3D<float> initializeBox3D(ScriptBlock& sb);
}