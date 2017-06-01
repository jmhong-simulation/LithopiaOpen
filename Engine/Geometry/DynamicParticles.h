// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include "GENERIC_DEFINITIONS.h"
#include "DataStructure/Array1D.h"
#include "DataStructure/LinkedArray.h"

class DynamicParticles
{
public:
	Array1D<TV> pos_, vel_;
};