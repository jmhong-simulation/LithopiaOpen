// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.

#pragma once

#include "GENERIC_DEFINITIONS.h"

#include <list>
#include <vector>
#include <iostream>

class DynamicVertex;
class DynamicTriangle;

////////////////////////////////////////////////////////////////////
// EDGE
////////////////////////////////////////////////////////////////////
class DynamicEdge
{
public:
	DynamicVertex* vertices_[2];
	TV3 cut_pos_;
	bool is_cutting_;

public:
	DynamicEdge();
	DynamicEdge(DynamicVertex* v0,DynamicVertex* v1);
	~DynamicEdge();

public:
	void	Initlize(DynamicVertex* v0, DynamicVertex* v1, TV3 cut_pus);
	bool	IsSame(DynamicEdge* edge);
	bool	IsEmpty();
	void	Draw();
	T		GetLength();
};