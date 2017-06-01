// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#pragma once

#include <list>

class DynamicVertex;
class DynamicEdge;

////////////////////////////////////////////////////////////////////
// HOLE
////////////////////////////////////////////////////////////////////
class DynamicHole
{
public:
	std::list <DynamicVertex* > vertices_;
	std::list <DynamicEdge* >	edges_;

public:
	DynamicEdge* AddEdge(DynamicVertex* v0,DynamicVertex* v1);
	DynamicEdge* AddEdge(DynamicEdge* edge);
	void Draw();	
};