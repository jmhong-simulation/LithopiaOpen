// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "CONVENTIONAL_MACROS.h"
#include "DynamicHole.h"
#include "DynamicEdge.h"

DynamicEdge* DynamicHole::AddEdge(DynamicVertex*  v0, DynamicVertex*  v1)
{
	DynamicEdge* edge = new DynamicEdge(v0, v1);
	DynamicHole::edges_.push_back(edge);
	return edge;
}

DynamicEdge* DynamicHole::AddEdge(DynamicEdge*  edge)
{
	DynamicHole::edges_.push_back(edge);
	return edge;
}

/*
void DynamicHole::Draw()
{
	TRAVERSE_EDGES
	{
		(*itr_edge)->Draw();	
	}
}
*/