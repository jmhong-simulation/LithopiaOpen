// The Lithopia project initiated by Jeong-Mo Hong for 3D Printing Community.
// Copyright (c) 2015 Jeong-Mo Hong - All Rights Reserved. 
// This file is subject to the terms and conditions defined in 
// file 'LICENSE.txt', which is part of this source code package.
#include "CONVENTIONAL_MACROS.h"
#include "DynamicEdge.h"
#include "DynamicTriangle.h"
#include "DynamicVertex.h"

DynamicEdge::DynamicEdge()
{
	vertices_[0] = NULL;
	vertices_[1] = NULL;
	is_cutting_ = false;
}

DynamicEdge::DynamicEdge(DynamicVertex* v0, DynamicVertex* v1)
{
	Initlize(v0, v1, TV3());
}

DynamicEdge::~DynamicEdge(){}

void DynamicEdge::Initlize(DynamicVertex* v0, DynamicVertex* v1, TV3 cut_pos)
{
	vertices_[0] = v0;
	vertices_[1] = v1;
	cut_pos_ = cut_pos;
	
	if (cut_pos != TV3())
	{
		is_cutting_ = true;
	}
}

bool DynamicEdge::IsSame(DynamicEdge* edge)
{
	DynamicVertex *v00,*v01,*v10,*v11;
	v00 = vertices_[0];
	v01 = vertices_[1];
	v10 = edge->vertices_[0];
	v11 = edge->vertices_[1];
	
	if((v00==v10) && (v01==v11))return true;
	else if((v00==v11) && (v01==v10))return true;
	else return false;
}

bool DynamicEdge::IsEmpty()
{
	if (vertices_[0] == NULL || vertices_[1] == NULL)
		return true;
	else
		return false;
}

T DynamicEdge::GetLength()
{
	if (IsEmpty()) return 0;

	TV deviation = vertices_[0]->x_ - vertices_[1]->x_;

	return deviation.getSqrMagnitude();
}

/*
void DynamicEdge::Draw()
{
	glBegin(GL_LINES);
#ifdef USE_FLOAT_T		
		glVertex3fv(vertices_[0]->GetPosition());
		glVertex3fv(vertices_[1]->GetPosition());
#else
		glVertex3dv(vertices_[0]->GetPosition());
		glVertex3dv(vertices_[1]->GetPosition());
#endif			
	glEnd();
}
*/