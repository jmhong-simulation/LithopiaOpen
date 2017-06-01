#pragma once

#include "GL_ELEMENT.h"

class GL_ELEMENT_LINE_SEGMENTS : public GL_ELEMENT
{
public:
	const int num_vbo_;

	int num_vertices_;

	void InitializeGeometry();
	void InitializeVBO(const LinkedArray<Vector3D<float> >& positions);
	void drawArrays();

	GL_ELEMENT_LINE_SEGMENTS() : num_vbo_(1)
	{}
};