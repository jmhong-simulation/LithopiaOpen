#pragma once

#include "GL_ELEMENT.h"

class GL_ELEMENT_COLOR_POINTS : public GL_ELEMENT
{
public:
	int num_vertices_;

	void InitializeVBO(LinkedArray<TV>& vertices, LinkedArray<TV>& colors);
	void drawArrays();
};