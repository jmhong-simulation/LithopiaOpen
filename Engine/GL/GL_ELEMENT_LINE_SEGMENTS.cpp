#include "GL_ELEMENT_LINE_SEGMENTS.h"
#include "../DataStructure/Vector3D.h"

void GL_ELEMENT_LINE_SEGMENTS::InitializeVBO(const LinkedArray<Vector3D<float> >& positions)
{
	releaseBuffers();

	// Generate vertex buffer objects
	vbo_ids_.initialize(num_vbo_);
	glGenBuffers(vbo_ids_.num_elements_, vbo_ids_.values_);

	int vbo_count = 0;
	bindBufferData(vbo_ids_[vbo_count++], GL_ARRAY_BUFFER, positions);

	num_vertices_ = positions.num_elements_;
}

void GL_ELEMENT_LINE_SEGMENTS::drawArrays()
{
	int vbo_count = 0;

	bindBufferAttribute(vbo_ids_[vbo_count++], "a_position", 3);				// position 3 floats converted to vec 4 in vshader

	glDrawArrays(GL_LINES, 0, num_vertices_);
}
