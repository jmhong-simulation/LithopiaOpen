#include "GL_ELEMENT_COLOR_POINTS.h"
#include "assert.h"

void GL_ELEMENT_COLOR_POINTS::InitializeVBO(LinkedArray<TV>& vertices, LinkedArray<TV>& colors)
{
	releaseBuffers();

	// Generate vertex buffer objects
	const int num_vbo = 2;
	vbo_ids_.initialize(num_vbo);
	glGenBuffers(vbo_ids_.num_elements_, vbo_ids_.values_);

	int vbo_count = 0;

	bindBufferData(vbo_ids_[vbo_count++], GL_ARRAY_BUFFER, vertices);
	bindBufferData(vbo_ids_[vbo_count++], GL_ARRAY_BUFFER, colors);

	num_vertices_ = vertices.num_elements_;
}

void GL_ELEMENT_COLOR_POINTS::drawArrays()
{
	int vbo_count = 0;

	bindBufferAttribute(vbo_ids_[vbo_count++], "a_position", 3);			// position 3 floats converted to vec 4 in vshader
	bindBufferAttribute(vbo_ids_[vbo_count++], "a_color", 3);				// normals has 3 floats converted to vec 4 in vshader

	glDrawArrays(GL_POINTS, 0, num_vertices_);
}

// NOTE: we don't use GL_ELEMENT_ARRAY_BUFFER and glDrawElements due to it's limited number of triangles.
// type must be GL_UNSIGNED_BYTE or GL_UNSIGNED_SHORT.
// glDrawElements(GL_TRIANGLES, surface_.v_ix_of_triangles_.num_elements_ * 3, GL_UNSIGNED_SHORT, 0);