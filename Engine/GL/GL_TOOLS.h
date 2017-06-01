#include "GENERIC_DEFINITIONS.h"
#include "DataStructure/GridUniform3D.h"
#include "DataStructure/LinkedArray.h"
#include "Geometry/BOX_3D.h"

namespace GL_TOOLS
{
	void AddCubeEdges(const BOX_3D<T>& cube, LinkedArray<Vector3D<T> >& v_pos_buffer);
	void AddUniformGridEdges(const GridUniform3D& grid, LinkedArray<Vector3D<T> >& v_pos_buffer);
	void addCircleLines(const TV& center, const T& radius, const int num_segments, LinkedArray<Vector3D<T> >& v_pos_buffer);
}