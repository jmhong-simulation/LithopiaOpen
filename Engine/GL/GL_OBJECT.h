#pragma once

#include <list>
#include "Framework/OBJECT.h"
#include "Geometry/BOX_3D.h"
#include "GL/GL_ELEMENT_LINE_SEGMENTS.h"
#include "GL/GL_ELEMENT_COLOR_POINTS.h"
#include "GL/GL_VIEW.h"

class GL_OBJECT : protected QGLFunctions
{
public:
	std::atomic<bool> is_locked_;

	Object *object_;

	bool draw_surface_;
	bool draw_edges_;
	bool draw_vertices_;
	bool draw_short_edges_;
	bool draw_high_curvature_vertices_;
	bool draw_surface_oobb_;
	bool draw_slicer_;

	GLenum polygon_mode_face_;
	GLenum polygon_mode_mode_;

	GL_LIGHTING_PARAMETERS lighting_parameters_;
	GL_MATERIAL_PARAMETERS surface_material_;

	GL_ELEMENT gl_triangles_;
	GL_ELEMENT gl_surface_oobb_;
	GL_ELEMENT gl_short_edge_triangles_;
	GL_ELEMENT gl_slicer_;
	GL_ELEMENT gl_phi_layers_;
	GL_ELEMENT_COLOR_POINTS gl_high_curvature_vertices_;
//	GL_ELEMENT_LINE_SEGMENTS gl_short_edges_;

	SinglyLinkedList<GL_ELEMENT*> gl_element_list_;

public:
//	GL_OBJECT();
	GL_OBJECT(Object* object_input);
	~GL_OBJECT();

	void Draw(const GL_VIEW& camera);
	
	void DrawSurface(const GL_VIEW& camera);
	void DrawEdges(const GL_VIEW& camera);
	void DrawSlicer(const GL_VIEW& camera);
	void DrawShortEdgeTriangles(const GL_VIEW& camera);
	void DrawHighCurvatureVertices(const GL_VIEW& camera);
	void drawSimAABB(const GL_VIEW& camera);
	void drawSimObjectAABB(const GL_VIEW& camera);

	void UpdateSurface();
	void UpdateSurface(StaticTriangularSurface& surface_);
	void UpdateAMR();
	void UpdateShortEdges(const T& min_length);
	void UpdateHighCurvatureVertices(const T& kappa);
	void UpdateSlicer();
	void updateSimAABB(const BOX_3D<T>& aabb);
	void updateSimObjectAABB(const Array1D<BOX_3D<T> >& aabb_list);

	void update();

	void resetElementList()
	{
		for (gl_element_list_.begin(); gl_element_list_.valid();)
		{
			GL_ELEMENT *temp = gl_element_list_.getItr();
			gl_element_list_.next();
			SAFE_DELETE(temp);		// is temp necessary?	SAFE_DELETE(gl_element_list_.getItr()); gl_element_list_.next(); would work fine!
		}

		gl_element_list_.reset();
	}
};
