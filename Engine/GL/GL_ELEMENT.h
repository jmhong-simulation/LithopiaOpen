#pragma once

#include <QGLFunctions>
#include <QGLShaderProgram>
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include "GL_VIEW.h"

#include "DataStructure/Array1D.h"
#include "DataStructure/LinkedArray.h"
#include "DataStructure/Vector3D.h"
#include "Geometry/StaticTriangularSurface.h"
#include "GL/GL_MATERIAL_PARAMETERS.h"
#include "GL/GL_LIGHTING_PARAMETERS.h"

enum VERTEX_TYPE
{
	USE_POSITION	= 0x01,
	USE_NORMAL		= 0x02,
	USE_COLOR		= 0x04,
	USE_TEXTURE		= 0x08
};

class GL_ELEMENT : protected QGLFunctions
{
public:
	std::string name_;

	bool draw_;
	bool update_;

	GLfloat line_width_;
	GLfloat point_size_;
	glm::vec4 color_;

	GLenum primitive_type_;			// GL_POINTS, GL_LINES, GL_TRIANGLES, etc. for glDrawArrays
	GLenum facing_mode_;			// GL_FRONT_AND_BACK, GL_FRONT, GL_BACK. for glPolygonMode
	GLenum rasterize_mode_;			// GL_POINT, GL_LINE, GL_FILL for glPolygonMode

	QGLShaderProgram *program_;
	Array1D<GLuint>	 vbo_ids_;

	int num_vertices_;
	unsigned int vertex_type_;		// USE_POSITION, USE_NORMAL, USE_COLOR, USE_TEXTURE

public:
	GL_ELEMENT();
	~GL_ELEMENT();

	void draw(const glm::mat4& m_matrix);
	void initializeQTGL();
	void setProgram(QGLShaderProgram* program_input);
	void releaseBuffers();

	template<class TT> void bindBufferData(const int vbo_id, const int target, const Array1D<TT>& data_array);
	template<class TT> void bindBufferData(const int vbo_id, const int target, const LinkedArray<TT>& data_array);
	void bindBufferAttribute(const int vbo_id, const char* attribute_name, const int num_components);

	void applyMaterial(const GL_MATERIAL_PARAMETERS& material);
	void applyLighting(const GL_LIGHTING_PARAMETERS& lighting_parameters);
	void applyLight(const int num_light, const GL_LIGHT_PARAMETERS& light_parameters);

	void bindShader(QGLShaderProgram* _program);
	void releaseShader();

	void setUniformMVPMatrix(const glm::mat4& mvp_matrix);
	void setUniformColor(const glm::vec4& color);
	void setUniformNormalOffset(const T& normal_offset);
	
	GLuint getVBOID(const VERTEX_TYPE vf);

public:
	void begin(const GLenum draw_mode, const unsigned int vertex_flag);
	void end();

	template<class TT> void bindArray(const VERTEX_TYPE vf, const Array1D<TT>& data_array);
	template<class TT> void bindArray(const VERTEX_TYPE vf, const LinkedArray<TT>& data_array);

	void drawArrays();
	void drawArrays(const int _num_vertices);	// drawing a part of arrays	
};
