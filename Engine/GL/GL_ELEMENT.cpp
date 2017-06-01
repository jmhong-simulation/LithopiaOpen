#include "GL_ELEMENT.h"

#include <iostream>
#include <fstream>
#include <glm/glm.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/transform.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <ctime>
#include "../DataStructure/LinkedArray.h"
#include "../Geometry/StaticTriangularSurface.h"
#include "../Geometry/OBJ_FILE_READER.h"
#include "GL_WORLD.h"

GL_ELEMENT::GL_ELEMENT()
	: program_(nullptr), draw_(true), update_(true), line_width_(1.0f), point_size_(1.0f),	
	primitive_type_(GL_TRIANGLES), facing_mode_(GL_FRONT), rasterize_mode_(GL_FILL)
{
	initializeQTGL();
}

GL_ELEMENT::~GL_ELEMENT()
{
	releaseBuffers();
}

void GL_ELEMENT::initializeQTGL()
{
	initializeGLFunctions();
}

void GL_ELEMENT::releaseBuffers()
{
	if (vbo_ids_.num_elements_ == 0) return;

	glDeleteBuffers(vbo_ids_.num_elements_, vbo_ids_.values_);	// https://www.opengl.org/sdk/docs/man/html/glDeleteBuffers.xhtml
}

void GL_ELEMENT::setProgram(QGLShaderProgram* program_input)
{
	program_ = program_input;
}

template<class TT>
void GL_ELEMENT::bindBufferData(const int vbo_id, const int target, const Array1D<TT>& data_array)
{
	// target is GL_ARRAY_BUFFER or GL_ELEMENT_ARRAY_BUFFER
	glBindBuffer(target, vbo_id);

	const int temp = data_array.getSizeOfData();

	glBufferData(target, data_array.getSizeOfData(), data_array.values_, GL_STATIC_DRAW);
}

template<class TT>
void GL_ELEMENT::bindBufferData(const int vbo_id, const int target, const LinkedArray<TT>& data_array)
{
	// target is GL_ARRAY_BUFFER or GL_ELEMENT_ARRAY_BUFFER
	glBindBuffer(target, vbo_id);
	glBufferData(target, data_array.GetSizeOfData(), NULL, GL_STATIC_DRAW);
	for (LinkedArray<TT>::T_BLOCK *itr_block = data_array.head_; itr_block != nullptr; itr_block = itr_block->next_)
		glBufferSubData(target, itr_block->i_start_*itr_block->GetSizeOfType(), itr_block->GetSizeOfData(), itr_block->values_);
}

void GL_ELEMENT::bindBufferAttribute(const int vbo_id, const char* attribute_name, const int num_components)
{
	glBindBuffer(GL_ARRAY_BUFFER, vbo_id);
	const int location = program_->attributeLocation(attribute_name);
	program_->enableAttributeArray(location);
	glVertexAttribPointer(location, num_components, GL_FLOAT, GL_FALSE, 0, 0);		// components floats per vertex (x, y, z, w), (r, g, b, a), etc.
}

void GL_ELEMENT::applyMaterial(const GL_MATERIAL_PARAMETERS& material)
{
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, &material.Ka_[0]);
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, &material.Kd_[0]);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, &material.Ks_[0]);
	glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, &material.Ke_[0]);
	glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, material.Se_);
}

void GL_ELEMENT::applyLight(const int num_light, const GL_LIGHT_PARAMETERS& param)
{
	glLightfv(GL_LIGHT0, GL_SPOT_DIRECTION, &param.spot_direction_[0]);
	glLighti(GL_LIGHT0, GL_SPOT_EXPONENT, param.spot_exponent_);
	glLighti(GL_LIGHT0, GL_SPOT_CUTOFF, param.spot_cutoff_);

	glLightf(GL_LIGHT0, GL_CONSTANT_ATTENUATION, param.Kc_);
	glLightf(GL_LIGHT0, GL_LINEAR_ATTENUATION, param.Kl_);
	glLightf(GL_LIGHT0, GL_QUADRATIC_ATTENUATION, param.Kq_);

	glLightfv(GL_LIGHT0, GL_POSITION, &param.light_pos_[0]);
	glLightfv(GL_LIGHT0, GL_AMBIENT, &param.light_Ka_[0]);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, &param.light_Kd_[0]);
	glLightfv(GL_LIGHT0, GL_SPECULAR, &param.light_Ks_[0]);
}

void GL_ELEMENT::applyLighting(const GL_LIGHTING_PARAMETERS& param)
{
	if (param.use_lighting_) glEnable(GL_LIGHTING);
	if (param.use_light_[0]) glEnable(GL_LIGHT0);
	if (param.normalize_normals_) glEnable(GL_NORMALIZE);

	// Light model parameters:
	// -------------------------------------------

	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, &param.lmKa_[0]);
	glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, 1.0);
	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, 0.0);

	applyLight(GL_LIGHT0, param.light_params_[0]);
}

void GL_ELEMENT::bindShader(QGLShaderProgram* _program)
{
	if (program_ != nullptr) std::cout << "Program was not released after previous use" << std::endl;

	program_ = _program;

	program_->bind();
}

void GL_ELEMENT::releaseShader()
{
	if (program_ == nullptr) std::cout << "Program was not binded before the use" << std::endl;

	program_->release();

	program_ = nullptr;
}

void GL_ELEMENT::setUniformMVPMatrix(const glm::mat4& mvp_matrix)
{
	glUniformMatrix4fv(program_->uniformLocation("mvp_matrix"), 1, GL_FALSE, &mvp_matrix[0][0]);
}

void GL_ELEMENT::setUniformColor(const glm::vec4& color)
{
	glUniform4fv(program_->uniformLocation("color"), 1, &color[0]);
}

void GL_ELEMENT::setUniformNormalOffset(const T& normal_offset)
{
	glUniform1f(program_->uniformLocation("normal_offset"), normal_offset);
}

GLuint GL_ELEMENT::getVBOID(const VERTEX_TYPE vf)
{
	int num_vbo = 0;
 
	if (vf == USE_POSITION) return vbo_ids_[num_vbo];
	if (vertex_type_ & USE_POSITION)	++num_vbo;

	if (vf == USE_COLOR) return vbo_ids_[num_vbo];
	if (vertex_type_ & USE_COLOR)	++num_vbo;

	if (vf == USE_NORMAL) return vbo_ids_[num_vbo];
	if (vertex_type_ & USE_NORMAL)	++num_vbo;

	if (vf == USE_TEXTURE) return vbo_ids_[num_vbo];
	if (vertex_type_ & USE_TEXTURE)	++num_vbo;

	std::cout << "Incorrect VERTEX_FLAG in GL_ELEMENT::GetVBOID" << std::endl;

	exit(1);

	return vbo_ids_[num_vbo];
}

void GL_ELEMENT::begin(const GLenum draw_mode, const unsigned int vertex_flag)
{
	releaseBuffers();

	primitive_type_ = draw_mode;

	vertex_type_ = vertex_flag;

	int num_vbo = 0;

	if (vertex_type_ & USE_POSITION)	++num_vbo;
	if (vertex_type_ & USE_COLOR)	++num_vbo;
	if (vertex_type_ & USE_NORMAL)	++num_vbo;
	if (vertex_type_ & USE_TEXTURE)	++num_vbo;

	vbo_ids_.initialize(num_vbo);
	
	glGenBuffers(vbo_ids_.num_elements_, vbo_ids_.values_);
}

template<class TT> void GL_ELEMENT::bindArray(const VERTEX_TYPE vf, const Array1D<TT>& data_array)
{
	bindBufferData(getVBOID(vf), GL_ARRAY_BUFFER, data_array);

	num_vertices_ = data_array.num_elements_;
}

template<class TT> void GL_ELEMENT::bindArray(const VERTEX_TYPE vf, const LinkedArray<TT>& data_array)
{
	bindBufferData(getVBOID(vf), GL_ARRAY_BUFFER, data_array);

	num_vertices_ = data_array.num_elements_;
}

void GL_ELEMENT::end()
{
}

void GL_ELEMENT::drawArrays()
{
	drawArrays(num_vertices_);
}

void GL_ELEMENT::drawArrays(const int _num_vertices)
{
	int vbo_count = 0;

	if (vertex_type_ & USE_POSITION) bindBufferAttribute(vbo_ids_[vbo_count++], "a_position", 3);			// position 3 floats converted to vec 4 in vshader
	if (vertex_type_ & USE_COLOR)	bindBufferAttribute(vbo_ids_[vbo_count++], "a_color", 4);					// color has 4 floats
	if (vertex_type_ & USE_NORMAL)	bindBufferAttribute(vbo_ids_[vbo_count++], "a_normal", 3);				// normals has 3 floats converted to vec 4 in vshader

	glDrawArrays(primitive_type_, 0, _num_vertices);
}

void GL_ELEMENT::draw(const glm::mat4& m_matrix)
{
	glPointSize(point_size_);
	glLineWidth(line_width_);
	glPolygonMode(facing_mode_, rasterize_mode_);

	if (primitive_type_ == GL_TRIANGLES && rasterize_mode_ == GL_FILL)
		bindShader(&GL_WORLD::GetInstance().shader_programs_.program_phong_);
	else if (vertex_type_ & USE_COLOR)
		bindShader(&GL_WORLD::GetInstance().shader_programs_.program_pc_);
	else
		bindShader(&GL_WORLD::GetInstance().shader_programs_.program_p_);

//	applyMaterial(surface_material_);
	setUniformMVPMatrix(m_matrix);
	setUniformColor(color_);
	drawArrays();
	releaseShader();
}

template void GL_ELEMENT::bindBufferData(const int vbo_ix, const int target, const Array1D<float>& data_array);
template void GL_ELEMENT::bindBufferData(const int vbo_ix, const int target, const Array1D<glm::vec4>& data_array);
template void GL_ELEMENT::bindBufferData(const int vbo_ix, const int target, const Array1D<Vector3D<float> >& data_array);
template void GL_ELEMENT::bindBufferData(const int vbo_ix, const int target, const Array1D<Vector3D<int> >& data_array);
template void GL_ELEMENT::bindBufferData(const int vbo_ix, const int target, const Array1D<Vector3D<unsigned short> >& data_array);
template void GL_ELEMENT::bindBufferData(const int vbo_ix, const int target, const LinkedArray<Vector3D<float> >& data_array);
template void GL_ELEMENT::bindBufferData(const int vbo_ix, const int target, const LinkedArray<Vector3D<int> >& data_array);
template void GL_ELEMENT::bindBufferData(const int vbo_ix, const int target, const LinkedArray<Vector3D<unsigned short> >& data_array);

template void GL_ELEMENT::bindArray(const VERTEX_TYPE vf, const Array1D<float>& data_array);
template void GL_ELEMENT::bindArray(const VERTEX_TYPE vf, const Array1D<glm::vec4>& data_array);
template void GL_ELEMENT::bindArray(const VERTEX_TYPE vf, const Array1D<Vector3D<float> >& data_array);
template void GL_ELEMENT::bindArray(const VERTEX_TYPE vf, const Array1D<Vector3D<int> >& data_array);
template void GL_ELEMENT::bindArray(const VERTEX_TYPE vf, const Array1D<Vector3D<unsigned short> >& data_array);
template void GL_ELEMENT::bindArray(const VERTEX_TYPE vf, const LinkedArray<Vector3D<float> >& data_array);
template void GL_ELEMENT::bindArray(const VERTEX_TYPE vf, const LinkedArray<Vector3D<int> >& data_array);
template void GL_ELEMENT::bindArray(const VERTEX_TYPE vf, const LinkedArray<Vector3D<unsigned short> >& data_array);
template void GL_ELEMENT::bindArray(const VERTEX_TYPE vf, const LinkedArray<glm::vec4>& data_array);