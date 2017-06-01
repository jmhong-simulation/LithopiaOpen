#pragma once

#include <QGLShaderProgram>

class GL_SHADER_MANAGER : protected QGLFunctions
{
public:
	QGLShaderProgram program_pc_;
	QGLShaderProgram program_p_;
	QGLShaderProgram program_phong_;
	QGLShaderProgram program_pn_;		// position is deviated by normal

	void InitializeShaders()
	{
		initializeGLFunctions();

//		const QString shaders_path = "../../Engine/GL/shaders/";	//TODO: options
		const QString shaders_path = "./shaders/";	//TODO: options

		InitializeShader(program_pc_, shaders_path, "pc");
		InitializeShader(program_p_, shaders_path, "p");
		InitializeShader(program_phong_, shaders_path, "phong");
		InitializeShader(program_pn_, shaders_path, "pn");
	}

	void InitializeShader(QGLShaderProgram& program_, const QString& shader_path, const QString& prefix)
	{
		// Compile vertex shader
		if (!program_.addShaderFromSourceFile(QGLShader::Vertex, shader_path + prefix + "_vshader.glsl")) qFatal("Cannot initialize vertex shader.");

		// Compile fragment shader
		if (!program_.addShaderFromSourceFile(QGLShader::Fragment, shader_path + prefix + "_fshader.glsl")) qFatal("Cannot initialize fragment shader.");

		// Link shader pipeline
		if (!program_.link()) qFatal("Cannot link program to shader pipeline.");

		// Bind shader pipeline for use
		if (!program_.bind()) qFatal("Cannot link bind shader with pipeline.");
	}
};