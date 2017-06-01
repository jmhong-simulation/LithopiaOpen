#version 120

#ifdef GL_ES
// Set default precision to medium
precision mediump int;
precision mediump float;
#endif

uniform mat4 mvp_matrix;

attribute vec4 a_position;
attribute vec4 a_normal;

varying vec4 N;
varying vec4 v;

void main(void)
{
//	a_position.w = 1.0f;		// position
//	a_normal.w = 0.0f;			// vector

	vec4 a_pos = a_position;
	a_pos.w = 1.0f;

	vec4 a_nor = a_normal;
	a_nor.w = 0.0f;

	v = mvp_matrix * a_pos;
	N = normalize(mvp_matrix * a_nor);

	gl_Position = v;

// v = vec3(gl_ModelViewMatrix * gl_Vertex);       
// N = normalize(gl_NormalMatrix * gl_Normal);

// gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
}

// http://www.opengl.org/sdk/docs/tutorials/ClockworkCoders/lighting.php