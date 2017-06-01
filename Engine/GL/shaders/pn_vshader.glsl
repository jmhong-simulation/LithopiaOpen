#version 120

#ifdef GL_ES
// Set default precision to medium
precision mediump int;
precision mediump float;
#endif

uniform mat4 mvp_matrix;
uniform vec4 color;
uniform float normal_offset;

attribute vec4 a_position;
attribute vec4 a_normal;

varying vec4 vColor;

void main()
{
//	a_position.w = 1.0f;
//	a_normal.w = 0.0f;

	vec4 a_pos = a_position;
	a_pos.w = 1.0f;

	vec4 a_nor = a_normal;
	a_nor.w = 0.0f;


    // Calculate vertex position in screen space
    gl_Position = mvp_matrix * (a_position + a_normal*normal_offset);

    vColor = color;
}
