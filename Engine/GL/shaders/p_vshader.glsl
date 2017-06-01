#version 120

#ifdef GL_ES
// Set default precision to medium
precision mediump int;
precision mediump float;
#endif

uniform mat4 mvp_matrix;
uniform vec4 color;

attribute vec4 a_position;

varying vec4 vColor;

void main()
{
//	a_position.w = 1.0f;
	vec4 pos = a_position;
	pos.w = 1.0f;

    // Calculate vertex position in screen space
    gl_Position = mvp_matrix * pos;

    vColor = color;
}
