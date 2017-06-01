#version 120

#ifdef GL_ES
// Set default precision to medium
precision mediump int;
precision mediump float;
#endif

uniform mat4 mvp_matrix;

attribute vec4 a_position;
attribute vec4 a_color;

varying vec4 vColor;

void main()
{
    gl_Position = mvp_matrix * a_position;

    vColor = a_color;
}
