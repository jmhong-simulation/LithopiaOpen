#version 120

#ifdef GL_ES
// Set default precision to medium
precision mediump int;
precision mediump float;
#endif

varying vec4 vColor;

void main()
{
    gl_FragColor = vColor;
}
