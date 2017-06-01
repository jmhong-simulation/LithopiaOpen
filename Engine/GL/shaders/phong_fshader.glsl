#version 120

#ifdef GL_ES
// Set default precision to medium
precision mediump int;
precision mediump float;
#endif

varying vec4 N;
varying vec4 v;    

void main (void)  
{  
   vec3 L = (normalize(gl_LightSource[0].position.xyz - v.xyz)).xyz;   
   vec3 E = (normalize(-v.xyz)).xyz;		// we are in Eye Coordinates, so EyePos is (0,0,0)  
   vec3 R = (normalize(-reflect(L.xyz,N.xyz))).xyz;  
 
   //calculate Ambient Term:  
   vec4 Iamb = gl_FrontLightProduct[0].ambient;    

   //calculate Diffuse Term:  
   vec4 Idiff = gl_FrontLightProduct[0].diffuse * max(dot(N.xyz,L.xyz), 0.0);    
   
   // calculate Specular Term:
   vec4 Ispec = gl_FrontLightProduct[0].specular 
                * pow(max(dot(R,E),0.0),0.3*gl_FrontMaterial.shininess);

   // write Total Color:  
   gl_FragColor = gl_FrontLightModelProduct.sceneColor + Iamb + Idiff + Ispec;   
}