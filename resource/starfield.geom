#version 120
#extension GL_EXT_geometry_shader4 : enable

uniform mat4 rotmat0; // t=0 rotation
uniform mat4 rotmat1; // t=1 rotation
uniform vec2 res;

varying in  vec3 g_intensity[1];
varying out vec3 f_intensity;

void main(void) {
    vec4 p0 = gl_ProjectionMatrix * rotmat0 * gl_PositionIn[0];
    vec4 p1 = gl_ProjectionMatrix * rotmat1 * gl_PositionIn[0];
    
    float l = max(length((p1 / p1.w - p0 / p0.w) * vec4(res * 0.5, 0, 0)) * 0.1, 1);
    
    gl_Position = p0;
    f_intensity = g_intensity[0] / vec3(l);
    EmitVertex();
    
    gl_Position = p1;
    f_intensity = g_intensity[0] / vec3(l);
    EmitVertex();
    
    EndPrimitive();
}
