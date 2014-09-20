#version 120

varying vec3 f_intensity;

// TODO: darken stars to prevent brightening with motion.

void main() {
    gl_FragColor = vec4(f_intensity, 1.0);
}
