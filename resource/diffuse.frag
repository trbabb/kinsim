#version 120

uniform vec3  Lcolor1;
uniform vec3  Lcolor2;
uniform vec3  Lpos1;
uniform vec3  Lpos2;
uniform float Ks;
uniform float Kd;
uniform float gloss;
uniform mat4  w2cam;

varying vec3 normal;
varying vec3 color;
varying vec3 pos;


vec3 light(vec3 L, vec3 P, vec3 N, vec3 Cl, vec3 Kd_, float Ks_, float gloss_) {
    L = normalize(L);
    vec3 c = max(dot(L, N), 0.0) * Cl * Kd_;
    if (Ks_ > 0) {
        c += pow(dot(reflect(normalize(P), N), L), gloss_) * Ks_ * Cl;
    }
    return c;
}


void main() {
    vec3 n = normalize(normal);
    vec3 c = light((w2cam * vec4(Lpos1, 0)).xyz, pos, n, Lcolor1, color * Kd, Ks, gloss);
    c     += light((w2cam * vec4(Lpos2, 0)).xyz, pos, n, Lcolor2, color * Kd, Ks, gloss);
    gl_FragColor = vec4(c,1);
}
