#version 150 core
in vec2 position;
in float trans;
uniform mat4 MVP;

out vec4 v_color;

void main(){
    float f = log2(trans * 1.0 + 1.0);
    v_color = vec4(1.5 * f, 
        1.5 * pow(f, 3), 
        1.5 * pow(f, 6), 
        1.0);
    gl_Position = MVP * vec4(position, 0.0, 1.0);
}