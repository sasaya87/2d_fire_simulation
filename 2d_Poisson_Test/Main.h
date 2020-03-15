#pragma once

#include <iostream>
#include <gl/glew.h>
#include <GLFW/glfw3.h>
#include <string>
#include <fstream>
#include <sstream>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <chrono>
#include "Object.h"

using glm::vec3;
using glm::vec4;
using glm::mat4;

GLint makeShader(std::string vertexFileName, std::string fragmentFileName);
Object createTriangle();
void createPositions();

void init();
void pressureByCG();
void velo(bool isU);
void omega();
void kawamuraByBiCGSTAB(int mode);

