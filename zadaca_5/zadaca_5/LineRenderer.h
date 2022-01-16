#pragma once

#include <vector>
#include <iostream>
#include "glad.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include "Shader.h"
#include <fstream>

void InitLineRendering();

void RenderLine(std::vector<glm::vec2>& points);