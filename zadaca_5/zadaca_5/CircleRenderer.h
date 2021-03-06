#pragma once

#include <iostream>
#include "glad.h"
#include <glm/glm.hpp>
#include "Shader.h"
#include <fstream>
#include <sstream>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

void InitCircleRendering(int in_vertexCountPerCicle);

void RenderCircle(glm::vec2 position, float scale);
