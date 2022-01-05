#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <model.h>
#include <shader.h>

#include <vector>
#include <string>

using namespace std;

class Object {
public:
	Model *model;
	Shader *shader;
	glm::vec3 position;
	glm::vec3 scale;

	Object(const char* path, const char* vertex_shader_path, const char* fragment_shader_path, glm::vec3 position, glm::vec3 scale)
	{
		//model = md;
		model = new Model(path);
		shader = new Shader(vertex_shader_path, fragment_shader_path);
		//shader = sh;
		shader->use();
		shader->setInt("material.diffuse", 0);
		shader->setInt("material.specular", 1);

		this->position = position;
		this->scale = scale;
	}

	void Draw(const Camera& camera, const glm::mat4& projection, const glm::mat4& view, const vector<Light>& lights, const float& currentFrame = 0.0, const float& speed = 0.0)
	{
		shader->use();
		shader->setVec3("viewPos", camera.Position);
		shader->setFloat("material.shininess", 32.0f);

		for (int i = 0; i < LIGHT_NUM; ++i) {
			string str = "lights[" + to_string(i) + "].";

			shader->setVec3(str + "position", lights[i].position);
			shader->setVec3(str + "ambient", lights[i].ambient);
			shader->setVec3(str + "diffuse", lights[i].diffuse);
			shader->setVec3(str + "specular", lights[i].specular);
			shader->setFloat(str + "intensity", lights[i].intensity);
		}

		shader->setMat4("projection", projection);
		shader->setMat4("view", view);

		glm::mat4 toDraw = glm::mat4(1.0f);
		toDraw = glm::translate(toDraw, position);
		toDraw = glm::scale(toDraw, scale);

		toDraw = glm::rotate(toDraw, (glm::f32)glm::radians(currentFrame * speed), glm::vec3(0, 1, 0));

		shader->setMat4("model", toDraw);

		model->Draw(*shader);
	}
};