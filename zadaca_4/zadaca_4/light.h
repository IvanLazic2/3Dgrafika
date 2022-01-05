#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <model.h>
#include <shader.h>

using namespace std;

class Light {
public:
	Model *model;
	Shader *shader;
	glm::vec3 position;
	glm::vec3 scale;

	glm::vec3 ambient;
	glm::vec3 diffuse;
	glm::vec3 specular;
	float intensity;

	Light(const glm::vec3& position, 
		  const glm::vec3& scale, 
		  const glm::vec3& ambient, 
		  const glm::vec3& diffuse, 
		  const glm::vec3& specular,
		  const float& intensity)
	{
		model = new Model("objs/sfera.obj");
		shader = new Shader("shaders/light_shader.vs", "shaders/light_shader.fs");
		this->position = position;
		this->scale = scale;
		this->ambient = ambient;
		this->diffuse = diffuse;
		this->specular = specular;
		this->intensity = intensity;
	}

	void Draw(glm::mat4 projection, glm::mat4 view) 
	{
		shader->use();
		shader->setMat4("projection", projection);
		shader->setMat4("view", view);

		glm::mat4 toDraw = glm::mat4(1.0f);
		toDraw = glm::translate(toDraw, position);
		toDraw = glm::scale(toDraw, scale);
		shader->setMat4("model", toDraw);

		model->Draw(*shader);
	}
};