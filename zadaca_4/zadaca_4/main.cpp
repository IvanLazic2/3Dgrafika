#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <shader.h>
#include <camera.h>
#include <mesh.h>
#include <model.h>

const unsigned int LIGHT_NUM = 3;

#include <light.h>
#include <object.h>

#include <iostream>
#include <vector>

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void mouse_callback(GLFWwindow* window, double xpos, double ypos);
void scroll_callback(GLFWwindow* window, double xoffset, double yoffset);
void processInput(GLFWwindow* window);

const unsigned int WIDTH = 1024;
const unsigned int HEIGHT = 720;

float RotationSpeed = 1.0;

Camera camera(glm::vec3(0.0f, 0.0f, 5.0f));
float lastX = WIDTH / 2.0f;
float lastY = HEIGHT / 2.0f;
bool firstMouse = true;
float deltaTime = 0.0f;
float lastFrame = 0.0f;

vector<Light> Lights;
vector<Object> Objects;

void initLights() {
    Light light1(glm::vec3(0.0, -2.0, 0.0), 
                 glm::vec3(0.2, 0.2, 0.2),
                 glm::vec3(0.05f, 0.05f, 0.05f),
                 glm::vec3(0.8f, 0.8f, 0.8f),
                 glm::vec3(1.0f, 1.0f, 1.0f),
                 2.0f);

    Light light2(glm::vec3(3.0, 3.0, 3.0),
                 glm::vec3(0.5, 0.5, 0.5),
                 glm::vec3(0.05f, 0.05f, 0.05f),
                 glm::vec3(0.8f, 0.8f, 0.8f),
                 glm::vec3(1.0f, 1.0f, 1.0f),
                 4.0f);

    Light light3(glm::vec3(-4.0, 3.0, -4.0),
                 glm::vec3(0.1, 0.1, 0.1),
                 glm::vec3(0.05f, 0.05f, 0.05f),
                 glm::vec3(0.8f, 0.8f, 0.8f),
                 glm::vec3(1.0f, 1.0f, 1.0f),
                 2.0f);

    Lights.push_back(light1);
    Lights.push_back(light2);
    Lights.push_back(light3);
}

void initObjects() {
    Object obj1("objs/piramida/piramida.obj",
                "shaders/obj_shader.vs", "shaders/obj_shader.fs",
                glm::vec3(0.0f, 0.0f, 0.0f),
                glm::vec3(2.0f, 2.0f, 2.0f));

    Object obj2("objs/piramida/piramida.obj",
                "shaders/obj_shader.vs", "shaders/obj_shader.fs",
                glm::vec3(-3.0f, 0.0f, 0.0f),
                glm::vec3(1.0f, 1.0f, 1.0f));

    Object obj3("objs/cube/cube.obj",
                "shaders/obj_shader.vs", "shaders/obj_shader.fs",
                glm::vec3(2.0f, 0.0f, 2.0f),
                glm::vec3(0.3f, 0.3f, 0.3f));


    Object obj4("objs/cube/cube.obj",
                "shaders/obj_shader.vs", "shaders/obj_shader.fs",
                glm::vec3(4.0f, 0.0f, 4.0f),
                glm::vec3(0.75f, 0.75f, 0.75f));


    Object obj5("objs/djiraf/djiraf.obj",
                "shaders/obj_shader.vs", "shaders/obj_shader.fs",
                glm::vec3(6.0f, 0.0f, 6.0f),
                glm::vec3(0.05f, 0.05f, 0.05f));

    Object obj6("objs/car/car.obj",
                "shaders/obj_shader.vs", "shaders/obj_shader.fs",
                glm::vec3(-6.0f, 0.0f, -6.0f),
                glm::vec3(1.0f, 1.0f, 1.0f));

    Objects.push_back(obj1);
    Objects.push_back(obj2);
    Objects.push_back(obj3);
    Objects.push_back(obj4);
    Objects.push_back(obj5);
    Objects.push_back(obj6);
}


int main()
{
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    GLFWwindow* window = glfwCreateWindow(WIDTH, HEIGHT, "zadaca_4", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetScrollCallback(window, scroll_callback);

    glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    stbi_set_flip_vertically_on_load(true);
    
    glEnable(GL_DEPTH_TEST);

    initLights();
    initObjects();

    while (!glfwWindowShouldClose(window))
    {
        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        processInput(window);

        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        glm::mat4 projection = glm::perspective(glm::radians(camera.Zoom), (float)WIDTH / (float)HEIGHT, 0.1f, 100.0f);
        glm::mat4 view = camera.GetViewMatrix();

        for (auto&& object : Objects)
        {
            object.Draw(camera, projection, view, Lights, currentFrame, RotationSpeed);
        }

        for (auto&& light : Lights)
        {
            light.Draw(projection, view);
        }

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwTerminate();
    return 0;
}

bool isClicked = false;
bool isClicked2 = false;

void processInput(GLFWwindow* window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS)
        camera.ProcessKeyboard(FORWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS)
        camera.ProcessKeyboard(BACKWARD, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS)
        camera.ProcessKeyboard(LEFT, deltaTime);
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS)
        camera.ProcessKeyboard(RIGHT, deltaTime);

    if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS && isClicked == false) {
        RotationSpeed += 50;
        isClicked = true;
    }
    else if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_RELEASE && isClicked == true) {
        isClicked = false;
    }

    if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS && isClicked2 == false) {
        RotationSpeed -= 50;
        if (RotationSpeed < 0)
            RotationSpeed = 0;
        isClicked2 = true;
    }
    else if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_RELEASE && isClicked2 == true) {
        isClicked2 = false;
    }
        
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos;

    lastX = xpos;
    lastY = ypos;

    camera.ProcessMouseMovement(xoffset, yoffset);
}

void scroll_callback(GLFWwindow* window, double xoffset, double yoffset)
{
    camera.ProcessMouseScroll(yoffset);
}