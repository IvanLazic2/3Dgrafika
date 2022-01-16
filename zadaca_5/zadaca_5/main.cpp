#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>
#include <limits>
#include "glad/glad.h"
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include "CircleRenderer.h"
#include "LineRenderer.h"

using namespace std;
using namespace glm;

vec2* selected_point = nullptr;

int Factorial(int x) {
    if (x > 1) return x * Factorial(x - 1);
    else return 1;
}

int BinomialCoefficiant(int n, int k) {
    return Factorial(n) / (Factorial(k) * Factorial(n - k));
}

struct BezierCurve {
    vec2 p0, p1, p2, p3;
    vector<vec2> points;
    vec2* selected_point = nullptr;

    BezierCurve(vec2 p0, vec2 p1, vec2 p2, vec2 p3) : p0(p0), p1(p1), p2(p2), p3(p3) {}

    vector<vec2>& GetCurve() {
        return points;
    }

    void CreateCurve(float res) {
        for (float i = 0; i <= 1.0; i += res) {
            points.push_back((float)pow(1 - i, 3) * p0 + 3 * (float)pow(1 - i, 2) * i * p1 + 3 * (1 - i) * (float)pow(i, 2) * p2 + (float)pow(i, 3) * p3);
        }
    }

    void ClearPoints() {
        points.clear();
    }
};
BezierCurve bezier(vec2(50.0, 50.0), vec2(800.0, 50.0), vec2(50.0, 300.0), vec2(900.0, 500));

void mouse_callback(GLFWwindow* window, int button, int action, int mods) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            double x;
            double y;

            glfwGetCursorPos(window, &x, &y);

            if (sqrt(pow(x - bezier.p0.x, 2) + pow(y - bezier.p0.y, 2)) < 5) 
                selected_point = &bezier.p0;
            else if (sqrt(pow(x - bezier.p1.x, 2) + pow(y - bezier.p1.y, 2)) < 5) 
                selected_point = &bezier.p1;
            else if (sqrt(pow(x - bezier.p2.x, 2) + pow(y - bezier.p2.y, 2)) < 5) 
                selected_point = &bezier.p2;
            else if (sqrt(pow(x - bezier.p3.x, 2) + pow(y - bezier.p3.y, 2)) < 5) 
                selected_point = &bezier.p3;
        }
        else if (action == GLFW_RELEASE)
            selected_point = nullptr;
    }
}

void cursor_callback(GLFWwindow* window, double x, double y) {
    if (selected_point != nullptr) {
        selected_point->x = x;
        selected_point->y = y;
        bezier.ClearPoints();
        bezier.CreateCurve(0.01);
    }
}

// minGW: g++ main.cpp glad\glad.c glfw\lib-mingw-w64\libglfw3.a LineRenderer.cpp CircleRenderer.cpp -l gdi32 -l opengl32 -std=c++17
int main() {
    GLFWwindow* window;

    if (!glfwInit()) cout << "Error : could not initilize GLFW";

    int width = 1000;
    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);
    glfwWindowHint(GLFW_SAMPLES, 4);
    window = glfwCreateWindow(width, width * 9 / 16, "Hello World", NULL, NULL);
    if (!window) {
        cout << "Error : could not create window";
        glfwTerminate();
    }

    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) cout << "Error : could not initilize Glad";

    glfwSwapInterval(1);

    InitCircleRendering(32);
    InitLineRendering();

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    bezier.CreateCurve(0.01);

    glfwSetMouseButtonCallback(window, mouse_callback);
    glfwSetCursorPosCallback(window, cursor_callback);

    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT);

        RenderCircle(bezier.p0, 5);
        RenderCircle(bezier.p1, 5);
        RenderCircle(bezier.p2, 5);
        RenderCircle(bezier.p3, 5);

        RenderLine(bezier.GetCurve());

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    return 0;
}