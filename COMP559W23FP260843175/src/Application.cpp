/*
 *	COMP559 Winter 2023 Final Project
 *	Submission by Edwin Pan 260843175
 *
 *	Simulator code reference is Matthias M�ller's online "Ten Minute Physics" resource, particularly their video on "How to write a FLIP water / fluid simulation" - youtube https://www.youtube.com/watch?v=XmzBREkK8kY and github https://github.com/matthias-research/pages/blob/master/tenMinutePhysics/18-flip.html .
 *	Instead of reverse engineering opengl code provided for assignments I figured I'd just take an opengl tutorial starting here https://www.youtube.com/watch?v=OR4fNpBjmq8&list=PLlrATfBNZ98foTJPJ_Ev03o2oq3-GGOS2&index=2
 *
 */


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <filesystem>

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>

#include "OpenGlHelpers.cpp"

#define ASSERT(x) if (!(x)) __debugbreak();



int main() {

	GLFWwindow* window;

	/* Initialize GLFW */
	if (!glfwInit()) {
		return -1;
	}

	/* Create a windowed mode window and its OpenGL context */
	window = glfwCreateWindow(1024, 768, "COMP559 Winter 2023 Final Project 260843175", NULL, NULL);
	if (!window) {
		glfwTerminate();
		return -1;
	}

	// Mark current flgw context and initialize glew
	glfwMakeContextCurrent(window);
	if (glewInit() != GLEW_OK) {
		std::cout << "Error! Failed to initialize glew!" << std::endl;
		return -2;
	}
	std::cout << glGetString(GL_VERSION) << std::endl;

	/* Provide Data about Stuff To Render */
	float positions[] = {
		-0.5f, -0.5f, //0
		 0.5f, -0.5f, //1
		 0.5f,  0.5f, //2
		-0.5f,  0.5f, //3
	};
	unsigned int indices[] = {
		0, 1, 2,
		2, 3, 0,
	};
	// Vertex Buffer
	unsigned int buffer;
	glGenBuffers(1, &buffer);
	glBindBuffer(GL_ARRAY_BUFFER, buffer); // Select this buffer for rendering next
	glNamedBufferData(buffer, 8 * sizeof(float), positions, GL_STATIC_DRAW);
	glEnableVertexAttribArray(0);
	glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, sizeof(float) * 2, (const void*)0);
	// Index Buffer
	unsigned int ibo; // Index Buffer Object
	glGenBuffers(1, &ibo);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo); // Select this buffer for rendering next
	glNamedBufferData(ibo, 6 * sizeof(unsigned int), indices, GL_STATIC_DRAW);
	// Shader 
	ShaderProgramSource source = ParseShader("res/shaders/basic.shader");
	std::cout << "Vertex Source:" << std::endl;
	std::cout << source.VertexSource << std::endl;
	std::cout << "Fragment Source:" << std::endl;
	std::cout << source.FragmentSource << std::endl;
	unsigned int shader = CreateShader(source.VertexSource, source.FragmentSource);
	glUseProgram(shader);

	float r = 0.0f;
	float rIncrement = 0.05f;
	GLCall(int location = glGetUniformLocation(shader, "u_Color"));
	ASSERT(location != -1);
	GLCall(glUniform4f(location, r, 0.3f, 0.8f, 1.0f));


	// Rendering loop
	while (!glfwWindowShouldClose(window)) {
		glClear(GL_COLOR_BUFFER_BIT);

		GLCall(glUniform4f(location, r, 0.3f, 0.8f, 1.0f));
		GLCall(glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, nullptr));

		r += rIncrement;
		if (r >= 1) {
			r = 1;
			rIncrement *= -1;
		}
		else if (r <= 0) {
			r = 0;
			rIncrement *= -1;
		}

		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	// Clean up
	glDeleteProgram(shader);

	glfwTerminate();
	return 0;

}