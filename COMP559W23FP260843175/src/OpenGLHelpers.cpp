/*
* A collection of OpenGL helper functions.
* 
* Based on Youtube Channel "The Cherno"'s tutorial series on OpenGL, one of the first episodes of which are found here: https://www.youtube.com/watch?v=OR4fNpBjmq8&list=PLlrATfBNZ98foTJPJ_Ev03o2oq3-GGOS2&index=2
* 
* @author Edwin Pan (260843175) for COMP559 Winter 2023 Final Project
*/

#include "includes.h"

#define ASSERT(x) if (!(x)) __debugbreak();
#define GLCall(x) GLClearError (); x; ASSERT(GLLogCall(#x, __FILE__, __LINE__))

void GLClearError() {
	while (glGetError() != GL_NO_ERROR);
}

bool GLLogCall(const char* function, const char* file, int line) {
	while (GLenum error = glGetError()) {
		std::cout << "[OpenGL Error] (" << error << "): " << function << " " << file << ":" << line << std::endl;
		return false;
	}
	return true;
}

struct ShaderProgramSource {
	std::string VertexSource;
	std::string FragmentSource;
};

ShaderProgramSource ParseShader(const std::string& filepath) {
	// Open file
	std::ifstream stream(filepath);
	// Set up shader source buffers
	enum class ShaderType {
		NONE = -1, VERTEX = 0, FRAGMENT = 1
	};
	std::stringstream ss[2];
	ShaderType type = ShaderType::NONE;
	// Extract shader sources from file
	std::string line;
	while (getline(stream, line)) {
		if (line.find("#shader") != std::string::npos) {
			if (line.find("vertex") != std::string::npos) {
				// set mode to vertex
				type = ShaderType::VERTEX;
			}
			else if (line.find("fragment") != std::string::npos) {
				// set mode to fragment
				type = ShaderType::FRAGMENT;

			}
		}
		else if (type != ShaderType::NONE) {
			ss[(int)type] << line << '\n';
		}
	}
	// return shader sources
	return { ss[0].str(), ss[1].str() };
}

unsigned int CompileShader(unsigned int type, const std::string& source) {
	// Initialize, prepare source, and compile shader
	unsigned int id = glCreateShader(type);
	const char* src = source.c_str();
	glShaderSource(id, 1, &src, nullptr);
	glCompileShader(id);
	// Query for errors
	int result;
	glGetShaderiv(id, GL_COMPILE_STATUS, &result);
	if (result == GL_FALSE) {
		int length;
		glGetShaderiv(id, GL_INFO_LOG_LENGTH, &length);
		char* message = (char*)malloc(length * sizeof(char));
		glGetShaderInfoLog(id, length, &length, message);
		std::cout << "Failed to compile " << (type == GL_VERTEX_SHADER ? "vertex" : "fragment") << " shader" << std::endl;
		std::cout << message << std::endl;
		free(message);
		glDeleteShader(id);
	}
	// Return shader
	return id;
}

unsigned int CreateShader(const std::string& vertexShader, const std::string& fragmentShader) {
	unsigned int program = glCreateProgram();
	unsigned int vs = CompileShader(GL_VERTEX_SHADER, vertexShader);
	unsigned int fs = CompileShader(GL_FRAGMENT_SHADER, fragmentShader);
	// Assemble
	glAttachShader(program, vs);
	glAttachShader(program, fs);
	glLinkProgram(program);
	glValidateProgram(program);
	// Delete shaders as parts, keeping program and keeping source code
	glDeleteShader(vs);
	glDeleteShader(fs);
	// Return the shader
	return program;
}

/*
* Because why didn't The Cherno on youtube just make it a single function?
*/
unsigned int ParseAndCreateShader(const std::string& filepath) {
	std::cout << "Reading in shader source code at " << filepath << std::endl;
	ShaderProgramSource source = ParseShader(filepath);
	std::cout << "Vertex Source:" << std::endl;
	std::cout << source.VertexSource << std::endl;
	std::cout << "Fragment Source:" << std::endl;
	std::cout << source.FragmentSource << std::endl;
	return CreateShader(source.VertexSource, source.FragmentSource);
}

