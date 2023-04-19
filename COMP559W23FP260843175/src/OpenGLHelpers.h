/*
* Header file for opengl helper functions
*/

#ifndef OPENGLHELPERS_H
#define OPENGLHELPERS_H

#define ASSERT(x) if (!(x)) __debugbreak();
#define GLCall(x) GLClearError (); x; ASSERT(GLLogCall(#x, __FILE__, __LINE__))
void GLClearError();
bool GLLogCall(const char* function, const char* file, int line);
struct ShaderProgramSource {
	std::string VertexSource;
	std::string FragmentSource;
};
ShaderProgramSource ParseShader(const std::string& filepath);
unsigned int CompileShader(unsigned int type, const std::string& source);
unsigned int CreateShader(const std::string& vertexShader, const std::string& fragmentShader);
unsigned int ParseAndCreateShader(const std::string& filepath);

#endif
