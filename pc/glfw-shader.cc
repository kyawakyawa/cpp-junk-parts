/*
MIT License

Copyright (c) 2019-2020 kyawakyawa

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

// reference http://marina.sys.wakayama-u.ac.jp/~tokoi/GLFWdraft.pdf

#include <GL/glew.h>
//
#include <GLFW/glfw3.h>

#include <fstream>
#include <iostream>

using namespace std;

const int g_width  = 640;
const int g_height = 480;

static GLuint crateShader() {
  //バーテックスシェーダのコンパイル
  GLuint vShaderId    = glCreateShader(GL_VERTEX_SHADER);
  string vertexShader = R"#(
    attribute vec3 position;
    attribute vec4 color;
    varying vec4 vColor;
    void main(void){
        gl_Position = vec4(position, 1.0);
        vColor = color;
    }
    )#";
  const char* vs      = vertexShader.c_str();
  glShaderSource(vShaderId, 1, &vs, nullptr);
  glCompileShader(vShaderId);

  //フラグメントシェーダのコンパイル
  GLuint fShaderId      = glCreateShader(GL_FRAGMENT_SHADER);
  string fragmentShader = R"#(
    varying vec4 vColor;
    void main(void){
        gl_FragColor = vColor;
    }
    )#";
  const char* fs        = fragmentShader.c_str();
  glShaderSource(fShaderId, 1, &fs, nullptr);
  glCompileShader(fShaderId);

  //プログラムオブジェクトの作成
  GLuint programId = glCreateProgram();
  glAttachShader(programId, vShaderId);
  glAttachShader(programId, fShaderId);

  // リンク
  glLinkProgram(programId);

  glUseProgram(programId);

  return programId;
}

int main() {
  if (!glfwInit()) {
    return -1;
  }

  GLFWwindow* window =
      glfwCreateWindow(g_width, g_height, "Simple", nullptr, nullptr);
  if (!window) {
    glfwTerminate();
    return -1;
  }

  // このwindowをターゲットにする
  glfwMakeContextCurrent(window);

  // Init GLEW after making glfw's window
  glewExperimental = GL_TRUE;
  if (glewInit() != GLEW_OK) {
    fprintf(stderr, "failed initializing GLEW.\n");
    exit(EXIT_FAILURE);
  }
  fprintf(stderr, "success initializing GLEW.\n");

  // モニタとの同期
  glfwSwapInterval(1);

  GLuint programId = crateShader();

  // ゲームループ
  while (!glfwWindowShouldClose(window)) {
    // 画面の初期化
    glClearColor(0.2f, 0.2f, 0.2f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearDepth(1.0);

    // 頂点データ
    float vertex_position[] = {1.f, 1.f, -1.f, 1.f, -1.f, -1.f, 1.f, -1.f};

    // 色情報データ
    float vertex_color[] = {
        1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0,
        0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0,
    };

    // 何番目のattribute変数か
    int positionLocation = glGetAttribLocation(programId, "position");
    int colorLocation    = glGetAttribLocation(programId, "color");

    // attribute属性を有効にする
    glEnableVertexAttribArray(positionLocation);
    glEnableVertexAttribArray(colorLocation);

    // attribute属性を登録
    glVertexAttribPointer(positionLocation, 2, GL_FLOAT, false, 0,
                          vertex_position);
    glVertexAttribPointer(colorLocation, 4, GL_FLOAT, false, 0, vertex_color);

    // モデルの描画
    glDrawArrays(GL_TRIANGLE_FAN, 0, 4);

    // バッファの入れ替え
    glfwSwapBuffers(window);

    // Poll for and process events
    glfwPollEvents();
  }

  glfwTerminate();

  return 0;
}
