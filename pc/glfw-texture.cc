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
#include <vector>

using namespace std;

const int g_width = 640;
const int g_height = 480;

static GLuint crateShader() {
  //バーテックスシェーダのコンパイル
  GLuint vShaderId = glCreateShader(GL_VERTEX_SHADER);
  string vertexShader = R"#(
    attribute vec3 position;
    attribute vec2 uv;
    varying vec2 vuv;
    void main(void){
        gl_Position = vec4(position, 1.0);
        vuv = uv;
    }
    )#";
  const char* vs = vertexShader.c_str();
  glShaderSource(vShaderId, 1, &vs, nullptr);
  glCompileShader(vShaderId);

  //フラグメントシェーダのコンパイル
  GLuint fShaderId = glCreateShader(GL_FRAGMENT_SHADER);
  string fragmentShader = R"#(
    varying vec2 vuv;
    uniform sampler2D texture;
    void main(void){
        gl_FragColor = texture2D(texture, vuv);
    }
    )#";
  const char* fs = fragmentShader.c_str();
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

static GLuint loadTexture(const string& filename) {
  (void)filename;
  // テクスチャIDの生成
  GLuint texID = 0;
  glGenTextures(1, &texID);

  // // ファイルの読み込み
  // std::ifstream fstr(filename, std::ios::binary);
  // const size_t fileSize = static_cast<size_t>(fstr.seekg(0,
  // fstr.end).tellg()); fstr.seekg(0, fstr.beg); char* textureBuffer = new
  // char[fileSize]; fstr.read(textureBuffer, fileSize);
  std::vector<unsigned char> textureBuffer(256 * 256 * 3, 255);

  // テクスチャをGPUに転送
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glBindTexture(GL_TEXTURE_2D, texID);
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 256, 256, 0, GL_RGB, GL_UNSIGNED_BYTE,
               textureBuffer.data());

  // テクスチャの設定
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

  // テクスチャのアンバインド
  glBindTexture(GL_TEXTURE_2D, 0);

  return texID;
}

int main() {
  if (glfwInit() == 0) {
    return -1;
  }

  GLFWwindow* window =
      glfwCreateWindow(g_width, g_height, "Simple", nullptr, nullptr);
  if (window == nullptr) {
    glfwTerminate();
    return -1;
  }

  // このwindowをターゲットにする
  glfwMakeContextCurrent(window);

  // Init GLEW after making glfw's window
  glewExperimental = GL_TRUE;
  if (glewInit() != GLEW_OK) {
    // NOLINTNEXTLINE
    fprintf(stderr, "failed initializing GLEW.\n");
    // NOLINTNEXTLINE
    exit(EXIT_FAILURE);
  }
  // NOLINTNEXTLINE
  fprintf(stderr, "success initializing GLEW.\n");

  // モニタとの同期
  glfwSwapInterval(1);

  GLuint programId = crateShader();

  GLuint texID = loadTexture("cat.raw");

  // ゲームループ
  while (glfwWindowShouldClose(window) == 0) {
    // 画面の初期化
    glClearColor(0.2f, 0.2f, 0.2f, 0.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearDepth(1.0);

    // 頂点データ
    float vertex_position[] = {0.5f,  0.5f,  -0.5f, 0.5f,
                               -0.5f, -0.5f, 0.5f,  -0.5f};

    const GLfloat vertex_uv[] = {
        1, 0, 0, 0, 0, 1, 1, 1,
    };

    // 何番目のattribute変数か
    int positionLocation = glGetAttribLocation(programId, "position");
    int uvLocation = glGetAttribLocation(programId, "uv");
    int textureLocation = glGetUniformLocation(programId, "texture");

    // attribute属性を有効にする
    glEnableVertexAttribArray(positionLocation);
    glEnableVertexAttribArray(uvLocation);

    // uniform属性を設定する
    glUniform1i(textureLocation, 0);

    // attribute属性を登録
    glVertexAttribPointer(positionLocation, 2, GL_FLOAT, 0u, 0,
                          vertex_position);
    glVertexAttribPointer(uvLocation, 2, GL_FLOAT, 0u, 0, vertex_uv);

    // モデルの描画
    glBindTexture(GL_TEXTURE_2D, texID);
    glDrawArrays(GL_TRIANGLE_FAN, 0, 4);

    // バッファの入れ替え
    glfwSwapBuffers(window);

    // イベント待ち
    glfwPollEvents();
  }

  glfwTerminate();

  return 0;
}
