#include <GL/glut.h> // 或者在 macOS 上是 <GLUT/glut.h>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>   // 用于 std::stringstream
#include <cmath>     // 用于 sin, cos, sqrt, M_PI
#include <cstdlib>   // For rand() and srand()
#include <ctime>     // For time() for srand()
#include "SOIL.h" 
// --- 常量定义 ---
const int RESOLUTION = 60; // 水面网格精细度 (可以适当增加以适应大水面)
const int MAX_WAVES = 10;
const int MIN_WAVES = 1;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- 全局变量 ---
float waterScale = 10.0f;    // 水面在每个方向上扩展到 -waterScale 到 +waterScale (例如 3.0f)

// 摄像机和交互
static float rotate_x = 20.0f;
static float rotate_y = 30.0f;      // 初始稍微俯视
static float translate_z = 4.0f * waterScale; // 根据水面大小调整初始摄像机距离

static int xold, yold;
static int left_click = GLUT_UP;
static int right_click = GLUT_UP;

// 显示模式
bool wire_frame = false;
bool normals = false;

// 震源相关
int numWaves = MIN_WAVES;
float centers[MAX_WAVES][2];    // 震源中心坐标 (x, z)
float amplitude[MAX_WAVES];     // 振幅
float wavelength[MAX_WAVES];    // 波长
float speed[MAX_WAVES];         // 波速

// 高亮相关
bool highlight_sources = false;
float highlight_start_time = 0.0f;

// 顶点数据数组
// 注意：RESOLUTION 较大时，这些数组会很大。对于非常大的 RESOLUTION，考虑动态分配。
static float surface[6 * RESOLUTION * (RESOLUTION + 1)]; // 2 triangles per quad * 3 vertices * (x,y,z)
static float normal[6 * RESOLUTION * (RESOLUTION + 1)];  // Matching normals
static float texCoords[4 * RESOLUTION * (RESOLUTION + 1)];// 2 triangles per quad * 2 texcoords * (s,t)

// 光照参数
GLfloat LightAmbient[] = {0.6f, 0.6f, 0.7f, 1.0f}; // 调整光照以配合水面纹理
GLfloat LightDiffuse[] = {0.8f, 0.8f, 0.9f, 1.0f};
GLfloat LightPosition[] = {1.5f * waterScale, 2.0f * waterScale, 1.0f * waterScale, 0.0f}; // 光源位置也随水面调整

// 纹理ID
GLuint waterTextureID = 0;

// --- 结构体定义 ---
struct Vector3 {
    float x;
    float y;
    float z;
};

// --- 函数声明 ---
void initWaveSources();
float dot(int i, float x, float y_coord);
float wave(int i, float x, float y_coord, float time);
float waveHeight(float x, float y_coord, float time);
float dWavedx(int i, float x, float y_coord, float time);
float dWavedy(int i, float x, float y_coord, float time);
Vector3 waveNormal(float x, float y_coord, float time);
void drawText(float x_screen, float y_screen, const std::string &text);
bool LoadGLTextures();
void renderScene(void);
void Keyboard(unsigned char key, int x_mouse, int y_mouse);
void Mouse(int button, int state, int x_mouse, int y_mouse);
void mouseMotion(int x_mouse, int y_mouse);
void InitGL();
void changeSize(int w, int h);
float randomFloat(float min, float max);

// --- 函数实现 ---

float randomFloat(float min, float max) {
    if (min >= max) return min;
    return min + static_cast<float>(rand()) / (static_cast<float>(RAND_MAX / (max - min)));
}

void initWaveSources() {
    for (int i = 0; i < numWaves; ++i) {
        centers[i][0] = randomFloat(-waterScale * 0.85f, waterScale * 0.85f);
        centers[i][1] = randomFloat(-waterScale * 0.85f, waterScale * 0.85f);

        amplitude[i] = randomFloat(0.01f * waterScale, 0.03f * waterScale);
        wavelength[i] = randomFloat(0.2f * waterScale, 0.4f * waterScale);
        // 确保波长不为0，或非常小
        if (wavelength[i] < 0.1f * waterScale) wavelength[i] = 0.1f * waterScale;
        speed[i] = randomFloat(-0.35f, -0.1f);
    }
    std::cout << "Initialized " << numWaves << " random wave sources. Water scale: " << waterScale << std::endl;
}

float dot(int i, float x, float y_coord) {
    float xc = x - centers[i][0];
    float zc = y_coord - centers[i][1];
    return sqrt(xc * xc + zc * zc);
}

float wave(int i, float x, float y_coord, float time) {
    if (wavelength[i] == 0) return 0.0f;
    float frequency = 2.0f * M_PI / wavelength[i];
    float phase = speed[i] * frequency;
    float dist = dot(i, x, y_coord);
    return amplitude[i] * sin(dist * frequency + time * phase);
}

float waveHeight(float x, float y_coord, float time) {
    float height = 0.0f;
    for (int i = 0; i < numWaves; ++i) {
        height += wave(i, x, y_coord, time);
    }
    return height;
}

float dWavedx(int i, float x, float y_coord, float time) {
    if (wavelength[i] == 0) return 0.0f;
    float frequency = 2.0f * M_PI / wavelength[i];
    float phase = speed[i] * frequency;
    float dist = dot(i, x, y_coord);
    if (dist == 0) return 0.0f;
    float cos_term = cos(dist * frequency + time * phase);
    float dx_dist = (x - centers[i][0]) / dist;
    return amplitude[i] * frequency * cos_term * dx_dist;
}

float dWavedy(int i, float x, float y_coord, float time) {
    if (wavelength[i] == 0) return 0.0f;
    float frequency = 2.0f * M_PI / wavelength[i];
    float phase = speed[i] * frequency;
    float dist = dot(i, x, y_coord);
    if (dist == 0) return 0.0f;
    float cos_term = cos(dist * frequency + time * phase);
    float dy_dist = (y_coord - centers[i][1]) / dist;
    return amplitude[i] * frequency * cos_term * dy_dist;
}

Vector3 waveNormal(float x, float y_coord, float time) {
    float total_dHdx = 0.0f;
    float total_dHdz = 0.0f;
    for (int i = 0; i < numWaves; ++i) {
        total_dHdx += dWavedx(i, x, y_coord, time);
        total_dHdz += dWavedy(i, x, y_coord, time);
    }
    Vector3 n;
    n.x = -total_dHdx; n.y = 1.0f; n.z = -total_dHdz;
    float l = sqrt(n.x * n.x + n.y * n.y + n.z * n.z);
    if (l != 0) {
        n.x /= l; n.y /= l; n.z /= l;
    } else {
        n.x = 0; n.y = 1; n.z = 0;
    }
    return n;
}

void drawText(float x_screen, float y_screen, const std::string &text) {
    glMatrixMode(GL_PROJECTION); glPushMatrix(); glLoadIdentity();
    gluOrtho2D(0, glutGet(GLUT_WINDOW_WIDTH), 0, glutGet(GLUT_WINDOW_HEIGHT));
    glMatrixMode(GL_MODELVIEW); glPushMatrix(); glLoadIdentity();
    glPushAttrib(GL_LIGHTING_BIT | GL_CURRENT_BIT);
    glDisable(GL_LIGHTING); glColor3f(1.0f, 1.0f, 1.0f);
    glRasterPos2f(x_screen, y_screen);
    for (char c : text) { glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, c); }
    glPopAttrib();
    glMatrixMode(GL_PROJECTION); glPopMatrix();
    glMatrixMode(GL_MODELVIEW); glPopMatrix();
}

bool LoadGLTextures() {
    waterTextureID = SOIL_load_OGL_texture(
        "reflection.jpg", SOIL_LOAD_AUTO, SOIL_CREATE_NEW_ID,
        SOIL_FLAG_MIPMAPS | SOIL_FLAG_INVERT_Y | SOIL_FLAG_NTSC_SAFE_RGB | SOIL_FLAG_COMPRESS_TO_DXT);
    if (waterTextureID == 0) {
        std::cerr << "SOIL loading error: '" << SOIL_last_result() << "' (reflection.jpg)" << std::endl;
        return false;
    }
    glBindTexture(GL_TEXTURE_2D, waterTextureID);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    std::cout << "Texture reflection.jpg loaded with ID: " << waterTextureID << std::endl;
    return true;
}

void renderScene(void) {
    const float t = glutGet(GLUT_ELAPSED_TIME) / 1000.0f;
    const float totalSpan = 2.0f * waterScale;
    const float delta = totalSpan / RESOLUTION;
    const unsigned int strip_length_vertices = 2 * (RESOLUTION + 1);

    unsigned int vert_idx, tex_idx;
    float current_x, current_z_coord;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glTranslatef(0.0f, 0.0f, -translate_z); // Camera distance controlled by translate_z
    glRotatef(rotate_y, 1.0f, 0.0f, 0.0f); // Pitch
    glRotatef(rotate_x, 0.0f, 1.0f, 0.0f); // Yaw

    // 计算水面顶点、法线和纹理坐标
    for (unsigned int j = 0; j < RESOLUTION; ++j) {
        current_z_coord = (j + 1) * delta - waterScale;
        for (unsigned int i = 0; i <= RESOLUTION; ++i) {
            vert_idx = 6 * (i + j * (RESOLUTION + 1));
            tex_idx = 4 * (i + j * (RESOLUTION + 1));
            current_x = i * delta - waterScale;

            // 下边点
            surface[vert_idx + 3] = current_x;
            surface[vert_idx + 4] = waveHeight(current_x, current_z_coord, t);
            surface[vert_idx + 5] = current_z_coord;
            texCoords[tex_idx + 2] = (current_x + waterScale) / totalSpan * waterScale * 0.5f; // 纹理重复次数调整 (e.g. * waterScale * 0.5f)
            texCoords[tex_idx + 3] = (current_z_coord + waterScale) / totalSpan * waterScale * 0.5f;

            // 上边点
            if (j != 0) {
                unsigned int prev_vert_idx = 6 * (i + (j - 1) * (RESOLUTION + 1));
                unsigned int prev_tex_idx = 4 * (i + (j - 1) * (RESOLUTION + 1));
                surface[vert_idx + 0] = surface[prev_vert_idx + 3];
                surface[vert_idx + 1] = surface[prev_vert_idx + 4];
                surface[vert_idx + 2] = surface[prev_vert_idx + 5];
                texCoords[tex_idx + 0] = texCoords[prev_tex_idx + 2];
                texCoords[tex_idx + 1] = texCoords[prev_tex_idx + 3];
            } else {
                surface[vert_idx + 0] = current_x;
                surface[vert_idx + 1] = waveHeight(current_x, -waterScale, t);
                surface[vert_idx + 2] = -waterScale;
                texCoords[tex_idx + 0] = (current_x + waterScale) / totalSpan * waterScale * 0.5f;
                texCoords[tex_idx + 1] = (-waterScale + waterScale) / totalSpan * waterScale * 0.5f; // Should be 0
            }
            Vector3 n_top = waveNormal(surface[vert_idx + 0], surface[vert_idx + 2], t);
            normal[vert_idx + 0]=n_top.x; normal[vert_idx + 1]=n_top.y; normal[vert_idx + 2]=n_top.z;
            Vector3 n_bottom = waveNormal(surface[vert_idx + 3], surface[vert_idx + 5], t);
            normal[vert_idx + 3]=n_bottom.x; normal[vert_idx + 4]=n_bottom.y; normal[vert_idx + 5]=n_bottom.z;
        }
    }

    // 高亮震源
    if (highlight_sources) {
        float current_time_highlight = glutGet(GLUT_ELAPSED_TIME) / 1000.0f;
        if (current_time_highlight - highlight_start_time > 1.0f) {
            highlight_sources = false;
        } else {
            glPushAttrib(GL_POINT_BIT | GL_LIGHTING_BIT | GL_CURRENT_BIT | GL_TEXTURE_BIT);
            glDisable(GL_LIGHTING); glDisable(GL_TEXTURE_2D);
            glColor3f(1.0f, 1.0f, 0.0f); glPointSize(10.0f);
            glBegin(GL_POINTS);
            for (int k = 0; k < numWaves; ++k) { // iterator k to avoid conflict
                glVertex3f(centers[k][0], waveHeight(centers[k][0], centers[k][1], t) + 0.05f * waterScale, centers[k][1]);
            }
            glEnd();
            glPopAttrib();
        }
    }

    // 绘制水面
    if (wire_frame) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDisable(GL_TEXTURE_2D);
    } else {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, waterTextureID);
    }
    glColor4f(1.0f, 1.0f, 1.0f, 0.85f); // 白色与纹理混合，调整alpha
    glEnable(GL_BLEND); glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LIGHTING);
    glEnableClientState(GL_VERTEX_ARRAY); glEnableClientState(GL_NORMAL_ARRAY);
    if (!wire_frame) {
        glEnableClientState(GL_TEXTURE_COORD_ARRAY);
        glTexCoordPointer(2, GL_FLOAT, 0, texCoords);
    }
    glVertexPointer(3, GL_FLOAT, 0, surface);
    glNormalPointer(GL_FLOAT, 0, normal);
    for (unsigned int j_draw = 0; j_draw < RESOLUTION; ++j_draw) { // iterator j_draw
        glDrawArrays(GL_TRIANGLE_STRIP, j_draw * strip_length_vertices, strip_length_vertices);
    }
    glDisableClientState(GL_VERTEX_ARRAY); glDisableClientState(GL_NORMAL_ARRAY);
    if (!wire_frame) {
        glDisableClientState(GL_TEXTURE_COORD_ARRAY);
        glDisable(GL_TEXTURE_2D);
    }
    glDisable(GL_BLEND);

    // 绘制法线
    if (normals) {
        glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_CURRENT_BIT);
        glDisable(GL_LIGHTING); glDisable(GL_TEXTURE_2D);
        glColor3f(1.0f, 0.0f, 0.0f);
        glBegin(GL_LINES);
        for (unsigned int j_norm = 0; j_norm < RESOLUTION; ++j_norm) { // iterator j_norm
            for (unsigned int i_norm = 0; i_norm <= RESOLUTION; ++i_norm) { // iterator i_norm
                vert_idx = 6 * (i_norm + j_norm * (RESOLUTION + 1));
                glVertex3fv(&(surface[vert_idx]));
                glVertex3f(surface[vert_idx] + normal[vert_idx] / 20.0f,
                           surface[vert_idx+1] + normal[vert_idx+1] / 20.0f,
                           surface[vert_idx+2] + normal[vert_idx+2] / 20.0f);
            }
        }
        glEnd();
        glPopAttrib();
    }

    // 绘制文本
    std::stringstream ss; ss << "Sources: " << numWaves;
    drawText(10.0f, glutGet(GLUT_WINDOW_HEIGHT) - 25.0f, ss.str());
    glutSwapBuffers(); glutPostRedisplay();
}

void Keyboard(unsigned char key, int x_mouse, int y_mouse) {
    switch (key) {
        case 'q': case 27: exit(0); break;
        case '1': wire_frame = !wire_frame; break;
        case 'n': normals = !normals; break;
        case '+': case '=': if (numWaves < MAX_WAVES) { numWaves++; initWaveSources(); } break;
        case '-': if (numWaves > MIN_WAVES) { numWaves--; initWaveSources(); } break;
        case 's': case 'S':
            highlight_sources = true;
            highlight_start_time = glutGet(GLUT_ELAPSED_TIME) / 1000.0f;
            break;
    }
    glutPostRedisplay();
}

void Mouse(int button, int state, int x_mouse, int y_mouse) {
    if (button == GLUT_LEFT_BUTTON) left_click = state;
    if (button == GLUT_RIGHT_BUTTON) right_click = state;
    xold = x_mouse; yold = y_mouse;
}

void mouseMotion(int x_mouse, int y_mouse) {
    if (left_click == GLUT_DOWN) {
        rotate_y += (y_mouse - yold) / 5.0f;
        rotate_x += (x_mouse - xold) / 5.0f;
        if (rotate_y > 89.0f) rotate_y = 89.0f; if (rotate_y < -89.0f) rotate_y = -89.0f; // Avoid gimbal lock issues
    }
    if (right_click == GLUT_DOWN) {
        translate_z += (yold - y_mouse) / (20.0f / waterScale) ; // Scale zoom speed
        float min_zoom = 0.5f * waterScale;
        float max_zoom = 20.0f * waterScale;
        if (translate_z < min_zoom) translate_z = min_zoom;
        if (translate_z > max_zoom) translate_z = max_zoom;
    }
    xold = x_mouse; yold = y_mouse;
    glutPostRedisplay();
}

void InitGL() {
    srand(static_cast<unsigned int>(time(0)));
    if (!LoadGLTextures()) { /* Handle error or proceed without textures */ }
    glClearColor(0.05f, 0.05f, 0.15f, 0.0f); // 更暗的深蓝背景
    glClearDepth(1.0f);
    glEnable(GL_DEPTH_TEST); glDepthFunc(GL_LEQUAL);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glLightfv(GL_LIGHT1, GL_AMBIENT, LightAmbient);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuse);
    glLightfv(GL_LIGHT1, GL_POSITION, LightPosition);
    glEnable(GL_LIGHT1); glEnable(GL_LIGHTING);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    // glEnable(GL_CULL_FACE); // 可以开启面剔除提高一点效率，如果法线方向正确
    // glCullFace(GL_BACK);
    initWaveSources();
}

void changeSize(int w, int h) {
    if (h == 0) h = 1;
    float ratio = 1.0f * w / h;
    glMatrixMode(GL_PROJECTION); glLoadIdentity();
    glViewport(0, 0, w, h);
    gluPerspective(45.0f, ratio, 0.1f * waterScale, 100.0f * waterScale); // 调整近远裁剪面
    glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH /*| GLUT_MULTISAMPLE*/); // 可以尝试开启多重采样抗锯齿
    glutInitWindowSize(1024, 768); // 更大的窗口
    glutCreateWindow("Water Surfaces");
    InitGL();
    glutDisplayFunc(renderScene);
    glutReshapeFunc(changeSize);
    glutKeyboardFunc(Keyboard);
    glutMouseFunc(Mouse);
    glutMotionFunc(mouseMotion);
    // glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_CONTINUE_EXECUTION); // 使得关闭窗口时程序继续执行（如果需要调试某些状态）
    glutMainLoop();
    return 0;
}