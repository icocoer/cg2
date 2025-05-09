#include <GL/glut.h>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include "SOIL.h" // 确保 SOIL.h 路径正确

// --- 常量定义 ---
const int RESOLUTION = 50;
const int MAX_WAVES = 8;
const int MIN_WAVES = 1;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// --- 全局变量 ---
float waterScale = 3.0f;
float gerstner_Q = 0.6f;

static float rotate_x = 15.0f;
static float rotate_y = 25.0f;
static float translate_z = 3.5f * waterScale;

static int xold_mouse, yold_mouse;
static int left_click_state = GLUT_UP;
static int right_click_state = GLUT_UP;

bool wire_frame_mode = false;
bool show_normals = false;

int numWaves = MIN_WAVES;

struct Vector3 { float x, y, z; };

Vector3 wave_directions[MAX_WAVES];
float wave_amplitudes[MAX_WAVES];
float wave_wavelengths[MAX_WAVES];
float wave_speeds[MAX_WAVES];

// 高亮功能相关
bool highlight_sources = false; // Re-enabled for Gerstner wave component visualization
float highlight_start_time = 0.0f;

static float surface_vertices[6 * RESOLUTION * (RESOLUTION + 1)];
static float surface_normals[6 * RESOLUTION * (RESOLUTION + 1)];
static float surface_texcoords[4 * RESOLUTION * (RESOLUTION + 1)];

GLfloat LightAmbient[] = {0.5f, 0.5f, 0.6f, 1.0f};
GLfloat LightDiffuse[] = {0.8f, 0.8f, 0.9f, 1.0f};
GLfloat LightPosition[] = {1.0f * waterScale, 1.5f * waterScale, 0.5f * waterScale, 0.0f};

GLuint waterTextureID = 0;

// --- 函数声明 ---
void initWaveSources_Gerstner();
Vector3 calculateGerstnerPosition(float ox, float oz, float time);
Vector3 calculateGerstnerNormal(float ox, float oz, float time);
void drawText(float x_screen, float y_screen, const std::string &text);
bool LoadGLTextures_Gerstner();
void renderScene_Gerstner(void);
void Keyboard_Gerstner(unsigned char key, int x, int y);
void Mouse_Gerstner(int button, int state, int x, int y);
void mouseMotion_Gerstner(int x, int y);
void InitGL_Gerstner();
void changeSize_Gerstner(int w, int h);
float randomFloat_Gerstner(float min, float max);

// --- 函数实现 ---

float randomFloat_Gerstner(float min, float max) {
    if (min >= max) return min;
    return min + (static_cast<float>(rand()) / static_cast<float>(RAND_MAX)) * (max - min);
}

void initWaveSources_Gerstner() {
    for (int i = 0; i < numWaves; ++i) {
        wave_amplitudes[i] = randomFloat_Gerstner(0.02f * waterScale, 0.05f * waterScale);
        wave_wavelengths[i] = randomFloat_Gerstner(0.5f * waterScale, 1.5f * waterScale);
        if (wave_wavelengths[i] < 0.1f * waterScale) wave_wavelengths[i] = 0.1f * waterScale;
        wave_speeds[i] = randomFloat_Gerstner(0.2f * waterScale, 0.5f * waterScale);
        float angle = randomFloat_Gerstner(0.0f, 2.0f * M_PI);
        wave_directions[i].x = cos(angle);
        wave_directions[i].z = sin(angle);
        wave_directions[i].y = 0;
    }
    std::cout << "Initialized " << numWaves << " Gerstner wave components. Q: " << gerstner_Q << std::endl;
}

Vector3 calculateGerstnerPosition(float ox, float oz, float time) {
    float new_x = ox;
    float new_z = oz;
    float new_y_height = 0.0f;
    for (int i = 0; i < numWaves; ++i) {
        if (wave_wavelengths[i] == 0) continue;
        float wi = 2.0f * M_PI / wave_wavelengths[i];
        float Ai = wave_amplitudes[i];
        float dot_D_xy = wave_directions[i].x * ox + wave_directions[i].z * oz;
        float phase_phi_t = wi * wave_speeds[i] * time;
        float arg = wi * dot_D_xy + phase_phi_t;
        float cos_arg = cos(arg);
        float sin_arg = sin(arg);
        new_x += gerstner_Q * Ai * wave_directions[i].x * cos_arg;
        new_z += gerstner_Q * Ai * wave_directions[i].z * cos_arg;
        new_y_height += Ai * sin_arg;
    }
    return {new_x, new_y_height, new_z};
}

Vector3 calculateGerstnerNormal(float ox, float oz, float time) {
    float sum_Nx = 0.0f;
    float sum_Nz = 0.0f;
    float sum_Ny_S_term = 0.0f;
    for (int i = 0; i < numWaves; ++i) {
        if (wave_wavelengths[i] == 0) continue;
        float wi = 2.0f * M_PI / wave_wavelengths[i];
        float Ai = wave_amplitudes[i];
        float dot_D_xy = wave_directions[i].x * ox + wave_directions[i].z * oz;
        float phase_phi_t = wi * wave_speeds[i] * time;
        float arg = wi * dot_D_xy + phase_phi_t;
        float cos_arg = cos(arg);
        float sin_arg = sin(arg);
        float WAi = wi * Ai;
        sum_Nx += wave_directions[i].x * WAi * cos_arg;
        sum_Nz += wave_directions[i].z * WAi * cos_arg;
        sum_Ny_S_term += gerstner_Q * WAi * sin_arg;
    }
    Vector3 n_vec;
    n_vec.x = -sum_Nx; n_vec.z = -sum_Nz; n_vec.y = 1.0f - sum_Ny_S_term;
    float l = sqrt(n_vec.x*n_vec.x + n_vec.y*n_vec.y + n_vec.z*n_vec.z);
    if (l > 1e-6) {
        n_vec.x /= l; n_vec.y /= l; n_vec.z /= l;
    } else {
        n_vec.x = 0; n_vec.y = 1; n_vec.z = 0;
    }
    return n_vec;
}

void renderScene_Gerstner(void) {
    const float t = glutGet(GLUT_ELAPSED_TIME) / 1000.0f;
    const float totalSpan = 2.0f * waterScale;
    const float delta = totalSpan / RESOLUTION;
    const unsigned int strip_length_verts = 2 * (RESOLUTION + 1);

    unsigned int vert_idx, tex_idx;
    float original_x, original_z_bottom, original_z_top;

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glTranslatef(0.0f, 0.0f, -translate_z);
    glRotatef(rotate_y, 1.0f, 0.0f, 0.0f);
    glRotatef(rotate_x, 0.0f, 1.0f, 0.0f);

    for (unsigned int j = 0; j < RESOLUTION; ++j) {
        original_z_bottom = (j + 1) * delta - waterScale;
        original_z_top = j * delta - waterScale;
        for (unsigned int i = 0; i <= RESOLUTION; ++i) {
            vert_idx = 6 * (i + j * (RESOLUTION + 1));
            tex_idx = 4 * (i + j * (RESOLUTION + 1));
            original_x = i * delta - waterScale;

            Vector3 p_bottom = calculateGerstnerPosition(original_x, original_z_bottom, t);
            surface_vertices[vert_idx + 3] = p_bottom.x;
            surface_vertices[vert_idx + 4] = p_bottom.y;
            surface_vertices[vert_idx + 5] = p_bottom.z;
            surface_texcoords[tex_idx + 2] = (original_x + waterScale) / totalSpan * waterScale * 0.4f;
            surface_texcoords[tex_idx + 3] = (original_z_bottom + waterScale) / totalSpan * waterScale * 0.4f;

            Vector3 p_top = calculateGerstnerPosition(original_x, original_z_top, t);
            surface_vertices[vert_idx + 0] = p_top.x;
            surface_vertices[vert_idx + 1] = p_top.y;
            surface_vertices[vert_idx + 2] = p_top.z;
            surface_texcoords[tex_idx + 0] = (original_x + waterScale) / totalSpan * waterScale * 0.4f;
            surface_texcoords[tex_idx + 1] = (original_z_top + waterScale) / totalSpan * waterScale * 0.4f;

            Vector3 n_top = calculateGerstnerNormal(original_x, original_z_top, t);
            surface_normals[vert_idx+0]=n_top.x; surface_normals[vert_idx+1]=n_top.y; surface_normals[vert_idx+2]=n_top.z;
            Vector3 n_bottom = calculateGerstnerNormal(original_x, original_z_bottom, t);
            surface_normals[vert_idx+3]=n_bottom.x; surface_normals[vert_idx+4]=n_bottom.y; surface_normals[vert_idx+5]=n_bottom.z;
        }
    }

    // --- 高亮 "震源" (波分量方向) ---
    if (highlight_sources) {
        float current_time_highlight = glutGet(GLUT_ELAPSED_TIME) / 1000.0f;
        if (current_time_highlight - highlight_start_time > 1.0f) {
            highlight_sources = false;
        } else {
            glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_CURRENT_BIT | GL_LINE_BIT);
            glDisable(GL_LIGHTING);
            glDisable(GL_TEXTURE_2D);
            glColor3f(1.0f, 1.0f, 0.0f); // Yellow color for direction vectors
            glLineWidth(2.0f); // Make lines a bit thicker

            glBegin(GL_LINES);
            for (int i = 0; i < numWaves; ++i) {
                float line_length = 0.5f * waterScale; // Length of the direction vector visualization
                float start_y = 0.2f * waterScale;   // Draw slightly above the average water plane for visibility

                // Start point of the line (e.g., from origin, slightly elevated)
                glVertex3f(0.0f, start_y, 0.0f);
                // End point of the line
                glVertex3f(wave_directions[i].x * line_length,
                           start_y, // Keep it in the horizontal plane for clarity of direction
                           wave_directions[i].z * line_length);
            }
            glEnd();
            glPopAttrib(); // Restore attributes
        }
    }


    if (wire_frame_mode) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); glDisable(GL_TEXTURE_2D);
    } else {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D, waterTextureID);
    }
    glColor4f(0.8f, 0.9f, 1.0f, 0.80f);
    glEnable(GL_BLEND); glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LIGHTING);
    glEnableClientState(GL_VERTEX_ARRAY); glEnableClientState(GL_NORMAL_ARRAY);
    if (!wire_frame_mode) {
        glEnableClientState(GL_TEXTURE_COORD_ARRAY);
        glTexCoordPointer(2, GL_FLOAT, 0, surface_texcoords);
    }
    glVertexPointer(3, GL_FLOAT, 0, surface_vertices);
    glNormalPointer(GL_FLOAT, 0, surface_normals);
    for (unsigned int jd=0; jd<RESOLUTION; ++jd)
        glDrawArrays(GL_TRIANGLE_STRIP, jd*strip_length_verts, strip_length_verts);
    glDisableClientState(GL_VERTEX_ARRAY); glDisableClientState(GL_NORMAL_ARRAY);
    if (!wire_frame_mode) {
        glDisableClientState(GL_TEXTURE_COORD_ARRAY); glDisable(GL_TEXTURE_2D);
    }
    glDisable(GL_BLEND);

    if (show_normals) {
        glPushAttrib(GL_LIGHTING_BIT | GL_TEXTURE_BIT | GL_CURRENT_BIT);
        glDisable(GL_LIGHTING); glDisable(GL_TEXTURE_2D);
        glColor3f(1.0f,0.0f,0.0f); glBegin(GL_LINES);
        for (unsigned int jn=0; jn<RESOLUTION; ++jn) {
            for (unsigned int in=0; in<=RESOLUTION; ++in) {
                vert_idx = 6*(in + jn*(RESOLUTION+1));
                for(int k_vert=0; k_vert<2; ++k_vert) {
                    int current_vert_offset = k_vert * 3;
                     glVertex3f(surface_vertices[vert_idx+current_vert_offset+0], surface_vertices[vert_idx+current_vert_offset+1], surface_vertices[vert_idx+current_vert_offset+2]);
                     glVertex3f(surface_vertices[vert_idx+current_vert_offset+0] + surface_normals[vert_idx+current_vert_offset+0]/15.f,
                                surface_vertices[vert_idx+current_vert_offset+1] + surface_normals[vert_idx+current_vert_offset+1]/15.f,
                                surface_vertices[vert_idx+current_vert_offset+2] + surface_normals[vert_idx+current_vert_offset+2]/15.f);
                }
            }
        }
        glEnd(); glPopAttrib();
    }

    std::stringstream ss; ss << "Gerstner: " << numWaves << " waves | Q: " << gerstner_Q;
    drawText(10.0f, glutGet(GLUT_WINDOW_HEIGHT)-25.0f, ss.str());
    glutSwapBuffers(); glutPostRedisplay();
}

void Keyboard_Gerstner(unsigned char key, int x, int y) {
    switch (key) {
        case 27: exit(0); break; // ESC for exit
        case '1': wire_frame_mode = !wire_frame_mode; break;
        case 'n': show_normals = !show_normals; break;
        case '+': case '=': if (numWaves < MAX_WAVES) { numWaves++; initWaveSources_Gerstner(); } break;
        case '-': if (numWaves > MIN_WAVES) { numWaves--; initWaveSources_Gerstner(); } break;
        case 's': case 'S': // Handle 's' for highlighting
            highlight_sources = true;
            highlight_start_time = glutGet(GLUT_ELAPSED_TIME) / 1000.0f;
            std::cout << "Highlighting wave component directions..." << std::endl;
            break;
        case 'Q':
        case 'q':
            gerstner_Q += 0.05f;
            if (gerstner_Q > 1.5f) gerstner_Q = 1.5f;
            std::cout << "Gerstner Q set to: " << gerstner_Q << std::endl;
            break;
        case 'E':
        case 'e':
            gerstner_Q -= 0.05f;
            if (gerstner_Q < 0.05f) gerstner_Q = 0.05f;
            std::cout << "Gerstner Q set to: " << gerstner_Q << std::endl;
            break;
    }
    glutPostRedisplay();
}

// Mouse, InitGL, changeSize, LoadGLTextures_Gerstner, drawText, main functions remain the same as the previous corrected Gerstner code
void Mouse_Gerstner(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON) left_click_state = state;
    if (button == GLUT_RIGHT_BUTTON) right_click_state = state;
    xold_mouse = x; yold_mouse = y;
}

void mouseMotion_Gerstner(int x, int y) {
    if (left_click_state == GLUT_DOWN) {
        rotate_y += (y - yold_mouse) / 5.0f;
        rotate_x += (x - xold_mouse) / 5.0f;
        if (rotate_y > 89.0f) rotate_y = 89.0f; if (rotate_y < -89.0f) rotate_y = -89.0f;
    }
    if (right_click_state == GLUT_DOWN) {
        translate_z += (yold_mouse - y) / (20.0f / waterScale);
        float min_zoom = 0.3f * waterScale; float max_zoom = 15.0f * waterScale;
        if (translate_z < min_zoom) translate_z = min_zoom;
        if (translate_z > max_zoom) translate_z = max_zoom;
    }
    xold_mouse = x; yold_mouse = y;
    glutPostRedisplay();
}

void InitGL_Gerstner() {
    srand(static_cast<unsigned int>(time(0)));
    if (!LoadGLTextures_Gerstner()) { /* Handle texture loading error */ }
    glClearColor(0.02f, 0.03f, 0.1f, 0.0f);
    glClearDepth(1.0f);
    glEnable(GL_DEPTH_TEST); glDepthFunc(GL_LEQUAL);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glLightfv(GL_LIGHT1, GL_AMBIENT, LightAmbient);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuse);
    glLightfv(GL_LIGHT1, GL_POSITION, LightPosition);
    glEnable(GL_LIGHT1); glEnable(GL_LIGHTING);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glEnable(GL_COLOR_MATERIAL);
    initWaveSources_Gerstner();
}

void changeSize_Gerstner(int w, int h) {
    if (h == 0) h = 1;
    float ratio = 1.0f * w / h;
    glMatrixMode(GL_PROJECTION); glLoadIdentity();
    glViewport(0, 0, w, h);
    gluPerspective(50.0f, ratio, 0.1f * waterScale, 150.0f * waterScale);
    glMatrixMode(GL_MODELVIEW);
}

bool LoadGLTextures_Gerstner() {
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

int main(int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_MULTISAMPLE);
    glutInitWindowSize(1280, 720);
    glutCreateWindow("Gerstner Wave Water Surface");
    InitGL_Gerstner();
    glutDisplayFunc(renderScene_Gerstner);
    glutReshapeFunc(changeSize_Gerstner);
    glutKeyboardFunc(Keyboard_Gerstner);
    glutMouseFunc(Mouse_Gerstner);
    glutMotionFunc(mouseMotion_Gerstner);
    glutMainLoop();
    return 0;
}