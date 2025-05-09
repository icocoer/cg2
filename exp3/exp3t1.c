#include <GL/glut.h>   
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>   
#include <time.h>
#include "SOIL.h" 

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MAX_PARTICLES 2000 

// 粒子结构定义
typedef struct {
    bool active;
    float life;
    float fade;
    float r, g, b;
    float x, y, z;
    float v_x, v_y, v_z;
    float a_x, a_y, a_z;
} particle;

particle particles[MAX_PARTICLES];

// 全局变量
float slowdown = 2.5f;
float zoom = -50.0f;
float portalRadius = 10.0f;
GLuint texture[1];
float particle_size = 1.8f;

// 函数原型
void InitGL(void);
int LoadGLTextures(void);
void renderScene(void);
void changeSize(int w, int h);
void TimerFunction(int value);
float random_float(float min, float max);

// 主函数
int main(int argc, char **argv) {
    srand(time(NULL));

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(800, 600);
    glutCreateWindow("Doctor Strange Portal");

    glutDisplayFunc(renderScene);
    glutReshapeFunc(changeSize);
    glutTimerFunc(20, TimerFunction, 1);

    InitGL();

    glutMainLoop();
    return 0;
}

// 生成指定范围的随机浮点数
float random_float(float min, float max) {
    if (min > max) {
        float temp = min;
        min = max;
        max = temp;
    }
    return min + ((float)rand() / (float)RAND_MAX) * (max - min);
}

// 加载纹理函数
int LoadGLTextures() {
    texture[0] = SOIL_load_OGL_texture(
        "Particle.bmp",
        SOIL_LOAD_AUTO,
        SOIL_CREATE_NEW_ID,
        SOIL_FLAG_INVERT_Y
    );

    if (texture[0] == 0) {
        printf("SOIL loading error for 'Particle.bmp': '%s'\n", SOIL_last_result());
        return false;
    }

    glBindTexture(GL_TEXTURE_2D, texture[0]);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    return true;
}

// OpenGL及粒子初始化
void InitGL(void) {
    if (!LoadGLTextures()) {
        printf("Failed to load textures! Exiting.\n");
        exit(1);
    }

    glEnable(GL_TEXTURE_2D);
    glShadeModel(GL_SMOOTH);
    glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
    glClearDepth(1.0f);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

    for (int i = 0; i < MAX_PARTICLES; i++) {
        particles[i].active = true;
        particles[i].life = random_float(0.4f, 1.2f);    
        particles[i].fade = random_float(0.015f, 0.05f); 

        float r_color_factor = random_float(0.0f, 1.0f);
        if (r_color_factor < 0.75f) {
            particles[i].r = 1.0f;
            particles[i].g = random_float(0.3f, 0.6f);
            particles[i].b = random_float(0.0f, 0.1f);
        } else {
            particles[i].r = 1.0f;
            particles[i].g = random_float(0.7f, 1.0f);
            particles[i].b = random_float(0.2f, 0.5f);
        }

        float emit_angle = random_float(0.0f, 2.0f * M_PI);
        particles[i].x = portalRadius * cosf(emit_angle);
        particles[i].y = portalRadius * sinf(emit_angle);
        particles[i].z = random_float(-0.5f, 0.5f);

        float outward_speed_base = 3.0f; // 向外速度基准
        float outward_speed_random_factor = random_float(0.8f, 2.5f); // 向外速度随机因子
        float outward_speed = outward_speed_base * outward_speed_random_factor;

        float tangent_speed_base = 3.0f; // 切向速度基准 (保持一定的旋转)
        float tangent_speed_random_factor = random_float(0.7f, 1.3f);
        float tangent_speed = tangent_speed_base * tangent_speed_random_factor;

        float z_spread_speed = random_float(-1.5f, 1.5f); // Z轴扩散

        float v_out_x = cosf(emit_angle);
        float v_out_y = sinf(emit_angle);
        float v_tan_x = -sinf(emit_angle);
        float v_tan_y = cosf(emit_angle);

        particles[i].v_x = v_out_x * outward_speed + v_tan_x * tangent_speed;
        particles[i].v_y = v_out_y * outward_speed + v_tan_y * tangent_speed;
        particles[i].v_z = z_spread_speed;

        particles[i].a_x = 0.0f;
        particles[i].a_y = 0.0f;
        particles[i].a_z = 0.0f;
    }
}

// 渲染场景
void renderScene(void) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    glBindTexture(GL_TEXTURE_2D, texture[0]);

    for (int i = 0; i < MAX_PARTICLES; i++) {
        if (particles[i].active) {
            float p_x = particles[i].x;
            float p_y = particles[i].y;
            float p_z = particles[i].z + zoom;

            glColor4f(particles[i].r, particles[i].g, particles[i].b, particles[i].life);

            float half_size = particle_size / 2.0f;

            glBegin(GL_TRIANGLE_STRIP);
                glTexCoord2f(1.0f, 1.0f); glVertex3f(p_x + half_size, p_y + half_size, p_z);
                glTexCoord2f(0.0f, 1.0f); glVertex3f(p_x - half_size, p_y + half_size, p_z);
                glTexCoord2f(1.0f, 0.0f); glVertex3f(p_x + half_size, p_y - half_size, p_z);
                glTexCoord2f(0.0f, 0.0f); glVertex3f(p_x - half_size, p_y - half_size, p_z);
            glEnd();

            particles[i].x += particles[i].v_x / (slowdown * 100.0f);
            particles[i].y += particles[i].v_y / (slowdown * 100.0f);
            particles[i].z += particles[i].v_z / (slowdown * 100.0f);

            particles[i].life -= particles[i].fade;

            if (particles[i].life < 0.0f) {

                particles[i].life = random_float(0.4f, 1.2f);
                particles[i].fade = random_float(0.015f, 0.05f);

                float r_color_factor = random_float(0.0f, 1.0f);
                 if (r_color_factor < 0.75f) {
                    particles[i].r = 1.0f; particles[i].g = random_float(0.3f, 0.6f); particles[i].b = random_float(0.0f, 0.1f);
                } else {
                    particles[i].r = 1.0f; particles[i].g = random_float(0.7f, 1.0f); particles[i].b = random_float(0.2f, 0.5f);
                }

                float emit_angle = random_float(0.0f, 2.0f * M_PI);
                particles[i].x = portalRadius * cosf(emit_angle);
                particles[i].y = portalRadius * sinf(emit_angle);
                particles[i].z = random_float(-0.5f, 0.5f);

                float outward_speed_base = 3.0f;
                float outward_speed_random_factor = random_float(0.8f, 2.5f);
                float outward_speed = outward_speed_base * outward_speed_random_factor;

                float tangent_speed_base = 3.0f;
                float tangent_speed_random_factor = random_float(0.7f, 1.3f);
                float tangent_speed = tangent_speed_base * tangent_speed_random_factor;

                float z_spread_speed = random_float(-1.5f, 1.5f);

                float v_out_x = cosf(emit_angle); float v_out_y = sinf(emit_angle);
                float v_tan_x = -sinf(emit_angle); float v_tan_y = cosf(emit_angle);

                particles[i].v_x = v_out_x * outward_speed + v_tan_x * tangent_speed;
                particles[i].v_y = v_out_y * outward_speed + v_tan_y * tangent_speed;
                particles[i].v_z = z_spread_speed;

                particles[i].a_x = 0.0f; particles[i].a_y = 0.0f; particles[i].a_z = 0.0f;
            }
        }
    }
    glutSwapBuffers();
}

// 窗口大小改变回调
void changeSize(int w, int h) {
    if (h == 0) h = 1;
    float ratio = w * 1.0 / h;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glViewport(0, 0, w, h);
    gluPerspective(45.0f, ratio, 1.0f, 200.0f);
    glMatrixMode(GL_MODELVIEW);
}

// 动画计时器回调
void TimerFunction(int value) {
    glutPostRedisplay();
    glutTimerFunc(20, TimerFunction, 1);
}