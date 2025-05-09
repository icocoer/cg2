#include <GL/glut.h> // 或者根据您的设置可能是 <glut.h>
#include <cstdio>    // 用于 printf
#include <cstdlib>   // 用于 exit
#include <cmath>     // 用于 sqrt (如果main中直接使用，或已被Physics头文件包含)
#include <vector>    // 如果在main.cpp中也直接使用vector (通常Physics头文件会包含)

#include "Physics2.h" // 假设 ElasticSolidSimulation 在 Physics2.h 中


ElasticSolidSimulation* solidSimulation = nullptr;
const float timeStep = 0.001f; // 全局时间步长

// 摄像机变量
float camAngleX = 20.0f, camAngleY = -30.0f;
float camDistance = 20.0f;
int lastMouseX, lastMouseY;
bool mouseLeftDown = false;
bool mouseRightDown = false;


void Deinitialize() {
    if (solidSimulation) {
        solidSimulation->release(); // 调用派生类的 release
        delete solidSimulation;
        solidSimulation = nullptr;
    }
}

void InitGL() {
    // 创建 ElasticSolidSimulation 实例时，传递 timeStep
    solidSimulation = new ElasticSolidSimulation(
        6, 6, 6,                        // dimX, dimY, dimZ
        0.1f,                           // massVal
        2000.0f, 0.5f,                  // springConstStructural, springLenStructural
        1500.0f, std::sqrt(2.0f),       // springConstShear, springLenShearFactor
        500.0f, 2.0f,                   // springConstBend, springLenBendFactor
        0.75f,                          // frictionConstSpring
        Vector3D(0, -9.81f, 0),         // gravitation
        0.0f,                           // groundLevel
        0.4f,                           // restitution coefficient
        0.6f, 0.4f,                     // staticFriction, dynamicFriction
        timeStep                        // ★★★ 将全局 timeStep 传递给构造函数 ★★★
    );

    glClearColor(0.1f, 0.1f, 0.2f, 1.0f); // 深蓝色背景
    glClearDepth(1.0f);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_DEPTH_TEST);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    GLfloat light_ambient[] = { 0.2f, 0.2f, 0.2f, 1.0f };
    GLfloat light_diffuse[] = { 0.8f, 0.8f, 0.8f, 1.0f };
    GLfloat light_specular[] = { 0.5f, 0.5f, 0.5f, 1.0f };
    GLfloat light_position[] = { 10.0f, 15.0f, 10.0f, 1.0f };
    glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, light_specular);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 50.0f);
}


void changeSize(int w, int h) {
    if (h == 0) h = 1;
    float ratio = w * 1.0f / h;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glViewport(0, 0, w, h);
    gluPerspective(45.0f, ratio, 0.1f, 200.0f);
    glMatrixMode(GL_MODELVIEW);
}

void UpdateSimulation() {
    if (solidSimulation) {
        solidSimulation->operate(timeStep); // 调用基类的 operate，它会调用派生类的 solve
    }
}

void renderScene() {
    UpdateSimulation(); // 先更新物理

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    // 摄像机变换
    glTranslatef(0.0f, -1.0f, -camDistance);
    glRotatef(camAngleX, 1.0f, 0.0f, 0.0f);
    glRotatef(camAngleY, 0.0f, 1.0f, 0.0f);


    if (solidSimulation) {
        // 1. 绘制质量点 (小球)
        glColor3f(0.8f, 0.8f, 0.2f); // 黄色
        for (int i = 0; i < solidSimulation->numOfMasses; ++i) {
            if (solidSimulation->masses[i]) { // 确保 mass 存在
                glPushMatrix();
                Vector3D* pos = &solidSimulation->masses[i]->pos;
                glTranslatef(pos->x, pos->y, pos->z);
                glutSolidSphere(0.1, 8, 8);
                glPopMatrix();
            }
        }

        // 2. 绘制弹簧 (线段)
        glColor3f(0.2f, 0.8f, 0.2f); // 绿色
        glLineWidth(1.5f);
        glBegin(GL_LINES);
        for (const Spring* s : solidSimulation->springs) { // 使用 ElasticSolidSimulation 的 springs 成员
            if (s && s->mass1 && s->mass2) {
                 glVertex3fv(&s->mass1->pos.x);
                 glVertex3fv(&s->mass2->pos.x);
            }
        }
        glEnd();

        // 3. 绘制地面
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glColor4f(0.5f, 0.5f, 0.5f, 0.5f); // 半透明灰色
        glBegin(GL_QUADS);
          glNormal3f(0,1,0);
          glVertex3f(-20.0f, solidSimulation->groundLevel, -20.0f);
          glVertex3f( 20.0f, solidSimulation->groundLevel, -20.0f);
          glVertex3f( 20.0f, solidSimulation->groundLevel,  20.0f);
          glVertex3f(-20.0f, solidSimulation->groundLevel,  20.0f);
        glEnd();
        glDisable(GL_BLEND);
    }

    glutSwapBuffers();
}

void processNormalKeys(unsigned char key, int x, int y) {
    if (key == 27) { // ESC 键
        Deinitialize();
        exit(0);
    }
    if (key == 'r' || key == 'R') { // 重置模拟
        printf("Resetting simulation...\n");
        Deinitialize(); // 先释放旧的
        InitGL();       // 再创建新的
    }
}

void processSpecialKeys(int key, int x, int y) {
    Vector3D forceToApply(0,0,0);
    float forceMagnitude = 2.0f; //力量

    switch(key) {
        case GLUT_KEY_UP:    forceToApply.z -= forceMagnitude; break;
        case GLUT_KEY_DOWN:  forceToApply.z += forceMagnitude; break;
        case GLUT_KEY_LEFT:  forceToApply.x -= forceMagnitude; break;
        case GLUT_KEY_RIGHT: forceToApply.x += forceMagnitude; break;
        default: return; // 如果不是方向键，则不做任何事
    }

    if (solidSimulation) {
        // 调用 ElasticSolidSimulation 的方法来设置键盘力
        // 需要将 solidSimulation 指针转换为 ElasticSolidSimulation* 类型才能调用其特有方法
        ((ElasticSolidSimulation*)solidSimulation)->setExternalForceForKey(forceToApply);
    }
}

// ★★★ 新增/修改：处理特殊按键释放的函数 ★★★
void processSpecialKeysUp(int key, int x, int y) {
    switch (key) {
        case GLUT_KEY_UP:
        case GLUT_KEY_DOWN:
        case GLUT_KEY_LEFT:
        case GLUT_KEY_RIGHT:
            if (solidSimulation) {
                // 当按键释放时，将键盘力重置为零
                ((ElasticSolidSimulation*)solidSimulation)->setExternalForceForKey(Vector3D(0,0,0));
            }
            break;
    }
}

void mouseButton(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON) {
        mouseLeftDown = (state == GLUT_DOWN);
    } else if (button == GLUT_RIGHT_BUTTON) {
        mouseRightDown = (state == GLUT_DOWN);
    }
    lastMouseX = x;
    lastMouseY = y;
}

void mouseMove(int x, int y) {
    if (mouseLeftDown) { // 旋转摄像机
        camAngleY += (x - lastMouseX) * 0.5f;
        camAngleX += (y - lastMouseY) * 0.5f;
        if (camAngleX < -90.0f) camAngleX = -90.0f;
        if (camAngleX > 90.0f) camAngleX = 90.0f;
    }
    if (mouseRightDown) { // 缩放摄像机
        camDistance -= (y - lastMouseY) * 0.1f;
        if (camDistance < 1.0f) camDistance = 1.0f;
        if (camDistance > 100.0f) camDistance = 100.0f; // 之前的 far clipping plane 是 100，这里也限制到100，或者调整gluPerspective
    }
    lastMouseX = x;
    lastMouseY = y;
    glutPostRedisplay();
}


int main(int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(800, 600);
    glutCreateWindow("3D Elastic Solid Simulation");

    InitGL(); // 调用初始化

    glutDisplayFunc(renderScene);
    glutReshapeFunc(changeSize);
    glutIdleFunc(renderScene); // 使用 renderScene 作为 idle 函数以保持更新

    glutKeyboardFunc(processNormalKeys);
    glutSpecialFunc(processSpecialKeys);
    glutSpecialUpFunc(processSpecialKeysUp); // ★★★ 注册按键释放回调函数 ★★★
    glutMouseFunc(mouseButton);
    glutMotionFunc(mouseMove);



    printf("Controls:\n");
    printf("ESC: Exit\n");
    printf("R: Reset Simulation\n");
    printf("Mouse Left Drag: Rotate Camera\n");
    printf("Mouse Right Drag: Zoom Camera\n");
    printf("Arrow Keys: Apply global force (hold to apply)\n"); // 更新提示


    glutMainLoop();

    return 0;
}