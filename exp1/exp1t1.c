// 配置glut并且显示一个窗口.

#include <GL/glut.h>

void renderScene() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glFlush();
}

void setupRC() {
    glClearColor(0.4f , 0.8f, 1.0f, 1.0f);
}
int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(500, 500);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE | GLUT_DEPTH);
    glutCreateWindow("FIRST GLUT PROGRAM");

    glutDisplayFunc(renderScene);
    setupRC();
    glutMainLoop();
    return 0;
}