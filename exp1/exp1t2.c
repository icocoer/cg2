// 显示一个三角形.

#include <GL/glut.h>

void renderScene() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glBegin(GL_TRIANGLES);
        glVertex3f(-1.0f, -1.0f, 0.0f);
        glVertex3f(1.0f, 0.0f, 0.0f);   
        glVertex3f(0.0f, 1.0f, 0.0f);
    glEnd();

    glFlush();
}

void setupRC() {
    glClearColor(1.0f , 0.0f, 0.0f, 1.0f);
    glColor3ub(66, 132, 255);
}
int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(500, 500);
    glutInitDisplayMode(GLUT_RGB | GLUT_SINGLE | GLUT_DEPTH);
    glutCreateWindow("TRIANGLE GLUT PROGRAM");

    glutDisplayFunc(renderScene);
    setupRC();
    glutMainLoop();
    return 0;
}