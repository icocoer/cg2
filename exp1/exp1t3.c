// 可变窗口维持不形变.

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

void reshapeSceneType1(int w, int h) {
    GLfloat aspectRatio;

    if(h == 0) h = 1;
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    aspectRatio = (GLfloat)w / (GLfloat)h;
    if(w <= h) {
        glOrtho(-1.0f, 1.0f, -1.0f / aspectRatio, 1.0f / aspectRatio, -1.0f, 1.0f);
    } else {
        glOrtho(-1.0f * aspectRatio, 1.0f * aspectRatio, -1.0f, 1.0f, -1.0f, 1.0f);
    }
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void reshapeSceneType2(int w, int h) {
    GLfloat aspectRatio;
    if(h == 0) h = 1;
    
    aspectRatio = (GLfloat)w / (GLfloat)h;
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glViewport(0, 0, w, h);

    gluPerspective(45.0f, aspectRatio, 1.0f, 1000.0f);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0.0f, 0.0f, 5.0f, 0.0f, 0.0f, -1.0f, 0.0f, 1.0f, 0.0f);
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
    glutCreateWindow("RESIZE GLUT PROGRAM");

    glutDisplayFunc(renderScene);
    glutReshapeFunc(reshapeSceneType2);

    setupRC();
    glutMainLoop();
    return 0;
}