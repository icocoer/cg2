// 通过键盘控制程序退出和三角形变色.

#include <GL/glut.h>

GLfloat angle = 0.0f;
GLfloat red = 1.0f,green = 0.0f,blue = 0.0f;
void renderScene() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glColor3f(red, green, blue);    

    glPushMatrix();
    glRotatef(angle, 0.0f, 1.0f, 0.0f);


    glBegin(GL_TRIANGLES);
        glVertex3f(-1.0f, -1.0f, 0.0f);
        glVertex3f(1.0f, 0.0f, 0.0f);   
        glVertex3f(0.0f, 1.0f, 0.0f);
    glEnd();

    glPopMatrix();
    glutSwapBuffers();
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

void timeFunction(int value) {
    angle += 2.0f;
    glutPostRedisplay();
    glutTimerFunc(23, timeFunction, 1);

}
void setupRC() {
    glClearColor(0.5f , 0.5f, 0.5f, 1.0f);
}

void processNormalKeys(unsigned char key, int x, int y){
    if (key == 27 ) exit(0);
}
void processSpecialKeys(int key, int x, int y) {
    switch (key)
    {
    case GLUT_KEY_F1:
        red = 1.0f;
        green = 0.0f;
        blue = 0.0f;break;
    case GLUT_KEY_F2:
        red = 0.0f;
        green = 1.0f;
        blue = 0.0f;break;
    case GLUT_KEY_F3:
        red = 0.0f;
        green = 0.0f;
        blue = 1.0f;break;
    
    default:
        break;
    }
}
int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(500, 500);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("KEY CONTROL GLUT PROGRAM");

    glutDisplayFunc(renderScene);
    glutReshapeFunc(reshapeSceneType2);
    glutTimerFunc(23, timeFunction, 1);

    glutKeyboardFunc(processNormalKeys);
    glutSpecialFunc(processSpecialKeys);

    setupRC();

    glutMainLoop();
    return 0;
}