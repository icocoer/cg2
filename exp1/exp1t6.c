// 正方形回弹程序.

#include <GL/glut.h>

GLfloat x=0.0f,y=0.0f,rsize=25.0f;

GLfloat xstep = 1.0f , ystep = 1.0f;

GLfloat windowHeight,windowWidth;
void renderScene() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glColor3f(1.0 , 0.0 , 0.0);
    
    glRectf(x, y, x + rsize, y - rsize);

    glutSwapBuffers();
}

void timeFunction(int value) {
    if (x >= windowWidth - rsize || x <= -windowWidth)
        xstep = -xstep;
    
    if (y >= windowHeight|| y <= -windowHeight+ rsize )
        ystep = -ystep;

    x += xstep;
    y += ystep;

    if (x>(windowWidth - rsize + xstep)) x = windowWidth - rsize - 1;
        else if (x<(-windowWidth + xstep)) x = -windowWidth + 1;

    if (y>(windowHeight + ystep)) y = windowHeight- 1;
        else if (y<(-windowHeight + ystep)) y = -windowHeight +rsize - 1;

    glutPostRedisplay();
    glutTimerFunc(23, timeFunction, 1);

}
void setupRC() {
    glClearColor(0.0f , 0.0f, 1.0f, 1.0f);
}

void reshapeSize(int w, int h) {
    GLfloat ratio;

    if (h == 0) h = 1;

    glViewport(0, 0, w, h);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    ratio = (GLfloat)w / (GLfloat)h;

    if (w<=h){
        windowWidth = 100;
        windowHeight = 100 / ratio;
        glOrtho(-100, 100, -windowHeight, windowHeight, 1, -1);
    }else{
        windowWidth = 100 * ratio;
        windowHeight = 100;
        glOrtho(-windowWidth, windowWidth, -100, 100, 1, -1);
    }
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}
int main(int argc, char** argv) {
    glutInit(&argc, argv);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(800, 600);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
    glutCreateWindow("BOUNCE GLUT PROGRAM");

    glutDisplayFunc(renderScene);
    glutReshapeFunc(reshapeSize);
    glutTimerFunc(23, timeFunction, 1);

    setupRC();

    glutMainLoop();
    return 0;
}