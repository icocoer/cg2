// light cube L button toggle
#include<GL/glut.h>
#include<math.h>
#include<stdbool.h>


#define glf GLfloat

bool light = 0;

GLfloat LightAmbient[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat LightDiffuse[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat LightPosition[] = { 0.0, 0.0, 2.0, 1.0 };


glf xrotate = 0.0;
glf yrotate = 0.0;
glf xspeed = 0.0;
glf yspeed = 0.0;
glf z = -5.0;

void changeSize(int w, int h)
{
    if (h == 0)
        h = 1;
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0, (GLfloat)w / (GLfloat)h, 1.0, 100.0);
    glMatrixMode(GL_MODELVIEW);
}

int initGL()
{
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    glClearColor(0.0, 0.0, 0.0, .5);
    glClearDepth(1.0);
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glLightfv(GL_LIGHT1, GL_AMBIENT, LightAmbient);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuse);
    glLightfv(GL_LIGHT1, GL_POSITION, LightPosition);
    glEnable(GL_LIGHT1);

    return 1;
}

void renderScene(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();

    glTranslatef(0.0, 0.0, z);

    glRotatef(xrotate, 1.0, 0.0, 0.0);
    glRotatef(yrotate, 0.0, 1.0, 0.0);

    glColor4f(1.0, .0, .0, 1.0);
    glutSolidCube(1.0);

    glutSwapBuffers();

    xrotate += xspeed;
    yrotate += yspeed;

}

void processNormalKeys(unsigned char key, int x, int y)
{   
    if(key == 27)
        exit(0);
    if(key == 'l'){
        light = !light;
        light ? glEnable(GL_LIGHTING) : glDisable(GL_LIGHTING);
    }    
    
}
void processSpecialKeys(int key, int x, int y)
{
    if(key == GLUT_KEY_PAGE_UP)
        z += 0.02;
    else if(key == GLUT_KEY_PAGE_DOWN)
        z -= 0.02;
    else if(key == GLUT_KEY_UP)
        xspeed += 0.01;
    else if(key == GLUT_KEY_DOWN)
        xspeed -= 0.01;
    else if(key == GLUT_KEY_LEFT)
        yspeed -= 0.01;
    else if(key == GLUT_KEY_RIGHT)
        yspeed += 0.01;
}

int main(int argc, char **argv)
{   
    glutInit(&argc, argv);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(500, 500);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("3D CUBE GLUT PROGRAM");
    initGL();
    glutDisplayFunc(renderScene);
    glutIdleFunc(renderScene);
    glutReshapeFunc(changeSize);
    glutSpecialFunc(processSpecialKeys);
    glutKeyboardFunc(processNormalKeys);
    glutMainLoop();
}