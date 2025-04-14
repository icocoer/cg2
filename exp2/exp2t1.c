//t1 画一堆雪人

#include<GL/glut.h>
#include<stdio.h>
#include<math.h>

GLfloat _angle = 0.0f;
GLfloat _lx =0.0f, _lz = -1.0f;
GLfloat _x = 0.0f, _z = 5.0f;


void changeSize(GLsizei w, GLsizei h) {
    GLfloat faspect;

    if(h == 0) h = 1;

    glViewport(0, 0, w, h);

    faspect = (GLfloat)w / (GLfloat)h;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    gluPerspective(60.0f, faspect, 1.0f, 400.0f);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}
void drawSnowman() {
    glColor3f(1.0f, 1.0f, 1.0f);

    //body
    glTranslatef(0.0f, 0.75f, 0.0f);
    glutSolidSphere(0.75f, 20, 20);

    //head
    glTranslatef(0.0f,1.0f ,0.0f);
    glutSolidSphere(0.25f, 20, 20);

    //eye
    glPushMatrix();
        glColor3f(0.0f, 0.0f, 0.0f);
        glTranslatef(0.1f, 0.03f, 0.23f);
        glutSolidSphere(0.03f, 10, 10);

        glTranslatef(-0.2f, 0.0f, 0.0f);
        glutSolidSphere(0.03f, 10, 10);
    glPopMatrix();

    //nose
    glColor3f(1.0f, 0.0f, 0.0f);
    glTranslatef(0.0f, 0.0f, 0.25f);
    glutSolidCone(0.05f, 0.05f, 10, 2);
}

void renderScene(void) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();

    gluLookAt(_x, 1.0f, _z, 
              _x +_lx, 1.0f,_z + _lz,
              0.0f, 1.0f, 0.0f);

    glColor3f(0.7f, 0.7f, 0.7f);
    
    glBegin(GL_QUADS);
        glVertex3f(-100.0f, 0.0f, -100.0f);
        glVertex3f(-100.0f, 0.0f, 100.0f);  
        glVertex3f(100.0f, 0.0f, 100.0f);
        glVertex3f(100.0f, 0.0f, -100.0f);
    glEnd();
    float zoom = 8.0f;
    for(int i = -3; i < 3; i++) {
        for(int j = -3; j < 3; j++) {
            glPushMatrix();
            glTranslatef(i * zoom, 0.0f, j * zoom);
            drawSnowman();
            glPopMatrix();
        }
    }

    glutSwapBuffers();
}

void processSpecialKeys(int key, int xx, int yy) {
    float fraction = 0.1f;

    switch (key) {
        case GLUT_KEY_UP : 
            _x += _lx * fraction;
            _z += _lz * fraction;
            break;
        case GLUT_KEY_DOWN :
            _x -= _lx * fraction;
            _z -= _lz * fraction;
            break;
        case GLUT_KEY_LEFT :
            _angle -= 0.01f;
            _lx = sin(_angle);
            _lz = -cos(_angle);
            break;
        case GLUT_KEY_RIGHT :
            _angle += 0.01f;
            _lx = sin(_angle);
            _lz = -cos(_angle);
            break;
    }
}

int main(int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(800, 600);
    glutCreateWindow("Snowman");

    glutDisplayFunc(renderScene);
    glutReshapeFunc(changeSize);
    glutIdleFunc(renderScene);
    glutSpecialFunc(processSpecialKeys);

    glEnable(GL_DEPTH_TEST);

    glutMainLoop();
}