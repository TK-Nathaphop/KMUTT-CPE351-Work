#include <GL/freeglut.h> 
#include <stdlib.h> 
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#define screenX 1366
#define screenY 768
#define FILENAME "heatMatrixGL.txt"

#define MAX_ITERATION 1000


float r[screenX][screenY],g[screenX][screenY],b[screenX][screenY];
float * panel, *panelNext;
int row, col;
int iteration = 0;

void idle()
{                      

	for(int i=0;i<screenX;i++)
	{
		for(int j=0;j<screenY;j++)
		{
			//r[i][j]=(float)((rand() % 9))/8;
			//g[i][j]=(float)((rand() % 9))/8;
			//b[i][j]=(float)((rand() % 9))/8;
			r[i][j]=(float)0;
			g[i][j]=(float)0;
			b[i][j]=(float)0;
		}
	
	}	
	usleep(100000); //sleep 0.1 second
	glutPostRedisplay();

}

void magic_dots(void)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, screenX, 0.0, screenY);

	

	for(int i=0;i<screenX;i++)
	{
		for(int j=0;j<screenY;j++)
		{
			glColor3f(r[i][j],g[i][j],b[i][j]); 
			glBegin(GL_POINTS);
			glVertex2i (i,j);
			glEnd();
		}
	
	}		

	

	glFlush();	
}


void create_heatdot(void)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0.0, screenX, 0.0, screenY);
	int index = 0;
	float r,g,b;


	/*compute center of panel*/

	if(iteration < MAX_ITERATION)
	{
		iteration++;
		for (int ansRow = row - 2; ansRow-- ;)
			{
			int row1 = ansRow * col;
			int row2 = (ansRow + 1) * col ;
			int row3 = (ansRow + 2) * col ;
			for( int ansCol = col - 2; ansCol-- ;)
				{
				int rowcol1 = row1 + ansCol;
				int rowcol2 = row2 + ansCol;
				int rowcol3 = row3 + ansCol;

				float value = ( panel[rowcol1] + panel[rowcol1 + 1] + panel[rowcol1 + 2] 
							 + panel[rowcol2] + panel[rowcol2 + 1] + panel[rowcol2 + 2] 
							 + panel[rowcol3] + panel[rowcol3 + 1] + panel[rowcol3 + 2] )/9;
				panelNext[rowcol2 + 1] = value;

				value = 2 * (value-0) / (255 - 0);
	    			b = (1 - value);
				if(b < 0)
					b = 0;

	    			r = (value - 1);
				if(r < 0)
					r = 0;

	    			g = 1 - b - r;
	    		

				//value = (float) j / screenY;
				glColor3f(r, g, b); 
				glBegin(GL_POINTS);
				glVertex2i (ansRow,ansCol);
				glEnd();
				}
			}
		/*move next iteration panel to current panel*/
		float * temp = panel;
		panel = panelNext;
		panelNext = temp;

		
	}

	glFlush();
}
int main(int argc,char** argv)
{
	FILE * file = NULL;
	
	

	file = fopen(FILENAME,"r"); 
	fscanf(file, "%d %d", &row, &col);
	panel = malloc(row * col * sizeof(float));

	int i;
	i = 0;
	while (fscanf(file, "%f", &panel[i++]) == 1);
	fclose(file);
	
	panelNext = malloc(row * col * sizeof(float));
	memcpy( panelNext, panel, row*col*sizeof(float));
	

	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_SINGLE);
	glutInitWindowSize(screenX, screenY);
	glutCreateWindow("Heat transfer assignment1 (60070503420-3461)");
	glClearColor(0,0,0,0);
	glClear(GL_COLOR_BUFFER_BIT);
	glutDisplayFunc(create_heatdot);
	glutIdleFunc(idle);
	glutMainLoop();
	
	return 0;
}
