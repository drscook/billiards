#ifdef __APPLE__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <random>
#include <cstring>

//opengl vis parameters
#define XWindowSize 1500
#define YWindowSize 1500
#define FAR 50.0
#define DT 0.005

//wall types for collision 
#define passive 0
#define heated  1
#define sliding 2

#define MIN(x, y) (x < y) ? x : y
#define MAX(x, y) (x > y) ? x : y

std:: random_device generator;

//physical parameters
int DIMENSION;
int N;
float MAX_CUBE_DIM;
float default_radius;
float default_mass;

//non-physical parameters
int STOP_TIME;
int steps_per_record;
int visualize;
int track_large_particle;
int ignore_particle_interaction;
int all_particles_diffused = 0;
int diff_number = 2;
int packing_scheme = 0;

//walls
float normal[6][3];    // normal is only ever 6 x 3, so ok to declare on stack
int WALL_TAG[6];       // wall tag also fixed at 6 long. Will not overflow stack
float WALL_TEMP[6][2]; // temperature for each wall
float alpha[6];        // prob of thermal effect for each wall
float max_temp;
float min_temp;

//physical constants
float BOLTZ_CONST = 1.0; //1.38064852e-23;
int max_walls = 3;  // number of walls that a particle can hit at a time

//file names for IO
char dir_name[256] = "\0./";
char *in_fname = NULL;

//particle
float **p     = NULL;
float **v     = NULL;
float *mass   = NULL;
float *radius = NULL;
float *max_square_displacement = NULL;
float *p_temp = NULL;
float default_p_temp;


//visualization
float EYE;

// memory management to make sure particles collide correctly
int *tag = NULL;
int *what_w_hit = NULL;
int *what_p_hit = NULL;
int *how_many_p = NULL;
int *how_many_w = NULL;
float *global_min_time = NULL;

//calculation variables
float total_entropy;
float total_heat_flow;
float gas_temp;
float wall_entropy[6];
float wall_collisions[6];
float *tot_dist = NULL;
float *p_collisions = NULL;
float *w_collisions_heated = NULL;
float *w_collisions_passive = NULL;
float slice_width;
int n_slices = 10;
int n_in_slice[10];
float t_in_slice[10];

// reads input file of parameters (temperature of walls, particle size, &c)
void set_macros()
{
  	FILE * fp = NULL;
	const int bdim = 132;
	char buff[bdim];
	int i;
	int d;
	double f, g;
	char s[256];

	//Default values
	//physical parameters
	DIMENSION = 3;
	N = 25;
	MAX_CUBE_DIM = 4.0;
	default_radius = 0.1;
	default_mass = 2.0;

	//non-physical parameters
	STOP_TIME = 1000;
	steps_per_record = 50;
	visualize = false;
	track_large_particle = false;
	ignore_particle_interaction = false;

	//walls
	WALL_TAG[0] = WALL_TAG[1] = heated;
	WALL_TEMP[0][0] = WALL_TEMP[0][1] = 45.0;
	WALL_TEMP[1][0] = WALL_TEMP[1][1] = 90.0;
	WALL_TAG[2] = WALL_TAG[3] = WALL_TAG[4] = WALL_TAG[5] = passive;
	WALL_TEMP[2][0] = WALL_TEMP[3][0] = WALL_TEMP[4][0] = WALL_TEMP[5][0] = 0.0;
	WALL_TEMP[2][1] = WALL_TEMP[3][1] = WALL_TEMP[4][1] = WALL_TEMP[5][1] = 0.0;

	
	//Read from input file & replace default values
	if( (fp = fopen(in_fname,"r")) == NULL)
	{
		printf("No input file.  Using default values.");
	}
	else
	{
		fgets(buff,bdim,fp);
		fgets(buff,bdim,fp);
		sscanf(buff, "%d", &d);
		DIMENSION = d;

		fgets(buff,bdim,fp);
		fgets(buff,bdim,fp);
		sscanf(buff, "%d", &d);
		visualize = d;
		
		fgets(buff, bdim, fp);
		fgets(buff, bdim, fp);
		sscanf(buff, "%d", &d);
		N = d;

		fgets(buff, bdim, fp);
		fgets(buff, bdim, fp);
		sscanf(buff, "%lf", &f);
		default_radius = f;
		if(default_radius > 0)
		{
		}
		else
		{
			ignore_particle_interaction = true;
		}

		fgets(buff,bdim,fp);
		fgets(buff,bdim,fp);
		sscanf(buff, "%d", &d);
		packing_scheme = d;

		fgets(buff,bdim,fp);
		fgets(buff,bdim,fp);
		sscanf(buff, "%d", &d);
		STOP_TIME = d;

		fgets(buff, bdim, fp);
		for(i = 0; i < 6; i++)
		{
			fgets(buff, bdim, fp);
			sscanf(buff, "%d %lf %lf", &d, &f, &g);
			WALL_TAG[i] = d;
			WALL_TEMP[i][0] = f;
			WALL_TEMP[i][1] = g;
		}

		fgets(buff, bdim, fp);
		fgets(buff, bdim, fp);
		sscanf(buff,"%d", &d);
		track_large_particle = d;

		fgets(buff, bdim, fp);
		fgets(buff, bdim, fp);
		sscanf(buff, "%s", &s);
		strcpy(dir_name, s);

		fgets(buff, bdim, fp);
		fgets(buff, bdim, fp);
		sscanf(buff, "%d", &d);
		steps_per_record = d;
		fclose(fp);
	}
}

//get overall temperature of gas based on individual particles' kinetic energy
void compute_T()
{
	int i;
	//float total_temp = 0.0;
	gas_temp = 0.0;

	for(i = 0; i < N; i++)
	{
		p_temp[i] = mass[i] * (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) / (1.0 * DIMENSION * BOLTZ_CONST);

		gas_temp += p_temp[i];
	}
	gas_temp /= (1.0*N);
	//KE /= (1.0 * N);

	//gas_temp = KE / (1.0 * DIMENSION * BOLTZ_CONST);
}

float wall_temperature(float x, int wall)
{
	return (x + MAX_CUBE_DIM) * (WALL_TEMP[wall][0] + (WALL_TEMP[wall][1] - WALL_TEMP[wall][0]) / (2.0 * MAX_CUBE_DIM));
}

float get_intersection_point(float time, int i)
{
	return p[i][0] + v[i][0] * time;
}


//initialize particles with position, velcoity, radius, etc. 
void pack_particles()
{
	int i, j, k, num, particles_per_side;
	float max_radius, hole_radius, temp;
	float T, tol = 1.0e-7;
	float x_length, x_start;
	float y_length, y_start;
	float z_length, z_start;

	//tag particles not to hit themselves
	for(i = 0;  i < max_walls * N; i++) tag[i] = i / max_walls;
	
	//set initial particle parameters
	for (i = 0; i < N; i++)
	{
		mass[i] = default_mass;
		radius[i] = default_radius;
		p_temp[i] = default_p_temp;
	}

	if (track_large_particle)
	{
		radius[0] = 3.0 * radius[0];
		mass[0] = 3.0 * mass[0];
	}
	max_radius = radius[0];

	if(packing_scheme < 1)
	{
		//square lattice
		hole_radius = max_radius * 1.1;
		temp = pow((float)N, 1.0 / (1.0 * DIMENSION)) + .99999;
		particles_per_side = (int)temp;
		x_length = hole_radius * (2.0 * particles_per_side);
		y_length = x_length;
		z_length = x_length;
		if (x_length > 2.0 * MAX_CUBE_DIM)
		{
			printf("Box not big enough to hold all the gas particles.");
			exit(1);
		}
		x_start = -x_length / 2.0 + hole_radius;
		y_start = -y_length / 2.0 + hole_radius;
		z_start = -z_length / 2.0 + hole_radius;
		num = 0;
		for (k = 0; k < particles_per_side; k++)
		{
			for (j = 0; j < particles_per_side; j++)
			{
				for (i = 0; i < particles_per_side; i++)
				{
					if (N <= num) break;
					p[num][0] = x_start + 2.0*hole_radius*i;
					p[num][1] = y_start + 2.0*hole_radius*j;
					p[num][2] = z_start + 2.0*hole_radius*k;
					num++;
				}
			}
		}
	}
	else
	{
		//tetrahedral lattice
		hole_radius = max_radius * 1.1;
		temp = pow((float)N, 1.0 / (1.0 * DIMENSION)) + .99999;
		particles_per_side = (int)temp;
		x_length = hole_radius * (2.0 * particles_per_side + 1.0);
		y_length = x_length;
		z_length = hole_radius * (sqrt(2) * (particles_per_side - 1) + 2.0);
		if (x_length > MAX_CUBE_DIM)
		{
			printf("Box not big enough to hold all the gas particles.");
			exit(1);
		}
		x_start = -x_length / 2.0 + hole_radius;
		y_start = -y_length / 2.0 + hole_radius;
		z_start = -z_length / 2.0 + hole_radius;
		num = 0;
		for (k = 0; k < particles_per_side; k++)
		{
			for (j = 0; j < particles_per_side; j++)
			{
				for (i = 0; i < particles_per_side; i++)
				{
					if (N <= num) break;
					p[num][0] = x_start + 2.0*hole_radius*i + (k % 2)*hole_radius;
					p[num][1] = y_start + 2.0*hole_radius*j + (k % 2)*hole_radius;
					p[num][2] = z_start + sqrt(2.0)*hole_radius*k;
					num++;
				}
			}
		}
	}

	if(DIMENSION < 3)
	{
		for(i = 0; i < N; i++)
		{
			p[i][2] = v[i][2] = 0.0;
		}
	}

	for (i = 0; i < N; i++)
	{
		max_square_displacement[i] = (MAX_CUBE_DIM - radius[i] + tol) * (MAX_CUBE_DIM - radius[i] + tol);
	}

	for (i = 0; i < N; i++)
	{
		T = sqrt(BOLTZ_CONST * p_temp[i] / mass[i]);
		std::normal_distribution<double> normal_dist(0, T);

		v[i][0] = normal_dist(generator);
		v[i][1] = normal_dist(generator);
		v[i][2] = normal_dist(generator);
		tot_dist[i] = 0.0;
		p_collisions[i] = w_collisions_heated[i] = w_collisions_passive[i] = 0.0;
	}
}

// initialize system with N particles in a box 
void set_initail_conditions()
{
	int i;

	set_macros();

	EYE = 2.0 * MAX_CUBE_DIM;

	// set up parameters for the box
	max_temp = min_temp = 0.0;
	for (i = 0; i < 6; i++)
	{
		max_temp = MAX(max_temp, WALL_TEMP[i][0]);
		min_temp = MAX(min_temp, WALL_TEMP[i][0]);
	}
	default_p_temp = (max_temp + min_temp) / 2.0;

	// set entropy to 0--nothing has accumulated yet
	for(i = 0; i < 6; i++){ wall_entropy[i] = wall_collisions[i] = 0; alpha[i] = 1.0;}

	// normal vector for each wall of the cube
	normal[0][0] = 1.0; normal[0][1] = 0.0; normal[0][2] = 0.0;
	normal[2][0] = 0.0; normal[2][1] = 1.0; normal[2][2] = 0.0;
	normal[4][0] = 0.0; normal[4][1] = 0.0; normal[4][2] = 1.0;
	normal[1][0] =-1.0; normal[1][1] = 0.0; normal[1][2] = 0.0;
	normal[3][0] = 0.0; normal[3][1] =-1.0; normal[3][2] = 0.0;
	normal[5][0] = 0.0; normal[5][1] = 0.0; normal[5][2] =-1.0;

	// allocate memory for particle position, velocity, mass, radius, and temperature
	p      = (float**)malloc(N * sizeof(float*));
	v      = (float**)malloc(N * sizeof(float*));
	mass   = (float* )malloc(N * sizeof(float ));
	radius = (float* )malloc(N * sizeof(float ));
	max_square_displacement = (float*)malloc(N * sizeof(float));
	p_temp = (float* )malloc(N * sizeof(float ));	 
	for(i = 0; i < N; i++)
	{
		p[i] = (float*)malloc(3 * sizeof(float));
		v[i] = (float*)malloc(3 * sizeof(float));
	}


	//variables for keeping track of mean free path
	tot_dist             = (float* )malloc(N * sizeof(float));
	p_collisions         = (float* )malloc(N * sizeof(float));
	w_collisions_heated  = (float* )malloc(N * sizeof(float));
	w_collisions_passive = (float* )malloc(N * sizeof(float));
	for(i = 0; i < N; i++)
	{
		w_collisions_heated[i]  = w_collisions_passive[i] = p_collisions[i] = tot_dist[i] = 0.0;
	}
	// memory management variables--make sure all particles collide correctly
	global_min_time = (float*)malloc(            N * sizeof(float));
	tag             = (int*  )malloc(max_walls * N * sizeof(int  ));
	what_w_hit      = (int*  )malloc(max_walls * N * sizeof(int  ));
	what_p_hit      = (int*  )malloc(            N * sizeof(int  ));
	how_many_p      = (int*  )malloc(            N * sizeof(int  ));
	how_many_w      = (int*  )malloc(            N * sizeof(int  ));

	pack_particles();
	compute_T();

	slice_width = (2.0 * MAX_CUBE_DIM) / (1.0 * n_slices );
}

//output OpenGL visualization
void draw_picture()
{
	int i;
	float red, green, blue, hi, mid;
	
	glClear(GL_COLOR_BUFFER_BIT);
	glClear(GL_DEPTH_BUFFER_BIT);
	
	glColor3d(1.0,1.0,0.5);
	for(i = 0; i < N; i++)
	{
		hi = MIN(p_temp[i] / max_temp, 1.7);
		mid = 0.75 - hi;
		red = hi * hi;
		green = 1.0 - mid * mid * mid * mid;
		blue = 1.5 - red;

		glColor3f(red, green, blue);
		glPushMatrix();
		glTranslatef(p[i][0], p[i][1], p[i][2]);
		glutSolidSphere(radius[i],20,20);
		glPopMatrix();
	}

	if(DIMENSION > 2)
	{
		glBegin(GL_QUADS);
		glColor3f(0.8, 0.8, 0.8);
		glVertex3f(-MAX_CUBE_DIM, -MAX_CUBE_DIM, -MAX_CUBE_DIM);
		glVertex3f(-MAX_CUBE_DIM,  MAX_CUBE_DIM, -MAX_CUBE_DIM);
		glVertex3f( MAX_CUBE_DIM,  MAX_CUBE_DIM, -MAX_CUBE_DIM);
		glVertex3f( MAX_CUBE_DIM, -MAX_CUBE_DIM, -MAX_CUBE_DIM);

		glColor3f(1.0, 1.0, 1.0);
		glVertex3f(-MAX_CUBE_DIM, -MAX_CUBE_DIM, -MAX_CUBE_DIM);
		glVertex3f( MAX_CUBE_DIM, -MAX_CUBE_DIM, -MAX_CUBE_DIM);
		glVertex3f( MAX_CUBE_DIM, -MAX_CUBE_DIM,  MAX_CUBE_DIM);
		glVertex3f(-MAX_CUBE_DIM, -MAX_CUBE_DIM,  MAX_CUBE_DIM);

		glColor3f(1.0, 1.0, 1.0);
		glVertex3f(-MAX_CUBE_DIM,  MAX_CUBE_DIM, -MAX_CUBE_DIM);
		glVertex3f( MAX_CUBE_DIM,  MAX_CUBE_DIM, -MAX_CUBE_DIM);
		glVertex3f( MAX_CUBE_DIM,  MAX_CUBE_DIM,  MAX_CUBE_DIM);
		glVertex3f(-MAX_CUBE_DIM,  MAX_CUBE_DIM,  MAX_CUBE_DIM);

		glColor3f(0.0, 1.0, 1.0);
		glVertex3f(-MAX_CUBE_DIM, -MAX_CUBE_DIM, -MAX_CUBE_DIM);
		glVertex3f(-MAX_CUBE_DIM, -MAX_CUBE_DIM,  MAX_CUBE_DIM);
		glVertex3f(-MAX_CUBE_DIM,  MAX_CUBE_DIM,  MAX_CUBE_DIM);
		glVertex3f(-MAX_CUBE_DIM,  MAX_CUBE_DIM, -MAX_CUBE_DIM);

		glColor3f(1.0, 0.15, 0.0);
		glVertex3f( MAX_CUBE_DIM, -MAX_CUBE_DIM, -MAX_CUBE_DIM);
		glVertex3f( MAX_CUBE_DIM, -MAX_CUBE_DIM,  MAX_CUBE_DIM);
		glVertex3f( MAX_CUBE_DIM,  MAX_CUBE_DIM,  MAX_CUBE_DIM);
		glVertex3f( MAX_CUBE_DIM,  MAX_CUBE_DIM, -MAX_CUBE_DIM);
		glEnd();
	}
	else
	{
		glBegin(GL_LINES);
		glColor3f(0.0, 1.0, 1.0);
		glVertex3f(-MAX_CUBE_DIM, -MAX_CUBE_DIM, 0.0);
		glVertex3f( MAX_CUBE_DIM, -MAX_CUBE_DIM, 0.0);

		glVertex3f( MAX_CUBE_DIM, -MAX_CUBE_DIM, 0.0);
		glVertex3f( MAX_CUBE_DIM,  MAX_CUBE_DIM, 0.0);

		glVertex3f( MAX_CUBE_DIM,  MAX_CUBE_DIM, 0.0);
		glVertex3f(-MAX_CUBE_DIM,  MAX_CUBE_DIM, 0.0);

		glVertex3f(-MAX_CUBE_DIM,  MAX_CUBE_DIM, 0.0);
		glVertex3f(-MAX_CUBE_DIM, -MAX_CUBE_DIM, 0.0);

		glEnd();
	}
	glutSwapBuffers();
}

// check if a given particle (i0) will hit a wall (i1), and if so at what time (t)
int particle_wall_collision(int i0, int i1, float * t)
{

	int collides = 0;
	float tt, max_cube_minus_radius;

	max_cube_minus_radius = MAX_CUBE_DIM - radius[i0];
	tt = -1.0;

	if     ( (i1 == 0) && (v[i0][0] * v[i0][0] > 0.0) )
	{
		tt = ( (-max_cube_minus_radius) - p[i0][0] ) / v[i0][0];
	}
	else if( (i1 == 2) && (v[i0][1] * v[i0][1] > 0.0) )
	{
		tt = ( (-max_cube_minus_radius) - p[i0][1] ) / v[i0][1];
	}
	else if( (i1 == 4) && (v[i0][2] * v[i0][2] > 0.0) )
	{
		tt = ( (-max_cube_minus_radius) - p[i0][2] ) / v[i0][2];
	}
	else if( (i1 == 1) && (v[i0][0] * v[i0][0] > 0.0) )
	{
		tt = ( max_cube_minus_radius - p[i0][0] ) / v[i0][0];
	}
	else if( (i1 == 3) && (v[i0][1] * v[i0][1] > 0.0) )
	{
		tt = ( max_cube_minus_radius - p[i0][1] ) / v[i0][1];
	}
	else if( (i1 == 5) && (v[i0][2] * v[i0][2] > 0.0) )
	{
		tt = ( max_cube_minus_radius - p[i0][2] ) / v[i0][2];
	}
	if( tt >= 0.0)
	{
	  	(*t) = tt;
		collides = 1;
	}

	return collides;
}

// Rayleigh distribution. Requires uniform[0,1] distributed variable and 
// a standard deviation
float gen_chi_sq(float u, float sigma)
{
	return sigma * sqrt(fabs(2.0 * log(1.0 - u)));
}

// implements cosine law to reflect particle off of heated wall 
void heated_wall(float * v_in, float T, float m, float * n, float * v_out)
{
	float t1[3], t2[3];
	float u, v, w, len, sigma = sqrt(BOLTZ_CONST * T / m);
	int ok = 0;

	std::normal_distribution<double> normal_dist(0, sigma);
	std::uniform_real_distribution<double> uniform_dist(0, 1.0);

	u = uniform_dist(generator);
	u = gen_chi_sq(u, sigma);
	v = normal_dist(generator);
	w = normal_dist(generator);

	//get vectors tangent to wall (perpendicular to n)
	// t1 is just a random vector perpendicular to n
	if( (n[0] * n[0] > 0) || (n[1] * n[1] > 0) )
	{
		t1[0] = -n[1];
		t1[1] =  n[0];
		t1[2] =  0.0;
		ok = 1;
	}
	else if(n[2] * n[2] > 0.0)
	{
		t1[0] = -n[2];
		t1[1] =  0.0;
		t1[2] =  n[0];
		ok = 1;
	}

	if(ok)
	{
	 	// normalize t1
	  	len = t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2];
		len = (len > 0.0) ? sqrt(len) : 1.0;
		t1[0] /= len; t1[1] /= len; t1[2] /= len;

		// now t2 is the cross product n x t1
		t2[0] = n[1] * t1[2] - n[2] * t1[1];
		t2[1] = n[2] * t1[0] - n[0] * t1[2];
		t2[2] = n[0] * t1[1] - n[1] * t1[0];

		// normalize t2
	  	len = t2[0] * t2[0] + t2[1] * t2[1] + t2[2] * t2[2];
		len = (len > 0.0) ? sqrt(len) : 1.0;
		t2[0] /= len; t2[1] /= len; t2[2] /= len;

		// compute outward going trajectory of velocity
		v_out[0] = (u * n[0] + v * t1[0] + w * t2[0]);
		v_out[1] = (u * n[1] + v * t1[1] + w * t2[1]);
		if(DIMENSION > 2) v_out[2] = (u * n[2] + v * t1[2] + w * t2[2]);
	}
	else
	{
		v_out[0] = v_out[1] = v_out[2] = 0.0;
	}
}

// particle reflection off of non-heated walls
void specular_reflect(float * v_in, float * n, float * v_out)
{
	float dotproduct;
  	dotproduct = -(v_in[0] * n[0] + v_in[1] * n[1] + v_in[2] * n[2]);
	v_out[0] = 2.0 * dotproduct * n[0] + v_in[0];
	v_out[1] = 2.0 * dotproduct * n[1] + v_in[1];
	v_out[2] = 2.0 * dotproduct * n[2] + v_in[2];
}

// checks if two particles (i0 and i1) will collide, and at what time (t)
int particle_particle_collision(int i0, int i1, float * t)
{
	float xc, yc, zc, xv, yv, zv;
	float discriminant, dd, a, b, c, dt = -1.0;
	int collides = 0;

	dd = (radius[i0] + radius[i1]) * (radius[i0] + radius[i1]);

	xc = p[i0][0] - p[i1][0];
	yc = p[i0][1] - p[i1][1];
	zc = p[i0][2] - p[i1][2];

	if( (xc * xc + yc * yc + zc * zc) > dd){

		xv = v[i0][0] - v[i1][0];
		yv = v[i0][1] - v[i1][1];
		zv = v[i0][2] - v[i1][2];

		a = xv * xv + yv * yv + zv * zv;
		b = 2.0 * (xc * xv + yc * yv + zc * zv);
		c = (xc * xc) + (yc * yc) + (zc * zc) - dd;
		discriminant = b * b - (4.0 * a * c);

		if(discriminant >= 0.0)
		{
			if(a * a > 0.0) // solve ax^2 + bx + c = 0
			{
				// choose the smallest positive root
				dt = (-b - sqrt(discriminant)) / (2.0 * a);
				if(dt < 0.0)
				{
					dt = (-b + sqrt(discriminant)) / (2.0 * a);
				}
			}
			else if(b * b > 0) // solve bx + c = 0 
			{
				dt = -c / b;
			}
			else if(c * c > 0)
			{
				dt = 0.0;
			}
			else
			{
				dt = -1.;
			}
			if(dt >= 0.0)
			{
				collides = 1;
				*t = dt;
			}
		}
	}
	return collides;
}


void resolve_particle_collision(int i0)
{
	float t[3], x_v1[3], y_v1[3], x_v2[3], y_v2[3], n[3], v0[3];
	float dt, len, m1m2;
	float ke_in, ke_out, entropy, u, temp_at_point, point;
	int wall, i, i1;
	std::uniform_real_distribution<double> u_dist(0, 1.0);

	if( how_many_w[i0] > 0) // if particle hits a wall
	{
	  	how_many_p[i0] = 0;
		v0[0] = v[i0][0]; v0[1] = v[i0][1]; v0[2] = v[i0][2];

		if(how_many_w[i0] > 1) // if hits multiple walls 
		{
		  	w_collisions_passive[i0] += 1.0;

			n[0] = n[1] = n[2] = 0.0;
			for(i = 0; i < how_many_w[i0]; i++)
			{
				wall = -(1 + what_w_hit[max_walls * i0 + i]);
				tag[max_walls * i0 + i] = what_w_hit[max_walls * i0 + i];
	
				n[0] += normal[wall][0];
				n[1] += normal[wall][1];
				n[2] += normal[wall][2];
			}

			len = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
			len = (len > 0) ? len : 1.0;
			n[0] /= len; n[1] /= len; n[2] /= len;
	
			specular_reflect(v0, n, v[i0]);
		}
		else 
		{
			tag[max_walls * i0] = what_w_hit[max_walls * i0];
			wall = -(1 + what_w_hit[max_walls * i0]);

			switch( WALL_TAG[wall] ){
				case passive:
		  			w_collisions_passive[i0] += 1.0;
					specular_reflect(v0, normal[wall], v[i0]);
					break;
				case heated: 
		  			w_collisions_heated[i0] += 1.0;
					if( (w_collisions_heated[i0] < (diff_number + 0.5)) && 
					    (w_collisions_heated[i0] > (diff_number - 0.5)) ) all_particles_diffused++;
					u = u_dist(generator);

					if(u > alpha[wall])
					{
						specular_reflect(v0, normal[wall], v[i0]);
					}
					else
					{
//					  	point = get_intersection_point(time, i0);
//						temp_at_point = wall_temperature(point, wall);
					  	temp_at_point = WALL_TEMP[wall][0];
						ke_in = 0.5 * mass[i0] * (v[i0][0] * v[i0][0] + v[i0][1] * v[i0][1] + v[i0][2] * v[i0][2]);

						heated_wall(v0, temp_at_point, mass[i0], normal[wall], v[i0]);

						ke_out = 0.5 * mass[i0] * (v[i0][0] * v[i0][0] + v[i0][1] * v[i0][1] + v[i0][2] * v[i0][2]);

						if(wall == 1) total_heat_flow += (ke_out - ke_in);
						wall_entropy[wall] += (ke_out - ke_in) / temp_at_point;
						wall_collisions[wall] += 1.0;
						total_entropy += (ke_in - ke_out) / temp_at_point;
					}

					break;
				case sliding:
					break;
				default:
					break;
			}
		}
		for(i = how_many_w[i0]; i < max_walls; i++) tag[max_walls * i0 + i] = N;
	}
	else if(how_many_p[i0] > 0) // if hits a particle
	{
		i1 = what_p_hit[i0];

		if(i1 < N)
		{
			tag[max_walls * i0] = i1;
		  	tag[max_walls * i1] = i0;
			for(i = 1; i < max_walls; i++) tag[max_walls * i1 + i] = tag[max_walls * i0 + i] = N;

			m1m2 = mass[i0] + mass[i1];
			m1m2 = (m1m2 * m1m2 > 0) ? m1m2 : 1.0;

			// t is the vector from center of i0 to center of i1
			t[0] = p[i0][0] - p[i1][0];
			t[1] = p[i0][1] - p[i1][1];
			t[2] = p[i0][2] - p[i1][2];

			// normalize t
			len = sqrt(t[0] * t[0] + t[1] * t[1] + t[2] * t[2]);
			len = (len > 0) ? len : 1.0;
			t[0] /= len;
			t[1] /= len;
			t[2] /= len;


			//--------------------Particle 1-----------------------//
			// get vector component of v[i0] parallel to t (x_v1)
			len = t[0] * v[i0][0] + t[1] * v[i0][1] + t[2] * v[i0][2];
			x_v1[0] = len * t[0];
			x_v1[1] = len * t[1];
			x_v1[2] = len * t[2];
			// get vector component of v[i0] perpendicular to t (y_v1)
			y_v1[0] = v[i0][0] - x_v1[0];
			y_v1[1] = v[i0][1] - x_v1[1];
			y_v1[2] = v[i0][2] - x_v1[2];
			//-----------------------------------------------------//


			//--------------------Particle 2-----------------------//
			// get vector component of v[i1] parallel to t (x_v2)
			t[0] = -t[0]; t[1] = -t[1]; t[2] = -t[2];
			len = t[0] * v[i1][0] + t[1] * v[i1][1] + t[2] * v[i1][2];
			x_v2[0] = len * t[0];
			x_v2[1] = len * t[1];
			x_v2[2] = len * t[2];
			// get vector component of v[i1] perpendicular to t (y_v2)
			y_v2[0] = v[i1][0] - x_v2[0];
			y_v2[1] = v[i1][1] - x_v2[1];
			y_v2[2] = v[i1][2] - x_v2[2];
			//-----------------------------------------------------//

			// update vidx0 and vidx1 to hold new velocities
			v[i0][0] = (x_v1[0] * (mass[i0] - mass[i1]) + x_v2[0] * 2.0 * mass[i1]) / m1m2 + y_v1[0];
			v[i0][1] = (x_v1[1] * (mass[i0] - mass[i1]) + x_v2[1] * 2.0 * mass[i1]) / m1m2 + y_v1[1];
			v[i0][2] = (x_v1[2] * (mass[i0] - mass[i1]) + x_v2[2] * 2.0 * mass[i1]) / m1m2 + y_v1[2];
                                                                                               
			v[i1][0] = (x_v2[0] * (mass[i1] - mass[i0]) + x_v1[0] * 2.0 * mass[i0]) / m1m2 + y_v2[0];
			v[i1][1] = (x_v2[1] * (mass[i1] - mass[i0]) + x_v1[1] * 2.0 * mass[i0]) / m1m2 + y_v2[1];
			v[i1][2] = (x_v2[2] * (mass[i1] - mass[i0]) + x_v1[2] * 2.0 * mass[i0]) / m1m2 + y_v2[2];

			how_many_p[i1] = how_many_w[i1] = 0;
		}
	}
}

void smooth_vis(float time)
{
	float dt;
	int i, j, num;
	float **new_pos = (float **)malloc(N * sizeof(float *));
	for (i = 0; i<N; i++)
		new_pos[i] = (float *)malloc(3 * sizeof(float));
	//float new_pos[N][3];

	dt = DT;
	num = (int)((time / dt) + 0.5);

	for(i = 0; i < N; i++)
	{
		new_pos[i][0] = p[i][0];
		new_pos[i][1] = p[i][1];
		new_pos[i][2] = p[i][2];
	}

	dt = time / (1.0 * (num - 1.0));
	if((num - 1) < 1){dt = time;}

	for(i = 0; i < num; i++)
	{
		for(j = 0; j < N; j++)
		{
			p[j][0] = new_pos[j][0] + (v[j][0] * dt * i);
			p[j][1] = new_pos[j][1] + (v[j][1] * dt * i);
			p[j][2] = new_pos[j][2] + (v[j][2] * dt * i);
		}
		draw_picture();
	}

	for(i = 0; i < N; i++)
	{
		p[i][0] = new_pos[i][0];
		p[i][1] = new_pos[i][1];
		p[i][2] = new_pos[i][2];
	}
}

int find_global_min_dt(float * minimum_time)
{
	float t, min, timesteptol = 1.0e-4;
  	int * wall_ct = NULL;
	int p1, p2, w, i;
	int p_first_collision, first_collision, collides;
	first_collision = 1;

	for(p1 = 0; p1 < N; p1++)
	{
		p_first_collision = 1;
		global_min_time[p1] = -1;

		how_many_p[p1] = how_many_w[p1] = 0;
		what_p_hit[p1] = N;
		for(i = 0; i < max_walls; i++) what_w_hit[max_walls * p1 + i] = N;

	 	// check current particle against all others for collision
		if(ignore_particle_interaction)
		{
		  	// do nothing. We're only looking at walls
		}
		else
		{
			for(p2 = 0; p2 < N; p2++)
			{
			  	//make sure this particle isn't current particle, and isn't one that was just hit.
				if( (p1 != p2) && ( (tag[max_walls * p1] != p2) || (tag[max_walls * p2] != p1) ) )
				{
					collides = particle_particle_collision(p1, p2, &t);
					if(collides > 0)
					{
						if( (t < global_min_time[p1]) || (p_first_collision > 0) )
						{
							p_first_collision = 0;
							global_min_time[p1] = t;
							how_many_p[p1] = 1;
							what_p_hit[p1] = p2;
							min = t;
						}
						else if( t <= global_min_time[p1] )
						{
							global_min_time[p1] = t;
							how_many_p[p1]++;
						}
					}
				}
			}
		}
	 	// check current particle against walls for collision
		for(w = 0; w < 6; w++)
		{
			p2 = 0; // make sure it's not a wall hit in last time step
			for(i = 0; i < max_walls; i++)
			{
				if( tag[max_walls * p1 + i] == -(w + 1) )
				{
					p2 = 1;
				}
			}
			if( p2 < 1 ) 
			{
				collides = particle_wall_collision(p1, w, &t);
				if(collides > 0)
				{
					if( (t < global_min_time[p1]) || (p_first_collision > 0) )
					{
					  	p_first_collision = 0;
						global_min_time[p1] = t;

						how_many_w[p1] = 1;
						what_w_hit[max_walls * p1] = -(w + 1);

						how_many_p[p1] = 0;
						what_p_hit[p1] = N;

						min = t;
					}
					else if(t <= global_min_time[p1] )
					{
						what_w_hit[max_walls * p1 + how_many_w[p1] ] = -(w + 1);
						how_many_w[p1]++;
					}
				}
			}
		}
		if(p_first_collision < 1) first_collision = 0;
	}
	// get the mininmum timestep computed amongst all particles that did hit something.
	for(p1 = 0; p1 < N; p1++) min = ( (how_many_p[p1] > 0) || (how_many_w[p1] > 0) ) ? MIN(min, global_min_time[p1]) : min;

	// set size of new timestep for all particles
	for(p1 = 0; p1 < N; p1++)
	{
		//if particle doesn't hit anything by minimum timestep, set its time to the minimum
		if( ((how_many_w[p1] < 1) && (how_many_p[p1] < 1)) || (global_min_time[p1] > (min + timesteptol)) )
		{
			global_min_time[p1] = min;
			how_many_w[p1] = 0;
			how_many_p[p1] = 0;
		}
	}

	// set global minimum time for all particles
	(*minimum_time) = min;

	return first_collision;
}

// check to see if multiple particles and/or walls are involved in a single collisions
// n.b. seems to be some dependence upon timing for last if statement and for loop in 
//      this routine: global minimum time step requires separate treatment. 
void check_complex_collisions(float * min_time)
{
	int i, type_1, type_2;
	float t = global_min_time[0];

	type_1 = N;
	type_2 = N;

	for(i = 0; i < N; i++)
	{
		if( (how_many_p[i] > 0) && (how_many_w[i] > 0) )
		{
			type_1 = i;
		}
		else if( how_many_p[i] > 0 )
		{
			type_2 = i;
		}
	}
	type_2 = (type_1 < N) ? N : type_2;
	if( type_1 < N ) how_many_p[type_1] = 0;

	if( (type_1 < N) || (type_2 < N) )
	{
	  	t = 0.95 * t;
		for(i = 0; i < N; i++)
		{
			if( (i != type_1) && (i != type_2) )
			{
				how_many_p[i] = how_many_w[i] = 0;
			}
		}
	}
	if( (type_1 < N) || (type_2 < N) ) (*min_time) = 0.95 * (*min_time);

	for(i = 0; i < N; i++) 
	{
		if( ( i != type_1 ) && ( i != type_2 ) )
		{
			global_min_time[i] = (*min_time);
		}
	}

}

void add_to_slices()
{
	int i, j;
	float x;
	for(i = 0; i < n_slices; i++) n_in_slice[i] = 0;
	for(i = 0; i < n_slices; i++) t_in_slice[i] = 0.0;

	for(i = 0; i < N; i++)
	{
		x = p[i][0];
		for(j = 0; j < n_slices; j++)
		{
			if (((-MAX_CUBE_DIM + j*slice_width) <= x) && (x <= (-MAX_CUBE_DIM + (j+1)*slice_width)))
			{
				t_in_slice[j] += p_temp[i];
				n_in_slice[j]++;
				break;
			}
		}
	}
}


void n_body()
{
	float fixed_dt = 0.0, time_elapsed;
	int burn_in_period = 1;
	FILE *param_file, *pos_file, *temp_file, *slice_file, *ent_file, *mfp_file;
	char dir[256];
	int i, j, time, first_collision;	
	float min_time;
	float tt = 0.0;
	total_entropy = 0.0;
	total_heat_flow = 0.0;
	time = 0;

	//*/ 
	// * * * * * * * output file stuff * * * * * * * //
	param_file = fopen(strcat(strcpy(dir, dir_name),"parameters.txt"), "w");
	fprintf(param_file, "dimension num_gas box_half_length radius temp1 temp2 temp3 temp4 temp5 temp6\n");
	fprintf(param_file, "%d %li %lf %lf %lf %lf %lf %lf %lf %lf\n", DIMENSION, N, MAX_CUBE_DIM, radius[0], WALL_TEMP[0][0], WALL_TEMP[1][0], WALL_TEMP[2][0], WALL_TEMP[3][0], WALL_TEMP[4][0], WALL_TEMP[5][0]);
	fclose(param_file);
	
	pos_file = fopen(strcat(strcpy(dir, dir_name),"positions.txt"), "w");	
	temp_file = fopen(strcat(strcpy(dir, dir_name),"temperatures.txt"), "w");	
	slice_file = fopen(strcat(strcpy(dir, dir_name),"slices.txt"), "w");
	ent_file = fopen(strcat(strcpy(dir, dir_name), "entropy.txt"),"w");
	mfp_file = fopen(strcat(strcpy(dir, dir_name), "mean_free_path.txt"),"w");
	// * * * * * * * * * * * * * * * * * * * * * * * //
	//*/ 

	while(time < STOP_TIME || burn_in_period > 0)
	{
		//Do the evolution
		// * * * * * * * * * * * * * * * * * * * * * * * //
		if (DIMENSION < 3)
		{
			for (i = 0; i < N; i++) p[i][2] = v[i][2] = 0.0;
		}
		//get minimum time step until next collision event
		first_collision = find_global_min_dt(&min_time);

		if (first_collision > 0)
		{
			printf("\nEND: NO MORE COLLISIONS\n");
			time = STOP_TIME;
		}

		// now check for complex collisions.
		// This is wen any collision has more than two colliding components. 
		// In this case, what happens is two of the components are chosen to collide, 
		// and all other particles are updated to a slightly smaller time step
		check_complex_collisions(&min_time);

		//*/
		if(burn_in_period > 0)
		{
		  	time++;
			fixed_dt += min_time;

			if( all_particles_diffused > (N - 1) )
			{
			  	//reset everything to 0 so we can 'begin' our simulation
				fixed_dt *= (1.0 * steps_per_record) / (1.0 * time);
				time_elapsed = 0.0;
				burn_in_period = time = 0;
				tt = total_entropy = total_heat_flow = 0.0;
				for(i = 0; i < N; i++) p_collisions[i] = w_collisions_heated[i] = w_collisions_passive[i] = 0;
			}
		}
		else
		{
			fprintf(ent_file, "%.5e, %.5e %.5e\n", tt, total_entropy, total_heat_flow);

			if( (time_elapsed < fixed_dt) && ( (time_elapsed + min_time) >= fixed_dt) )
			{
			  	time++;
				for(i = 0; i < N; i++)
				{
					p[i][0] += (fixed_dt - time_elapsed) * v[i][0];
					p[i][1] += (fixed_dt - time_elapsed) * v[i][1];
					p[i][2] += (fixed_dt - time_elapsed) * v[i][2];
				}

				fprintf(pos_file, "%d %lf ", time, tt);
				for(i = 0; i < N; i++)
				{
					fprintf(pos_file, "%lf %lf %lf ", p[i][0], p[i][1], p[i][2]);
				}
				fprintf(pos_file, "\n");

				compute_T();
				fprintf(temp_file, "%d %lf ", time, tt);
				for(i = 0; i < N; i++)
				{
					fprintf(temp_file, "%lf ", p_temp[i]);
				}
				fprintf(temp_file, "\n");

				add_to_slices();
				fprintf(slice_file, "%d ", time);
				for(i = 0; i < n_slices; i++){ fprintf(slice_file, "%d ", n_in_slice[i]);}
				for(i = 0; i < n_slices; i++){ fprintf(slice_file, "%d ", t_in_slice[i]);}
				fprintf(slice_file, "\n");

				for(i = 0; i < N; i++)
				{
					p[i][0] -= (fixed_dt - time_elapsed) * v[i][0];
					p[i][1] -= (fixed_dt - time_elapsed) * v[i][1];
					p[i][2] -= (fixed_dt - time_elapsed) * v[i][2];
				}
				time_elapsed = 0.0;
			}
			else
			{
				time_elapsed += min_time;
			}
		}
		//*/

		if (visualize && min_time > 0) smooth_vis(min_time);

		// update position of all particles to new time 
		for (i = 0; i < N; i++)
		{
			p[i][0] += global_min_time[i] * v[i][0];
			p[i][1] += global_min_time[i] * v[i][1];
			p[i][2] += global_min_time[i] * v[i][2];

			tot_dist[i] += sqrt(v[i][0] * v[i][0] +
				v[i][1] * v[i][1] +
				v[i][2] * v[i][2]) *
				global_min_time[i];
		}
		//find new velocity for particle(s) that did collide
		for (i = 0; i < N; i++)
		{
			if (how_many_p[i] > 0 || how_many_w[i] > 0)
			{
				resolve_particle_collision(i);
			}
		}

		if (DIMENSION < 3)
		{
			for (i = 0; i < N; i++) p[i][2] = v[i][2] = 0.0;
		}

		//fireproofing: check at each time step that no particles escaped. 
		for (i = 0; i < N; i++)
		{
			for (j = 0; j < 3; j++)
			{
				if (p[i][j] * p[i][j] > max_square_displacement[i])
				{
					printf("\n----------------------------------------\n");
					printf("Error at time %6d: particle %2d moved outside of the box! ", time, i);
					printf("\n----------------------------------------\n");
					time = STOP_TIME;
					exit(1);
				}
			}
		}
		tt += min_time;
		//end of this time step
	}

	printf("%i gas particles, %.1f steps, %.4f seconds in time\n", N, 1.0 * time, tt);

	for(i = 0; i < N; i++)
	{
		fprintf(mfp_file, "%lf, %lf, %lf, %lf\n", tot_dist[i], p_collisions[i], w_collisions_heated[i], w_collisions_passive[i]);
	}
	fclose(temp_file);
	fclose(pos_file);
	fclose(slice_file);
	fclose(ent_file);
	fclose(mfp_file);
}


void control()
{
	clock_t time_0, time_1;

	if (visualize) draw_picture();

	time_0 = clock();
    	n_body();
	time_1 = clock();

	printf("\n DONE \n");
	printf("\n Runtime %.5f seconds\n", (float)(time_1 - time_0) / CLOCKS_PER_SEC );
}

void Display(void)
{
  	if(DIMENSION == 3){
		gluLookAt(0.0, 0.0, 1.5 * EYE, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	}
	else
	{
		gluLookAt(0.0, 0.0, EYE, 0.0, 0.0, -EYE, 0.0, 1.0, 0.0);
	}
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	control();
}

void reshape(int w, int h)
{
	glViewport(0, 0, (GLsizei) w, (GLsizei) h);

	glMatrixMode(GL_PROJECTION);

	glLoadIdentity();

	glFrustum(-0.2, 0.2, -0.2, 0.2, 0.2, FAR);

	glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char** argv)
{
	glutInit(&argc,argv);

	//now that glut commands have been stripped from the 
	if(--argc < 1)
	{
		printf("Without input file, reverting to default parameters\n");
	}
	else
	{
		in_fname = argv[1];
	}
	set_initail_conditions();
	if(visualize)
	{
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGB);
		glutInitWindowSize(XWindowSize, YWindowSize);
		glutInitWindowPosition(0, 0);
		glutCreateWindow("Gas Particles");
		GLfloat light_position[] = { 1.0, 1.0, 1.0, 0.0 };
		GLfloat light_ambient[] = { 0.0, 0.0, 0.0, 1.0 };
		GLfloat light_diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
		GLfloat light_specular[] = { 1.0, 1.0, 1.0, 1.0 };
		GLfloat lmodel_ambient[] = { 0.2, 0.2, 0.2, 1.0 };
		GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
		GLfloat mat_shininess[] = { 10.0 };
		glClearColor(0.0, 0.0, 0.0, 0.0);
		glShadeModel(GL_SMOOTH);
		glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
		glLightfv(GL_LIGHT0, GL_POSITION, light_position);
		glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
		glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
		glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
		glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
		glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
		glEnable(GL_LIGHTING);
		glEnable(GL_LIGHT0);
		glEnable(GL_COLOR_MATERIAL);
		glEnable(GL_DEPTH_TEST);
		glutDisplayFunc(Display);
		glutReshapeFunc(reshape);
		glutMainLoop();
	}
	else
	{
		control();
	}
	return 0;
}
