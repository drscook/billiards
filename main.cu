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

//openGL vis parameters
#define XWindowSize 2500
#define YWindowSize 2500
#define FAR 50.0
#define DT 0.002

//wall types for collision 
#define passive 0
#define heated  1
#define sliding 2

#define MIN(x, y) (x < y) ? x : y
#define MAX(x, y) (x > y) ? x : y

std:: random_device generator;

dim3 block, grid;

//physical parameters
int DIMENSION = 3;
int N = 10;
float MAX_CUBE_DIM = 2.0;
float default_radius = 0.1;
float default_mass = 2.0;

//non-physical parameters
int STOP_TIME = 1000;
int steps_per_record = 50;
int visualize = 1;
int track_large_particle = 0;
int ignore_particle_interaction;
int all_particles_diffused = 0;
int diff_number = 2;
int packing_scheme = 0;

//physical constants
float BOLTZ_CONST = 1.0; //1.38064852e-23;
int max_walls = 3;  // number of walls that a particle can hit at a time

//file names for IO
char *in_fname = NULL;
char dir_name[256] = "\0./";
int n_events;
int collides[2];
int collider[2];

//particle
float3 *p_CPU, *p_GPU; // position
float3 *v_CPU, *v_GPU; // velocity
float *mass_CPU, *mass_GPU; // mass
float *radius_CPU, *radius_GPU; // radius
float *max_square_displacement = NULL; // dist from wall
float *p_temp_CPU; // kinetic energy
float *t_CPU, *t_GPU; // time particle will take to collide
float default_p_temp = 1.0; // default kinetic energy
float *p_collisions = NULL;//mean free path
float *w_collisions_heated = NULL;//mean free path
float *w_collisions_passive= NULL;//mean free path
float3 collision_normal;

//walls
float3 normal[6];
float alpha[6];
int WALL_TAG[6];
float WALL_TEMP[6][2];
float max_temp = 90.0;
float min_temp = 45.0;

//visualization
float EYE;

//memory management so particle collisions happen in order
int *tag_CPU, *tag_GPU;
int *what_w_CPU, *what_w_GPU;
int *what_p_CPU, *what_p_GPU;
int *how_many_p_CPU, *how_many_p_GPU;
int *how_many_w_CPU, *how_many_w_GPU;

//calculation variables
float gas_temp;

// reads input file of parameters (temperature of walls, particle size, &c)
void set_macros()
{
  	FILE * fp = NULL;
	const int bdim = 132;
	char buff[bdim];
	int i, d;
	double f, g;
	char s[256];

	//Default values
	//physical parmaeters
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
	WALL_TAG[2] = WALL_TAG[3] = WALL_TAG[4] = WALL_TAG[5] = passive;
	WALL_TEMP[0][0] = WALL_TEMP[0][1] = 45.0;
	WALL_TEMP[1][0] = WALL_TEMP[1][1] = 90.0;
	WALL_TEMP[2][0] = WALL_TEMP[2][1] = WALL_TEMP[3][0] = WALL_TEMP[3][1] = 0.0;
	WALL_TEMP[4][0] = WALL_TEMP[4][1] = WALL_TEMP[5][0] = WALL_TEMP[5][1] = 0.0;

	if( (fp = fopen(in_fname,"r")) == NULL)
	{
		printf("No input file. Using default values.\n");
	}
	else
	{
		fgets(buff,bdim,fp);
		fgets(buff,bdim,fp);
		sscanf(buff, "%d", &d);
		DIMENSION = d;

		//fgets(buff,bdim,fp);
		//fgets(buff,bdim,fp);
		//sscanf(buff, "%d", &d);
		//visualize = d;

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
		sscanf(buff, "%s", s);
		strcpy(dir_name, s);

		fgets(buff, bdim, fp);
		fgets(buff, bdim, fp);
		sscanf(buff, "%d", &d);
		steps_per_record = d;
		fclose(fp);
	}
}

// compute gas temperature based on kinetic energy of particles
void compute_t()
{
	int i;
	gas_temp = 0.0;

	for(i = 0; i < N; i++)
	{
		p_temp_CPU[i] = mass_CPU[i] * (v_CPU[i].x * v_CPU[i].x + v_CPU[i].y * v_CPU[i].y + v_CPU[i].z * v_CPU[i].z)
				/ (1.0 * DIMENSION * BOLTZ_CONST);
		gas_temp += p_temp_CPU[i];
	}
	gas_temp /= (1.0 * N);
}

// evaluates temperature of wall at a point assuming T varies linearly between ends
float wall_temperature(float x, int wall)
{
	return (x + MAX_CUBE_DIM) * (WALL_TEMP[wall][0] + (WALL_TEMP[wall][1] - WALL_TEMP[wall][0]) / (2.0 * MAX_CUBE_DIM));
}

float get_intersection_point(float time, int p)
{
	return p_CPU[p].x + v_CPU[p].x * time;
}

//initialize particles with position, velcoity, radius, etc. 
void pack_particles()
{
	int i, j, k, num, particles_per_side;
	float max_radius, hole_radius, temp, rad_sep = 1.0 + 1.0 / 20.0;
	float T, tol = 1.0e-7;
	float x_length, x_start;
	float y_length, y_start;
	float z_length, z_start;

	//tag particles not to hit themselves
	for(i = 0;  i < max_walls * N; i++) tag_CPU[i] = i / max_walls;
	
	//set initial particle parameters
	for (i = 0; i < N; i++)
	{
		mass_CPU[i] = default_mass;
		radius_CPU[i] = default_radius;
		p_temp_CPU[i] = default_p_temp;
	}

	if (track_large_particle)
	{
		radius_CPU[0] = 3.0 * radius_CPU[0];
		mass_CPU[0] = 3.0 * mass_CPU[0];
	}
	max_radius = radius_CPU[0];

	if(packing_scheme < 1)
	{
		//square lattice
		hole_radius = max_radius * rad_sep;
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
					p_CPU[num].x = x_start + 2.0*hole_radius*i;
					p_CPU[num].y = y_start + 2.0*hole_radius*j;
					p_CPU[num].z = z_start + 2.0*hole_radius*k;
					num++;
				}
			}
		}
	}
	else
	{
		//tetrahedral lattice
		hole_radius = max_radius * rad_sep;
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
					p_CPU[num].x = x_start + 2.0*hole_radius*i + (k % 2)*hole_radius;
					p_CPU[num].y = y_start + 2.0*hole_radius*j + (k % 2)*hole_radius;
					p_CPU[num].z = z_start + sqrt(2.0)*hole_radius*k;
					num++;
				}
			}
		}
	}

	if(DIMENSION < 3)
	{
		for(i = 0; i < N; i++)
		{
			p_CPU[i].z = v_CPU[i].z = 0.0;
		}
	}

	for (i = 0; i < N; i++)
	{
		max_square_displacement[i] = (MAX_CUBE_DIM - radius_CPU[i] + tol) * (MAX_CUBE_DIM - radius_CPU[i] + tol);
	}

	for (i = 0; i < N; i++)
	{
		T = sqrt(BOLTZ_CONST * p_temp_CPU[i] / mass_CPU[i]);
		std::normal_distribution<double> normal_dist(0, T);

		v_CPU[i].x = normal_dist(generator);
		v_CPU[i].y = normal_dist(generator);
		v_CPU[i].z = normal_dist(generator);
		p_collisions[i] = w_collisions_heated[i] = w_collisions_passive[i] = 0.0;
	}
}

void set_initial_conditions()
{
	int i;
	set_macros();

	EYE = 2.0 * MAX_CUBE_DIM; // set visualization vantage point

	//CPU MEMORY ALLOCATION
	p_CPU          = (float3*)malloc(             N * sizeof(float3) );
	v_CPU          = (float3*)malloc(             N * sizeof(float3) );
	radius_CPU     = (float* )malloc(             N * sizeof(float ) );
	max_square_displacement = (float*)malloc(N * sizeof(float));
	mass_CPU       = (float* )malloc(             N * sizeof(float ) );
	t_CPU          = (float* )malloc(             N * sizeof(float ) );
	p_temp_CPU     = (float* )malloc(             N * sizeof(float ) );
	tag_CPU        = (int*   )malloc( max_walls * N * sizeof(int   ) );
	how_many_p_CPU = (int*   )malloc(             N * sizeof(int   ) );
	how_many_w_CPU = (int*   )malloc(             N * sizeof(int   ) );
	what_p_CPU     = (int*   )malloc(             N * sizeof(int   ) );
	what_w_CPU     = (int*   )malloc( max_walls * N * sizeof(int   ) );
	p_collisions   = (float* )malloc(             N * sizeof(float ) );
	w_collisions_heated = (float* )malloc(        N * sizeof(float ) );
	w_collisions_passive= (float* )malloc(        N * sizeof(float ) );

	// GPU MEMORY ALLOCATION
	block.x = 512;
	block.y = 1;
	block.z = 1;

	grid.x = (N - 1) / block.x + 1;
	grid.y = 1;
	grid.z = 1;

	cudaMalloc( (void**)&p_GPU,       N *sizeof(float3) );
	cudaMalloc( (void**)&v_GPU,       N *sizeof(float3) );
	cudaMalloc( (void**)&radius_GPU,  N *sizeof(float ) );
	cudaMalloc( (void**)&mass_GPU,    N *sizeof(float ) );

	cudaMalloc( (void**)&tag_GPU,    max_walls * N *sizeof(int  ) );
	cudaMalloc( (void**)&t_GPU,                  N *sizeof(float) );
	cudaMalloc( (void**)&how_many_p_GPU,         N *sizeof(int  ) );
	cudaMalloc( (void**)&how_many_w_GPU,         N *sizeof(int  ) );
	cudaMalloc( (void**)&what_p_GPU,             N *sizeof(int  ) );
	cudaMalloc( (void**)&what_w_GPU, max_walls * N *sizeof(int  ) );



	// set up wall parameters for box
	max_temp = min_temp = 0.0; // average temperature of walls
	for(i = 0; i < 6; i++)
	{
		max_temp = MAX(max_temp, WALL_TEMP[i][0]);
		min_temp = MIN(min_temp, WALL_TEMP[i][0]);
	}
	default_p_temp = (max_temp + min_temp) / 2.0;

	//set entropy and collisions to 0
	for (i = 0; i < 6; i++){ alpha[i] = 1.0;}

	//normal vectors to walls 
	normal[0].x = 1.0; normal[0].y = 0.0; normal[0].z = 0.0;
	normal[2].x = 0.0; normal[2].y = 1.0; normal[2].z = 0.0;
	normal[4].x = 0.0; normal[4].y = 0.0; normal[4].z = 1.0;
	normal[1].x =-1.0; normal[1].y = 0.0; normal[1].z = 0.0;
	normal[3].x = 0.0; normal[3].y =-1.0; normal[3].z = 0.0;
	normal[5].x = 0.0; normal[5].y = 0.0; normal[5].z =-1.0;

	// set up particle parameters
	pack_particles();
	compute_t();

	// copy CPU initialization to GPU
	cudaMemcpy( p_GPU,      p_CPU,                  N *sizeof(float3), cudaMemcpyHostToDevice );
	cudaMemcpy( v_GPU,      v_CPU,                  N *sizeof(float3), cudaMemcpyHostToDevice );
	cudaMemcpy( mass_GPU,   mass_CPU,               N *sizeof(float ), cudaMemcpyHostToDevice );
	cudaMemcpy( radius_GPU, radius_CPU,             N *sizeof(float ), cudaMemcpyHostToDevice );
	cudaMemcpy( tag_GPU,    tag_CPU,    max_walls * N *sizeof(int   ), cudaMemcpyHostToDevice );

}

void draw_picture()
{
	int i;
	float red, green, blue, hi, mid;
	
	glClear(GL_COLOR_BUFFER_BIT);
	glClear(GL_DEPTH_BUFFER_BIT);
	
	glColor3d(1.0,1.0,0.5);
	for(i=0; i<N; i++)
	{
		hi = MIN(p_temp_CPU[i] / max_temp, 1.7);
		mid = 0.75 - hi;
		red = hi * hi;
		green = 1.0 - mid * mid * mid * mid;
		blue = 1.5 - red;

		glColor3f(red, green, blue);
		glPushMatrix();
		glTranslatef(p_CPU[i].x, p_CPU[i].y, p_CPU[i].z);
		glutSolidSphere(radius_CPU[i],20,20);
		glPopMatrix();
	}

	if(DIMENSION > 2)//if we're in 3D
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
	else // if we're in 2D
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

__device__ int particle_wall_collision(float3 * p, float3 * v, float * radius, float max_cube_dim, int i0, int i1, float * t)
{

	int collides = 0;
	float tt, max_cube_minus_radius;
	max_cube_minus_radius = max_cube_dim - radius[i0];
	tt = -2;

	if( (i1 == 0) && (v[i0].x * v[i0].x > 0.0) )
	{
		tt = ( (-max_cube_minus_radius) - p[i0].x ) / v[i0].x;
	}
	else if( (i1 == 2) && (v[i0].y * v[i0].y > 0.0) )
	{
		tt = ( (-max_cube_minus_radius) - p[i0].y ) / v[i0].y;
	}
	else if( (i1 == 4) && v[i0].z * v[i0].z > 0.0)  
	{
		tt = ( (-max_cube_minus_radius) - p[i0].z ) / v[i0].z;
	}
	else if( (i1 == 1) && (v[i0].x * v[i0].x > 0.0) )
	{
		tt = ( max_cube_minus_radius - p[i0].x ) / v[i0].x;
	}
	else if( (i1 == 3) && (v[i0].y * v[i0].y > 0.0) )
	{
		tt = ( max_cube_minus_radius - p[i0].y ) / v[i0].y;
	}
	else if( (i1 == 5) && (v[i0].z * v[i0].z > 0.0) )
	{
		tt = ( max_cube_minus_radius - p[i0].z ) / v[i0].z;
	}
	if( tt >= 0.0)
	{
	  	(*t) = tt;
		collides = 1;
	}

	return collides;
}

__device__ int particle_particle_collision(float3 * p, float3 * v, float * radius, int i0, int i1, float max_cube_dim, float * t)
{
	float xc, yc, zc, xv, yv, zv;
	float discriminant, dd, a, b, c, dt;
	int collides = 0;

	dd = (radius[i0] + radius[i1]) * (radius[i0] + radius[i1]);

	xc = p[i0].x - p[i1].x;
	yc = p[i0].y - p[i1].y;
	zc = p[i0].z - p[i1].z;

	if( (xc * xc + yc * yc + zc * zc) > dd){

		xv = v[i0].x - v[i1].x;
		yv = v[i0].y - v[i1].y;
		zv = v[i0].z - v[i1].z;

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

__global__ void find(float3 * p, float3 * v, float * radius, float * mass,  // particle data--position, velocity, radius 
			int * tag, int * how_many_p, int * how_many_w, int * what_p_hit, int * what_w_hit, // memory management--what each particle hits and has hit
			int n, int max_walls, float max_cube_dim, int ignore_particle_interaction, float * min_dt) // macros--number of particles, shape of geometry, &c.
{
	float t, global_min_dt = 0.0;
	int j, k, ok, collides, first_collision, this_particle;

	first_collision = 1;
	global_min_dt = -1.0;

	this_particle = blockDim.x * blockIdx.x + threadIdx.x;

	if(this_particle < n)
	{
		how_many_p[this_particle] = 0;
		how_many_w[this_particle] = 0;

		first_collision = 1;

		if(ignore_particle_interaction)
		{
			// do nothing-only interested in particle-wall collisions
		}
		else
		{
			// check current particle against all particles for collision
			for(j = 0; j < n; j++)
			{
				if( (this_particle != j) && ( (tag[max_walls * this_particle] != j) || (tag[max_walls * j] != this_particle) ) )
				{
				  	collides = particle_particle_collision(p, v, radius, this_particle, j, max_cube_dim, &t);
					if( (collides > 0) && ( (t < global_min_dt) || (first_collision > 0) ) )
					{
					  	first_collision = 0;
						what_p_hit[this_particle] = j;
						how_many_p[this_particle] = 1;
						global_min_dt = t;
					}
					else if( (collides > 0) && (t <= global_min_dt) )
					{
						global_min_dt = t;
						how_many_p[this_particle]++;
					}
				}
			}
		}

		// check current particle against walls for collision
		for(j = 0; j < 6; j++)
		{
			ok = 1;
			for(k = 0; k < max_walls; k++) 
			{
				if(tag[max_walls * this_particle + k] == -(j + 1) )
				{
					ok = 0;
				}
			}

			if(ok)
			{
				collides = particle_wall_collision(p, v, radius, max_cube_dim, this_particle, j, &t);

				if( (collides > 0) && ( (t < global_min_dt) || (first_collision > 0) ) )
				{
					first_collision = 0;
					global_min_dt = t;

					how_many_p[this_particle] = 0; // since we're re-setting the min time, it doesn't hit any other particles
					what_p_hit[this_particle] = n; // reset the number of particles this has hit to be 0 for new min timestep

					how_many_w[this_particle] = 1;
					what_w_hit[max_walls * this_particle] = -(j + 1);
				}
				else if( (collides > 0) && (t <= global_min_dt) )
				{
					what_w_hit[max_walls * this_particle + how_many_w[this_particle]] = -(j + 1);
					how_many_w[this_particle]++;
					global_min_dt = t;
				}
			}
		}
		min_dt[this_particle] = (first_collision < 1) ? global_min_dt : -10;
	}
}


float3 specular_reflect(float3 v_in, float3 n)
{
	float dotproduct;
	float3 v_out;

  	dotproduct = - (v_in.x * n.x + v_in.y * n.y + v_in.z * n.z);
	v_out.x = v_in.x + 2.0 * dotproduct * n.x;
	v_out.y = v_in.y + 2.0 * dotproduct * n.y;
	v_out.z = v_in.z + 2.0 * dotproduct * n.z;

	return (v_out);
}


float chi_sq(float u, float sigma)
{
	return sigma * sqrt(fabs(2.0 * log(1.0 - u)));
}

float3 heated_wall(float3 v_in, float T, float m, float3 n)
{
	float3 t1, t2, v_out;
	float u, v, w, len, sigma = sqrt(BOLTZ_CONST * T / m);
	int ok = 0;

	std:: normal_distribution<double> norm_dist(0, sigma);
	std:: uniform_real_distribution<double> uniform_dist(0, 1.0);
	u = uniform_dist(generator);
	u = chi_sq(u, sigma);
	v = norm_dist(generator);
	w = norm_dist(generator);

	if( (n.x * n.x > 0.0) || (n.y * n.y > 0.0) )
	{
		t1.x = -n.y;
		t1.y =  n.x;
		t1.z =  0.0;

		ok = 1;
	}
	else if( n.z * n.z > 0.0 )
	{
		t1.x = -n.z;
		t1.y =  0.0;
		t1.z =  n.x;
		ok = 1;
	}
	if(ok)
	{
		len = t1.x * t1.x + t1.y * t1.y + t1.z * t1.z;
		len = (len > 0) ? sqrt(len) : 1.0;
		t1.x /= len; t1.y /= len; t1.z /= len;

		t2.x = n.y * t1.z - n.z * t1.y;
		t2.y = n.z * t1.x - n.x * t1.z;
		t2.z = n.x * t1.y - n.y * t1.x;

		len = t2.x * t2.x + t2.y * t2.y + t2.z * t2.z;
		len = (len > 0) ? sqrt(len) : 1.0;
		t2.x /= len; t2.y /= len; t2.z /= len;

		v_out.x = (u * n.x + v * t1.x + w * t2.x);
		v_out.y = (u * n.y + v * t1.y + w * t2.y);
		if(DIMENSION > 2) v_out.z = (u * n.z + v * t1.z + w * t2.z);
	}
	return v_out;
}
void resolve_particle_collision(int i0, float time)
{
	float3 t, x_v1, y_v1, x_v2, y_v2, n;
	float len, m1m2, u;
	int wall;
	int i, i1;
	float temp_at_point, point;
	std:: uniform_real_distribution<double> u_dist(0, 1.0);

	n_events = 0;

	if( how_many_w_CPU[i0] > 0) // if particle hits a wall
	{
		t = v_CPU[i0];
		n_events = 1;
		collides[0] = i0;
		collider[0] = what_w_CPU[i0];

		if(how_many_w_CPU[i0] > 1) // if it hits mutliple walls
		{
			// if it hits multiple walls, then it hits a 
			// corner (either of 2 or 3 walls). 
			// Solution is to find average of all 
			// normal vectors at that corner, and 
			// have particle rebound specularly 
			// about that averaged normal
			w_collisions_passive[i0] += 1.0;
			n.x = n.y = n.z = 0.0;
			for(i = 0; i < how_many_w_CPU[i0]; i++)
			{
				wall = -(1 + what_w_CPU[max_walls * i0 + i]);
				tag_CPU[max_walls * i0 + i] = what_w_CPU[max_walls * i0 + i];

				n.x += normal[wall].x;
				n.y += normal[wall].y;
				n.z += normal[wall].z;
			}
			len = sqrt(n.x * n.x + n.y * n.y + n.z * n.z);
			len = (len > 0) ? len : 1.0;
			n.x /= len; n.y /= len; n.z /= len;

			v_CPU[i0] = specular_reflect(t, n);
			collision_normal = n;
		}
		else // if hits single wall
		{
			tag_CPU[max_walls * i0] = what_w_CPU[max_walls * i0];
			wall = -(1 + what_w_CPU[max_walls * i0]);

			if( WALL_TAG[wall] == passive )
			{
				v_CPU[i0] = specular_reflect(t, normal[wall]);
				w_collisions_passive[i0] += 1.0;
			}
			else if( WALL_TAG[wall] == heated )
			{
				w_collisions_heated[i0] += 1.0;
				if( (w_collisions_heated[i0] < (diff_number + 0.5)) && (w_collisions_heated[i0] > (diff_number - 0.5)) ) all_particles_diffused++;

				u = u_dist(generator);
				// alpha is percentage of time particles are affected by wall temperature 
				if(u > alpha[wall])
				{
					v_CPU[i0] = specular_reflect(t, normal[wall]);
				}
				else
				{
					point = get_intersection_point(time, i0);
//					temp_at_point = wall_temperature(point, wall);
					temp_at_point = WALL_TEMP[wall][0];

					v_CPU[i0] = heated_wall(t, temp_at_point, mass_CPU[i0], normal[wall]);

				}
			}
		collision_normal = normal[wall];
		}

		for(i = how_many_w_CPU[i0]; i < max_walls; i++) tag_CPU[max_walls * i0 + i] = N;
	}
	else if(how_many_p_CPU[i0] > 0) // if it hits a particle
	{
		i1 = what_p_CPU[i0];
		p_collisions[i0] += 1.0;
		p_collisions[i1] += 1.0;
		n_events = 2;
		collides[0] = i0;
		collides[1] = i1;
		collider[0] = i1;
		collider[1] = i0;

		if(i1 < N)
		{
			for(i = 0; i < max_walls; i++) tag_CPU[max_walls * i1 + i] = tag_CPU[max_walls * i0 + i] = N;

			tag_CPU[max_walls * i0] = i1;
			tag_CPU[max_walls * i1] = i0;

			how_many_w_CPU[i0] = how_many_w_CPU[i1] = how_many_p_CPU[i1] = how_many_p_CPU[i0] = 0;

			m1m2 = mass_CPU[i0] + mass_CPU[i1];
			m1m2 = (m1m2 * m1m2 > 0) ? m1m2 : 1.0;

			// t is the vector from center of i0 to center of i1
			t.x = p_CPU[i0].x - p_CPU[i1].x;
			t.y = p_CPU[i0].y - p_CPU[i1].y;
			t.z = p_CPU[i0].z - p_CPU[i1].z;

			// normalize t
			len = sqrt(t.x * t.x + t.y * t.y + t.z * t.z);
			len = (len > 0) ? len : 1.0;
			t.x /= len;
			t.y /= len;
			t.z /= len;

			//--------------------Particle 1-----------------------//
			// get vector component of v_CPU[i0] parallel to t (x_v1)
			len = t.x * v_CPU[i0].x + t.y * v_CPU[i0].y + t.z * v_CPU[i0].z;
			x_v1.x = len * t.x;
			x_v1.y = len * t.y;
			x_v1.z = len * t.z;
			// get vector component of v_CPU[i0] perpendicular to t (y_v1)
			y_v1.x = v_CPU[i0].x - x_v1.x;
			y_v1.y = v_CPU[i0].y - x_v1.y;
			y_v1.z = v_CPU[i0].z - x_v1.z;
			//-----------------------------------------------------//


			//--------------------Particle 2-----------------------//
			// get vector component of v_CPU[i1] parallel to t (x_v2)
			t.x = -t.x; t.y = -t.y; t.z = -t.z;
			len = t.x * v_CPU[i1].x + t.y * v_CPU[i1].y + t.z * v_CPU[i1].z;
			x_v2.x = len * t.x;
			x_v2.y = len * t.y;
			x_v2.z = len * t.z;
			// get vector component of v_CPU[i1] perpendicular to t (y_v2)
			y_v2.x = v_CPU[i1].x - x_v2.x;
			y_v2.y = v_CPU[i1].y - x_v2.y;
			y_v2.z = v_CPU[i1].z - x_v2.z;
			//-----------------------------------------------------//

			// update vidx0 and vidx1 to hold new velocities
			v_CPU[i0].x = (x_v1.x * (mass_CPU[i0] - mass_CPU[i1]) + x_v2.x * 2.0 * mass_CPU[i1]) / m1m2 + y_v1.x;
			v_CPU[i0].y = (x_v1.y * (mass_CPU[i0] - mass_CPU[i1]) + x_v2.y * 2.0 * mass_CPU[i1]) / m1m2 + y_v1.y;
			v_CPU[i0].z = (x_v1.z * (mass_CPU[i0] - mass_CPU[i1]) + x_v2.z * 2.0 * mass_CPU[i1]) / m1m2 + y_v1.z;

			v_CPU[i1].x = (x_v2.x * (mass_CPU[i1] - mass_CPU[i0]) + x_v1.x * 2.0 * mass_CPU[i0]) / m1m2 + y_v2.x;
			v_CPU[i1].y = (x_v2.y * (mass_CPU[i1] - mass_CPU[i0]) + x_v1.y * 2.0 * mass_CPU[i0]) / m1m2 + y_v2.y;
			v_CPU[i1].z = (x_v2.z * (mass_CPU[i1] - mass_CPU[i0]) + x_v1.z * 2.0 * mass_CPU[i0]) / m1m2 + y_v2.z;
		}
	}
}

void smooth_vis(float time)
{
	//float3 new_pos[N];
	float3 * new_pos = (float3*)malloc(N * sizeof(float3));
	float dt;
	int i, j, num;

	dt = DT;
	num = (int)((time / dt) + 0.5);

	compute_t();

	for(i = 0; i < N; i++)
	{
		new_pos[i].x = p_CPU[i].x;
		new_pos[i].y = p_CPU[i].y;
		new_pos[i].z = p_CPU[i].z;
	}

	dt = time / (1.0 * (num - 1.0));
	if((num - 1) < 1){dt = time;}

	for(i = 0; i < num; i++)
	{
		for(j = 0; j < N; j++)
		{
			p_CPU[j].x = new_pos[j].x + (v_CPU[j].x * dt * i);
			p_CPU[j].y = new_pos[j].y + (v_CPU[j].y * dt * i);
			p_CPU[j].z = new_pos[j].z + (v_CPU[j].z * dt * i);
		}
		draw_picture();
	}

	for(i = 0; i < N; i++)
	{
		p_CPU[i].x = new_pos[i].x;
		p_CPU[i].y = new_pos[i].y;
		p_CPU[i].z = new_pos[i].z;
	}
}

//float find_min_dt_and_update(float * times, int * idx, int * i0, int * i1)
void find_min_dt(float * times, float * min_dt)
{
	int i;
	float t = -1., timetol = 1.0e-4;

	// find minimum value in times, and the corresponding value in idx
	//              both at array index (*i0), and save:
	//
	//          (*i1) = idx[(*i0)],    min(t) = times[(*i0)]
	for(i = 0; i < N; i++) 
	{
		if( (how_many_p_CPU[i] > 0) || (how_many_w_CPU[i] > 0) )
		{
			t = times[i];
		}
	}
	for(i = 0; i < N; i++)
	{
		if( (t >= times[i]) && ( (how_many_p_CPU[i] > 0) || (how_many_w_CPU[i] > 0) ) )
		{
			t = times[i];
		}
	}

	for(i = 0; i < N; i++)
	{
		if( ( (how_many_p_CPU[i] < 1) && (how_many_w_CPU[i] < 1) ) || (times[i] > (t + timetol) ) )
		{
			how_many_p_CPU[i] = 0;
			how_many_w_CPU[i] = 0;
		}
		times[i] = 0.99 * t;
	}

	(*min_dt) = 0.99 * t;
}

void errorCheck(int num, const char * message)
{
	cudaError_t error;
	error = cudaGetLastError();
	if(error != cudaSuccess)
	{
		printf("Cuda Error at time %d: %s = %s\n", num, message, cudaGetErrorString(error));
	}
}

void check_complex_collisions(float * t, float * particle_t)
{
	int i, type_1, type_2;

	type_1 = N; 
	type_2 = N;

	for(i = 0; i < N; i++)
	{
		if( (how_many_p_CPU[i] > 0) && (how_many_w_CPU[i] > 0) )
		{
			type_1 = i;
		}
		else if( how_many_p_CPU[i] > 1)
		{
			type_2 = i;
		}
	}
	type_2 = (type_1 < N) ? N : type_2;

	if(type_1 < N) how_many_p_CPU[type_1] = 0;

	if( (type_1 < N) || (type_2 < N) )
	{
		for(i = 0; i < N; i++)
		{
			if( (i != type_1) && (i != type_2) )
			{
				how_many_p_CPU[i] = how_many_w_CPU[i] = 0;
			}
		}
	}
	if( (type_1 < N) || (type_2 < N) ) (*t) = 0.95 * (*t);

	for(i = 0; i < N; i++)
	{
		if( ( i != type_1 ) && (i != type_2 ) )
		{
			particle_t[i] = (*t);
		}
	}
}

void n_body()
{
	float t, tt = 0.0;
	FILE * vis_file;
	char dir[256];
	int burn_in_period = 0, i, j, time, n;

	/*/		OUTPUT FILE STUFF		 /*/
	vis_file = fopen(strcat(strcpy(dir, dir_name), "visualization.csv"), "w");
	fprintf(vis_file, "#box dimension\n box, %lf\n", MAX_CUBE_DIM);
	fprintf(vis_file, "#particle radii\n");
	for(i = 0; i < N; i++)
	{
		fprintf(vis_file, "r, %d, %lf, %lf\n", i, radius_CPU[i], mass_CPU[i]);
	}
	for(i = 0; i < N; i++)
	{
		fprintf(vis_file, "c, %d, %d, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", 
					i, 0, tt, 
					p_CPU[i].x, p_CPU[i].y, p_CPU[i].z, 
					0.0, 0.0, 0.0,
					v_CPU[i].x, v_CPU[i].y, v_CPU[i].z,
					collision_normal.x, collision_normal.y, collision_normal.z
			);
	}

	time = 0;
	n = N;
	while(time++ < STOP_TIME || burn_in_period > 0)
	{
	  	// on GPU - find smallest time step s.t. any particle(s) collide either 
	  	// with each other or a wall and update all particles to that time step

		find<<<grid, block>>>(p_GPU, v_GPU, radius_GPU, mass_GPU, 
					tag_GPU, how_many_p_GPU, how_many_w_GPU, what_p_GPU, what_w_GPU, 
					n, max_walls, MAX_CUBE_DIM, ignore_particle_interaction, t_GPU);
		errorCheck(time, "find");

		//copy minimum time step and index of corresponding colliding element onto CPU 
		cudaMemcpy( how_many_p_CPU, how_many_p_GPU,             N * sizeof(int  ), cudaMemcpyDeviceToHost);
		cudaMemcpy( how_many_w_CPU, how_many_w_GPU,             N * sizeof(int  ), cudaMemcpyDeviceToHost);
		cudaMemcpy(     what_p_CPU,     what_p_GPU,             N * sizeof(int  ), cudaMemcpyDeviceToHost);
		cudaMemcpy(     what_w_CPU,     what_w_GPU, max_walls * N * sizeof(int  ), cudaMemcpyDeviceToHost);
		cudaMemcpy(          t_CPU,          t_GPU,             N * sizeof(float), cudaMemcpyDeviceToHost);

		find_min_dt(t_CPU, &t);

		// if no collisions were detected, we are done. 
		if(t < 0.0)
		{
			printf("\nEND: NO MORE COLLISIONS\n");
			break;
		}

		check_complex_collisions(&t, t_CPU);

		if(visualize) smooth_vis(t);
		for (i = 0; i < N; i++)
		{
			p_CPU[i].x += t_CPU[i] * v_CPU[i].x;
			p_CPU[i].y += t_CPU[i] * v_CPU[i].y;
			p_CPU[i].z += t_CPU[i] * v_CPU[i].z;
		}

		for(i = 0; i < N; i++)
		{
			if( how_many_p_CPU[i] > 0 || how_many_w_CPU[i] > 0)
			{
				n_events = 0;

				resolve_particle_collision(i, t_CPU[i]);

				for(j = 0; j < n_events; j++)
				{
					fprintf(vis_file, "c, %d, %d, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", 
						collides[j], collider[j], tt, 
						p_CPU[collides[j]].x, p_CPU[collides[j]].y, p_CPU[collides[j]].z, 
						0.0, 0.0, 0.0,
						v_CPU[collides[j]].x, v_CPU[collides[j]].y, v_CPU[collides[j]].z,
						collision_normal.x, collision_normal.y, collision_normal.z
					);
				}
			}
			if(DIMENSION < 3) p_CPU[i].z = v_CPU[i].z = 0.0;
		}

		// update position on GPU to new time step
		// update velocity on GPU to match CPU 
		// (and also tag which keeps track of most recent collision for each particle)
		cudaMemcpy(   p_GPU,   p_CPU,             N * sizeof(float3), cudaMemcpyHostToDevice );
		cudaMemcpy(   v_GPU,   v_CPU,             N * sizeof(float3), cudaMemcpyHostToDevice );
		cudaMemcpy( tag_GPU, tag_CPU, max_walls * N * sizeof(int   ), cudaMemcpyHostToDevice );

		tt += t;

		//fireproofing: check at each time step that no particles escaped.
		for (i = 0; i < N; i++)
		{
			if  ((p_CPU[i].x * p_CPU[i].x > max_square_displacement[i]) ||
				 (p_CPU[i].y * p_CPU[i].y > max_square_displacement[i]) ||
				 (p_CPU[i].z * p_CPU[i].z > max_square_displacement[i]))
			{
				printf("\nError in time %6d: particle %2d escaped!!! last hit %2d %2d %2d\n",
					time, i, tag_CPU[max_walls * i], tag_CPU[max_walls * i + 1], tag_CPU[max_walls * i + 2]);
				time = STOP_TIME;
				break;
			}
		}
		//end of this time step
	}
	printf("%i gas particles, %.1f timesteps, %.4f seconds in time\n", N, 1.0 * time, tt);

	for(i = 0; i < N; i++)
	{
		fprintf(vis_file, "c, %d, %d, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", 
					i, 0, tt, 
					p_CPU[i].x, p_CPU[i].y, p_CPU[i].z, 
					0.0, 0.0, 0.0,
					v_CPU[i].x, v_CPU[i].y, v_CPU[i].z,
					collision_normal.x, collision_normal.y, collision_normal.z
			);
	}
	fclose(vis_file);
}



void control()
{	
	clock_t time_0, time_1;
	
	if(visualize) draw_picture();

	time_0 = clock();
    	n_body();
	time_1 = clock();
	
	printf("\n Runtime %.5f seconds\n", (float)(time_1 - time_0) / CLOCKS_PER_SEC);
	printf("\n DONE \n");
}

void Display(void)
{
  	if(DIMENSION == 3){
		gluLookAt(0.0, 0.0, 1.5 * EYE, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
	}
	else
	{
		gluLookAt(0.0, 0.0, EYE, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
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
	//glutInit(&argc,argv);

	//now that glut commands have been stripped from the 
	if(--argc < 1)
	{
		printf("Without input file, reverting to default parameters\n");
	}
	else
	{
		in_fname = argv[1];
	}
	set_initial_conditions();

  visualize = false;
	if (visualize)
	{
		glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGB);
		glutInitWindowSize(XWindowSize, YWindowSize);
		glutInitWindowPosition(500, 0);
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
