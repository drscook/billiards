# What is this repository for? 

This project is a multi-particle code for non-continuum gas dynamics simulations. 

* A user-input number of particles are initialized in a cube. 
They have elastic collisions with other particles, and wall interactions that 
follow knudsen's cosine law for thermal effects. 

* Can produce visualizations (using the GLUT library for OpenGL) and data for analysis. 
(data includes mean free path, velocities, particles positions, etc.)

* There is a serial version of the code for single processor machines (serial branch), 
and a parallel version that uses Nvidia's CUDA platform (master branch). 

# Hardware/Software requirements for running 

* C++ compiler that supports c++11 libraries (specifically random number generation)

* Cuda-enabled device (for parallel code, not required for serial)

# To run simulations:
* Required files: 

    * main.cu (main code)

    * input_file (options for customizing each simulation)

* Command line arguments: 
```
#!bash
    > make
    > ./a.out file
```
