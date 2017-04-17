# Kidnapped Vehicle Project
Udacity Self-Driving Car Engineer Nanodegree Program

---

## Project Introduction
This is the 8th project of the Udacity self-driving car engineering nanodegree program. A robot has been kidnapped and transported to a new location! Luckily it has a map of this location, a (noisy) GPS estimate of its initial location, and lots of (noisy) sensor and control data.

In this project a 2-D particle filter is implemented in C++. 

## Running the Code
Run the following in the repository folder:  

```
> ./clean.sh
> ./build.sh
> ./run.sh
```

Note: ./build-rev.sh is created for windowd environment.

Time step: 2444
Cumulative mean weighted error: x .1 y .1 yaw .02
Runtime (sec): 38.187226
Success! Your particle filter passed!


## Implementing the Particle Filter
The directory structure of this repository is as follows:

```
root
|   build.sh
|   clean.sh
|   CMakeLists.txt
|   README.md
|   run.sh
|
|___data
|   |   control_data.txt
|   |   gt_data.txt
|   |   map_data.txt
|   |
|   |___observation
|       |   observations_000001.txt
|       |   ... 
|       |   observations_002444.txt
|   
|___src
    |   helper_functions.h
    |   main.cpp
    |   map.h
    |   particle_filter.cpp
    |   particle_filter.h
```

The only file that is modified for the proejct is `particle_filter.cpp` in the `src` directory. The file contains the entire code for a `ParticleFilter` class and some associated methods. 

## Inputs to the Particle Filter
You can find the inputs to the particle filter in the `data` directory. 

## The Map*
`map_data.txt` includes the position of landmarks (in meters) on an arbitrary Cartesian coordinate system. Each row has three columns
1. x position
2. y position
3. landmark id

> * Map data provided by 3D Mapping Solutions GmbH.


## Control Data
`control_data.txt` contains rows of control data. Each row corresponds to the control data for the corresponding time step. The two columns represent
1. vehicle speed (in meters per second)
2. vehicle yaw rate (in radians per second)

## Observation Data
The `observation` directory includes around 2000 files. Each file is numbered according to the timestep in which that observation takes place. 

These files contain observation data for all "observable" landmarks. Here observable means the landmark is sufficiently close to the vehicle. Each row in these files corresponds to a single landmark. The two columns represent:
1. x distance to the landmark in meters (right is positive) RELATIVE TO THE VEHICLE. 
2. y distance to the landmark in meters (forward is positive) RELATIVE TO THE VEHICLE.


