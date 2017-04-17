/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#define _USE_MATH_DEFINES

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <list>

#include "particle_filter.h"

using namespace std;


/*
* Calculates the bivariate normal pdf of a point given a mean and std and assuming zero correlation
*/
inline double bivariate_normal(double x, double y, double mu_x, double mu_y, double sig_x, double sig_y) {
  return exp(-((x - mu_x)*(x - mu_x) / (2 * sig_x*sig_x) + (y - mu_y)*(y - mu_y) / (2 * sig_y*sig_y))) / (2 * M_PI*sig_x*sig_y);
}

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  // Set  the number of particles to draw
  num_particles = 25; 

  // noise generation
  default_random_engine gen;
  normal_distribution<double> N_x_init(0, std[0]);
  normal_distribution<double> N_y_init(0, std[1]);
  normal_distribution<double> N_theta_init(0, std[2]);

  // add particles at iniital position + noise with weight 1.0 for all
  for (int i = 0; i < num_particles; ++i) {
    Particle p_tmp;
    p_tmp.id = i;
    p_tmp.x = x + N_x_init(gen);
    p_tmp.y = y + N_y_init(gen);
    p_tmp.theta = theta + N_theta_init(gen);
    p_tmp.weight = 1.0;
    weights.push_back(p_tmp.weight);
    particles.push_back(p_tmp);
  }

  // Change the flag to inidicate it is initialized now
  is_initialized = true;
}



void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

  // noise generation
  default_random_engine gen;
  normal_distribution<double> N_x_init(0, std_pos[0]);
  normal_distribution<double> N_y_init(0, std_pos[1]);
  normal_distribution<double> N_theta_init(0, std_pos[2]);

  for (unsigned int i=0; i < particles.size(); ++i) {
    
    // define temporary variables
    double x_new;
    double y_new;
    double theta_new;

    // update next state
    if (yaw_rate < 0.001) {
      x_new = particles[i].x + velocity*cos(particles[i].theta)*delta_t;
      y_new = particles[i].y + velocity*sin(particles[i].theta)*delta_t;
      theta_new = particles[i].theta + yaw_rate*delta_t;
    } else {
      x_new = particles[i].x + velocity / yaw_rate*(sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
      y_new = particles[i].y + velocity / yaw_rate*(-cos(particles[i].theta + yaw_rate*delta_t) + cos(particles[i].theta));
      theta_new = particles[i].theta + yaw_rate*delta_t;
    }

    // add noise and update particle
    particles[i].x = x_new + N_x_init(gen);
    particles[i].y = y_new + N_y_init(gen);
    particles[i].theta = theta_new + N_theta_init(gen);

  }

}


void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

  // iterate through observations
  for (unsigned i = 0; i < observations.size(); ++i) {

    // define and initialize temporary variables
    double cur_dist = 1e6;
    int closest_j = -1;

    // iterate through predicted landmarks to find the closest
    for (unsigned j = 0; j < predicted.size(); ++j) {

      // calculate Euclidian distance between 2 points
      double eval_dist = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);

      // assign min
      if (eval_dist < cur_dist) {
        cur_dist = eval_dist;
        closest_j = j;
      }
    }
    // assign the closest id to the obeservation
    observations[i].id = predicted[closest_j].id;
  }
}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html


  // emptying the list of weights for all particles
  weights.clear();

  // loop through all particles
  for (unsigned i = 0; i < particles.size(); ++i) {

    // 1- Transform the coordinates of the observations from local vehicle coordinates to map coordinates
    std::vector<LandmarkObs> obs_map_coord;
    for (unsigned int j = 0; j < observations.size(); ++j) {
      if (dist(observations[j].x, observations[j].y, 0, 0) <= sensor_range) {
        LandmarkObs obs_tmp;
        obs_tmp.x = particles[i].x + observations[j].x*cos(particles[i].theta) - observations[j].y*sin(particles[i].theta);
        obs_tmp.y = particles[i].y + observations[j].x*sin(particles[i].theta) + observations[j].y*cos(particles[i].theta);
        obs_tmp.id = -1;
        obs_map_coord.push_back(obs_tmp);
      }
    }
    

    // 2- create a list of close-by landmarks in map coordintes
    std::vector<LandmarkObs> close_landmarks;
    for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); ++j) {
      if (dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f) <= sensor_range) {
        LandmarkObs obs_tmp;
        obs_tmp.x = map_landmarks.landmark_list[j].x_f;
        obs_tmp.y = map_landmarks.landmark_list[j].y_f;
        obs_tmp.id = map_landmarks.landmark_list[j].id_i;
        close_landmarks.push_back(obs_tmp);
      }
    }

    
    // 3- find the closest landmark id for each observaton by calling the dataAssociation function (first predicted_obs should be calculated)
    dataAssociation(close_landmarks, obs_map_coord);

    // 4- calculate and assign the wieghts by multiplying the probability pdf values - i.e. bivarite probability of each relevant
    //    observation given the mean of map and std_landmark
    double weight = 1;
    for (unsigned int j = 0; j < close_landmarks.size(); j++) {
      double min_dist = 1e6;
      int min_k = -1;
      // find the map coordinates for the landmark
      for (unsigned int k = 0; k < obs_map_coord.size(); ++k) {
        // find the min_distance observation only for those observations that are assigned to the landmark - i.e. closest to the landmark
        if (obs_map_coord[k].id == close_landmarks[j].id) {
          double eval_dist = dist(close_landmarks[j].x, close_landmarks[j].y, obs_map_coord[k].x, obs_map_coord[k].y);
          if (eval_dist < min_dist) {
            min_dist = eval_dist;
            min_k = k;
          }
        }
      }
      if (min_k != -1) {
        weight *= bivariate_normal(obs_map_coord[min_k].x, obs_map_coord[min_k].y, close_landmarks[j].x, close_landmarks[j].y, std_landmark[0], std_landmark[1]);
      }
    }

    // update the weight of the particle
    weights.push_back(weight);
    particles[i].weight = weight;

  } 

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  // random sampler setup
  default_random_engine gen;
  auto first = weights.cbegin();
  auto last = weights.cend();
  auto count = distance(first, last);
  discrete_distribution<int> dist(
    count,
    -0.5,
    -0.5 + count,
    [&first](size_t i)
  {
    return *std::next(first, i);
  });

  // re-sample
  std::vector<Particle> resampled_particles;
  for (int i = 0; i<num_particles; ++i) {
    resampled_particles.push_back(particles[dist(gen)]);
  }
  particles = resampled_particles;

}

void ParticleFilter::write(std::string filename) {
	// You don't need to modify this file.
	std::ofstream dataFile;
	dataFile.open(filename, std::ios::app);
	for (int i = 0; i < num_particles; ++i) {
		dataFile << particles[i].x << " " << particles[i].y << " " << particles[i].theta << "\n";
	}
	dataFile.close();
}
