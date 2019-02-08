/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	if(is_initialized){
		return;	
	}
	num_particles = 100;
	default_random_engine gen;

	normal_distribution<double> init_noise_x(x,std[0]);
	normal_distribution<double> init_noise_y(y,std[1]);
	normal_distribution<double> init_noise_theta(theta,std[2]);

	for(int i=0;i<num_particles;i++){
		Particle p;
		p.id = i;
		p.weight = 1.0;
		p.x = init_noise_x(gen);
		p.y = init_noise_y(gen);
		p.theta = init_noise_theta(gen);		
		particles.push_back(p);
	}

	is_initialized = true;
	cout <<"Particles Initialized"<<endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
	normal_distribution<double> pred_noise_x(0,std_pos[0]);
	normal_distribution<double> pred_noise_y(0,std_pos[1]);
	normal_distribution<double> pred_noise_theta(0,std_pos[2]);

	for(int i=0; i<num_particles; i++){
		if( fabs(yaw_rate) < 0.0001){ 
			particles[i].x += velocity * delta_t * cos(particles[i].theta);
			particles[i].y += velocity * delta_t * sin(particles[i].theta);
		} 
		else{
			particles[i].x += velocity / yaw_rate * ( sin( particles[i].theta + yaw_rate*delta_t ) - sin(particles[i].theta) );
			particles[i].y += velocity / yaw_rate * ( cos( particles[i].theta ) - cos( particles[i].theta + yaw_rate*delta_t ) );
			particles[i].theta += yaw_rate * delta_t;
		}
		particles[i].x += pred_noise_x(gen);
    	particles[i].y += pred_noise_y(gen);
		particles[i].theta += pred_noise_theta(gen);
	}
}


void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
  //   implement this method and use it as a helper during the updateWeights phase.

  
}


void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
    const std::vector<LandmarkObs> &observations,const Map &map_landmarks) {
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

  for (int i = 0; i < num_particles; i++) {
    Particle p = particles[i];
   
	// 1. Store landmarks within sensor range
	vector<LandmarkObs> preds;
	for(unsigned int j=0;j<map_landmarks.landmark_list.size();j++){
		auto lm = map_landmarks.landmark_list[j];
		double e_dist = dist(p.x,p.y,lm.x_f,lm.y_f);
		if(e_dist<sensor_range){
			preds.push_back(LandmarkObs{lm.id_i,lm.x_f,lm.y_f});
		}
	}
	
	// 2. Translation and rotation
	double x_comp = cos(p.theta);
	double y_comp = sin(p.theta);
	vector<LandmarkObs> obs_trans;
	for(unsigned int j=0;j<observations.size();j++){
		auto obs = observations[j];
		LandmarkObs temp;
		temp.x = obs.x*x_comp - obs.y*y_comp + p.x;
		temp.y = obs.x*y_comp + obs.y*x_comp + p.y;
		obs_trans.push_back(temp);
	}

	// 3. Data Association 
    for (unsigned int i = 0; i < obs_trans.size(); i++) {
		auto obs = obs_trans[i];
		double min_dist = 1000000.0;    
		for (unsigned int j = 0; j < preds.size(); j++) {
			auto pred = preds[j];
			double e_dist = dist(obs.x, obs.y, pred.x, pred.y);
			if (e_dist < min_dist) {
				min_dist = e_dist;
				obs_trans[i].id = pred.id;
			}
		}
  	}

	// 4. Update weights
	particles[i].weight = 1.0;
	for(unsigned int j=0;j<obs_trans.size();j++){
		LandmarkObs obs_ = obs_trans[j];
		LandmarkObs pre_;
		for (unsigned int k = 0; k < preds.size(); k++) {
			if (preds[k].id == obs_.id) {
				pre_ = preds[k];
			}
		}
		double x_comp = pow(obs_.x - pre_.x, 2) / (2 * pow(std_landmark[0], 2));
		double y_comp = pow(obs_.y - pre_.y, 2) / (2 * pow(std_landmark[1], 2));
		double w_comp = exp(-(x_comp + y_comp)) / (2 * M_PI * std_landmark[0] * std_landmark[1]);
		particles[i].weight*=w_comp;
	}
  }
}

void ParticleFilter::resample() {
  // TODO: Resample particles with replacement with probability proportional to their weight. 
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  
  vector<Particle> resampled;
  vector<double> weights;
  default_random_engine gen;

  for (int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
  }

  double max_w = *max_element(weights.begin(), weights.end());
  uniform_real_distribution<double> real_dist(0.0, max_w);
  uniform_int_distribution<int> int_dist(0, num_particles-1);
  int idx = int_dist(gen);

  double beta = 0.0;
  for (int i = 0; i < num_particles; i++) {
    beta += real_dist(gen) * 2.0;
    while (beta > weights[idx]) {
      beta -= weights[idx];
      idx = (idx + 1) % num_particles;
    }
	resampled.push_back(particles[idx]);
  }
  particles = resampled;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

	particle.sense_x.clear();
	particle.sense_y.clear();
	particle.associations.clear();

    particle.associations= associations;
    particle.sense_x = sense_x;
	particle.sense_y = sense_y;
	
	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
