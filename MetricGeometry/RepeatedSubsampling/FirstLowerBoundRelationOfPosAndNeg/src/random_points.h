/*
 * random_points.h
 *
 *  Created on: 26.03.2018
 *      Author: Berens
 *	
 *	This function gets two numbers, the first is the number of points in a file,
 *		the second is the number of points for that the bound should be calculated
 */

#ifndef RANDOM_POINTS_H_
#define RANDOM_POINTS_H_

#include <random>
#include <algorithm>

//#include <gsl/gsl_sf_bessel.h>
// #include <gsl/gsl_randist.h>

using namespace std;

vector<unsigned int> random_points_old(unsigned int number_of_all_points, unsigned int number_to_selected);

/*
 * Claimed to be a very fast random generator
 * https://stackoverflow.com/questions/1640258/need-a-fast-random-generator-for-c
 */
//static unsigned long x=123456789, y=362436069, z=521288629;

unsigned long xorshf96(void);

extern std::random_device rd;
extern std::mt19937 g;



void random_points(const unsigned int &number_of_all_points, const unsigned int &number_to_select,
//							vector<unsigned int> &list_of_numbers,
							vector<unsigned int> &indices_selected);

extern std::default_random_engine generator;

extern std::discrete_distribution<unsigned int> distribution_X_pos;
void set_distribution_X_pos(vector<double> &mu);

extern std::discrete_distribution<unsigned int> distribution_X_neg;
void set_distribution_X_neg(vector<double> &mu);

extern std::discrete_distribution<unsigned int> distribution_Y_pos;
void set_distribution_Y_pos(vector<double> &mu);

extern std::discrete_distribution<unsigned int> distribution_Y_neg;
void set_distribution_Y_neg(vector<double> &mu);



void random_points_with_measure(std::discrete_distribution<unsigned int> &distribution,
								const unsigned int &number_to_select,
								vector<unsigned int> &indices_selected);


#endif /* RANDOM_POINTS_H_ */
