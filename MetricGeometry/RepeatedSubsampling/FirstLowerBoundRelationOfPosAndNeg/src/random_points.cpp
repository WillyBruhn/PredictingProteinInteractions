/*
 * random_points.h
 *
 *  Created on: 26.03.2018
 *      Author: Berens
 *
 *	This function gets two numbers, the first is the number of points in a file,
 *		the second is the number of points for that the bound should be calculated
 */

//#ifndef RANDOM_POINTS_H_
//#define RANDOM_POINTS_H_

#include <random>
#include <algorithm>
#include "random_points.h"
#include <iostream>


//#include <gsl/gsl_sf_bessel.h>
//#include <gsl/gsl_randist.h>

using namespace std;

vector<unsigned int> random_points_old(unsigned int number_of_all_points, unsigned int number_to_selected)
{
	/*
		random_points:
		Samples randomly number_to_selected integers from range 0 to number_of_all_points.

		Input
		-----
		number_of_all_points - an int, which stands for the number of points on the isosurface
		number_to_selected - an int, which stands for the number of points that should be selected

		Output
		------
		Is a vector of size number_to_selected, filled with the randomly drawn integers

		Complexity
    		----------
		O(number_of_all_points)

		Example
		-------
		vector <int> A = random_points(10, 3);
		cout << A[0] << " " << A[1] << " " << A[2];

		=> 2 0 4
	*/
	if(number_of_all_points < number_to_selected)
	{
		cerr << "Can't select more points, then there are. Number of all points: " << number_of_all_points << " Number of selected points: " << number_to_selected;
		exit(-1);
	}

	vector<unsigned int> list_of_numbers;
	vector<unsigned int> out;
	unsigned int random_number;

	for(unsigned int i = 0; i < number_of_all_points; i++)
	{
		list_of_numbers.push_back(i);
	}

	for(unsigned int i = 0; i < number_to_selected; i++)
	{
		random_number = rand() % list_of_numbers.size();
		out.push_back(list_of_numbers[random_number]);
		list_of_numbers.erase(list_of_numbers.begin() + random_number);
	}

	return(out);
}

/*
 * Claimed to be a very fast random generator
 * https://stackoverflow.com/questions/1640258/need-a-fast-random-generator-for-c
 */
static unsigned long x=123456789, y=362436069, z=521288629;

unsigned long xorshf96(void) {          //period 2^96-1
unsigned long t;
    x ^= x << 16;
    x ^= x >> 5;
    x ^= x << 1;

   t = x;
   x = y;
   y = z;
   z = t ^ x ^ y;

  return z;
}

std::random_device rd;
std::mt19937 g(rd());



void random_points(const unsigned int &number_of_all_points, const unsigned int &number_to_select,
//							vector<unsigned int> &list_of_numbers,
							vector<unsigned int> &indices_selected)
{
	/*
		random_points:
		Samples randomly number_to_selected integers from range 0 to number_of_all_points.

		Input
		-----
		number_of_all_points - an int, which stands for the number of points on the isosurface
		number_to_selected - an int, which stands for the number of points that should be selected

		Output
		------
		Is a vector of size number_to_selected, filled with the randomly drawn integers

		Complexity
    		----------
		O(number_of_all_points)

		Example
		-------
		vector <int> A = random_points(10, 3);
		cout << A[0] << " " << A[1] << " " << A[2];

		=> 2 0 4
	*/

//	if(number_of_all_points < number_to_selected)
//	{
//		cerr << "Can't select more points, then there are. Number of all points: " << number_of_all_points << " Number of selected points: " << number_to_selected;
//		exit(-1);
//	}


	vector<unsigned int> list_of_numbers(number_of_all_points,0);

	for(unsigned int i = 0; i < number_of_all_points; i++)
	{
		list_of_numbers[i] = i;
	}

	//--------------------------------------------------
	double random_number = 0;

	for(unsigned int i = 0; i < number_to_select; i++)
	{
//		indices_selected[i] = rand() % number_of_all_points;

		random_number = rand() % list_of_numbers.size();

		indices_selected[i] = list_of_numbers[random_number];
		list_of_numbers.erase(list_of_numbers.begin() + random_number);

//		if(indices_selected[i] >= number_of_all_points) cerr << "ERROR" << endl;
	}


//	std::shuffle(indices.begin(), indices.end(), g);
//	std::copy(indices.begin(), indices.begin() + number_to_select, indices_selected.begin());

}

std::default_random_engine generator;

std::discrete_distribution<unsigned int> distribution_X_pos;
void set_distribution_X_pos(vector<double> &mu){
//	distribution_X_pos(mu.begin(), mu.end());

	std::discrete_distribution<unsigned int> d(mu.begin(), mu.end());
	distribution_X_pos = d;
}

std::discrete_distribution<unsigned int> distribution_X_neg;
void set_distribution_X_neg(vector<double> &mu){
	std::discrete_distribution<unsigned int> d(mu.begin(), mu.end());
//	distribution_X_neg.param(d.param());

	distribution_X_neg = d;
}

std::discrete_distribution<unsigned int> distribution_Y_pos;
void set_distribution_Y_pos(vector<double> &mu){
	std::discrete_distribution<unsigned int> d(mu.begin(), mu.end());
//	distribution_Y_pos.param(d.param());

	distribution_Y_pos = d;
}

std::discrete_distribution<unsigned int> distribution_Y_neg;
void set_distribution_Y_neg(vector<double> &mu){
	std::discrete_distribution<unsigned int> d(mu.begin(), mu.end());
//	distribution_Y_neg.param(d.param());

	distribution_Y_neg = d;
}



void random_points_with_measure(std::discrete_distribution<unsigned int> &distribution,
								const unsigned int &number_to_select,
								vector<unsigned int> &indices_selected)
{
	/*
	 * Theorem 5.1 e suggests that points should be sampled with replacement.
	 * Points with higher probability can be sampled multiple times.
	 */

//	  std::discrete_distribution<unsigned int> distribution(mu.begin(), mu.end());

	  for(unsigned int i=0; i<number_to_select; ++i) {
	    unsigned int number = distribution(generator);

	    indices_selected[i] = number;
	  }
}


//#endif /* RANDOM_POINTS_H_ */
