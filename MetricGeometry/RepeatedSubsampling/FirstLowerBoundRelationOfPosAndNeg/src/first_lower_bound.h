/*
 * first_lower_bound.h
 *
 *  Created on: 17.10.2018
 *      Author: Berens
 *       Function to compute a lower bound of the Gromov Wasserstein Metric for discrete points,
 *          see section 3.4.1 and 4.2.1 "First Lower Bound" in "Quantitative comparison of protein isosurfaces 
 *			with approximated Gromov-Wasserstein-distance" from Felix Berens.
 */

#ifndef FIRST_LOWER_BOUND_H_
#define FIRST_LOWER_BOUND_H_

#include "mathFunctions.h"
#include <vector>
#include <algorithm>
#include <assert.h>
#include <sstream>
#include <fstream>

#include <boost/algorithm/string.hpp>
//#include <boost/filesystem.hpp>

using namespace std;

void list_of_unique_eccent(vector<double> &ecc1, vector<double> &ecc2, vector<double> &listEcc);
double eccentricities(vector< vector <double> > &d, vector<double> &distribution, int i);
void all_eccentricities(vector< vector <double> > &d, vector<double> &distribution,vector<double> &ecc);
void which_index(vector<double> &ecc, vector<double> &list, double u);
double distr_of_eccentricities(vector<double> &ecc, vector<double> &distribution, double u);
double flb(vector< vector <double> > &points1, vector< vector <double> > &points2, vector<double> &distribution1, 
	   vector<double> &distribution2, int number_of_selected_points);



struct thread_data {
	int  thread_id;
	vector< vector <double> > points1_pos_all;
	vector< vector <double> > points2_pos_all;

	vector< vector <double> > points1_neg_all;
	vector< vector <double> > points2_neg_all;

    vector<double> distribution1_pos_all;
    vector<double> distribution2_pos_all;

    vector<double> distribution1_neg_all;
    vector<double> distribution2_neg_all;



    // everything for flb below
	vector< vector <double> > points1_pos;
	vector< vector <double> > points2_pos;

	vector< vector <double> > points1_neg;
	vector< vector <double> > points2_neg;

    vector<double> distribution1_pos;
    vector<double> distribution2_pos;

    vector<double> distribution1_neg;
    vector<double> distribution2_neg;

	// positive and negative
	vector < vector<double> > d1_pos;
	vector < vector<double> > d2_pos;

	vector < vector<double> > d1_neg;
	vector < vector<double> > d2_neg;

	// merging the positive and negative of one protein
	vector< vector <double> > points1_pos_and_neg;
	vector< vector <double> > points2_pos_and_neg;


	vector <double> distribution1_pos_and_neg;
	vector <double> distribution2_pos_and_neg;

	// and for the relation between positive and negative
	vector < vector<double> > d1_pos_and_neg;
	vector < vector<double> > d2_pos_and_neg;

	vector<double> ecc1_pos;
	vector<double> ecc2_pos;

	vector<double> ecc1_neg;
	vector<double> ecc2_neg;

	// and for the relation between the potentials within one protein
	vector<double> ecc1_pos_and_neg;
	vector<double> ecc2_pos_and_neg;

	vector<double> listEcc_pos;
	vector<double> listEcc_neg;

	// and for the relation between the potentials within one protein
	vector<double> listEcc_pos_and_neg;

	double c1,c2,c3;

	unsigned int number_of_selected_points;

	vector<unsigned int> random_indices_list;

	void resize(const unsigned int number_of_selected_points){
		points1_pos.resize(3,vector<double>(number_of_selected_points));
		points2_pos.resize(3,vector<double>(number_of_selected_points));

		points1_neg.resize(3,vector<double>(number_of_selected_points));
		points2_neg.resize(3,vector<double>(number_of_selected_points));

	    distribution1_pos.resize(number_of_selected_points);
	    distribution2_pos.resize(number_of_selected_points);

	    distribution1_neg.resize(number_of_selected_points);
	    distribution2_neg.resize(number_of_selected_points);

		// positive and negative
		d1_pos.resize(number_of_selected_points, vector<double>(number_of_selected_points,0));
		d2_pos.resize(number_of_selected_points, vector<double>(number_of_selected_points,0));

		d1_neg.resize(number_of_selected_points, vector<double>(number_of_selected_points,0));
		d2_neg.resize(number_of_selected_points, vector<double>(number_of_selected_points,0));

		// merging the positive and negative of one protein
		points1_pos_and_neg.resize(3,vector<double>(number_of_selected_points*2));
		points2_pos_and_neg.resize(3,vector<double>(number_of_selected_points*2));

		distribution1_pos_and_neg.resize(number_of_selected_points*2,1/(number_of_selected_points*2));
		distribution2_pos_and_neg.resize(number_of_selected_points*2,1/(number_of_selected_points*2));

		// and for the relation between positive and negative
		d1_pos_and_neg.resize(2*number_of_selected_points, vector<double>(2*number_of_selected_points,0));
		d2_pos_and_neg.resize(2*number_of_selected_points, vector<double>(2*number_of_selected_points,0));

		ecc1_pos.resize(number_of_selected_points);
		ecc2_pos.resize(number_of_selected_points);

		ecc1_neg.resize(number_of_selected_points);
		ecc2_neg.resize(number_of_selected_points);

		// and for the relation between the potentials within one protein
		ecc1_pos_and_neg.resize(2*number_of_selected_points);
		ecc2_pos_and_neg.resize(2*number_of_selected_points);

		listEcc_pos.resize(number_of_selected_points+number_of_selected_points);
		listEcc_neg.resize(number_of_selected_points+number_of_selected_points);

		// and for the relation between the potentials within one protein
		listEcc_pos_and_neg.resize(number_of_selected_points*4);

		this->number_of_selected_points = number_of_selected_points;

		random_indices_list.resize(number_of_selected_points);
	}
};



void list_of_unique_eccent(vector<double> &ecc1, vector<double> &ecc2, vector<double> &listEcc)
{
	/*
		list_of_unique_eccent: 
		Gets two eccentricities vectors and gives a vector (listEcc) of all occurring eccentricities out.
		This list is sorted in ascending order.

		Input 
		-----
		ecc1 - vector containing the eccentricities for points of the first object
		ecc2 - vector containing the eccentricities for points of the second object
		listEcc - here all eccentricities occurring eccentricities will be stored

		Output
		------
		No output

		Complexity
    		----------
		O(N·log(N)), where N = std::distance(listEcc.begin(), listEcc.end())

		Example
		-------
		vector<double> listEcc(12);
		vector <double > ecc1{ 5, 10, 13, 2, 2, 2};
		vector <double > ecc2{ 3, 34, 1, 2 ,5, 1};
		
		list_of_unique_eccent(ecc1, ecc2, listEcc);
		
		=> listEcc = {1, 2, 3, 5, 10, 13, 34}
	*/
	listEcc.resize(ecc1.size() + ecc2.size());
	
	for(unsigned int i = 0; i < ecc1.size(); i++)
	{
		listEcc[i] = ecc1[i];
	}
	
	for(unsigned int i = 0; i < ecc2.size(); i++)
	{
		listEcc[i+ecc1.size()] = ecc2[i];
	}

	sort(listEcc.begin(), listEcc.end());
  	vector<double>::iterator last = std::unique(listEcc.begin(), listEcc.end());
	listEcc.erase(last, listEcc.end());
}

double eccentricities(vector< vector <double> > &d, vector<double> &distribution, int i)
{
	/*
		eccentricities: 
		See "3.4.1 First Lower Bound" in "Quantitative comparison of protein isosurfaces 
 			with approximated Gromov-Wasserstein-distance" Definition 3.21

		Input 
		-----
		d - is from type vector, which represent a matrix, with dimension
			number_of_selected_points x number_of_selected_points
		distribution - vector containing the weight for the points
		i - integer in range 0 to number_of_selected_points, which is the
			index of the point for that the eccentricity is to calculate

		Output
		------
		Is a double that is the eccentricity of the i-th point

		Complexity
    		----------
		O(n), where n = number_of_selected_points

		Example
		-------
		vector < vector<double> > d{{0, 2, 2},
						{2, 0, 2},	
						{2, 2, 0}};
		vector <double> distribution{0.25,0.5,0.25};
		int i = 1;
		
		eccentricities(d, distribution, i);

		=> 1
	*/
	double s = 0;

	for(unsigned int k = 0; k < distribution.size(); k++)
	{
		s += distribution[k] * d[i][k];
	}
	return(s);
}

void all_eccentricities(vector< vector <double> > &d, vector<double> &distribution,vector<double> &ecc)
{	
	/*
		all_eccentricities: 
		calculates all eccentricities for a given set of points

		Input 
		-----
		d - is from type vector, which represent a matrix, with dimension
			number_of_selected_points x number_of_selected_points
		distribution - vector containing the weight for the points
		ecc - is a vector, here all eccentrisities will be stored

		Output
		------
		No return

		Complexity
    		----------
		O(n^2), where n = number_of_selected_points

		Example
		-------
		vector < vector<double> > d{{0, 2, 2},
						{2, 0, 2},
						{2, 2, 0}};
		vector <double> distribution{0.25,0.5,0.25};
		vector <double> ecc(3);
		
		all_eccentricities(d, distribution, ecc);

		=> 1.5 1 1.5
	*/
	for(unsigned int k = 0; k < distribution.size(); k++)
	{
		ecc[k] = eccentricities(d, distribution, k);
	}
}

void which_index(vector<double> &ecc, vector<double> &list, double u)
{
	/*
		which_index: 
		determines all indices of points that have an eccentricity lower than u

		Input 
		-----
		ecc - vector that contains all eccentricities of the points in an objects
		list - here the indices that have an eccentricity lower than u will be stored
		u - a double, that specifies how low the eccentricities must be

		Output
		------
		No return

		Complexity
    		----------
		O(n^2), where n = number_of_selected_points

		Example
		-------
		vector <double> list;
		vector <double> ecc{1.5, 1, 1.5};
		double u = 1;

		which_index(ecc, list, u);
		
		=> 1
	*/
	int count = 0;
	list.resize(ecc.size());

	for(unsigned int i = 0; i < ecc.size(); i++)
	{
		if(ecc[i] <= u)
		{
			list[count] = i;
			count++;
		}
	}
	list.resize(count);
}

double distr_of_eccentricities(vector<double> &ecc, vector<double> &distribution, double u)
{
	/*
		distr_of_eccentricities: 
		Calculates the probability if one choose a point randomly, with the given distribution, 
		that the eccentricity is lower than u.

		Input 
		-----
		ecc - vector that contains all eccentricities of the points in an objects
		distribution - vector containing the weight for the points
		u - a double, that specifies how low the eccentricities must be

		Output
		------
		Is a double that is the distribution of eccentricities.

		Complexity
    		----------
		O(n^2), where n = number_of_selected_points

		Example
		-------
		vector <double> distribution{0.25,0.5,0.25};
		vector <double> ecc{1.5, 1, 1.5};
		double u = 1;
		
		distr_of_eccentricities(ecc, distribution, u);

		=> 0.5
	*/
	
	vector<double>list_index_of_points_ecc_le_u;
	double mu = 0;

	which_index(ecc, list_index_of_points_ecc_le_u, u);

	for(unsigned int i = 0; i < list_index_of_points_ecc_le_u.size(); i++)
	{
		mu = mu + distribution[list_index_of_points_ecc_le_u[i]];
	}
	return(mu);
}

double distr_of_eccentricities_sorted(const vector<double> &ecc_sorted, const vector<double> &distribution, const double &u, unsigned int &last_index, double &last_integral){
	/*
	 * Assuming the eccentricities list ist sorted.
	 * And the distribution is permutated in the same way.
	 */

	for(unsigned int i = last_index; i < ecc_sorted.size(); i++)
	{

		if(ecc_sorted[i] > u)
		{
//			i--;
			last_index = i;

			return last_integral;
		}

		last_integral += distribution[i];
	}

	last_index = ecc_sorted.size();

	return last_integral;
}

void sort_ecc_and_dist(vector<double> &ecc, vector<double> &distribution){
	/*
	 * sorts the eccentricities.
	 * permutates distribution in the same way.
	 */

	vector<unsigned int> index(ecc.size(), 0);
	for (unsigned int i = 0 ; i < index.size() ; i++) {
	    index[i] = i;
	}

	sort(index.begin(), index.end(),
	    [&](const int& a, const int& b) {
	        return (ecc[a] < ecc[b]);
	    }
	);

	vector<double> distribution_copy = distribution;
	vector<double> ecc_copy = ecc;
	for (int i = 0 ; i < index.size() ; i++) {
	    distribution[i] = distribution_copy[index[i]];
	    ecc[i] = ecc_copy[index[i]];
	}
}



double flb(vector< vector <double> > &points1, vector< vector <double> > &points2, vector<double> &distribution1, 
	   vector<double> &distribution2, int number_of_selected_points)
{
	/*
		flb: 
		Calculates the first lower bound. See section 3.4.1 and 4.2.1 "First Lower Bound" in 
			"Quantitative comparison of protein isosurfaces 
 			with approximated Gromov-Wasserstein-distance" from Felix Berens.

		Input 
		-----
		points1 - vector with dimension number_of_selected_points x 3. Points of the first object as 3D coordinates
		points2 - vector with dimension number_of_selected_points x 3. Points of the first object as 3D coordinates
		distribution1 - vector containing the weight for the points of the first object
		distribution2 - vector containing the weight for the points of the second object
		number_of_selected_points - a integer, for how many points the first lower bound will be calculated
		
		Output
		------
		Is a double that is the first lower bound.

		Complexity
    		----------
		O(n^2), where n = number_of_selected_points

		Example
		-------
		vector < vector<double> > points1{{0,0,1},
						{1,0,0},
						{0,1,0}};
		vector < vector<double> > points2{{2,0,1},
						{1,0,2},
						{0,1,3}};
		vector <double> distribution1{0.4,0.2,0.2};
		vector <double> distribution2{0.25,0.5,0.25};
		int number_of_selected_points = 3;
		flb(points1, points2, distribution1, distribution2, number_of_selected_points);

		=> 0.372717
	*/
	vector < vector<double> > d1(number_of_selected_points, vector<double>(number_of_selected_points,0));
	vector < vector<double> > d2(number_of_selected_points, vector<double>(number_of_selected_points,0));
	vector <double> eccentricities1(number_of_selected_points);
	vector <double> eccentricities2(number_of_selected_points);
	vector<double> ecc1(number_of_selected_points);
	vector<double> ecc2(number_of_selected_points);
	vector<double> listEcc(number_of_selected_points+number_of_selected_points);
	double u;
	double out = 0;
	
	// calculation of the distance matrices
	dist_matrix(d1, points1);
	dist_matrix(d2, points2);

	// calculation of all eccentricities
	all_eccentricities(d1, distribution1,ecc1);
	all_eccentricities(d2, distribution2,ecc2);

	// mergin and sorting the list of eccentricities
	list_of_unique_eccent(ecc1, ecc2, listEcc);

	// calculating the lower bound
	for(unsigned int i = 0; i < listEcc.size() - 1; i++)
	{
		u = listEcc[i];
		out += (listEcc[i+1]-listEcc[i])*abs(distr_of_eccentricities(ecc1, distribution1, u) 
						     - distr_of_eccentricities(ecc2, distribution2, u));
	}
	return(0.5 * out);
}


void appendPoints(const vector< vector <double> > &points_pos, const vector< vector <double> > &points_neg, vector< vector <double> > &points){
/*
*	Input: points_pos, points_neg
*	Output: points <- is a merge of points_pos and points_neg
*/
	/*cout << "merging" << endl;
	cout << points_pos[0].size() << endl;
	cout << points_pos.size() << endl;
	cout << points_neg[0].size() << endl;
	cout << points_neg.size() << endl;
	cout << points[0].size() << endl;
	cout << points.size() << endl;
	*/

	for(unsigned int i= 0; i < points_pos[0].size(); i++){
		points[0][i] = points_pos[0][i];
		points[1][i] = points_pos[1][i];
		points[2][i] = points_pos[2][i];
	}

	for(unsigned int i= 0; i < points_neg[0].size(); i++){
		points[0][i + points_pos[0].size()] = points_neg[0][i];
		points[1][i + points_pos[0].size()] = points_neg[1][i];
		points[2][i + points_pos[0].size()] = points_neg[2][i];

	}
}

void appendDistributions(const vector <double> &distribution_pos, const vector <double> &distribution_neg, vector <double> &distribution){
/*
*	Input: points_pos, points_neg
*	Output: points <- is a merge of points_pos and points_neg
*/
	
	for(unsigned int i= 0; i < distribution_pos.size(); i++){
		distribution[i] = distribution_pos[i];
	}

	for(unsigned int i= 0; i < distribution_neg.size(); i++){
		distribution[i + distribution_pos.size()] = distribution_neg[i];
	}
}

double flb_with_relationOLD(	vector< vector <double> > &points1_pos, vector<double> &distribution1_pos,
			 	vector< vector <double> > &points1_neg, vector<double> &distribution1_neg,

				vector< vector <double> > &points2_pos, vector<double> &distribution2_pos,
				vector< vector <double> > &points2_neg, vector<double> &distribution2_neg, 
   				int number_of_selected_points)
{
	/*
		flb: 
		Calculates the first lower bound. See section 3.4.1 and 4.2.1 "First Lower Bound" in 
			"Quantitative comparison of protein isosurfaces 
 			with approximated Gromov-Wasserstein-distance" from Felix Berens.

		Input 
		-----
		points1 - vector with dimension number_of_selected_points x 3. Points of the first object as 3D coordinates
		points2 - vector with dimension number_of_selected_points x 3. Points of the first object as 3D coordinates
		distribution1 - vector containing the weight for the points of the first object
		distribution2 - vector containing the weight for the points of the second object
		number_of_selected_points - a integer, for how many points the first lower bound will be calculated
		
		Output
		------
		Is a double that is the first lower bound.

		Complexity
    		----------
		O(n^2), where n = number_of_selected_points

		Example
		-------
		vector < vector<double> > points1{{0,0,1},
						{1,0,0},
						{0,1,0}};
		vector < vector<double> > points2{{2,0,1},
						{1,0,2},
						{0,1,3}};
		vector <double> distribution1{0.4,0.2,0.2};
		vector <double> distribution2{0.25,0.5,0.25};
		int number_of_selected_points = 3;
		flb(points1, points2, distribution1, distribution2, number_of_selected_points);

		=> 0.372717
	*/

	// positive and negative
	vector < vector<double> > d1_pos(number_of_selected_points, vector<double>(number_of_selected_points,0));
	vector < vector<double> > d2_pos(number_of_selected_points, vector<double>(number_of_selected_points,0));

	vector < vector<double> > d1_neg(number_of_selected_points, vector<double>(number_of_selected_points,0));
	vector < vector<double> > d2_neg(number_of_selected_points, vector<double>(number_of_selected_points,0));

	// merging the positive and negative of one protein
	vector< vector <double> > points1_pos_and_neg(3,vector<double>(number_of_selected_points*2));

	appendPoints(points1_pos, points1_neg, points1_pos_and_neg);

	vector< vector <double> > points2_pos_and_neg(3,vector<double>(number_of_selected_points*2));
	appendPoints(points2_pos, points2_neg, points2_pos_and_neg);

	vector <double> distribution1_pos_and_neg(number_of_selected_points*2,1/(number_of_selected_points*2));
	appendDistributions(distribution1_pos, distribution1_neg, distribution1_pos_and_neg);

	vector <double> distribution2_pos_and_neg(number_of_selected_points*2,1/(number_of_selected_points*2));
	appendDistributions(distribution2_pos, distribution2_neg, distribution2_pos_and_neg);


	// and for the relation between positive and negative
	vector < vector<double> > d1_pos_and_neg(2*number_of_selected_points, vector<double>(2*number_of_selected_points,0));
	vector < vector<double> > d2_pos_and_neg(2*number_of_selected_points, vector<double>(2*number_of_selected_points,0));


//	vector <double> eccentricities1(number_of_selected_points);
//	vector <double> eccentricities2(number_of_selected_points);

	vector<double> ecc1_pos(number_of_selected_points);
	vector<double> ecc2_pos(number_of_selected_points);

	vector<double> ecc1_neg(number_of_selected_points);
	vector<double> ecc2_neg(number_of_selected_points);

	// and for the relation between the potentials within one protein
	vector<double> ecc1_pos_and_neg(2*number_of_selected_points);
	vector<double> ecc2_pos_and_neg(2*number_of_selected_points);

	vector<double> listEcc_pos(number_of_selected_points+number_of_selected_points);
	vector<double> listEcc_neg(number_of_selected_points+number_of_selected_points);

	// and for the relation between the potentials within one protein
	vector<double> listEcc_pos_and_neg(number_of_selected_points*4);



	double u_pos;
	double u_neg;
	double u_pos_and_neg;

	double out_pos = 0;
	double out_neg = 0;
	double out_pos_and_neg = 0;

	double c1 = 1.0;
	double c2 = 1.0;
	double c3 = 1.0;
	
	// calculation of the distance matrices
	dist_matrix(d1_pos, points1_pos);
	dist_matrix(d2_pos, points2_pos);

	dist_matrix(d1_neg, points1_neg);
	dist_matrix(d2_neg, points2_neg);

	dist_matrix(d1_pos_and_neg, points1_pos_and_neg);
	dist_matrix(d2_pos_and_neg, points2_pos_and_neg);




	// calculation of all eccentricities
	all_eccentricities(d1_pos, distribution1_pos,ecc1_pos);
	all_eccentricities(d2_pos, distribution2_pos,ecc2_pos);

	all_eccentricities(d1_neg, distribution1_neg,ecc1_neg);
	all_eccentricities(d2_neg, distribution2_neg,ecc2_neg);

	all_eccentricities(d1_pos_and_neg, distribution1_pos_and_neg, ecc1_pos_and_neg);
	all_eccentricities(d2_pos_and_neg, distribution2_pos_and_neg, ecc2_pos_and_neg);


	// mergin and sorting the list of eccentricities
	list_of_unique_eccent(ecc1_pos, ecc2_pos, listEcc_pos);
	list_of_unique_eccent(ecc1_neg, ecc2_neg, listEcc_neg);

	list_of_unique_eccent(ecc1_pos_and_neg, ecc2_pos_and_neg, listEcc_pos_and_neg);


	/*
	printPoints(points1_pos,"points1_pos");
	printMatrix(d1_pos,"d1_pos");

	printPoints(points1_neg,"points1_neg");
	printMatrix(d1_neg,"d1_neg");

	printPoints(points1_pos_and_neg,"points1_pos_and_neg");
	printMatrix(d1_pos_and_neg,"d1_pos_and_neg");
	*/


	// calculating the lower bound POS
	for(unsigned int i = 0; i < listEcc_pos.size() - 1; i++)
	{
		u_pos = listEcc_pos[i];
		out_pos += (listEcc_pos[i+1]-listEcc_pos[i])*abs(distr_of_eccentricities(ecc1_pos, distribution1_pos, u_pos) 
						     - distr_of_eccentricities(ecc2_pos, distribution2_pos, u_pos));
	}

	// calculating the lower bound NEG
	for(unsigned int i = 0; i < listEcc_neg.size() - 1; i++)
	{
		u_neg = listEcc_neg[i];
		out_neg += (listEcc_neg[i+1]-listEcc_neg[i])*abs(distr_of_eccentricities(ecc1_neg, distribution1_neg, u_neg) 
						     - distr_of_eccentricities(ecc2_neg, distribution2_neg, u_neg));
	}



	// calculating the lower bound POS and NEG
	for(unsigned int i = 0; i < listEcc_pos_and_neg.size() - 1; i++)
	{
		u_pos_and_neg = listEcc_pos_and_neg[i];
		out_pos_and_neg += (listEcc_pos_and_neg[i+1]-listEcc_pos_and_neg[i])
					*abs(
						distr_of_eccentricities(ecc1_pos_and_neg, distribution1_pos_and_neg, u_pos_and_neg) 
						- distr_of_eccentricities(ecc2_pos_and_neg, distribution2_pos_and_neg, u_pos_and_neg)
					);
	}

	/*
	cout << out_pos << endl;
	cout << out_neg << endl;
	cout << out_pos_and_neg << endl;
	*/


	return(0.5 * (c1*out_pos + c2*out_neg + c3*out_pos_and_neg));
}


double calc_flb_sorted(	vector<double> &ecc1, vector<double> &distribution1,
						vector<double> &ecc2, vector<double> &distribution2,
						vector<double> &listEcc){
	/*
	 * first sorts the eccentricities.
	 * permutates distribution in the same way.
	 *
	 * then calculates the integral efficiently
	 */


	sort_ecc_and_dist(ecc1,distribution1);
	sort_ecc_and_dist(ecc2,distribution2);

	double out = 0.0;

	unsigned int last_index1 = 0;
	double last_integral1 = 0.0;

	unsigned int last_index2 = 0;
	double last_integral2 = 0.0;

	double u = 0;
	// calculating the lower bound POS and NEG sorted
	for(unsigned int i = 0; i < listEcc.size() - 1; i++)
	{
		u = listEcc[i];

		out += (listEcc[i+1]-listEcc[i])
					*abs(
						distr_of_eccentricities_sorted(ecc1, distribution1, u, last_index1, last_integral1)
						- distr_of_eccentricities_sorted(ecc2, distribution2, u, last_index2, last_integral2)
					);
	}

	return out;
}

//double flb_with_relation(	vector< vector <double> > &points1_pos, vector<double> &distribution1_pos,
//							vector< vector <double> > &points1_neg, vector<double> &distribution1_neg,
//							vector< vector <double> > &points2_pos, vector<double> &distribution2_pos,
//							vector< vector <double> > &points2_neg, vector<double> &distribution2_neg,
//							int number_of_selected_points,
//							double c1, double c2, double c3,
//							vector<double> &solution)
//{
//	/*
//		flb:
//		Calculates the first lower bound. See section 3.4.1 and 4.2.1 "First Lower Bound" in
//			"Quantitative comparison of protein isosurfaces
// 			with approximated Gromov-Wasserstein-distance" from Felix Berens.
//
//		Input
//		-----
//		points1 - vector with dimension number_of_selected_points x 3. Points of the first object as 3D coordinates
//		points2 - vector with dimension number_of_selected_points x 3. Points of the first object as 3D coordinates
//		distribution1 - vector containing the weight for the points of the first object
//		distribution2 - vector containing the weight for the points of the second object
//		number_of_selected_points - a integer, for how many points the first lower bound will be calculated
//
//		Output
//		------
//		Is a double that is the first lower bound.
//
//		Complexity
//    		----------
//		O(n^2), where n = number_of_selected_points
//
//		Example
//		-------
//		vector < vector<double> > points1{{0,0,1},
//						{1,0,0},
//						{0,1,0}};
//		vector < vector<double> > points2{{2,0,1},
//						{1,0,2},
//						{0,1,3}};
//		vector <double> distribution1{0.4,0.2,0.2};
//		vector <double> distribution2{0.25,0.5,0.25};
//		int number_of_selected_points = 3;
//		flb(points1, points2, distribution1, distribution2, number_of_selected_points);
//
//		=> 0.372717
//	*/
//
////	// positive and negative
////	vector < vector<double> > d1_pos(number_of_selected_points, vector<double>(number_of_selected_points,0));
////	vector < vector<double> > d2_pos(number_of_selected_points, vector<double>(number_of_selected_points,0));
////
////	vector < vector<double> > d1_neg(number_of_selected_points, vector<double>(number_of_selected_points,0));
////	vector < vector<double> > d2_neg(number_of_selected_points, vector<double>(number_of_selected_points,0));
////
////	// merging the positive and negative of one protein
////	vector< vector <double> > points1_pos_and_neg(3,vector<double>(number_of_selected_points*2));
////
////	appendPoints(points1_pos, points1_neg, points1_pos_and_neg);
////
////	vector< vector <double> > points2_pos_and_neg(3,vector<double>(number_of_selected_points*2));
////	appendPoints(points2_pos, points2_neg, points2_pos_and_neg);
////
////	vector <double> distribution1_pos_and_neg(number_of_selected_points*2,1/(number_of_selected_points*2));
////	appendDistributions(distribution1_pos, distribution1_neg, distribution1_pos_and_neg);
////
////	vector <double> distribution2_pos_and_neg(number_of_selected_points*2,1/(number_of_selected_points*2));
////	appendDistributions(distribution2_pos, distribution2_neg, distribution2_pos_and_neg);
////
////
////	// and for the relation between positive and negative
////	vector < vector<double> > d1_pos_and_neg(2*number_of_selected_points, vector<double>(2*number_of_selected_points,0));
////	vector < vector<double> > d2_pos_and_neg(2*number_of_selected_points, vector<double>(2*number_of_selected_points,0));
////
////
////	vector<double> ecc1_pos(number_of_selected_points);
////	vector<double> ecc2_pos(number_of_selected_points);
////
////	vector<double> ecc1_neg(number_of_selected_points);
////	vector<double> ecc2_neg(number_of_selected_points);
////
////	// and for the relation between the potentials within one protein
////	vector<double> ecc1_pos_and_neg(2*number_of_selected_points);
////	vector<double> ecc2_pos_and_neg(2*number_of_selected_points);
////
////	vector<double> listEcc_pos(number_of_selected_points+number_of_selected_points);
////	vector<double> listEcc_neg(number_of_selected_points+number_of_selected_points);
////
////	// and for the relation between the potentials within one protein
////	vector<double> listEcc_pos_and_neg(number_of_selected_points*4);
//
//	//------------------------------------------same as above
//
////	// positive and negative
////	vector < vector<double> > d1_pos(number_of_selected_points, vector<double>(number_of_selected_points,0));
////	vector < vector<double> > d2_pos(number_of_selected_points, vector<double>(number_of_selected_points,0));
////
////	vector < vector<double> > d1_neg(number_of_selected_points, vector<double>(number_of_selected_points,0));
////	vector < vector<double> > d2_neg(number_of_selected_points, vector<double>(number_of_selected_points,0));
////
////	// merging the positive and negative of one protein
////	vector< vector <double> > points1_pos_and_neg(3,vector<double>(number_of_selected_points*2));
////
////	vector< vector <double> > points2_pos_and_neg(3,vector<double>(number_of_selected_points*2));
////	appendPoints(points2_pos, points2_neg, points2_pos_and_neg);
////
////	vector <double> distribution1_pos_and_neg(number_of_selected_points*2,1/(number_of_selected_points*2));
////	vector <double> distribution2_pos_and_neg(number_of_selected_points*2,1/(number_of_selected_points*2));
////
////	// and for the relation between positive and negative
////	vector < vector<double> > d1_pos_and_neg(2*number_of_selected_points, vector<double>(2*number_of_selected_points,0));
////	vector < vector<double> > d2_pos_and_neg(2*number_of_selected_points, vector<double>(2*number_of_selected_points,0));
////
////	vector<double> ecc1_pos(number_of_selected_points);
////	vector<double> ecc2_pos(number_of_selected_points);
////
////	vector<double> ecc1_neg(number_of_selected_points);
////	vector<double> ecc2_neg(number_of_selected_points);
////
////	// and for the relation between the potentials within one protein
////	vector<double> ecc1_pos_and_neg(2*number_of_selected_points);
////	vector<double> ecc2_pos_and_neg(2*number_of_selected_points);
////
////	vector<double> listEcc_pos(number_of_selected_points+number_of_selected_points);
////	vector<double> listEcc_neg(number_of_selected_points+number_of_selected_points);
////
////	// and for the relation between the potentials within one protein
////	vector<double> listEcc_pos_and_neg(number_of_selected_points*4);
//
//	appendPoints(points1_pos, points1_neg, points1_pos_and_neg);
//	appendPoints(points2_pos, points2_neg, points2_pos_and_neg);
//	appendDistributions(distribution1_pos, distribution1_neg, distribution1_pos_and_neg);
//	appendDistributions(distribution2_pos, distribution2_neg, distribution2_pos_and_neg);
//
//
//	dist_matrix(d1_pos_and_neg, points1_pos_and_neg);
//	dist_matrix(d2_pos_and_neg, points2_pos_and_neg);
//
//
//    submatrix(d1_pos_and_neg, 0, number_of_selected_points, 0, number_of_selected_points, d1_pos);
//    submatrix(d2_pos_and_neg, 0, number_of_selected_points, 0, number_of_selected_points, d2_pos);
//
//    submatrix(d1_pos_and_neg, number_of_selected_points, 2*number_of_selected_points, number_of_selected_points, 2*number_of_selected_points, d1_neg);
//    submatrix(d2_pos_and_neg, number_of_selected_points, 2*number_of_selected_points, number_of_selected_points, 2*number_of_selected_points, d2_neg);
//
//
//	// calculation of all eccentricities
//	all_eccentricities(d1_pos, distribution1_pos,ecc1_pos);
//	all_eccentricities(d2_pos, distribution2_pos,ecc2_pos);
//
//	all_eccentricities(d1_neg, distribution1_neg,ecc1_neg);
//	all_eccentricities(d2_neg, distribution2_neg,ecc2_neg);
//
//	all_eccentricities(d1_pos_and_neg, distribution1_pos_and_neg, ecc1_pos_and_neg);
//	all_eccentricities(d2_pos_and_neg, distribution2_pos_and_neg, ecc2_pos_and_neg);
//
//
//	// mergin and sorting the list of eccentricities
//	list_of_unique_eccent(ecc1_pos, ecc2_pos, listEcc_pos);
//	list_of_unique_eccent(ecc1_neg, ecc2_neg, listEcc_neg);
//
//	list_of_unique_eccent(ecc1_pos_and_neg, ecc2_pos_and_neg, listEcc_pos_and_neg);
//
//
//	double out_pos = calc_flb_sorted(ecc1_pos, distribution1_pos, ecc2_pos, distribution2_pos, listEcc_pos);
//	double out_neg = calc_flb_sorted(ecc1_neg, distribution1_neg, ecc2_neg, distribution2_neg, listEcc_neg);
//	double out_pos_and_neg = calc_flb_sorted(ecc1_pos_and_neg, distribution1_pos_and_neg, ecc2_pos_and_neg, distribution2_pos_and_neg, listEcc_pos_and_neg);
//
//
//
//	solution[0] = out_pos;
//	solution[1] = out_neg;
//	solution[2] = out_pos_and_neg;
//
//	//return(0.5 * (c1*out_pos + c2*out_neg + c3*out_neg_vs_pos + c4*out_pos_vs_neg));
//	return(0.5 * (c1*out_pos + c2*out_neg + c3*out_pos_and_neg));
//}


double flb_with_relation(thread_data *m, vector<double> &solution)
{
	/*
		flb:
		Calculates the first lower bound. See section 3.4.1 and 4.2.1 "First Lower Bound" in
			"Quantitative comparison of protein isosurfaces
 			with approximated Gromov-Wasserstein-distance" from Felix Berens.

		Input
		-----
		points1 - vector with dimension number_of_selected_points x 3. Points of the first object as 3D coordinates
		points2 - vector with dimension number_of_selected_points x 3. Points of the first object as 3D coordinates
		distribution1 - vector containing the weight for the points of the first object
		distribution2 - vector containing the weight for the points of the second object
		number_of_selected_points - a integer, for how many points the first lower bound will be calculated

		Output
		------
		Is a double that is the first lower bound.

		Complexity
    		----------
		O(n^2), where n = number_of_selected_points

		Example
		-------
		vector < vector<double> > points1{{0,0,1},
						{1,0,0},
						{0,1,0}};
		vector < vector<double> > points2{{2,0,1},
						{1,0,2},
						{0,1,3}};
		vector <double> distribution1{0.4,0.2,0.2};
		vector <double> distribution2{0.25,0.5,0.25};
		int number_of_selected_points = 3;
		flb(points1, points2, distribution1, distribution2, number_of_selected_points);

		=> 0.372717
	*/


	appendPoints(m->points1_pos, m->points1_neg, m->points1_pos_and_neg);
	appendPoints(m->points2_pos, m->points2_neg, m->points2_pos_and_neg);
	appendDistributions(m->distribution1_pos, m->distribution1_neg, m->distribution1_pos_and_neg);
	appendDistributions(m->distribution2_pos, m->distribution2_neg, m->distribution2_pos_and_neg);


	dist_matrix(m->d1_pos_and_neg, m->points1_pos_and_neg);
	dist_matrix(m->d2_pos_and_neg, m->points2_pos_and_neg);


    submatrix(m->d1_pos_and_neg, 0, m->number_of_selected_points, 0, m->number_of_selected_points, m->d1_pos);
    submatrix(m->d2_pos_and_neg, 0, m->number_of_selected_points, 0, m->number_of_selected_points, m->d2_pos);

    submatrix(m->d1_pos_and_neg, m->number_of_selected_points, 2*m->number_of_selected_points, m->number_of_selected_points, 2*m->number_of_selected_points, m->d1_neg);
    submatrix(m->d2_pos_and_neg, m->number_of_selected_points, 2*m->number_of_selected_points, m->number_of_selected_points, 2*m->number_of_selected_points, m->d2_neg);


	// calculation of all eccentricities
	all_eccentricities(m->d1_pos, m->distribution1_pos,m->ecc1_pos);
	all_eccentricities(m->d2_pos, m->distribution2_pos,m->ecc2_pos);

	all_eccentricities(m->d1_neg, m->distribution1_neg,m->ecc1_neg);
	all_eccentricities(m->d2_neg, m->distribution2_neg,m->ecc2_neg);

	all_eccentricities(m->d1_pos_and_neg, m->distribution1_pos_and_neg, m->ecc1_pos_and_neg);
	all_eccentricities(m->d2_pos_and_neg, m->distribution2_pos_and_neg, m->ecc2_pos_and_neg);


	// mergin and sorting the list of eccentricities
	list_of_unique_eccent(m->ecc1_pos, m->ecc2_pos, m->listEcc_pos);
	list_of_unique_eccent(m->ecc1_neg, m->ecc2_neg, m->listEcc_neg);

	list_of_unique_eccent(m->ecc1_pos_and_neg, m->ecc2_pos_and_neg, m->listEcc_pos_and_neg);


	double out_pos = calc_flb_sorted(m->ecc1_pos, m->distribution1_pos, m->ecc2_pos, m->distribution2_pos, m->listEcc_pos);
	double out_neg = calc_flb_sorted(m->ecc1_neg, m->distribution1_neg, m->ecc2_neg, m->distribution2_neg, m->listEcc_neg);
	double out_pos_and_neg = calc_flb_sorted(m->ecc1_pos_and_neg, m->distribution1_pos_and_neg, m->ecc2_pos_and_neg, m->distribution2_pos_and_neg, m->listEcc_pos_and_neg);



	solution[0] = out_pos;
	solution[1] = out_neg;
	solution[2] = out_pos_and_neg;

	//return(0.5 * (c1*out_pos + c2*out_neg + c3*out_neg_vs_pos + c4*out_pos_vs_neg));
	return(0.5 * (m->c1*out_pos + m->c2*out_neg + m->c3*out_pos_and_neg));
}


double flb_with_relation_only_pos(thread_data *m, vector<double> &solution)
{
	/*
		flb:
		Calculates the first lower bound. See section 3.4.1 and 4.2.1 "First Lower Bound" in
			"Quantitative comparison of protein isosurfaces
 			with approximated Gromov-Wasserstein-distance" from Felix Berens.

		Input
		-----
		points1 - vector with dimension number_of_selected_points x 3. Points of the first object as 3D coordinates
		points2 - vector with dimension number_of_selected_points x 3. Points of the first object as 3D coordinates
		distribution1 - vector containing the weight for the points of the first object
		distribution2 - vector containing the weight for the points of the second object
		number_of_selected_points - a integer, for how many points the first lower bound will be calculated

		Output
		------
		Is a double that is the first lower bound.

		Complexity
    		----------
		O(n^2), where n = number_of_selected_points

		Example
		-------
		vector < vector<double> > points1{{0,0,1},
						{1,0,0},
						{0,1,0}};
		vector < vector<double> > points2{{2,0,1},
						{1,0,2},
						{0,1,3}};
		vector <double> distribution1{0.4,0.2,0.2};
		vector <double> distribution2{0.25,0.5,0.25};
		int number_of_selected_points = 3;
		flb(points1, points2, distribution1, distribution2, number_of_selected_points);

		=> 0.372717
	*/



	dist_matrix(m->d1_pos, m->points1_pos);
	dist_matrix(m->d2_pos, m->points2_pos);

	// calculation of all eccentricities
	all_eccentricities(m->d1_pos, m->distribution1_pos,m->ecc1_pos);
	all_eccentricities(m->d2_pos, m->distribution2_pos,m->ecc2_pos);


	// mergin and sorting the list of eccentricities
	list_of_unique_eccent(m->ecc1_pos, m->ecc2_pos, m->listEcc_pos);

	double out_pos = calc_flb_sorted(m->ecc1_pos, m->distribution1_pos, m->ecc2_pos, m->distribution2_pos, m->listEcc_pos);

	solution[0] = out_pos;
	solution[1] = 0;
	solution[2] = 0;

	//return(0.5 * (c1*out_pos + c2*out_neg + c3*out_neg_vs_pos + c4*out_pos_vs_neg));
	return(0.5 * (m->c1*out_pos));
}



void calculateDistancesToActiveCenter(const vector< vector <double> > &points, vector<double> &distances, 
                                        const vector<double> &actCenter){
    /*
    *   Calculate the distances to the active center of all points
    */

    for(unsigned int i = 0; i < points[0].size(); i++){
        distances[i] = euclidean_distance(points[0][i], points[1][i], points[2][i], actCenter[0], actCenter[1], actCenter[2]);
    }
}


void writeMeasureToFile(const vector<double> &distribution, const string fName){

    ofstream myfile;
    myfile.open (fName);
	for(unsigned int i=0; i < distribution.size(); i++)
	{
		myfile << distribution[i] << "\n";
	}

    myfile.close();
}

void getMeasureWithDistToAC(  const vector< vector <double> > &points_pos, vector<double> &distribution_pos, 
                                    const vector< vector <double> > &points_neg, vector<double> &distribution_neg, 
                                    const vector<double> &actCenter,
                                    const double closest){
    /*
    *   Calculate the weights according to the distance to the active center.
    *   Points closer to the active center get a higher weight.
    *   closest ... defines how much more times the weight is in the closest point
    *               compared to the farthest point
    */

    double furthest = 1.0;

    distribution_pos.resize(points_pos[0].size());
    distribution_neg.resize(points_neg[0].size());

    vector <double> distances_pos(distribution_pos.size());
    calculateDistancesToActiveCenter(points_pos, distances_pos, actCenter);

    vector <double> distances_neg(distribution_neg.size());
    calculateDistancesToActiveCenter(points_neg, distances_neg, actCenter);

    double max_distance = max(distances_pos);
    if(max(distances_neg) > max_distance) max_distance = max(distances_neg);

    double min_distance = min(distances_pos);
    if(min(distances_neg) > min_distance) min_distance = min(distances_neg);

    double dx =  max_distance - min_distance;
    double dy =  furthest - closest;
    double m = dy/dx;
    
    for(unsigned int i = 0; i < distribution_pos.size(); i++){
//        distribution_pos[i] = 1.1*max_distance - distances_pos[i];
        distribution_pos[i] = closest + m*(distances_pos[i]-min_distance);
//        distribution_pos[i] = 1.0/distances_pos[i];

//    	distribution_pos[i] = distances_pos[i];
    }

    for(unsigned int i = 0; i < distribution_neg.size(); i++){
//        distribution_neg[i] = 1.1*max_distance - distances_neg[i];
        distribution_neg[i] = closest + m*(distances_neg[i]-min_distance);
//        distribution_neg[i] = 1.0/distances_neg[i];

//    	distribution_neg[i] = distances_neg[i];
    }

    double sum_distribution = sum(distribution_pos) + sum(distribution_neg);

    scale(distribution_pos,1.0/sum_distribution);
    scale(distribution_neg,1.0/sum_distribution);

}

bool readActiveCenter(vector<double> &actCenter, const string active_centre_fileName){
//	active_centre_fileName

	int lc = 0;
	  string line;
	  ifstream myfile (active_centre_fileName);
	  if (myfile.is_open())
	  {
	    while ( getline (myfile,line) )
	    {

//	    	cout << line << endl;

	    	if(lc > 0){
	  	      // tokenize
	  	    	
                vector<string> strs;
	  	    	boost::split(strs,line,boost::is_any_of(";"));

	  	    	//setActiveCenter(Vec(std::stod(strs[0]), std::stod(strs[1]), std::stod(strs[2])));

                actCenter[0] = std::stod(strs[0]);
                actCenter[1] = std::stod(strs[1]);
                actCenter[2] = std::stod(strs[2]);

	    	}

	    	lc++;

	    }
	    myfile.close();

        return true;
	  }

	  else cerr << "Unable to open file " << active_centre_fileName << endl;

  return false;
}

bool readMeasure(const string Dataname_mu, vector<double> &distribution){
    /* reads in the measure from the file
    * make sure that distribution.size() == points[0].size()
    */

    string line;
    ifstream myfile (Dataname_mu);

    unsigned int expected_size = distribution.size();

    unsigned int ind = 0;
    if (myfile.is_open())
    {
        while ( getline (myfile,line) )
        {
        	if(expected_size <= ind) return false;

            distribution[ind] = std::stod(line);
            ind++;
        }
        myfile.close();

        if(expected_size != ind) return false;
        return true;
    }

    else cerr << "Unable to open file " << Dataname_mu << endl; 

    return false;
}

void getMeasureWithDistToAC(  const vector< vector <double> > &points_pos, vector<double> &distribution_pos, 
                                    const vector< vector <double> > &points_neg, vector<double> &distribution_neg, 
                                    const string Dataname1_mu_pos, 
                                    const string Dataname1_mu_neg, 
                                    const string Dataname1_actCent,
                                    const double closest){

//	cout << points_pos[0].size() << endl;
//	cout << points_neg[0].size() << endl;
//	cout << distribution_pos.size() << endl;
//	cout << distribution_neg.size() << endl;

    // check if measurefile allready exists, if it does read it
    bool pos_ok = readMeasure(Dataname1_mu_pos,distribution_pos);
    bool neg_ok = readMeasure(Dataname1_mu_neg,distribution_neg);

    //bool sizes_ok = (points_pos[0].size() == distribution_pos.size()) && (points_neg[0].size() == distribution_neg.size());

    //cout << (!pos_exists || !neg_exists || !sizes_ok) << " " << pos_exists << " " << neg_exists << " " << sizes_ok << endl;

    // else calculate measure and create the file
//    if(!pos_ok || !neg_ok){
//        cout << "reading active center" << endl;
        vector<double> actCenter(3,0);
        readActiveCenter(actCenter, Dataname1_actCent);
        
//        cout << "creating measure" << endl;
        getMeasureWithDistToAC(points_pos, distribution_pos, points_neg, distribution_neg, actCenter, closest);

        writeMeasureToFile(distribution_pos,Dataname1_mu_pos);
        writeMeasureToFile(distribution_neg,Dataname1_mu_neg);
//    }
}

#endif /* FIRST_LOWER_BOUND_H_ */
