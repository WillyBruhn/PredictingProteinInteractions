/*
 * mathFunctions.h
 *
 *  Created on: 26.03.2018
 *      Author: Berens
 *
 *		Simple mathematics functions that are needed
 */

#ifndef MATHFUNCTIONS_H_
#define MATHFUNCTIONS_H_

#include <cmath>
#include <iomanip>

#include <assert.h>

#include <algorithm>

#include <iostream>

#include <vector>
#include "mathFunctions.h"

using namespace std;


double euclidean_distance(double x1,double y1, double z1,double x2,double y2, double z2)
{
	/*
		euclidean_distance:
		Calculates the euclidean distance for 2 3D points

		Input
		-----
		x1, y1, z1 - are the coordinates of the first point, all types are double
		x2, y2, z2 - are the coordinates of the second point, all types are double

		Output
		------
		The output is a double, which is the euclidean distance between (x1,y1,z1) and (x2,y2,z2)

		Complexity
    		----------
		O(1)

		Example
		-------
		double x1 = 0;
		double y1 = 0;
		double z1 = 0;
		double x2 = 1;
		double y2 = 1;
		double z2 = 1;
		euclidean_distance(x1,y1,z1,x2,y2,z2);

		=> 1.73205
	*/
//	return(sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2)));
	return(sqrt((x1-x2)*(x1-x2)+ (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)));
}

void dist_matrix(vector < vector<double> > &d, const vector< vector <double> > &points )
{
	/*
		dist_matrix:
		Calculates a distance matrix for given 3D points

		Input
		-----
		d - is from type vector, which represent a matrix, with dimension
			number_of_selected_points x number_of_selected_points.
			The distance matrix will be stored there.
		points - is from type vector, with dimenson number_of_selected_points x 3.
			The points for which the distance matrix will be calculated are stored here.

		Output
		------
		No return

		Complexity
    		----------
		O(number_of_selected_points)

		Example
		-------
		vector < vector<double> > d(3, vector<double>(3,0));
		vector < vector<double> > points{ { 1, 0, 0 },
						  { 0, 0, 1 },
						  { 0, 1, 0 } };
		dist_matrix(d,points);

			(0	 1.41421 1.41421)
		=> d =	(1.41421 0	 1.41421)
			(1.41421 1.41421 0	)
	*/
	for(unsigned int i=0; i < d.size(); i++)
	{
		for(unsigned int j=i;j< d.size(); j++)
		{
			d[i][j] = euclidean_distance(points[0][i],points[1][i],points[2][i],points[0][j],points[1][j],points[2][j]);
			d[j][i] = d[i][j];
		}
	}
}


void submatrix(	const vector < vector<double> > &d,
					const unsigned int i1, const unsigned int i2,
					const unsigned int j1, const unsigned int j2,
                    vector < vector<double> > &d_sub)
{
	/* select a matrix with the rows i1 to i2 and columns j1 to j2 of d
	*/

	unsigned int rows = i2-i1;
	unsigned int cols = j2-j1;

	//vector < vector<double> > ret(rows,vector<double>(cols));


	for(unsigned int i=0; i < rows; i++)
	{
		for(unsigned int j=0;j< cols; j++)
		{
			d_sub[i][j] = d[i1+i][j1+j];
		}
	}

	//return(ret);
}


void printMatrix(const vector < vector<double> > &d, const string name = "matrix"){
	/*
	* Willy 30.1.2019
	*/
	cout << "--------------------------------------" << endl;
	cout << name << endl;
	for(unsigned int i=0; i < d.size(); i++)
	{
		for(unsigned int j=0; j < d[i].size(); j++)
		{
			cout  <<  std::setw(8) << d[i][j] << std::setfill(' ') << " " ;
		}
		cout << "\n";
	}
	cout << "--------------------------------------" << endl;
}

void printPoints(const vector < vector<double> > &points, const string name = "points"){
	cout << "--------------------------------------" << endl;
	cout << name << endl;
	cout << points.size() << endl;
	cout << points[0].size() << endl;
	for(unsigned int i=0; i < points[0].size(); i++)
	{
		cout << points[0][i] << " " << points[1][i] << " " << points[2][i] << endl;
	}

	cout << "--------------------------------------" << endl;
}

//std::ostream& operator <<(std::ostream &os, vector<double> &point){
//	os << point[0] << " " << point[1] << " " << point[2] << endl;
//}

std::ostream& operator <<(std::ostream &os, vector<unsigned int> &hist){
	for(unsigned int i = 0; i < hist.size(); i++){
		os << hist[i] << "\n";
	}
	return os;
}

std::ostream& operator <<(std::ostream &os, vector<double> &hist){
	for(unsigned int i = 0; i < hist.size(); i++){
		os << hist[i] << "\n";
	}
	return os;
}

double sum(const vector<double> &v){
    double s = 0.0;
	for(unsigned int j=0; j < v.size(); j++)
	{
		s += v[j];
	}

    return s;
}

double max(const vector<double> &v){
    double max = v[0];
	for(unsigned int j=1; j < v.size(); j++)
	{
		if(v[j] > max) max = v[j];
	}

    return max;
}

double min(const vector<double> &v){
    double min = v[0];
	for(unsigned int j=1; j < v.size(); j++)
	{
		if(v[j] < min) min = v[j];
	}

    return min;
}

void scale(vector<double> &v, const double sc){
	for(unsigned int j=0; j < v.size(); j++)
	{
		v[j] = v[j]*sc;
	}
}


void normalize(vector<double> &v){
	scale(v,1.0/sum(v));
}

void hist(const vector<double> &Q, vector<unsigned int> &h, const unsigned int size, const double bin_width, const double min){
	/*
	 * create a histogram
	 */
	h.resize(size,0);
	unsigned int bin_idx = 0;
	for(unsigned int i = 0; i < Q.size(); i++){
//		bin_idx = (unsigned int)((Q[i] - min) / bin_width);
//		bin_idx = (unsigned int)floor(((Q[i] - min) / (bin_width+1)));

		bin_idx = (unsigned int)floor(((Q[i] - min)/bin_width));

		// this looks weird, but we need the compact intervals
		// we have to decide between < and > at the breaks
		// if we have one bin and two values then either of the values would be lost
		// therefore we use val < break_i at the breaks and at the last we use val <= break_last
		if(bin_idx == size) bin_idx = size-1;

//		cout << Q[i] << " " << min << " " << bin_width << " " << bin_idx << endl;

		h[bin_idx]++;
	}
}

void hist(	const vector<double> &P, vector<unsigned int> &P_hist,
			const vector<double> &Q, vector<unsigned int> &Q_hist,
			double &bin_width,
			const unsigned int num_bins){
	/*
	 * create two histograms
	 */

	double max = std::max(*max_element(std::begin(P), std::end(P)), *max_element(std::begin(Q), std::end(Q)));
	double min = std::min(*min_element(std::begin(P), std::end(P)), *min_element(std::begin(Q), std::end(Q)));

	bin_width = (max - min)/(num_bins);

	double dy = (num_bins-1)/(max-min);

//	for(unsigned int i = 0; i < num_bins; i++){
//		cout << bin_width*i + min << " | ";
//	}
//	cout << endl;
//
//
//	cout << max << " " << min << " " << num_bins << endl;
//	cout << "dy " << dy << endl;

	hist(P, P_hist, num_bins, bin_width, min);
	hist(Q, Q_hist, num_bins, bin_width, min);
}

double emd(const vector<double> &Q, const vector<double> &P, const unsigned int num_bins){

	// create histograms
	assert(Q.size() == P.size());
	vector<unsigned int> P_hist;
	vector<unsigned int> Q_hist;

	double bin_width;
	hist(P,P_hist, Q, Q_hist, bin_width, num_bins);

	// calculate emd
	vector<double> emd(P_hist.size()+1,0);
	for(unsigned int i = 0; i < emd.size()-1; i++){
		emd[i+1] = P_hist[i] + emd[i] - Q_hist[i];
	}

	double emd_val = 0.0;
	for(unsigned int i = 0; i < emd.size(); i++){
		emd_val += abs(emd[i]);
	}

	return emd_val/Q.size();
}

double emd2hist(vector<unsigned int> &P_hist,
				vector<unsigned int> &Q_hist,
				const unsigned int num_bins){
	/*
	 * calculate emd given two histograms
	 */

	// create histograms
	assert(Q_hist.size() == P_hist.size());

	// calculate emd
	vector<double> emd(P_hist.size()+1,0);
	for(unsigned int i = 0; i < emd.size()-1; i++){
		emd[i+1] = P_hist[i] + emd[i] - Q_hist[i];
	}

	double emd_val = 0.0;
	for(unsigned int i = 0; i < emd.size(); i++){
		emd_val += abs(emd[i]);
	}

	return emd_val/Q_hist.size();
}

#endif /* MATHFUNCTIONS_H_ */
