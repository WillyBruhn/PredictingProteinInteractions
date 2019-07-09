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

using namespace std;

double euclidean_distance(double x1,double y1, double z1,double x2,double y2, double z2);

void dist_matrix(vector < vector<double> > &d, const vector< vector <double> > &points );


void submatrix(	const vector < vector<double> > &d, 
					const unsigned int i1, const unsigned int i2,
					const unsigned int j1, const unsigned int j2,
                    vector < vector<double> > &d_sub);


void printMatrix(const vector < vector<double> > &d, const string name = "matrix");

void printPoints(const vector < vector<double> > &points, const string name = "points");

//std::ostream& operator <<(std::ostream &os, vector<double> &point){
//	os << point[0] << " " << point[1] << " " << point[2] << endl;
//}

std::ostream& operator <<(std::ostream &os, vector<unsigned int> &hist);

std::ostream& operator <<(std::ostream &os, vector<double> &hist);

double sum(const vector<double> &v);

double max(const vector<double> &v);

double min(const vector<double> &v);

void scale(vector<double> &v, const double sc);


void normalize(vector<double> &v);

void hist(const vector<double> &Q, vector<unsigned int> &h, const unsigned int size, const double bin_width, const double min);

void hist(	const vector<double> &P, vector<unsigned int> &P_hist,
			const vector<double> &Q, vector<unsigned int> &Q_hist,
			double &bin_width,
			const unsigned int num_bins);

double emd(const vector<double> &Q, const vector<double> &P, const unsigned int num_bins = 100);

double emd2hist(vector<unsigned int> &P_hist,
				vector<unsigned int> &Q_hist,
				const unsigned int num_bins = 100);

#endif /* MATHFUNCTIONS_H_ */
