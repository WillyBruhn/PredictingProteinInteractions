#include <iostream>

#include "gtest/gtest.h"
#include "../../src/mathFunctions.h"
#include <bits/stdc++.h>

#include <ostream>

#include "../../src/random_points.h"



#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <array>

#include "../../src/sampling.h"

#include "boost/date_time/posix_time/posix_time.hpp"
//namespace pt = boost::posix_time;

using namespace std;


std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}


void writeToFile(const string fName, vector<unsigned int > &hist){
	ofstream myfile (fName);
	if (myfile.is_open())
	{
		myfile << hist;

		myfile.close();
	}
	else cout << "Unable to open file " << fName << endl;
}

void writeToFile(const string fName, vector<double > &hist){
	ofstream myfile (fName);
	if (myfile.is_open())
	{
		myfile << hist;

		myfile.close();
	}
	else cout << "Unable to open file " << fName << endl;
}

void testHist(vector<double> &Q, vector<double> &P, unsigned int num_bins){
//	vector<double> Q{12,3,4,5,1,2,5,2,2,2,5,5,5,5,52,23,3,4,4,5,3,53,4,4,34,35,3,46,4,5,3,4,23};
	vector<unsigned int> Q_hist;

//	vector<double> P{14,5,2,2,2,2,23,4,5,53,5,5,35,3,43,6,41,2,2,3,4};
	vector<unsigned int> P_hist;

//	unsigned int num_bins = 2;
	double bin_width;

	hist(Q,Q_hist, P, P_hist, bin_width, num_bins);

	string path = "/home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/ComparingProteins/LowerBounds/FirstLowerBoundRelationOfPosAndNeg/Tests/TestSrc/";

	string Q_fName_hist = path+"Q_hist.csv";
	string Q_fName_vector = path+"Q_vector.csv";
	writeToFile(Q_fName_vector, Q);
	writeToFile(Q_fName_hist, Q_hist);

	string P_fName_hist = path+"P_hist.csv";
	string P_fName_vector = path+"P_vector.csv";
	writeToFile(P_fName_vector, P);
	writeToFile(P_fName_hist, P_hist);

	// compare with R histogram

	double emd_val = emd(Q,P,num_bins);


	string call = path + "./histTest.R "
						+ Q_fName_vector + " " + P_fName_vector + " "
						+ Q_fName_hist + " " + P_fName_hist + " "
						+ std::to_string(num_bins);

	string result = exec(call.c_str());

	cout << result << " " << endl;

	ASSERT_EQ(result,"[1] TRUE\n");



}


void testEmd(vector<double> &Q, vector<double> &P, unsigned int num_bins){
	vector<unsigned int> Q_hist;
	vector<unsigned int> P_hist;
	double bin_width;

	hist(Q,Q_hist, P, P_hist, bin_width, num_bins);

	string path = "/home/willy/RedoxChallenges/MasterThesis/IsoSurfSimilarity/ComparingProteins/LowerBounds/FirstLowerBoundRelationOfPosAndNeg/Tests/TestSrc/";

	string Q_fName_hist = path+"Q_hist.csv";
	string Q_fName_vector = path+"Q_vector.csv";
	writeToFile(Q_fName_vector, Q);
	writeToFile(Q_fName_hist, Q_hist);

	string P_fName_hist = path+"P_hist.csv";
	string P_fName_vector = path+"P_vector.csv";
	writeToFile(P_fName_vector, P);
	writeToFile(P_fName_hist, P_hist);

	// compare with R histogram

	double emd_val = emd(Q,P,num_bins);

	stringstream s;
	s << std::fixed << std::setprecision(4);

	s << emd_val;


	string call = path + "./emdTest.R "
						+ Q_fName_hist + " " + P_fName_hist + " "
						+ s.str();

//	cout << call << endl;

	string result = exec(call.c_str());

//	cout << "------------------------------------" << endl;
//	cout << result << endl;




//	if(result != "[1] TRUE\n"){
//			cout << call << endl;
//			cout << "------------------------------------" << endl;
//			cout << s.str() << " != " << result << endl;
//	}


	ASSERT_EQ(result,"[1] \"" + s.str() + "\"\n");

//	ASSERT_EQ(0,1);



}

TEST(hist, calculation)
{
	vector<double> Q{12,3,4,5,1,2,5,2,2,2,5,5,5,5,52,23,3,4,4,5,3,53,4,4,34,35,3,46,4,5,3,4,23};
	vector<double> P{14,5,2,2,2,2,23,4,5,53,5,5,35,3,43,6,41,2,2,3,4,5,4,16,12,23,13,14,25,18,12,10,11};
	testHist(Q,P,2);
	testHist(Q,P,3);
	testHist(Q,P,4);

	vector<double> Q2{1,2,2,3,6,7};

	testHist(Q2,Q2,3);

}

TEST(emd, calculation)
{
	vector<double> Q{12,3,4,5,1,2,5,2,2,2,5,5,5,5,52,23,3,4,4,5,3,53,4,4,34,35,3,46,4,5,3,4,23};
	vector<double> P{14,5,2,2,2,2,23,4,5,53,5,5,35,3,43,6,41,2,2,3,4,5,4,16,12,23,13,14,25,18,12,10,11};
	testEmd(Q,P,2);
	testEmd(Q,P,3);
	testEmd(Q,P,4);
	testEmd(Q,P,5);
	testEmd(Q,P,6);
	testEmd(Q,P,7);
	testEmd(Q,P,8);


}

TEST(emd, scaleInvariant)
{
	vector<double> Q{12,3,4,5,1,2,5,2,2,2,5,5,5,5,52,23,3,4,4,5,3,53,4,4,34,35,3,46,4,5,3,4,23};
	vector<double> P{14,5,2,2,2,2,23,4,5,53,5,5,35,3,43,6,41,2,2,3,4,5,4,16,12,23,13,14,25,18,12,10,11};

	unsigned int num_bins = 10;
	double sc = 5;
	double v = emd(Q,P,num_bins);

	vector<double> Q_scaled = Q;
	scale(Q_scaled,sc);
	vector<double> P_scaled = P;
	scale(P_scaled,sc);

	double v_scaled = emd(Q_scaled,P_scaled,num_bins);

	EXPECT_NEAR(v,v_scaled, 0.00000001);
}


void testHistScale(vector<double> &Q, vector<double> P, const unsigned int num_bins, const double sc){

//	unsigned int num_bins = 5;
//	double sc = 5;
	double v = emd(Q,P,num_bins);

	vector<unsigned int> P_hist;
	vector<unsigned int> Q_hist;

	double bin_width;
	hist(P,P_hist, Q, Q_hist, bin_width, num_bins);

	vector<double> Q_scaled = Q;
	scale(Q_scaled,sc);
	vector<double> P_scaled = P;
	scale(P_scaled,sc);


	vector<unsigned int> P_hist_scaled;
	vector<unsigned int> Q_hist_scaled;
	double bin_width_scaled;
	hist(P_scaled,P_hist_scaled, Q_scaled, Q_hist_scaled, bin_width_scaled, num_bins);


	for(unsigned int i = 0; i < P_hist_scaled.size(); i++){
//		cout << P_hist_scaled[i] << " =?= " << P_hist[i] << endl;

		EXPECT_EQ(P_hist_scaled[i],P_hist[i]);
	}

//	cout << "--------------------------------------------\n";
	for(unsigned int i = 0; i < Q_hist_scaled.size(); i++){
//		cout << Q_hist_scaled[i] << " =?= " << Q_hist[i] << endl;

		EXPECT_EQ(Q_hist_scaled[i],Q_hist[i]);
	}
}

TEST(hist, scaleInvariant)
{
	vector<double> Q{12,3,4,5,1,2,5,2,2,2,5,5,5,5,52,23,3,4,4,5,3,53,4,4,34,35,3,46,4,5,3,4,23};
	vector<double> P{14,5,2,2,2,2,23,4,5,53,5,5,35,3,43,6,41,2,2,3,4,5,4,16,12,23,13,14,25,18,12,10,11};

	testHistScale(Q,P,5,5);
	testHistScale(Q,P,5,50);

	testHistScale(Q,P,10,50);
	testHistScale(Q,P,100,50);

//	double v_scaled = emd(Q_scaled,P_scaled,num_bins);
}


TEST(random, samplePointsWithMeasure)
{

	vector<double> mu {10,1,1,2};
	unsigned int n = 2;
	vector<unsigned int> indices(n);

	vector<double> counts(mu.size());

	set_distribution_X_pos(mu);

	unsigned int m = 10000;
	for(unsigned int j = 0; j < m; j++){
		random_points_with_measure(distribution_X_pos,2,indices);
		for(unsigned int i = 0; i < indices.size(); i++){
			++counts[indices[i]];
		}
	}

	normalize(mu);
	normalize(counts);

	cout << "c      mu\n";
	for(unsigned int i = 0; i < counts.size(); i++){
			cout << counts[i]<< " " << mu[i] << "\n";
	}
	cout << endl;



	vector<double> mu2 {10,1};
	set_distribution_X_pos(mu2);


}

TEST(sampling, drawOnePoint)
{

	int n = 10;
	vector<vector<double>> points(3,vector<double>(n,0));
	for(unsigned int i = 0; i < n; i++){
		points[0][i] = i;
		points[1][i] = 0;
		points[2][i] = 0;
	}

	std::vector< std::vector<double> > d(points[0].size(), vector<double>(points[0].size(),0));
//	dist_matrix(d,points);

	vector<vector<double>> selected_points(3,vector<double>());
	vector<unsigned int> selected_points_indices;
	drawOnePoint(points,selected_points,selected_points_indices, 4,d);

	drawOnePoint(points,selected_points,selected_points_indices, 2,d);

	ASSERT_EQ(selected_points_indices[0],4);
	ASSERT_EQ(selected_points_indices[1],2);

	ASSERT_EQ(d[2][2],10000000000);
//	printMatrix(d);


}

TEST(min_element, euclidean)
{
	vector<double> points(10,0);
	for(unsigned int i = 0; i < 10; i++){
		points[i] = 10-i;
	}

	auto mm = std::min_element(points.begin(),points.end());
	unsigned int max_min_ind = std::distance(points.begin(), mm);

//	cout << "index is: " << max_min_ind << endl;

	ASSERT_EQ(max_min_ind,9);
}

void testEuclideanSampling(unsigned int K, vector<vector<double>> &points){
	vector<vector<double>> selected_points(3,vector<double>());
	vector<unsigned int> selected_points_indices = euclidean_farthest_point_sampling(points,selected_points,K);

	set<unsigned int> s( selected_points_indices.begin(), selected_points_indices.end() );
	selected_points_indices.assign( s.begin(), s.end() );

	ASSERT_EQ(s.size(),K);
}

void testEuclideanSamplingNPoints(unsigned int N, unsigned int K){
	int n = N;
	vector<vector<double>> points(3,vector<double>(n,0));
	for(unsigned int i = 0; i < n; i++){
		points[0][i] = i;
		points[1][i] = 0;
		points[2][i] = 0;
	}

	testEuclideanSampling(K,points);
}

TEST(sampling, euclidean)
{

	testEuclideanSamplingNPoints(3,1);
	testEuclideanSamplingNPoints(3,2);
	testEuclideanSamplingNPoints(3,3);

	testEuclideanSamplingNPoints(100,3);
	testEuclideanSamplingNPoints(100,50);

}


TEST(sampling, euclideanLargePc)
{

	unsigned int n = 10000;
	unsigned int k = 4000;

	// roughly takes 50 seconds

//	boost::posix_time::ptime tick = boost::posix_time::second_clock::local_time();
//	testEuclideanSamplingNPoints(n,k);

//	boost::posix_time::ptime now  = boost::posix_time::second_clock::local_time();
//
//	boost::posix_time::time_duration diff = tick - now;
//	cout << "sampling " << k << " points from " << n << " took: " << diff.total_milliseconds() << endl;

}

