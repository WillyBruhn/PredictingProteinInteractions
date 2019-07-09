/*
* main.cpp
*
*  Created on: 07.05.2018
*      Author: Berens
*
* modified by Willy Bruhn on 30.1.2019
*/

#include <iostream>
#include <vector>
#include <time.h>
#include <pthread.h>
#include "first_lower_bound.h"
#include "reading.h"
#include "random_points.h"
#include "sampling.h"

#include <math.h>
#include <set>

#include <queue>

//#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

using namespace std;

#define NUM_THREADS 20

vector< vector<double> > out;
//vector<double> distribution1_pos;
//vector<double> distribution1_neg;
//
//vector<double> distribution2_pos;
//vector<double> distribution2_neg;
int rounds;
int number_of_selected_points;
int number_of_selected_points_half;

bool normalize_measure;

double c1;
double c2;
double c3;

double measure;

unsigned int num_bins;

string emd_list_id;		// if you want to run the program with the same settings multiple times
						// the resulting file will have an _id<emd_list_id>.csv ending

vector< vector<double> > outX_v_X_training;
vector< vector< vector<double> > > outX_v_Y_training;


void *Thread(void *threadarg); // Threadfunction for multithreading


//void hist(const vector<double> &Q, vector<unsigned int> &h, const unsigned int size, const double bin_width, const double min){
//	/*
//	 * create a histogram
//	 */
//	h.resize(size,0);
//	unsigned int bin_idx = 0;
//	for(unsigned int i = 0; i < Q.size(); i++){
////		bin_idx = (unsigned int)((Q[i] - min) / bin_width);
//		bin_idx = (unsigned int)floor(((Q[i] - min) / bin_width));
//
//		h[bin_idx]++;
//	}
//}
//
//double emd(const vector<double> &Q, const vector<double> &P, const unsigned int num_bins = 100){
//
//	// create histograms
//	assert(Q.size() == P.size());
//	vector<unsigned int> P_hist;
//	vector<unsigned int> Q_hist;
//
////	cout << "emd" << endl;
//
//	double max = std::max(*max_element(std::begin(P), std::end(P)), *max_element(std::begin(Q), std::end(Q)));
//	double min = std::min(*min_element(std::begin(P), std::end(P)), *min_element(std::begin(Q), std::end(Q)));
//
////	cout << max << " " << min << endl;
//
////	unsigned int num_bins = (int)(max - min+ 1);
//	double bin_width = (max - min)/(num_bins-1);
//
//	hist(P, P_hist, num_bins, bin_width, min);
//	hist(Q, Q_hist, num_bins, bin_width, min);
//
////	cout << "ok" << endl;
//
//	// calculate emd
//	vector<double> emd(P_hist.size(),0);
//	for(unsigned int i = 1; i < emd.size(); i++){
//		emd[i] = P_hist[i] + emd[i-1] - Q_hist[i];
//	}
//
//	double emd_val = 0.0;
//	for(unsigned int i = 0; i < emd.size(); i++){
//		emd_val += abs(emd[i]);
//	}
//
//	return emd_val*bin_width;
//}


void getKnnIndices(const vector<double> &distances, vector< vector<double> > &points, unsigned int const k){
//	cout << "getKnnIndices" << endl;

	std::priority_queue<std::pair<double, int>> q;
	for (int i = 0; i < distances.size(); ++i) {
	  q.push(std::pair<double, int>(-distances[i], i));
	}

//	vector< vector<double> > points_copy = points;
	vector< vector<double> > points_copy(3,vector<double>(k,0));

	for (int i = 0; i < k; ++i) {
	  int ki = q.top().second;
//	  std::cout << "index[" << i << "] = " << ki << std::endl;

//	  indices[i] = ki;

	  points_copy[0][i] = points[0][ki];
	  points_copy[1][i] = points[1][ki];
	  points_copy[2][i] = points[2][ki];

	  q.pop();
	}

//	cout << "hi" << endl;

	points = points_copy;

//	return indices;
}

void keepOnlyClosestPointsToActCent(const int NNtoActCenter,
									const string &Dataname1_actCent,
									vector< vector<double> > &points_pos,
									vector< vector<double> > &points_neg){

    vector<double> actCenter(3,0);
    readActiveCenter(actCenter, Dataname1_actCent);


    if(NNtoActCenter < points_pos[0].size()){
        vector <double> distances_pos(points_pos[0].size());
        calculateDistancesToActiveCenter(points_pos, distances_pos, actCenter);

        // points_pos are reduced to the closest points to the active center
        getKnnIndices(distances_pos, points_pos, NNtoActCenter);
    }

    if(NNtoActCenter < points_neg[0].size()){
		vector <double> distances_neg(points_neg[0].size());
		calculateDistancesToActiveCenter(points_neg, distances_neg, actCenter);

		// points_neg are reduced to the closest points to the active center
		getKnnIndices(distances_neg, points_neg, NNtoActCenter);
    }

}

void readPointsAndMeasue(	const string &Dataname1_pos,
							vector< vector<double> > &points1_pos_all,
							const string &Dataname1_neg,
							vector< vector<double> > &points1_neg_all,
							const string &Dataname1_mu_pos,
							vector<double> &distribution1_pos_all,
							const string &Dataname1_mu_neg,
							vector<double> &distribution1_neg_all,
							const string &Dataname1_actCent,
							const int measure,
							const int NNtoActCenter = 0
							){

	points1_pos_all.clear();
	points1_pos_all.resize(3);

	points1_neg_all.clear();
	points1_neg_all.resize(3);


	reading(Dataname1_pos, points1_pos_all); // reads the first file containing the 3D points
	if(points1_pos_all[0].size() == 0){exit(-1);}

	reading(Dataname1_neg, points1_neg_all); // reads the first file containing the 3D points
	if(points1_neg_all[0].size() == 0){exit(-1);}

	if(NNtoActCenter != 0){
		// we only keep the closest points to the active center
		keepOnlyClosestPointsToActCent(NNtoActCenter, Dataname1_actCent, points1_pos_all, points1_neg_all);
	}


	distribution1_pos_all.resize(points1_pos_all[0].size(), 1.0/(points1_pos_all[0].size()+points1_neg_all[0].size()));
	distribution1_neg_all.resize(points1_neg_all[0].size(), 1.0/(points1_pos_all[0].size()+points1_neg_all[0].size()));


	// in the case of (measure == 1) we are assuming a uniform-distribution. Hence we are done here.
	if(measure != 1){
		getMeasureWithDistToAC(points1_pos_all, distribution1_pos_all, points1_neg_all, distribution1_neg_all,
								Dataname1_mu_pos, Dataname1_mu_neg, Dataname1_actCent, measure);

	}
    //        getMeasureWithDistToAC(points1_pos_all, distribution1_pos_all, points1_neg_all, distribution1_neg_all,
    //                                Dataname1_mu_pos, Dataname1_mu_neg, Dataname1_actCent, measure);

}

/*
 * Creates multiple threads.
 * In each thread the flb is calculated multiple times.
 * vals_X_v_X ... holds the resulting distances
 */
void sampleDistances(	vector<double> &vals_X_v_X,
						thread_data *td,
						pthread_t *threads,
						unsigned int rounds,

						vector< vector <double> > &points1_pos_all,
						vector< vector <double> > &points1_neg_all,
						vector<double> &distribution1_pos_all,
						vector<double> distribution1_neg_all,

						vector< vector <double> > &points2_pos_all,
						vector< vector <double> > &points2_neg_all,
						vector<double> &distribution2_pos_all,
						vector<double> distribution2_neg_all,

						const double c1, const double c2, const double c3){

	out.resize(rounds, vector<double>(3,-1));

	for(int i = 0; i < NUM_THREADS; i++ )
	{
		td[i].thread_id = i;
		td[i].points1_pos_all = points1_pos_all;
		td[i].points2_pos_all = points2_pos_all;

		td[i].points1_neg_all = points1_neg_all;
		td[i].points2_neg_all = points2_neg_all;

	    td[i].distribution1_pos_all = distribution1_pos_all;
	    td[i].distribution1_neg_all = distribution1_neg_all;

	    td[i].distribution2_pos_all = distribution2_pos_all;
	    td[i].distribution2_neg_all = distribution2_neg_all;

		int rc = pthread_create(&threads[i], NULL, Thread, (void *)&td[i]);

		if (rc) {
			cerr << "Error:unable to create thread," << rc << endl;
			exit(-1);
		}
	}

	for (int i = 0; i < NUM_THREADS; i++) {
		pthread_join(threads[i], NULL);
	}

	for (int i = 0; i < out.size(); i++) {
		vals_X_v_X[i] = out[i][0]*c1 + out[i][1]*c2 + out[i][2]*c3;
	}

}

bool readListOfProteinsToCompareWith(const string fName, vector<string> &protein_names){
    string line;
    ifstream myfile (fName);

    unsigned int ind = 0;
    if (myfile.is_open())
    {
        while ( getline (myfile,line) )
        {
        	protein_names.push_back(line);
        }
        myfile.close();

    } else {
    	cerr << "Unable to open file " << fName << endl;
    	return false;
    }

    return true;
}

void createNames(	const string &protein_name,
					const string &path,
					double &measure,
					string &pos_potential_points,
					string &neg_potential_points,
					string &pos_potential_mu,
					string &neg_potential_mu,
					string &act_center
					){

	string pos = "positive";
	string neg = "negative";

	string folder  = path + "/" + protein_name + "/";

	pos_potential_points = folder + protein_name + "_pot_" + pos + ".pts";
	neg_potential_points = folder + protein_name + "_pot_" + neg + ".pts";

	stringstream s;
	s << std::fixed << std::setprecision(4);
	s << measure;

	pos_potential_mu = folder + protein_name + "_mu_" + s.str() + pos + ".pts";
	neg_potential_mu = folder + protein_name + "_mu_" + s.str() + neg + ".pts";

	act_center = folder + protein_name + "_active_center.pts";
}


void writeEmdDistancesToFile(	const string &outpath,
								const unsigned int n,
								const unsigned int rounds,
								const double distribution,
								bool normalize_measure,
								const double c1,
								const double c2,
								const double c3,
								const vector<string> &protein_names_target,
								const vector<string> &protein_names,
								const vector<vector<double> > &emd_distances,
								const string &id,
								const int NNToActCent = 0){
	/*
	 * id is used in case of multiple iterations of the program
	 */

	stringstream s;
	s << std::fixed << std::setprecision(4);

	s << outpath << "/EMD_" << n << "_" << rounds << "_" << distribution
	<< "_" << c1 << "_" << c2 << "_" << c3 << "_id_" << id << "_NNact_" << NNToActCent << ".csv";



//	string fName = outpath + "/EMD_" + std::to_string(n) + "_" + std::to_string(rounds) + "_" + std::to_string(distribution)
//	  	  	  + "_" + std::to_string(c1) + "_" + std::to_string(c2) + "_" + std::to_string(c3) + ".csv";

	const char* path = outpath.c_str();
	boost::filesystem::path dir(path);
	if(boost::filesystem::create_directory(dir))
	{
	    std::cerr<< "Directory Created: "<<outpath<<std::endl;
	}


	string fName = s.str();
	ofstream myfile (fName);
	if (myfile.is_open())
	{
		myfile << "protein_name_target" << ";" << "protein_name" << ";" << "emd_distance" << "\n";
		for(unsigned int i = 0; i < emd_distances.size(); i++){
			for(unsigned int j = 0; j < emd_distances[i].size(); j++){
				myfile << protein_names_target[i] << ";" << protein_names[j] << ";" << emd_distances[i][j] << "\n";
			}
		}
	}
	else cout << "Unable to open file " << fName;
}

vector<string> readFunctionalProteinNames(const string &functional_proteins_path){

	vector<string> functional_protein_names;

    string line;
    ifstream myfile (functional_proteins_path);

    unsigned int ind = 0;
    if (myfile.is_open())
    {
        while ( getline (myfile,line) )
        {
        	functional_protein_names.push_back(line);
        }
        myfile.close();

    } else {
    	cerr << "Unable to open file " << functional_proteins_path << endl;
    }


	return functional_protein_names;
}

void sort_emd_and_protNames(	vector<string> &proteinNames,
								vector<double> &emd_distances){
	vector<unsigned int> index(emd_distances.size(), 0);
	for (unsigned int i = 0 ; i < index.size() ; i++) {
	    index[i] = i;
	}

	sort(index.begin(), index.end(),
	    [&](const int& a, const int& b) {
	        return (emd_distances[a] < emd_distances[b]);
	    }
	);

	vector<string> proteinNames_copy = proteinNames;
	vector<double> emd_distances_copy = emd_distances;
	for (int i = 0 ; i < index.size() ; i++) {
		proteinNames[i] = proteinNames_copy[index[i]];
		emd_distances[i] = emd_distances_copy[index[i]];
	}
}

unsigned int gaussSum(unsigned int n){
	return(0.5*n*(n+1));
}

double pseudoRoc(	const vector<string> &proteinNames_ordered,
					const set<string> &functional_protein_names_copy){

	set<string> functional_protein_names = functional_protein_names_copy;

	double s = 0.0;

	// 0.5*n*(n+1)
	double maxS = 	gaussSum(functional_protein_names.size())
					+ functional_protein_names.size()*(proteinNames_ordered.size()-functional_protein_names.size());

//	cout << "maxS: " << maxS << endl;

	unsigned int numOfFoundProteins = 0;
	for(unsigned int i = 0; i < proteinNames_ordered.size(); i++){

		if(functional_protein_names.find(proteinNames_ordered[i]) != functional_protein_names.end()){
			numOfFoundProteins++;
			functional_protein_names.erase(proteinNames_ordered[i]);
		}

		s += numOfFoundProteins;
	}

	return(s/maxS);
}

double testParameters(	const double c1, const double c2, const double c3,
						const vector<string> &proteinNames_target,
						const vector<string> &proteinNames,
						const std::set<std::string> &functional_protein_names,
						const unsigned int num_bins,
						const string &outPath){

	vector<double> vals_X_v_X(rounds,-1);
	vector<double> vals_X_v_Y(rounds,-1);


	vector<vector<double>> emd_distances(1, vector<double>(outX_v_Y_training.size(),0));
	for(unsigned int i = 0; i < outX_v_Y_training.size(); i++){

		for(unsigned int j = 0; j < rounds; j++){
			vals_X_v_X[j] = c1*outX_v_X_training[j][0] + c2*outX_v_X_training[j][1] + c3*outX_v_X_training[j][2];
			vals_X_v_Y[j] = c1*outX_v_Y_training[i][j][0] + c2*outX_v_Y_training[i][j][1] + c3*outX_v_Y_training[i][j][2];
		}

		emd_distances[0][i] = emd(vals_X_v_X, vals_X_v_Y, num_bins);


//		emd2hist(hist_X_v_X, hist_X_v_Y, rounds);
	}

	writeEmdDistancesToFile(outPath, number_of_selected_points, rounds,measure, normalize_measure, c1,c2,c3,proteinNames_target,proteinNames,emd_distances,emd_list_id);

	vector<string> proteinNames_sorted = proteinNames;
	// sorted in ascending order
	sort_emd_and_protNames(proteinNames_sorted,emd_distances[0]);

	double psR = pseudoRoc(proteinNames_sorted,functional_protein_names);

	return(psR);
}

vector<double> seq(double start, double end, unsigned int num){

	vector<double> sequence(num,0);

	double step = (end - start)/num;

	for(unsigned int i = 0; i < num; i++){
		sequence[i] = step*i + start;
	}
	return sequence;
}

vector<double> seqAndInverse(vector<double> v){

	vector<double> ret(v.size()*2);
	for(unsigned int i = 0; i < v.size(); i++){
		ret[i] = v[i];
		ret[i+v.size()] = 1/v[i];
	}

//	for(unsigned int i = 0; i < v.size(); i++){
//
//	}
	return ret;
}


std::set<double> c1_v_c2_ratio;
std::set<double> c1_v_c3_ratio;
std::set<double> c2_v_c3_ratio;

std::set<double> ratios;

struct parameterTripel{
	double c1;
	double c2;
	double c3;

	unsigned int num_bins;

	double result;

	parameterTripel(){
		this->c1 = 0;
		this->c2 = 0;
		this->c3 = 0;

		num_bins = 10;
		result = 0;
	}

	parameterTripel(double c1, double c2, double c3, unsigned int num_bins){
		this->c1 = c1;
		this->c2 = c2;
		this->c3 = c3;
		this->num_bins = num_bins;

		result = 0;
	}

	bool operator < (parameterTripel &other){
		return(other.result > result);
	}

};

ostream& operator<<(ostream& os, const parameterTripel& dt)
{
	os << "c1 = " << dt.c1 << ", c2 = " << dt.c2 << ", c3 = " << dt.c3 << ", bins = " << dt.num_bins << ", res = " << dt.result;
	return(os);
}

vector<parameterTripel> createParameterTripels(	vector<double> &c1_r_c2_vals,
												vector<double> &c1_r_c3_vals,
												vector<unsigned int> num_bins_vals){
	/*
	 * c1_r_c2_vals ... different values specifying c_1/c2
	 * c1_r_c3_vals ... different values specifying c_1/c3
	 */

	vector<parameterTripel> params((3 + c1_r_c2_vals.size()*c1_r_c3_vals.size())*num_bins_vals.size());

	for(unsigned int k = 0; k < num_bins_vals.size(); k++){
		unsigned int k_ind = k*(3 + c1_r_c2_vals.size()*c1_r_c3_vals.size());

		params[0 + k_ind] = parameterTripel(1,0,0, num_bins_vals[k]);
		params[1 + k_ind] = parameterTripel(0,1,0, num_bins_vals[k]);
		params[2 + k_ind] = parameterTripel(0,0,1, num_bins_vals[k]);

		for(unsigned int i = 0; i < c1_r_c2_vals.size(); i++){
			for(unsigned int j = 0; j < c1_r_c3_vals.size(); j++){
				params[i*c1_r_c3_vals.size() + j + 3 + k_ind] = parameterTripel(1,c1_r_c2_vals[i],c1_r_c3_vals[j], num_bins_vals[k]);
			}
		}
	}


	return(params);
}


void optimize_parameters(	const string &functional_proteins_path,
							const vector<string> &proteinNames_target,
							const vector<string> &proteinNames,
							const string &outPath){

	cout << "-----------------------------------------\n";
	cout << "Optimizing parameters ..." << endl;

	//read in functional proteins
	vector<string> f = readFunctionalProteinNames(functional_proteins_path);
	std::set<std::string> functional_protein_names(f.begin(), f.end());


	// only use values in range of 0 to 1
//	vector<double> c1_vals = seq(0,1,50);
//	vector<double> c2_seq = seq(1.1,9.9,10);
//	vector<double> c3_seq = seq(1.1,9.9,10);
//
//	vector<double> c2_vals = seqAndInverse(c2_seq);
//	vector<double> c3_vals = seqAndInverse(c3_seq);

	vector<double> c2_vals{0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9,10};
	vector<double> c3_vals{0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,2,3,4,5,6,7,8,9,10};


	vector<unsigned int> n_bins{num_bins};


	vector<parameterTripel> parameterTriples = createParameterTripels(c2_vals, c3_vals, n_bins);


	for(unsigned int i = 0; i < parameterTriples.size(); i++){
		double psR = testParameters(parameterTriples[i].c1,parameterTriples[i].c2,parameterTriples[i].c3,proteinNames_target, proteinNames,functional_protein_names, parameterTriples[i].num_bins, outPath);

		parameterTriples[i].result = psR;
	}




//	vector<double>::iterator maxpsR_it = std::max_element(std::begin(psR_results), std::end(psR_results));
//	unsigned int maxpsR_ind = std::distance(psR_results.begin(), maxpsR_it);
//	double maxpsR = psR_results[maxpsR_ind];


	vector<parameterTripel>::iterator maxpsR_it = std::max_element(std::begin(parameterTriples), std::end(parameterTriples));
	double maxpsR = maxpsR_it->result;

	cout << "tested " << parameterTriples.size() << " parameter-combinations " << endl;
//	cout << maxpsR << endl;

//	cout << "parameters: c1 = " << maxpsR_it->c1 << " c2 = " << maxpsR_it->c2 << " c3 = " << maxpsR_it->c3 << endl;

	cout << *maxpsR_it << endl;


	string fName = outPath + "/allParametersWithRoc.txt";

	bool existed = true;
	if ( !boost::filesystem::exists( fName ) )
	{
		std::ofstream outfile (fName);

		outfile << "n m measure c1 c2 c3 bins ROC\n";

		outfile.close();
	}


	ofstream myfile (fName , std::ios_base::app);
	if (myfile.is_open())
	{

		for(unsigned int i = 0; i < parameterTriples.size(); i++){
			myfile << number_of_selected_points << " " << rounds << " " << measure << " " <<
					parameterTriples[i].c1 << " " << parameterTriples[i].c2 << " " << parameterTriples[i].c3
					<< " " << parameterTriples[i].num_bins << " " << parameterTriples[i].result << " " << "\n";

		}
//		myfile << number_of_selected_points << " " << rounds << " " << measure << " " <<
//				maxpsR_it->c1 << " " << maxpsR_it->c2 << " " << maxpsR_it->c3 << " " << maxpsR_it->num_bins << " " << maxpsR_it->result << " " << endl;

	}
	else cout << "Unable to open file " << fName;
}

int readHelpPage(string fileName){
	  string line;
	  ifstream myfile (fileName);
	  if (myfile.is_open())
	  {
	    while ( getline (myfile,line) )
	    {
	      cout << line << '\n';
	    }
	    myfile.close();

	    return(0);
	  }

	  else cerr << "Unable to open the help-page " << fileName << endl;

	  return 1;
}

int main(int argc, char *argv[])
{
    if(argc < 14) {
        cerr << "need 13 arguments. exiting ..." << endl;

        readHelpPage("help.txt");

        return 1;
    }

    string path = argv[1];							// path to the folder with the proteins and pts-files
    string outPath = argv[2];						// path to the folder where the output of this programm should be stored
//    string Dataname1 = argv[3];
    string proteinsToCompareFile_target = argv[3];	// proteins to compare to, usually only Trx or all proteins
    string proteinsToCompareFile = argv[4];			// all other proteins

    measure = stod(argv[5]);						// model parameter, measure == 1 equals a uniform-distribution

//    normalize_measure = (bool)atoi(argv[6]);

//    cout << "----------------------------" << endl;
//    cout << normalize_measure << endl;
//    cout << "----------------------------" << endl;


	number_of_selected_points = atoi(argv[6]);		// model parameters n
	rounds = atoi(argv[7]);							// model parameters m

	 c1 = stod(argv[8]);		// model parameters
	 c2 = stod(argv[9]);		// model parameters
	 c3 = stod(argv[10]);		// model parameters

//	 num_bins = log2(rounds)+1;	// model parameters
     num_bins = 100;

	 emd_list_id = argv[11];	// for running multiple times with same parameters


	bool allParameterCombinations = (bool)atoi(argv[12]); // for optimizing the modelparameters

	int NNtoActCent = atoi(argv[13]); // = 0 then all points are taken



	string Dataname1_pos;
    string Dataname1_mu_pos;

	string Dataname1_neg;
	string Dataname1_mu_neg;

    string Dataname1_actCent;


	string Dataname2_pos;
    string Dataname2_mu_pos;

	string Dataname2_neg;
	string Dataname2_mu_neg;

    string Dataname2_actCent;


	pthread_t threads[NUM_THREADS];
	struct thread_data td[NUM_THREADS];

	for(unsigned int i = 0; i < NUM_THREADS; i++){
		td[i].resize(number_of_selected_points);
	}

	int rc;

	vector< vector <double> > points1_pos_all(3);
	vector< vector <double> > points1_neg_all(3);

	vector< vector <double> > points2_pos_all(3);
	vector< vector <double> > points2_neg_all(3);

	vector<double> distribution1_pos_all;
	vector<double> distribution1_neg_all;

	vector<double> distribution2_pos_all;
	vector<double> distribution2_neg_all;

	out.resize(rounds, vector<double>(3,-1));


	vector<string> protein_names;
	readListOfProteinsToCompareWith(proteinsToCompareFile, protein_names);


	vector<string> protein_names_target;
	readListOfProteinsToCompareWith(proteinsToCompareFile_target, protein_names_target);

	for(unsigned int i = 0; i < protein_names_target.size(); i++){
		cout << protein_names_target[i] << endl;
	}


	vector< vector <double> > emd_distances(protein_names_target.size(),vector<double>(protein_names.size(),0));
	for(unsigned int j = 0; j < protein_names_target.size(); j++){

		cout << "target protein: " << protein_names_target[j] << "\n";

		createNames(protein_names_target[j], path, measure, Dataname1_pos, Dataname1_neg, Dataname1_mu_pos, Dataname1_mu_neg, Dataname1_actCent);

		readPointsAndMeasue(Dataname1_pos, points1_pos_all, Dataname1_neg, points1_neg_all,
							Dataname1_mu_pos, distribution1_pos_all, Dataname1_mu_neg, distribution1_neg_all, Dataname1_actCent, measure, NNtoActCent);

		// Compare Protein1 with itself
		vector<double> vals_X_v_X(rounds,-1);
		sampleDistances(vals_X_v_X, td, threads, rounds,
						points1_pos_all, points1_neg_all,
						distribution1_pos_all, distribution1_neg_all,

						points1_pos_all, points1_neg_all,
						distribution1_pos_all, distribution1_neg_all,
						c1, c2, c3);

		outX_v_X_training = out;
		outX_v_Y_training.resize(protein_names.size());

//		vector<double> emd_distances(protein_names.size(),0);
		for(unsigned int i = 0; i < protein_names.size(); i++){

			createNames(protein_names[i], path, measure, Dataname2_pos, Dataname2_neg, Dataname2_mu_pos, Dataname2_mu_neg, Dataname2_actCent);

			readPointsAndMeasue(Dataname2_pos, points2_pos_all, Dataname2_neg, points2_neg_all,
								Dataname2_mu_pos, distribution2_pos_all, Dataname2_mu_neg, distribution2_neg_all, Dataname2_actCent, measure, NNtoActCent);

			// Compare Protein1 with Protein2
			vector<double> vals_X_v_Y(rounds,-1);
			sampleDistances(vals_X_v_Y, td, threads, rounds,
							points1_pos_all, points1_neg_all,
							distribution1_pos_all, distribution1_neg_all,

							points2_pos_all, points2_neg_all,
							distribution2_pos_all, distribution2_neg_all,
							c1, c2, c3);

			outX_v_Y_training[i] = out;

			double emd_val = emd(vals_X_v_X, vals_X_v_Y, num_bins);

	//	    cout << vals_X_v_X.size() << endl;
	//	    cout << out.size() << endl;
	//
	//	    cout << distribution2_neg_all.size() << endl;
	//
	//	    cout << emd_val << endl;
			emd_distances[j][i] = emd_val;

		}
	}

	if(allParameterCombinations){
		optimize_parameters("/home/willy/RedoxChallenges/MasterThesis/ExtrinsicDistances/data/106Redoxins/functionalProteins.txt",
							protein_names_target,
							protein_names,
							outPath);
	} else {
		writeEmdDistancesToFile(outPath, number_of_selected_points, rounds,measure, normalize_measure, c1,c2,c3,protein_names_target,protein_names,emd_distances,emd_list_id, NNtoActCent);
	}


	return(0);
}

//int main_old(int argc, char *argv[])
//{
//    if(argc < 18) {
//        cerr << "need 17 arguments. exiting ..." << endl;
//        return 1;
//    }
//
//	string Dataname1_pos = argv[1];
//    string Dataname1_mu_pos = argv[2];
//
//	string Dataname1_neg = argv[3];
//	string Dataname1_mu_neg = argv[4];
//
//    string Dataname1_actCent = argv[5];
//
//	string Dataname2_pos = argv[6];
//    string Dataname2_mu_pos = argv[7];
//
//	string Dataname2_neg = argv[8];
//	string Dataname2_mu_neg = argv[9];
//
//    string Dataname2_actCent = argv[10];
//
//
//	number_of_selected_points = atoi(argv[11]);
////	number_of_selected_points_doubled = number_of_selected_points*2;
//	rounds = atoi(argv[12]);
//
//	 c1 = stod(argv[13]);
//	 c2 = stod(argv[14]);
//	 c3 = stod(argv[15]);
//
//    // factor of weight that goes to the closest point
//    int measure = atoi(argv[16]);
//
//    string fileName = argv[17];
//
///*
//	cout << Dataname1_pos << endl;
//	cout << Dataname1_neg << endl;
//
//	cout << Dataname2_pos << endl;
//	cout << Dataname2_neg << endl;
//
//	cout << number_of_selected_points << endl;
//	cout << rounds << endl;
//
//*/
//
//    //cout << "parsing parameters done" << endl;
//
//	pthread_t threads[NUM_THREADS];
//	struct thread_data td[NUM_THREADS];
//
//	for(unsigned int i = 0; i < NUM_THREADS; i++){
//		td[i].resize(number_of_selected_points);
//	}
//
//	int rc;
//
//	vector< vector <double> > points1_pos_all(3);
//	vector< vector <double> > points1_neg_all(3);
//
//	vector< vector <double> > points2_pos_all(3);
//	vector< vector <double> > points2_neg_all(3);
//
//
//	out.resize(rounds, vector<double>(3,-1));
//
//	reading(Dataname1_pos, points1_pos_all); // reads the first file containing the 3D points
//	if(points1_pos_all[0].size() == 0){exit(-1);}
//
//	reading(Dataname1_neg, points1_neg_all); // reads the first file containing the 3D points
//	if(points1_neg_all[0].size() == 0){exit(-1);}
//
//	reading(Dataname2_pos, points2_pos_all); // reads the first file containing the 3D points
//	if(points2_pos_all[0].size() == 0){exit(-1);}
//
//	reading(Dataname2_neg, points2_neg_all); // reads the first file containing the 3D points
//	if(points2_neg_all[0].size() == 0){exit(-1);}
//
//	vector<double> distribution1_pos_all(points1_pos_all[0].size(), 1.0/(points1_pos_all.size()+points1_neg_all.size()));
//	vector<double> distribution1_neg_all(points1_neg_all[0].size(), 1.0/(points1_pos_all.size()+points1_neg_all.size()));
//
//	vector<double> distribution2_pos_all(points2_pos_all[0].size(), 1.0/(points2_pos_all.size()+points2_neg_all.size()));
//	vector<double> distribution2_neg_all(points2_neg_all[0].size(), 1.0/(points2_pos_all.size()+points2_neg_all.size()));
//
//
//    //if(measure != 1){
//        getMeasureWithDistToAC(points1_pos_all, distribution1_pos_all, points1_neg_all, distribution1_neg_all,
//                                Dataname1_mu_pos, Dataname1_mu_neg, Dataname1_actCent, measure);
//
//        getMeasureWithDistToAC(points2_pos_all, distribution2_pos_all, points2_neg_all, distribution2_neg_all,
//                                Dataname2_mu_pos, Dataname2_mu_neg, Dataname2_actCent, measure);
//    //}
//
//
//
//
//	// compare protein1 with itself
//	for(int i = 0; i < NUM_THREADS; i++ )
//	{
//		td[i].thread_id = i;
//		td[i].points1_pos_all = points1_pos_all;
//		td[i].points2_pos_all = points2_pos_all;
//
//		td[i].points1_neg_all = points1_neg_all;
//		td[i].points2_neg_all = points2_neg_all;
//
//        td[i].distribution1_pos_all = distribution1_pos_all;
//        td[i].distribution1_neg_all = distribution1_neg_all;
//
//        td[i].distribution2_pos_all = distribution2_pos_all;
//        td[i].distribution2_neg_all = distribution2_neg_all;
//
//		rc = pthread_create(&threads[i], NULL, Thread, (void *)&td[i]);
//
//		if (rc) {
//			cerr << "Error:unable to create thread," << rc << endl;
//			exit(-1);
//		}
//	}
//
//	for (int i = 0; i < NUM_THREADS; i++) {
//		pthread_join(threads[i], NULL);
//	}
//
//
//    ofstream myfile;
//    myfile.open (fileName);
//
//	myfile << "pos;" << "neg;" << "pos_and_neg" << "\n";
//	for(int i = 0; i < out.size();i++)
//	{
//
//		for(int j = 0; j < out[i].size();j++)
//		{
//			myfile << out[i][j];
//
//            if(j < out[i].size()-1){
//                myfile << ";";
//            }
//		}
//        myfile << "\n";
//	}
//    myfile.close();
//
//	return(0);
//}

void selectPointsRandomly(   const vector< vector <double> > &points,
                             vector< vector <double> > &selected_points,
                             const vector <double> &distribution,
                             vector <double> &selected_distribution,
                             const unsigned int &number_of_selected_points,
							 vector<unsigned int> &random_indices_list,
							 discrete_distribution<unsigned int> &dist){
	// Select number_of_selected_points indices randomly
//	vector<unsigned int> random_indices_list(number_of_selected_points);

//	vector<unsigned int> indices_list(points[0].size());
//	for(unsigned int i = 0; i < points[0].size(); i++){
//		indices_list[i] = i;
//	}


//	random_points(points[0].size(),number_of_selected_points, random_indices_list);


	// changed on 20.2.19, now points can be sampled multiple times
//	random_points_with_measure(distribution,number_of_selected_points,random_indices_list);

	random_points_with_measure(dist,number_of_selected_points,random_indices_list); // <-----------


	// changed this back on 15.4.2019. I think for uniform distribution this does not make a difference
//	random_points(points[0].size(),number_of_selected_points, random_indices_list);


	for(int i = 0; i < number_of_selected_points; i++)
	{
		selected_points[0][i] = points[0][random_indices_list[i]];
		selected_points[1][i] = points[1][random_indices_list[i]];
		selected_points[2][i] = points[2][random_indices_list[i]];

		// changed on 20.2.19
//        selected_distribution[i] = distribution[random_indices_list[i]];
		selected_distribution[i] = 1.0/number_of_selected_points;

	}

}

void *Thread(void *threadarg) {

	struct thread_data *my_data;
	my_data = (struct thread_data *) threadarg;

//	vector< vector <double> > selected_points1_pos(3, vector<double>(number_of_selected_points));
//	vector< vector <double> > selected_points2_pos(3, vector<double>(number_of_selected_points));
//
//	vector< vector <double> > selected_points1_neg(3, vector<double>(number_of_selected_points));
//	vector< vector <double> > selected_points2_neg(3, vector<double>(number_of_selected_points));
//
//    vector<double> selected_distribution1_pos(number_of_selected_points);
//    vector<double> selected_distribution1_neg(number_of_selected_points);
//    vector<double> selected_distribution2_pos(number_of_selected_points);
//    vector<double> selected_distribution2_neg(number_of_selected_points);

	vector<unsigned int> random_indices_list;
	double SolutionOfBound;

	// /* initialize random seed: */
	srand(time(NULL) + my_data->thread_id);

	set_distribution_X_pos(my_data->distribution1_pos_all);
	set_distribution_X_neg(my_data->distribution1_neg_all);
	set_distribution_Y_pos(my_data->distribution2_pos_all);
	set_distribution_Y_neg(my_data->distribution2_neg_all);

	for(int j = 0; j < rounds/NUM_THREADS; j++)
	{		
//		selectPointsRandomly(my_data->points1_pos, selected_points1_pos,
//                             my_data->distribution1_pos, selected_distribution1_pos,
//                             number_of_selected_points);

		selectPointsRandomly(my_data->points1_pos_all, my_data->points1_pos,
                             my_data->distribution1_pos_all, my_data->distribution1_pos,
                             number_of_selected_points, my_data->random_indices_list, distribution_X_pos);

		selectPointsRandomly(my_data->points1_neg_all, my_data->points1_neg,
                             my_data->distribution1_neg_all, my_data->distribution1_neg,
                             number_of_selected_points, my_data->random_indices_list, distribution_X_neg);

		selectPointsRandomly(my_data->points2_pos_all, my_data->points2_pos,
                             my_data->distribution2_pos_all, my_data->distribution2_pos,
                             number_of_selected_points, my_data->random_indices_list, distribution_Y_pos);

		selectPointsRandomly(my_data->points2_neg_all, my_data->points2_neg,
                             my_data->distribution2_neg_all, my_data->distribution2_neg,
                             number_of_selected_points, my_data->random_indices_list, distribution_Y_neg);


		vector<double > solution;
		solution.resize(3);
		// Calculates the first lower bound
//		SolutionOfBound = flb_with_relation(	selected_points1_pos, selected_distribution1_pos,
//						                        selected_points1_neg, selected_distribution1_neg,
//						                        selected_points2_pos, selected_distribution2_pos,
//						                        selected_points2_neg, selected_distribution2_neg,
//						                        number_of_selected_points,
//						                        c1,
//						                        c2,
//						                        c3,
//						                        solution);

		if(c1 > 0 && c2 == 0 &&  c3 == 0){

			SolutionOfBound = flb_with_relation_only_pos(my_data, solution);

//			SolutionOfBound = SolutionOfBound*sum(my_data->distribution1_pos)*sum(my_data->distribution2_pos);
		}else {
			SolutionOfBound = flb_with_relation(my_data, solution);

//			SolutionOfBound = SolutionOfBound*()*();
		}

//
		


		//out[j + (my_data->thread_id)*(rounds)/NUM_THREADS] = SolutionOfBound;
		out[j + (my_data->thread_id)*(rounds)/NUM_THREADS] = solution;
	}
	pthread_exit(NULL);

	return(NULL);
}
