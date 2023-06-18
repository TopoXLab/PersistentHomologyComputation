/**************************************************************************
* Copyright 2017 Rutgers University and CUNY Queens College.
* Code by Pengxiang Wu based on the implementations of Hubert Wagner and Chao Chen.
*
* License: GPLv3. (https://www.gnu.org/licenses/gpl.html)
*
* The SOFTWARE PACKAGE provided in this page is provided "as is", without any guarantee
* made as to its suitability or fitness for any particular use. It may contain bugs, so use of
* this tool is at your own risk. We take no responsibility for any damage of any sort that may
* unintentionally be caused through its use.
*
* If you have any questions, please contact Pengxiang Wu (pxiangwu@gmail.com)
* The class is wrapped by Fan Wang
**************************************************************************/

#include <cmath>
#include <map>
#include <set>
#include <list>
#include <cassert>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <deque>
#include <ctime>
#include <blitz/array.h>
#include <blitz/tinyvec-et.h>

#include "PersistenceIO.h"
#include "Debugging.h"
#include "Filtration/CubicalFiltration.h"

#include "InputRunner.h"
#include "PersistentPair.h"

#include "PersistenceCalculator.h"
#include "PersistenceCalcRunner.h"
#include "Globals.h"
#include "ParseCommandLine.h"
#include "PersistenceComputer.h"

using namespace std;

#ifdef enable_opencv
using namespace cv;
#endif

double Persistence_Computer::run() {
	time_t startTime, endTime;
	if (debug_enabled) debugStart(debug_path);
	time(&startTime);
	runPersistenceHomology(file_info, final_red_list_grand, final_boundary_list_grand, pers_V, pers_BD);
	time(&endTime);
	double ellapsed1 = difftime(endTime, startTime);
	if (debug_enabled) debugEnd();
	return ellapsed1;
}

void Persistence_Computer::source_from_file(const string& input_file, const string& output_file) {
	Globals::inputFileName = input_file;
	file_info.source_from_file(input_file, output_file);
	output_name = file_info.output_path;
}

#ifdef enable_opencv
void Persistence_Computer::source_from_mat(const string& output_file, const Mat& t) {
	file_info.source_from_mat(output_file, t);
	output_name = file_info.output_path;
}

void Persistence_Computer::source_from_mat_from_int(const string& output_file, vector<int>& t, const int height, const int width) {
	const int size = t.size();
	assert(size == height * width);
	Mat res(height, width, CV_32SC1, &t.begin()[0]);
	file_info.source_from_mat(output_file, res);
	output_name = file_info.output_path;
}

void Persistence_Computer::source_from_mat_from_double(const string& output_file, vector<double>& t, const int height, const int width) {
	const int size = t.size();
	assert(size == height * width);
	Mat res(height, width, CV_64FC1, &t.begin()[0]);
	file_info.source_from_mat(output_file, res);
	output_name = file_info.output_path;
}
#endif

void Persistence_Computer::set_pers_thd(double t) {
	Globals::reduction_threshold = t;
	Globals::use_optimal_alg = true;
}
void Persistence_Computer::set_algorithm(int t) {
	if (t >= 2 || t < 0) {
		cout << "0 for A* search; 1 for exhaustive search ..."; exit(1);
	}
	Globals::which_alg = t; Globals::use_optimal_alg = true;
}

void Persistence_Computer::set_output_file(const string& t) { file_info.output_path = t; }
void Persistence_Computer::set_max_dim(int t) { Globals::max_dim = t; }
void Persistence_Computer::set_num_threads(int t) { Globals::num_threads = t; }
void Persistence_Computer::set_verbose(bool t) { file_info.verbose = t; }
void Persistence_Computer::set_debug(bool t, const string& debug_path_) { debug_enabled = t; debug_path = debug_path_; }

void Persistence_Computer::write_bnd(const vector<vector<vector<vector<int>>>>& t) {
	string output_bnd_file = output_name + ".bnd";
	BinaryPersistentPairsSaver<1>::write_BNDorRED(t, output_bnd_file.c_str());
}

void Persistence_Computer::write_red(const vector<vector<vector<vector<int>>>>& t) {
	string output_red_file = output_name + ".red";
	BinaryPersistentPairsSaver<1>::write_BNDorRED(t, output_red_file.c_str());
}

void Persistence_Computer::write_pers_V(const vector<vector<vector<int>>>& pers_V) {
	string output_persV_file = output_name + ".pers";
	BinaryPersistentPairsSaver<1>::write_pers_V(pers_V, output_persV_file.c_str());
}

void Persistence_Computer::write_pers_BD(const vector<vector<vector<double>>>& pers_BD) {
	string output_persBD_file = output_name + ".pers.txt";
	BinaryPersistentPairsSaver<1>::write_pers_BD(pers_BD, output_persBD_file.c_str());
}

void Persistence_Computer::write_pers_BD(const vector<vector<double>>& pers_BD) {
	string output_persBD_file = output_name + ".pers.txt";
	BinaryPersistentPairsSaver<1>::write_pers_BD(pers_BD, output_persBD_file.c_str());
}

void Persistence_Computer::return_bnd(vector<vector<vector<vector<int>>>>& t) { t = final_boundary_list_grand; };
void Persistence_Computer::return_red(vector<vector<vector<vector<int>>>>& t) { t = final_red_list_grand; };
void Persistence_Computer::return_pers_V(vector<vector<vector<int>>>& t) { t = pers_V; };
void Persistence_Computer::return_pers_BD(vector<vector<vector<double>>>& t) { t = pers_BD; };

void Persistence_Computer::write_output() {
	write_bnd(final_boundary_list_grand);
	write_red(final_red_list_grand);
	write_pers_V(pers_V);
	write_pers_BD(pers_BD);
}
void Persistence_Computer::clear() {
	final_red_list_grand.clear();
	final_boundary_list_grand.clear();
	pers_V.clear();
	pers_BD.clear();
}

void Persistence_Computer::debugStart(const string& debug_path) {
	string logFile = debug_path + "/" + "log.txt";
	string errorFile = debug_path + "/" + "error.txt";
	DebuggerClass::init(false, logFile, errorFile);

	int max_pers_pts = Globals::BIG_INT; // the maximal number of persistence pairs allowed

	fstream filestr;
	filestr.open(logFile.c_str(), fstream::out | fstream::trunc);

	filestr << "################################################" << endl;
	filestr << "Start computing persistence" << endl;
	filestr << "Input file = \'" << Globals::inputFileName << "\'" << endl;
	filestr << "Maximal number of persistence points allowed = " << max_pers_pts << endl;
	filestr.close();
}

void Persistence_Computer::debugEnd() {
	DebuggerClass::finish();
}

#ifdef command_line
int main(int argc, const char* argv[])
{
	time_t startTime, endTime;

	parseCommandLine(argc, argv);

	vector<vector<vector<vector<int>>>> final_red_list_grand;
	vector<vector<vector<vector<int>>>> final_boundary_list_grand;
	vector<vector<vector<int>>> pers_V;
	vector<vector<vector<double>>> pers_BD;

	InputFileInfo input_file_info;
	input_file_info.source_from_file(Globals::inputFileName);

	Persistence_Computer::debugStart("."); // debug settings
	time(&startTime);

	// Run persistence homology algorithm
	runPersistenceHomology(input_file_info, final_red_list_grand, final_boundary_list_grand, pers_V, pers_BD);

	time(&endTime);

	double ellapsed1 = difftime(endTime, startTime);
	cout << "All calculations computed in " << setprecision(3) << ellapsed1 / 60.0 << " Min" << endl;

	Persistence_Computer::debugEnd();

	return 0;
}
#endif

void main() {
	//std::string file_name = "D:/Data/cremi_new/c/gen_00046.png";
	std::string file_name = "E:/Data2/AI_0007_vol2_sup.dat";
	//cv::Mat img = cv::imread(file_name, 0);

	vector<vector<vector<vector<int>>>> final_red_list_grand;
	vector<vector<vector<vector<int>>>> final_boundary_list_grand;
	vector<vector<vector<int>>> pers_V;
	vector<vector<vector<double>>> pers_BD;

	Persistence_Computer pc;
	//pc.source_from_mat(file_name, img);
	pc.source_from_file(file_name);
	pc.run();
	pc.write_output();

	pc.return_bnd(final_boundary_list_grand);
	pc.return_pers_V(pers_V);


	pc.clear();
	//system("pause");
}