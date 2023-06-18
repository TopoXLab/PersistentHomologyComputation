#include "InputFileInfo.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>

using namespace std;

#ifdef enable_opencv
using namespace cv;
#endif

InputFileInfo::InputFileInfo() {
	// common info
	binary          = true;			// if the input file is binary
	source_mat		= false;        // input from mat instead of file path
	verbose			= false;        // turn on/off text outputs
	input_path		= "";		    // input file path
	output_path		= "";           // output file path
	file_type		= 0;			// input file type (0: Image data; 1: Dense distance matrix; 2: sparse distance matrix)

	// info for cubical image data
	dimension = 0;						// the data dimension, exclusively for cubical image data and general simplicial complex

	// info for dense distance matrix
	numPoints = 0;						// number of points
	dimPoints = 0;						// dimension of each feature point
}

void InputFileInfo::source_from_file(const string &input_file_, const string &output_path_) {
	input_path = input_file_;
	source_mat = false;
	if (output_path_ == "") output_path = input_path;
	else {
		vector<string> split_res = split(input_path, '/');
		string file_name = split_res[split_res.size() - 1];
		output_path = output_path_ + "/" + file_name;
	}

	if (input_path.find(".txt") != string::npos) // plain text file
	{
		binary = false;
		fstream f(input_path.c_str());

		if (!f.is_open())
			cout << "Cannot find the file " << input_path << endl;

		f >> file_type;
		if (file_type == 0) // Image data
		{
			f >> dimension;
			assert(dimension <= 8);
		}
		else if (file_type == 1) // Dense distance matrix
		{
			f >> numPoints;
			f >> dimPoints;
		}
		else if (file_type == 2) // General simplicial complex
		{
			f >> dimension;
			f >> numPoints;
			f >> dimPoints;
		}

		f.close();
	}
	else // binary data file
	{
		binary = true;
		fstream f(input_path.c_str(), ios::in | ios::binary);

		if (!f.is_open())
			cout << "Cannot find the file " << input_path << endl;

		f.read(reinterpret_cast<char*>(&file_type), sizeof(int));
		if (file_type == 0) // Image data
		{
			f.read(reinterpret_cast<char*>(&dimension), sizeof(int));
			assert(dimension <= 8);

			cout << "The dimension of input data: " << dimension << endl;
		}
		else if (file_type == 1) // Dense distance matrix
		{
			f.read(reinterpret_cast<char*>(&numPoints), sizeof(int));
			f.read(reinterpret_cast<char*>(&dimPoints), sizeof(int));
		}
		else if (file_type == 2) // General simplicial complex
		{
			f.read(reinterpret_cast<char*>(&dimension), sizeof(int));
			f.read(reinterpret_cast<char*>(&numPoints), sizeof(int));
			f.read(reinterpret_cast<char*>(&dimPoints), sizeof(int));
		}

		f.close();
	}
}

#ifdef enable_opencv
void InputFileInfo::source_from_mat(const string& output_path_, const Mat& t) {
	// only image data is supported if pass from mat!!
	file_type = 0;
	source_mat = true;
	output_path = output_path_;
	t.convertTo(mat, CV_64F);
	dimension = mat.dims;
	assert(dimension <= 8);
}
#endif

vector<string> InputFileInfo::split(string strToSplit, char delimeter)
{
	stringstream ss(strToSplit);
	string item;
	vector<string> splittedStrings;
	while (getline(ss, item, delimeter))
	{
		splittedStrings.push_back(item);
	}
	return splittedStrings;
}
