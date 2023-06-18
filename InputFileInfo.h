#ifndef INPUT_FILE_INFO
#define INPUT_FILE_INFO

#include <cstring>
#include <vector>

#ifdef enable_opencv
#include <opencv2/opencv.hpp>
#endif

class InputFileInfo
{
public:
	bool binary;
	bool source_mat;
	bool verbose;
	std::string input_path;
	std::string output_path;

	int file_type;
	int dimension;

	int numPoints;
	int dimPoints;

#ifdef enable_opencv
	cv::Mat mat;
#endif

public:
	InputFileInfo();
	~InputFileInfo() {}

	void source_from_file(const std::string &input_file_, const std::string &output_path_ = "");

#ifdef enable_opencv
	void source_from_mat(const std::string& output_path_, const cv::Mat& t);
#endif

	std::vector<std::string> split(std::string strToSplit, char delimeter);
};


#endif
