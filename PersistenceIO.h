#ifndef PERSISTENCEIO_H
#define PERSISTENCEIO_H

#include <blitz/array.h>
#include <blitz/tinyvec-et.h>
#include <functional>
#include <vector>
#include <malloc.h>
#include <ctime>
#include "Debugging.h"

#define SWAP(a, b)  do { a ^= b; b ^= a; a ^= b; } while ( 0 )

typedef vector<int> MatrixListType;

/************************************************************************/
/* Save persistence values and reduction results and reduction results to binary file			 */
/************************************************************************/
void saveResults(const char *fileName,unsigned int *data, const vector<int> &header, int nrPairs);


template<int d>
struct TextPersistentPairsSaver
{
	template<typename ContT>
	void savePers(vector<ContT> &res, const char * output_fname)
	{
		OUTPUT_MSG( endl<< "---------- writing result ---------------" );

		cout << "writing text output to:  " << output_fname << endl;

		fstream output_filestr(output_fname, fstream::out | fstream::trunc);

		for (int i = 0; i < d; i++)
		{
			//sort(res[i].begin(), res[i].end());
			output_filestr << "[" << std::to_string(i) << "D-Cell" << ", " << std::to_string(i+1) << "D-Cell] Pairs = " << res[i].size() << endl;
			typename ContT::iterator veiter;
			for(veiter=res[i].begin(); veiter!=res[i].end(); veiter++)
				output_filestr << veiter->birth << "\t" << veiter->death << endl;
			output_filestr << endl;
		}	

		OUTPUT_MSG( endl<< "---------- writing result finished ---------------" );
	}
};


template<int d, int arrayDim = d, int vertexDim = d>
struct BinaryPersistentPairsSaver
{
	typedef blitz::TinyVector<int, vertexDim> Vertex;

	template<typename ContT>
	void saveOutput(vector<ContT> &res, vector< Vertex > & vList, const char * output_fname)
	{
		assert(res.size() == d);

		int count = 0;
		for (int i = 0; i < res.size(); i++) {
			count += res[i].size();
		}

		for (int i = 0; i < d; i++)
			cout << "[" << std::to_string(i) << "D-Cell" << ", " << std::to_string(i + 1) << "D-Cell] Pairs = " << res[i].size() << endl;

		unsigned int *pairArray = new unsigned int[count * 2 * vertexDim];

		int index = 0;

		vector<int> header(d);

		for (int i = 0; i < d; i++)
		{
			typename ContT::iterator veiter;
			for (veiter = res[i].begin(); veiter != res[i].end(); veiter++)
			{
				Vertex birthcoord = veiter->birthV;
				for (size_t ii = 0; ii < vertexDim; ii++)
					pairArray[index++] = birthcoord[vertexDim - 1 - ii] + 1;

				Vertex deathcoord = veiter->deathV;
				for (size_t ii = 0; ii < vertexDim; ii++)
					pairArray[index++] = deathcoord[vertexDim - 1 - ii] + 1;
			}
			header[i] = res[i].size();
		}

		saveResults(output_fname, pairArray, header, index);

		delete[] pairArray;
	}

	void saveOneDimReduction(const vector< MatrixListType > & final_red_list, const vector< Vertex > & vList, const char * output_fname, int whichOne)
	{
		switch (whichOne)
		{
		case 0: // Cubical saver
			saveReduction(final_red_list, vList, output_fname);
			break;
		case 1: // Full Rips saver
			saveReduction(final_red_list, output_fname);
			break;
		case 2: // General Simplicial Complex saver
			saveReduction(final_red_list, output_fname);
			break;
		default:
			break;
		}
	}

	void saveReduction(const vector< MatrixListType > & final_red_list, const vector< Vertex > & vList, const char * output_fname)
	{
		int count = 0;
		for (int j = 0; j < final_red_list.size(); j++)
		{
			assert(!final_red_list[j].empty());
			count += 1;
			count += final_red_list[j].size();
		}

		unsigned int *redArray = new unsigned int[count * d];
		int index = 0;
		vector<int> header(d);
		header[0] = final_red_list.size();

		for (vector< MatrixListType >::const_iterator red_list_iter = final_red_list.begin(); red_list_iter != final_red_list.end(); red_list_iter++)
		{
			// size of the reduction list of this dot
			assert(red_list_iter->size() != 0);
			redArray[index++] = red_list_iter->size();
			for (int j = 1; j < d; j++)
				redArray[index++] = 0;

			// for each vertex in the red list, write its coordinates
			for (MatrixListType::const_iterator cellid_iter = red_list_iter->begin(); cellid_iter != red_list_iter->end(); cellid_iter++)
			{
				Vertex coord = vList[*cellid_iter];

				for (size_t ii = 0; ii < d; ii++)
				{
					redArray[index++] = coord[d - 1 - ii] + 1;
				}
			}
		}

		saveResults(output_fname, redArray, header, index);

		delete[] redArray;
	}


	void saveReduction(const vector< MatrixListType > & final_red_list, const char * output_fname)
	{
		int count = 0;
		for (int j = 0; j < final_red_list.size(); j++)
		{
			assert(!final_red_list[j].empty());
			count += 1;
			count += final_red_list[j].size();
		}

		unsigned int *redArray = new unsigned int[count];
		int index = 0;
		vector<int> header(1);
		header[0] = final_red_list.size();

		for (vector< MatrixListType >::const_iterator red_list_iter = final_red_list.begin(); red_list_iter != final_red_list.end(); red_list_iter++)
		{
			assert(red_list_iter->size() != 0);
			redArray[index++] = red_list_iter->size(); // save the current reduction list size

			for (MatrixListType::const_iterator cellid_iter = red_list_iter->begin(); cellid_iter != red_list_iter->end(); cellid_iter++)
			{
				redArray[index++] = *cellid_iter + 1; // save the point index
			}
		}

		saveResults(output_fname, redArray, header, index);

		delete[] redArray;
	}

	// Added by Fan Wang
	void index2coord(const vector<MatrixListType>& final_list, vector<Vertex>& vList, vector<vector<vector<int>>>& res) {
		res.resize(final_list.size());
		for (int i = 0; i < final_list.size(); i++) {
			res[i].resize(final_list[i].size());
			for (int j = 0; j < final_list[i].size(); j++) {
				res[i][j].resize(d);
				Vertex coord = vList[final_list[i][j]];
				for (size_t k = 0; k < d; k++) {
#ifdef MATLAB_USE
					res[i][j][k] = coord[d - 1 - k] + 1;
#else
					res[i][j][k] = coord[d - 1 - k];
#endif
				}
			}
		}
	}

	// Added by Fan Wang
	template<typename ContT>
	void pers2vector(const vector<ContT>& pers, vector<vector<vector<int>>>& pers_V, vector<vector<vector<double>>>& pers_BD) {
		assert(pers.size() == d);
		vector<int> header(d);
		pers_V.resize(d);
		pers_BD.resize(d);

		for (size_t i = 0; i < d; i++) {
			header[i] = pers[i].size();
			pers_V[i].resize(pers[i].size());
			pers_BD[i].resize(pers[i].size());
			for (int j = 0; j < pers[i].size(); j++) {
				Vertex birthcoord = pers[i][j].birthV;
				Vertex deathcoord = pers[i][j].deathV;

				int index = 0;
				vector<double> tmp_BD = { pers[i][j].birth, pers[i][j].death };
				vector<int> tmp_V(2 * d);
				for (size_t ii = 0; ii < d; ii++) {
#ifdef MATLAB_USE
					tmp_V[index++] = birthcoord[d - 1 - ii] + 1;
#else
					tmp_V[index++] = birthcoord[d - 1 - ii];
#endif
				}
				for (size_t ii = 0; ii < d; ii++) {
#ifdef MATLAB_USE
					tmp_V[index++] = deathcoord[d - 1 - ii] + 1;
#else
					tmp_V[index++] = deathcoord[d - 1 - ii];
#endif
				}
				assert(index == 2 * d);
				pers_V[i][j] = tmp_V;
				pers_BD[i][j] = tmp_BD;
			}
		}
	}

	// Added by Fan Wang
	template<typename ContT>
	void pers2vector(const ContT& pers, vector<vector<int>>& pers_V, vector<vector<double>>& pers_BD) {
		pers_V.resize(pers.size());
		pers_BD.resize(pers.size());
		for (int i = 0; i < pers.size(); i++) {
			Vertex birthcoord = pers[i].birthV;
			Vertex deathcoord = pers[i].deathV;

			int index = 0;
			vector<double> tmp_BD = { pers[i].birth, pers[i].death };
			vector<int> tmp_V(2 * d);
			for (size_t j = 0; j < d; j++) {
#ifdef MATLAB_USE
				tmp_V[index++] = birthcoord[d - 1 - j] + 1;
#else
				tmp_V[index++] = birthcoord[d - 1 - j];
#endif
			}
			for (size_t j = 0; j < d; j++) {
#ifdef MATLAB_USE
				tmp_V[index++] = deathcoord[d - 1 - j] + 1;
#else
				tmp_V[index++] = deathcoord[d - 1 - j];
#endif
			}
			assert(index == 2 * d);
			pers_V[i]  = tmp_V;
			pers_BD[i] = tmp_BD;
		}
	}

	// Added by Fan Wang
	static void write_BNDorRED(const vector<vector<vector<vector<int>>>>& t, const char *output_name) {

		const int dim = t.size();

		// Determine the size of the data buffer
		int count = 0;
		vector<int> header(dim);
		for (size_t i = 0; i < dim; i++) {
			header[i] = t[i].size();
			count += t[i].size();
			for (int j = 0; j < t[i].size(); j++) {
				assert(t[i][j].size() != 0);
				count += t[i][j].size();
			}
		}
		unsigned int *data = new unsigned int[count * dim];
		int index = 0;
		for (size_t i = 0; i < dim; i++) {
			for (int j = 0; j < t[i].size(); j++) {
				data[index++] = t[i][j].size();
				for (int dum = 1; dum < dim; dum++) {
					data[index++] = 0;
				}

				for (int k = 0; k < t[i][j].size(); k++) {
					for (size_t l = 0; l < dim; l++) {
						data[index++] = t[i][j][k][l];
					}
				}
			}
		}
		assert(index == count * dim);
		saveResults(output_name, data, header, index);

		delete data;
	}

	// Added by Fan Wang
	static void write_pers_V(
		const vector<vector<vector<int>>>& pers_V,
		const char* output_fname_V
	) {
		// ===== write out dat file =====
		const int dim = pers_V.size();
		int count = 0;
		vector<int> header(dim);
		for (size_t i = 0; i < dim; i++) {
			count += pers_V[i].size();
			header[i] = pers_V[i].size();
		}

		int index = 0;
		unsigned int *data = new unsigned int[count * 2 * dim];
		for (size_t i = 0; i < dim; i++) {
			for (int j = 0; j < pers_V[i].size(); j++) {
				for (size_t ii = 0; ii < dim * 2; ii++) {
					data[index++] = pers_V[i][j][ii];
				}
			}
		}
		assert(index == count * 2 * dim);
		saveResults(output_fname_V, data, header, index);
		delete data;
	}

	// Added by Fan Wang
	static void write_pers_BD(
		const vector<vector<vector<double>>>& pers_BD,
		const char* output_fname_BD
	) {
		// ===== Write out txt file =====
		const int dim = pers_BD.size();
		string names[] = { "Vertex", "Edge", "Face", "Cube", "4D-Cell", "5D-Cell" };
		fstream output_filestr(output_fname_BD, fstream::out | fstream::trunc);
		for (size_t i = 0; i < dim; i++) {
			output_filestr << names[i] << " " << names[i + 1] << " Pairs, Number = " << pers_BD[i].size() << endl;
			for (int j = 0; j < pers_BD[i].size(); j++) {
				output_filestr << pers_BD[i][j][0] << "\t" << pers_BD[i][j][1] << endl;
			}
		}
		output_filestr.close();
	}

	// Added by Fan Wang
	static void write_pers_BD(
		const vector<vector<double>>& pers_BD,
		const char* output_fname_BD
	) {
		// ===== Write out txt file =====
		const int dim = pers_BD.size();
		string names[] = { "Vertex", "Edge", "Face", "Cube", "4D-Cell", "5D-Cell" };
		fstream output_filestr(output_fname_BD, fstream::out | fstream::trunc);
		for (size_t i = 0; i < dim; i++) {
			output_filestr << names[i] << " " << names[i + 1] << " Pairs, Number = " << pers_BD[i].size() / 2 << endl;
			for (int j = 0; j < pers_BD[i].size() / 2; j++) {
				output_filestr << pers_BD[i][j * 2] << "\t" << pers_BD[i][j * 2 + 1] << endl;
			}
		}
		output_filestr.close();
	}
};


// Save persistence values and reduction results to binary file
void saveResults(const char *fileName, unsigned int *data, const vector<int> &header, int count)
{
	FILE *outFile;

	outFile = fopen(fileName, "wb");
	unsigned int dim = header.size();

	fwrite((void*)&dim, sizeof(unsigned int), 1, outFile);

	fwrite((void*)&header[0], sizeof(unsigned int), dim, outFile);

	int writeCount = fwrite((void *)data, sizeof(unsigned int), count, outFile);

	//printf("Writing %d Data Reduction Results in File %s\n", writeCount, fileName);

	fclose(outFile);
}
#endif
