#ifndef INCLUDED_PERSISTENCE_CALC_RUNNER_H
#define INCLUDED_PERSISTENCE_CALC_RUNNER_H

#include "PersistentPair.h"
#include "PersistenceCalculator.h"
#include "Filtration/CubicalFiltration.h"
#include "Filtration/FullRipsFiltration.h"
#include "Filtration/SimComplexFiltration.h"

template<int dim>
struct PersistenceCalcRunnerCubical
{
	typedef blitz::TinyVector<int, dim> Vertex; 
	typedef vector<PersPair<Vertex>> PersResultContainer;	

	void go(
		blitz::Array<double, dim> *phi,
		double pers_thd,
		const InputFileInfo &info,
		vector<vector<vector<vector<int>>>>& final_red_list_grand,
		vector<vector<vector<vector<int>>>>& final_boundary_list_grand,
		vector<vector<vector<int>>>&         pers_V,
		vector<vector<vector<double>>>&      pers_BD
		)
	{	
		PersistenceCalculator<dim> calc;
		vector<PersResultContainer> res(dim);		

		vector<Vertex> vList;

		calc.calcPersistence(phi, pers_thd, res, vList, info, final_red_list_grand, final_boundary_list_grand, pers_V, pers_BD);

		//// local scope
		//{
		//	TextPersistentPairsSaver<dim> textSaver;

		//	stringstream output_file;
		//	output_file << info.input_path << ".pers.txt";
		//	
		//	textSaver.savePers(res, output_file.str().c_str());
		//}

		//{
		//	BinaryPersistentPairsSaver<dim> binSaver;

		//	stringstream output_file;
		//	output_file << info.input_path;
		//	output_file << ".pers";

		//	binSaver.saveOutput(res, vList, output_file.str().c_str());
		//}
	}
};



template<int dim>
struct PersistenceCalcRunnerFullRips
{
	typedef blitz::TinyVector<int, 2> Vertex;
	typedef vector<PersPair<Vertex>> PersResultContainer;

	void go(
		blitz::Array<double, 2>*				distMatrix,
		double									pers_thd,
		const InputFileInfo&					info,
		vector<vector<vector<vector<int>>>>&	final_red_list_grand,
		vector<vector<vector<vector<int>>>>&	final_boundary_list_grand,
		vector<vector<vector<int>>>&			pers_V,
		vector<vector<vector<double>>>&			pers_BD
	)
	{
		PersistenceCalculator<dim, 2, 2, FullRipsFiltration<dim>, 1> calc;
		vector<PersResultContainer> res(dim);

		vector<Vertex> vList;

		calc.calcPersistence(distMatrix, pers_thd, res, vList, info, final_red_list_grand, final_boundary_list_grand, pers_V, pers_BD);

		//// local scope
		//{
		//	TextPersistentPairsSaver<dim> textSaver;

		//	stringstream output_file;
		//	output_file << info.input_path << ".pers.txt";

		//	textSaver.savePers(res, output_file.str().c_str());
		//}

		//{
		//	BinaryPersistentPairsSaver<dim, 2, 2> binSaver;

		//	stringstream output_file;
		//	output_file << info.input_path;
		//	output_file << ".pers";

		//	binSaver.saveOutput(res, vList, output_file.str().c_str());
		//}
	}
};


template<int dim>
struct PersistenceCalcRunnerSimComplex
{
	typedef blitz::TinyVector<int, 1> Vertex;
	typedef vector<PersPair<Vertex>> PersResultContainer;

	void go(
		blitz::Array<double, 1>*				pointsVal,
		double									pers_thd,
		const InputFileInfo&					info,
		vector<vector<vector<vector<int>>>>&	final_red_list_grand,
		vector<vector<vector<vector<int>>>>&	final_boundary_list_grand,
		vector<vector<vector<int>>>&			pers_V,
		vector<vector<vector<double>>>&			pers_BD
	)
	{
		PersistenceCalculator<dim, 1, 1, SimComplexFiltration<dim>, 2> calc;
		vector<PersResultContainer> res(dim);

		vector<Vertex> vList;

		calc.calcPersistence(pointsVal, pers_thd, res, vList, info, final_red_list_grand, final_boundary_list_grand, pers_V, pers_BD);

		//// local scope
		//{
		//	TextPersistentPairsSaver<dim> textSaver;

		//	stringstream output_file;
		//	output_file << info.input_path << ".pers.txt";

		//	textSaver.savePers(res, output_file.str().c_str());
		//}

		//{
		//	BinaryPersistentPairsSaver<dim, 1, 1> binSaver;

		//	stringstream output_file;
		//	output_file << info.input_path;
		//	output_file << ".pers";

		//	binSaver.saveOutput(res, vList, output_file.str().c_str());
		//}
	}
};

#endif
