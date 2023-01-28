#ifndef A_STAR_H
#define A_STAR_H

#include <iostream>
#include <map>
#include <vector>
#include <queue> // for priority queue
#include <tuple>
#include <fstream>
#include <ctime>

#include <unordered_map>
#include <unordered_set>

#include "../PersistenceIO.h"
#include "../Algorithms/DijkstraShortestPath.h"
#include "../Algorithms/AnnotatingEdges.h"
#include "../BitSet.h"
#include "../Globals.h"
#include "../External/Mem_usage.h"
#include "../PriorityQueue.h"


/********************************************************************
* Description:	Given backtracing information, this function computes its correspoinding cycle
* Parameters:
* - source:						the source of search
* - target:						the target of search
* - previous:					the input backtracing information
* - edgeMap:					a map, mapping two endpoints to an edge
* - resShortestCycle:		the result shortest representative cycle
********************************************************************/
void backtraceShortestCycle(int source, int target, const vector<pair<int, int>> & previous,
	const map<pair<int, int>, int> & edgeMap, MatrixListType & resShortestCycle)
{
	resShortestCycle.clear();

	int mappedEdge;
	int ptr_1, ptr_2;
	ptr_1 = target;

	int pathSize = previous.size();
	for (int i = 0; i < pathSize; ++i)
	{
		ptr_1 = previous[i].first;
		ptr_2 = previous[i].second;
		if (ptr_1 > ptr_2)
			mappedEdge = edgeMap.at(std::make_pair(ptr_2, ptr_1));
		else
			mappedEdge = edgeMap.at(std::make_pair(ptr_1, ptr_2));

		resShortestCycle.push_back(mappedEdge);
	}

	// add the pivot edge
	ptr_1 = target; 
	ptr_2 = source;
	if (ptr_1 > ptr_2)
		SWAP(ptr_1, ptr_2);

	mappedEdge = edgeMap.at(std::make_pair(ptr_1, ptr_2));
	resShortestCycle.push_back(mappedEdge);

	// finally, adjust the cycle
	mysort(resShortestCycle);
}


/********************************************************************
* Description:	Given a homology class, this function employs A* algorithm to compute
						its shortest representative cycle.
* Parameters:
* - inputCycle:				the input homology class
* - cell2v_list:				a converter which projects cells to their corresponding constituent vertices
* - edgeAnnotations:		the annotations for all sentinel edges
* - edgeMap:					a map, mapping two endpoints to an edge
* - vertexNum:				the number of vertices in the whole topological space
* - resShortestCycle:		the result shortest representative cycle
********************************************************************/
void AStar_Optimal_Cycle(const MatrixListType & inputCycle, const vector<MatrixListType> & cell2v_list,
	const map<pair<int, int>, BitSet> & edgeAnnotations, const std::map<std::pair<int, int>, int> & edgeMap,
	int vertexNum, MatrixListType & resShortestCycle)
{
	resShortestCycle.clear();

	assert(!edgeAnnotations.empty());
	int BettiNum = (edgeAnnotations.cbegin()->second).getBitSize();

	BitSet targetAnnotation(BettiNum); // the target annotation we should reach finally
	computeCycleAnnotation(inputCycle, cell2v_list, edgeAnnotations, targetAnnotation);

	int low = inputCycle.back();
	MatrixListType edge, pivot = cell2v_list[low];
	int source = pivot[0];
	int target = pivot[1];
	if (source > target)
		SWAP(source, target);

	pair<int, int> key(source, target);
	map<pair<int, int>, BitSet>::const_iterator it_annotation;
	it_annotation = edgeAnnotations.find(key);
	if (it_annotation != edgeAnnotations.end()) // the pivot edge is a sentinel edge
		targetAnnotation ^= it_annotation->second; // exclude the pivot edge

	// -- Construct undirected graph for A*
	adjacency_list_t graph(vertexNum);
	for (size_t i = 0; i < low; i++)
	{
		edge = cell2v_list[i];
		graph[edge[0]].push_back(neighbor(edge[1]));
		graph[edge[1]].push_back(neighbor(edge[0]));
	}

	// -- Construct covering graphs for computing heuristics
	vector<adjacency_list_t> coveringGraphs(BettiNum);
	for (size_t i = 0; i < BettiNum; i++)
	{
		constructCoveringGraph1D(edgeAnnotations, cell2v_list, low, i, vertexNum, coveringGraphs[i]);
	}

	// -- Perform A* algorithm
	cgNode currentNode, neighborNode;
	cgNode sourceNode(BettiNum);
	sourceNode.vertex = source;

	priorityQueue searchQ; // priority queue for A* algorithm
	searchQ.push(sourceNode);

	unordered_map<pair<int, BitSet>, int, KeyHasher> computedHeuristics; // given a vertex, we store its already computed shortest path to the target for future usage
	unordered_map<pair<int, BitSet>, int, KeyHasher>::const_iterator it_heuristc;
	unordered_set<pair<int, BitSet>, KeyHasher> hasVisited; // a set which records nodes that have been visited
	unordered_set<pair<int, BitSet>, KeyHasher>::const_iterator ret;

	vector<unordered_map<int, int>> heuristic_database; // database of dynamically precomputed heuristic
	for (int i = 0; i < BettiNum; ++i)
		heuristic_database.push_back(unordered_map<int, int>());

	int cntExpanedNode = 0; // count the number of expaned nodes
	bool isInOpenSet = false;
	double lenHeuristicPath = 0;
	BitSet uvTargetAnnotation(BettiNum);
	pair<int, BitSet> keyH(-1, BitSet(BettiNum)); // search key for heuristics in computedHeuristics

	std::size_t physMemUsed; // for monitoring the memory footprint
	ofstream memoryFile(Globals::memoryFileName_HeuristicAlg, ios::out | ios::app);
	clock_t currTime;

	while (true)
	{
		currentNode = searchQ.top();
		searchQ.pop();
		hasVisited.insert(std::make_pair(currentNode.vertex, currentNode.sumAnnotation));
		++cntExpanedNode;

		if (currentNode.vertex == target && currentNode.sumAnnotation == targetAnnotation) // we have reached the target
		{
			backtraceShortestCycle(source, target, currentNode.previous, edgeMap, resShortestCycle); // find the shortest cycle
			break;
		}

		const std::vector<neighbor> & neighbors = graph[currentNode.vertex];
		for (const auto & nb : neighbors)
		{
			neighborNode = cgNode(BettiNum);
			neighborNode.vertex = nb.target;


			// -- update sumAnnotation
			neighborNode.sumAnnotation = currentNode.sumAnnotation;

			if (currentNode.vertex > nb.target) // still, we need to sort the endpoints
			{
				key.first = nb.target;
				key.second = currentNode.vertex;
			}
			else
			{
				key.first = currentNode.vertex;
				key.second = nb.target;
			}

			it_annotation = edgeAnnotations.find(key);
			if (it_annotation != edgeAnnotations.end()) // it is a sentinel edge
				neighborNode.sumAnnotation ^= it_annotation->second;


			// -- update fScore and gScore
			keyH.first = neighborNode.vertex;
			keyH.second = neighborNode.sumAnnotation;

			it_heuristc = computedHeuristics.find(keyH);
			if (it_heuristc != computedHeuristics.end()) // we have already computed that
			{
				lenHeuristicPath = it_heuristc->second;
			}
			else
			{
				uvTargetAnnotation = neighborNode.sumAnnotation ^ targetAnnotation;
				lenHeuristicPath = computeHeuristic(coveringGraphs, uvTargetAnnotation, 
					neighborNode.vertex, target, BettiNum, vertexNum, heuristic_database);
				computedHeuristics[keyH] = lenHeuristicPath;
			}
			neighborNode.gScore = currentNode.gScore + nb.weight;
			neighborNode.fScore = neighborNode.gScore + lenHeuristicPath;


			// -- update previous vector
			neighborNode.previous = currentNode.previous; // inherits its parent's data
			neighborNode.previous.push_back(std::make_pair(nb.target, currentNode.vertex));


			// -- insert it into priority queue
			// -- first, check if it is in close set
			ret = hasVisited.find(make_pair(neighborNode.vertex, neighborNode.sumAnnotation));
			if (ret != hasVisited.end()) // we have expanded this node before
				continue;

			// -- next, check if it is in open set
			isInOpenSet = false;
			isInOpenSet = searchQ.updateNodeInQueue(neighborNode);
			if (isInOpenSet == false) // this neighbor node is not in open set, so add it
				searchQ.push(neighborNode);
		}
	}// end for

	 // -- monitor the memory footprint
	physMemUsed = getPeakRSS() >> 20;
	currTime = clock();
	memoryFile << "Betti Number: " << BettiNum << "\t " << "Memory Footprint: " << physMemUsed << " (MB)" << endl;

	cout << "Number of expanded nodes: " << cntExpanedNode << endl;
	memoryFile.close();
}


/********************************************************************
* Description:	Interface for running A star algorithm.
* Parameters:
* - phi:							the input filter function
* - vList							a list storing the coordinates of vertices
* - lowerCellList:			a list storing the maxValues of cells
* - upperCellList:			a list stroing the maxValues of cells, whose dimension is greater
									than that of lowerCellList by 1
* - boundaryMatrix			the input boundary matrix, which will updated as the computed results
* - edgeMap:					a map, mapping two endpoints to an edge
* - cell2v_list:				a converter which projects cells to their corresponding constituent vertices
* - vertexNum:				the number of vertices in the whole topological space
* - low_array:				an array storing the pivot information
********************************************************************/
template<int arrayDim, int vertexDim = arrayDim>
void reduceND_AStar(blitz::Array<double, arrayDim> * phi, const vector<blitz::TinyVector<int, vertexDim>> & vList,
	const vector<int> & lowerCellList, const vector<CellNrType> & upperCellList, vector<MatrixListType> & boundaryMatrix,
	const std::map<std::pair<int, int>, int> & edgeMap, const vector<MatrixListType> &cell2v_list, int vertexNum,
	vector<int> &low_array)
{
	cout << "--- Using Heuristic-based Algorithm ---" << endl;
	ofstream memoryFile(Globals::memoryFileName_HeuristicAlg, ios::out | ios::trunc);
	memoryFile.close();

	clock_t startClock, endClock;

	std::map<int, MatrixListType> resCycles;
	for (int test = 0; test < boundaryMatrix.size(); test++)
	{
		if (boundaryMatrix[test].empty())
			continue;

		double birthTime, deathTime;
		double pers = computePersistence<arrayDim, vertexDim>(phi, vList, lowerCellList, upperCellList, boundaryMatrix, test, birthTime, deathTime);
		if (pers <= Globals::reduction_threshold)
			continue;

		// We first compute the annotations of all edges
		cout << "---------------------------------------" << endl;
		cout << "Compute edge annotations ..." << endl;
		map<pair<int, int>, BitSet> edgeAnnotations;
		computeAnnotations(boundaryMatrix, edgeMap, low_array, cell2v_list, test, vertexNum, edgeAnnotations);

		startClock = clock();
		cout << "Apply A* algorithm ..." << endl;
		MatrixListType resCycle;
		AStar_Optimal_Cycle(boundaryMatrix[test], cell2v_list, edgeAnnotations, edgeMap, vertexNum, resCycle);
		cout << "Size before: " << boundaryMatrix[test].size() << endl;
		cout << "Size after: " << resCycle.size() << endl;

		endClock = clock();
		cout << "Time consumed (s): " << (float)(endClock - startClock) / CLOCKS_PER_SEC << endl;
		cout << "---------------------------------------" << endl;

		resCycles.insert({ test, resCycle });
	}
	for (const auto & cycle : resCycles)
	{
		boundaryMatrix[cycle.first] = cycle.second;
	}
}

#endif // !A_STAR_H

