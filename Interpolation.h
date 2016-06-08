#ifndef __INTERPOLATION_H_INCLUDED__
#define __INTERPOLATION_H_INCLUDED__

#include <iostream>
#include <cmath>
#include <vector>


using namespace std;

class Edge;

typedef struct Point {
	double x;
	double y;
} Point;

class Interpolation {
 	//parameters 
 	//dont know them yet
	double mass;

	public:
	  Interpolation (double m_mass);

	  //returns vector of points of intersection between two edges
	  vector<Point> edge_edge_Intersection(Edge *P, Edge *Q);
	  	  
};

#endif