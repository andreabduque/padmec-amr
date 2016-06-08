#ifndef __EDGE_H_INCLUDED__
#define __EDGE_H_INCLUDED__

#include <iostream>
#include <Vector.h>
#include <stdio.h>

using namespace std;

class Edge: public Vector{

	Point P;
	Point Q;
	
	public:
		Vector* normal_to_edge;

		Edge (Point _P, Point _Q);		
		//signed distance calculation edge-point
		double power(Point other);
		double showMe(Point other);
		Point getP();
		Point getQ();
};

#endif


