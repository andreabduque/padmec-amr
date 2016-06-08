#include <Edge.h>
#include <Interpolation.h>
#include <Vector.h>
#include <iostream>
#include <stdio.h>
#include <vector>
#include <algorithm>

using namespace std;

int main(){
	
	Point P1, P2, Q1, Q2;
	double proj1, proj2, proj3, proj4;
	int it = 0;

	//inicializacao dos pontos correspondentes a extremidade das arestas
	P1.x = 0;
	P1.y = 0;
	P2.x = 2;
	P2.y = 2;
	Q1.x = 1;
	Q1.y = 1;
	Q2.x = 4;
	Q2.y = 4;

	/*
	
	vector<int> arg1;
	arg1.push_back(0);
	arg1.push_back(2);
	vector<int> arg2;
	arg2.push_back(1);
	arg2.push_back(3);

	int* intersection = new int[2];
	intersection[0] = max(min(arg1[0], arg1[1]), min(arg2[0], arg2[1]));
	intersection[1] = min(max(arg1[0], arg1[1]), max(arg2[0], arg2[1]));

	if (intersection[1] < intersection[0]) {
	  printf("empty");
	}
	else{
		printf("intersection interval: %d %d", intersection[0], intersection[1]);
	}
	

	*/
	Edge* edge_P = new Edge(P1, P2);
	Edge* edge_Q = new Edge(Q1, Q2);
	Interpolation* god = new Interpolation(123);
	vector<Point> intersects = god->edge_edge_Intersection(edge_P, edge_Q);

	proj1 = edge_P->showMe(Q1);
	proj2 = edge_P->showMe(Q2);
	proj3 = edge_Q->showMe(P1);
	proj4 = edge_Q->showMe(P2);

	for(it = 0; it < intersects.size(); it++){
		printf("\n Intersection Point: (%lf,%lf)\n", intersects[it].x, intersects[it].y);
	}
	if(!it)
		printf("There are no intersection points\n");

	return 0;
}