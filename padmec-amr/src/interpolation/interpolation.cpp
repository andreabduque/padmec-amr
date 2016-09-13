/*
 * interpolation.cpp
 *
 *  Created on: Dec 6, 2013
 *      Author: rogerio
 */

 /*Added Conservative Interpolation by Andrea*/

#include "interpolation.h"

void interpolation(InterpolationDataStruct* pIData, INTERPOLATION_OPTIONS opt){
	initialize(pIData);
	int dim = pIData->m1->getDim();
	switch ( opt ){
	case h_REFINEMENT:
		throw Exception(__LINE__,__FILE__,"Underconstruction!");
		break;
	case LINEAR:
		#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: Interpolation (Linear Method)\tIN\n";
		#endif
		calculate_GeometricCoefficients(pIData,dim);
		calculate_LinearInterpolation(pIData,dim);
		#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: Interpolation (Linear Method)\tOUT\n";
		#endif
		break;
	case QUADRATIC:
		#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: Interpolation (Quadratic Method)\tIN\n";
		#endif
		calculate_GeometricCoefficients(pIData,dim);
		calculate_LinearInterpolation(pIData,dim);
		for (int field=0; field<pIData->numFields; field++){
			calculate_Gradients(pIData,field);
			calculate_DerivativesError(pIData);
			calculate_QuadraticInterpolation(pIData,field);
		}
		#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: Interpolation (Quadratic Method)\tIN\n";
		#endif
		break;

	case CONSERVATIVE:
		#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: Interpolation (Conservative Method)\tIN\n";
		#endif
		//calculate_GeometricCoefficients(pIData,dim);
		//calculate_LinearInterpolation(pIData,dim);

		calculate_ConservativeInterpolation(pIData, dim);
	
		#ifdef TRACKING_PROGRAM_STEPS
		cout << "TRACKING_PROGRAM_STEPS: Interpolation (Conservative Method)\tIN\n";
		#endif
		cout << "parei motherfucker\n";
		STOP();

		break;
	default:
		throw Exception(__LINE__,__FILE__,"Interpolation method unknown. Exiting....");
	}
	finalize(pIData);
}


void initialize(InterpolationDataStruct* pIData){
	int dim = pIData->m1->getDim();
	pIData->theOctree = OctreeCreate2<iterall>(pIData->m2->beginall(dim),pIData->m2->endall(dim),dim);
	int nrows1 = M_numVertices(pIData->m1) + 1;	// New mesh (to)
	int nrows2 = M_numVertices(pIData->m2) + 1; // Old mesh (from)
	pIData->pNodeValue.allocateMemory(nrows2);
	pIData->pGrad.allocateMemory(nrows2,5);
	pIData->pGeomCoeff.allocateMemory(nrows1,dim+1);
	pIData->pInterpolatedVals.allocateMemory(nrows1);
}

void finalize(InterpolationDataStruct* pIData){
	pIData->pNodeValue.freeMemory();
	pIData->pInterpolatedVals.freeMemory();
	pIData->pGeomCoeff.freeMemory();
	pIData->pGrad.freeMemory();
	pIData->theOctree = 0;
}

double power(pPoint other, pPList  vertices_edge){

	pPoint P, Q;
	Vector* v;
	Vector* normal;
	//bool debug = true;

	pEntity vertex_P = (pEntity)PList_item(vertices_edge, 0);
	pEntity vertex_Q = (pEntity) PList_item(vertices_edge, 1);

	P =  V_point(vertex_P);
	Q =  V_point(vertex_Q);

	v = new Vector(P, Q);
	normal = v->normal_unit_pos_vector();

	double proj = (P_x(P) - P_x(other))*normal->getX() + (P_y(P) - P_y(other))*normal->getY();

	delete v;
	delete normal;

	return proj;

}
//Im very sorry  
bool its_me(double x1, double y1, double x2, double y2){
	return (abs(x1 - x2) < EPSILON) && (abs(y1 - y2) < EPSILON);
}
int who_is(double x, double y, vector<pPoint> vec){
	int who = 0;

	for(unsigned int i = 0; i < vec.size(); i++){
		if(its_me(x, y, P_x(vec[i]), P_y(vec[i]))){
			who = i + 1;
			break;
		}
	}
	return who;
}

bool tolerance(double value){
	return abs(value - 0) < pow(10,-6);

}

pPoint min_point(pPoint P, pPoint Q, int eixo){

	pPoint which;
	
	switch(eixo){
		case 0: if(P_x(P) < P_x(Q))
					which = P;
				else
					which = Q;

				break;

		case 1:
				if(P_y(P) < P_y(Q))
					which = P;					
				else
					which = Q;

				break;
	}
	return which;
}

pPoint max_point(pPoint P, pPoint Q, int eixo){
	pPoint which;
	
	switch(eixo){
		case 0: if(P_x(P) > P_x(Q))
					which = P;
				else
					which = Q;

				break;

		case 1:
				if(P_y(P) > P_y(Q))
					which = P;					
				else
					which = Q;

				break;
	}
	return which;
}

vector<pPoint> edge_intersection(pEntity edge1, pEntity edge2){

	pPList  vertices_edge2, vertices_edge1;
	double Power_P1, Power_P2, Power_Q1, Power_Q2, x, y;
	int numPoints = 0;
	pPoint A, B;
	vector<pPoint> intersect_points;
	bool _P1, _P2, _Q1, _Q2 = false;
	vector<int> id_points;

	//list of vertices for each edge
	vertices_edge1 = E_vertices(edge1);
	vertices_edge2 = E_vertices(edge2);

	pPoint P1 = V_point((pEntity)PList_item(vertices_edge1, 0));
	pPoint P2 = V_point((pEntity)PList_item(vertices_edge1, 1));
	pPoint Q1 = V_point((pEntity)PList_item(vertices_edge2, 0));
	pPoint Q2 = V_point((pEntity)PList_item(vertices_edge2, 1));

	Power_P1 = power(Q1, vertices_edge1);
	Power_P2 = power(Q2, vertices_edge1);
	Power_Q1 = power(P1, vertices_edge2);
	Power_Q2 = power(P2, vertices_edge2);

	//number of zero powers 
	int zero_powers = 0;

	//vector to evaluate intersection intervals
	pPoint* aux_x = new pPoint[2];
	pPoint* aux_y = new pPoint[2];

	if(_P1 = tolerance(Power_P1))
		zero_powers++;
	if(_P1 = tolerance(Power_P1))
		zero_powers++;
	if(_Q1 = tolerance(Power_Q1))
		zero_powers++;
	if(_Q2 = tolerance(Power_Q2))
		zero_powers++;

		

	switch (zero_powers){
		case 0:
				if(Power_P1*Power_P2 < 0 && Power_Q1*Power_Q2 < 0){
					x = P_x(P1) +(Power_Q1/(Power_Q1-Power_Q2))*(P_x(P2) - P_x(P1));
					y = P_y(P1) +(Power_Q1/(Power_Q1-Power_Q2))*(P_y(P2) - P_y(P1));


					if(abs(x - P_x(P1)) < pow(10,-6) && abs(y - P_y(P1)) < pow(10,-6))
						A = P1;
					else if (abs(x - P_x(P2)) < pow(10,-6) && abs(y - P_y(P2)) < pow(10,-6))
						A = P2;
					else if (abs(x - P_x(Q1)) < pow(10,-6) && abs(y - P_y(Q1)) < pow(10,-6))
						A = Q1;
					else if (abs(x - P_x(Q2)) < pow(10,-6) && abs(y - P_y(Q2)) < pow(10,-6))
						A = Q2;
					else{

						A = P_new();
						P_setPos(A, x, y, 0);
					}
					
						numPoints = 1;
				}

				break;
		case 1:	
				numPoints = 1;

				if (_P1 && Power_Q1*Power_Q2 < 0)
					A = Q1;		
				else if (_P2 && Power_Q1*Power_Q2 < 0)
					A = Q2;
				else if (_Q1 && Power_P1*Power_P2 < 0)
					A = P1;
				else if (_Q2 &&  Power_P1*Power_P2 < 0) 
					A = P2;
				else 
					numPoints = 0;

				break;
		case 2:	
				numPoints = 1;

				if((_P1 & _Q1) || (_P1 & _Q2))
					A = Q1;			
				else if ((_P2 & _Q1) || (_P2 & _Q2))
					A = Q2;
				else
					numPoints = 0;				
				break;
		case 4:		
				//intersection intervals in x				
				aux_x[0] = max_point(min_point(P1, P2, 0), min_point(Q1, Q2, 0), 0);
				aux_x[1] = min_point(max_point(P1, P2, 0), max_point(Q1, Q2, 0), 0);

				
				if (P_x(aux_x[1]) < P_x(aux_x[0])) {
				  numPoints = 0;
				  break;
				}
				else{
					//intersection intervals in y
					aux_y[0] = max_point(min_point(P1, P2, 0), min_point(Q1, Q2, 0), 0);
					aux_y[1] = min_point(max_point(P1, P2, 0), max_point(Q1, Q2, 0), 0);
					
					if (P_y(aux_y[1]) < P_y(aux_y[0])) {
				  		numPoints = 0;
				  		break;
					}
					else{
							if(EN_id(aux_x[0]) !=  EN_id(aux_y[0]))
								throw Exception(__LINE__,__FILE__,"shit happens...\n");

							if(EN_id(aux_x[1]) !=  EN_id(aux_y[1]))
								throw Exception(__LINE__,__FILE__,"shit always happens...\n");

							if(EN_id(aux_x[0]) != EN_id(aux_x[1])){
								numPoints = 2;
								A = aux_x[0];
								B = aux_x[1];
							}
							else{
								numPoints = 1;
								A = aux_x[0];
							}
					}
				}
				break;
	}

	//printf("\nPoints P1 e P2: (%lf, %lf) e (%lf, %lf)\n", P_x(P1),P_y(P1),  P_x(P2),P_y(P2));	
	//printf("\nPoints Q1 e Q2: (%lf, %lf) e (%lf, %lf)\n", P_x(Q1), P_y(Q1), P_x(Q2), P_y(Q2));	
	//printf("\nPowers: %lf %lf %lf %lf\n", Power_P1, Power_P2, Power_Q1, Power_Q2);
	
	if(numPoints == 1){
		intersect_points.push_back(A);
		//printf("uma intersecao: (%lf, %lf)\n", P_x(A), P_y(A));
	}
	else if (numPoints == 2){
		intersect_points.push_back(A);
		intersect_points.push_back(B);
		//printf("duas intersecoes: (%lf, %lf) e (%lf, %lf)\n", P_x(A), P_y(A),  P_x(B), P_y(B));		
	}
	return intersect_points;
	
	

}

template <typename T> int sgn(T val) {
	
 return ((T(0) + EPSILON) < val) - (val < (T(0) - EPSILON));
}
/*
template <typename T> int sgn(T val) {
	return (T(0) < val) - (val < T(0));
} */

bool point_insideTriangle(pVertex vertex, pFace triangle) {
	double a, b, c;
	double xyz[3];
	double p1x, p1y, p2x, p2y, p3x, p3y, px, py;

	//point coordinates
	V_coord(vertex, xyz);
	px = xyz[0];
	py = xyz[1];

	//printf("vertice\n" );
	//printf("%lf %lf\n", px, py);
	
	//triangle vertexes coordinates
	V_coord(F_vertex(triangle, 0),xyz);
	p1x = xyz[0];
	p1y = xyz[1];
	V_coord(F_vertex(triangle, 1), xyz);
	p2x = xyz[0];
	p2y = xyz[1];
	V_coord(F_vertex(triangle, 2), xyz);	
	p3x = xyz[0];
	p3y = xyz[1];


	//printf("triangulozinho\n");
	//printf("%lf %lf\n", p1x, p1y);
	//printf("%lf %lf\n", p2x, p2y);
	//printf("%lf %lf\n", p3x, p3y);

	a = (p1x - px)*(p2y - py) - (p2x - px)*(p1y - py);
	//printf("valor de a = %lf", a);
	b = (p2x - px)*(p3y - py) - (p3x - px)*(p2y - py);
	c = (p3x - px)*(p1y - py) - (p1x - px)*(p3y - py);
	
//	printf("sinais\n");
//	printf("a = %d  b = %d  c = %d ", sgn(a), sgn(b), sgn(c));

	if(a == 0 || b == 0 || c == 0)
		return true;

	return (sgn(a) == sgn(b) && sgn(b) == sgn(c));
}
/*
double triangle_area(cv::Point p1, cv::Point p2, cv::Point p3) {
	return (p1.x*p2.y + p2.x*p3.y + p3.x*p1.y - p2.x*p1.y - p3.x*p2.y - p1.x*p3.y) / 2000.0;
}

*/
/*
bool point_insideTriangle(pVertex vertex, pFace triangle){
	double xyz[3];
	double p0x, p0y, p1x, p1y, p2x, p2y, px, py, area, s, t;
	
	//point coordinates
	V_coord(vertex, xyz);
	px = xyz[0];
	py = xyz[1];

	//printf("vertice\n" );
	//printf("%lf %lf\n", px, py);
	
	//triangle vertexes coordinates
	V_coord(F_vertex(triangle, 0),xyz);
	p0x = xyz[0];
	p0y = xyz[1];
	V_coord(F_vertex(triangle, 1), xyz);
	p1x = xyz[0];
	p1y = xyz[1];
	V_coord(F_vertex(triangle, 2), xyz);	
	p2x = xyz[0];
	p2y = xyz[1];

	//printf("triangulozinho\n");
	//printf("%lf %lf\n", p0x, p0y);
	//printf("%lf %lf\n", p1x, p1y);
	//printf("%lf %lf\n", p2x, p2y);



	//signed area of the triangle
	area = (1/2)*(-p1y*p2x + p0y*(-p1x + p2x) + p0x*(p1y - p2y) + p1x*p2y);

	//barycentrics coordinates
	s = (1/(2*area))*(p0y*p2x - p0x*p2y + (p2y - p0y)*px + (p0x - p2x)*py);
	t = (1/(2*area))*(p0x*p1y - p0y*p1x + (p0y - p1y)*px + (p1x - p0x)*py);

	return ((s>=0) && (t>=0) && (1-s-t>=0));	
}
*/
void triangulate_cloud(vector<pPoint> cloud_points){
	


}
