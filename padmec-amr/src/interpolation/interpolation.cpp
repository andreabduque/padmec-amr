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

double signed_area(double p1x, double p1y, double p2x, double p2y, double p3x, double p3y){
	return (p1x*p2y + p2x*p3y + p3x*p1y - p2x*p1y - p3x*p2y - p1x*p3y);
}/*
bool point_insideTriangle(pVertex vertex, pFace triangle){
	double xyz[3];
	double p0x, p0y, p1x, p1y, p2x, p2y, px, py, area, s, t;
	
	//point coordinates
	V_coord(vertex, xyz);
	px = xyz[0];
	py = xyz[1];

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

	//signed area of the triangle
	area = (1/2)*(-p1y*p2x + p0y*(-p1x + p2x) + p0x*(p1y - p2y) + p1x*p2y);

	//barycentrics coordinates
	s = (1/(2*area))*(p0y*p2x - p0x*p2y + (p2y - p0y)*px + (p0x - p2x)*py);
	t = (1/(2*area))*(p0x*p1y - p0y*p1x + (p0y - p1y)*px + (p1x - p0x)*py);



	return ((1.0f <= s >= 0.0f) && (1.0f <= t>= 0.0f) && (1.0f <=(1-s-t)>= 0.0f));	
}
*/


bool comp_sgn(int a, int b){
		return a == b || a == 0 || b == 0;
}

bool point_insideTriangle(pVertex vertex, pFace triangle) {
	double a, b, c;
	double xyz[3];
	double p1x, p1y, p2x, p2y, p3x, p3y, px, py;

	//point coordinates
	V_coord(vertex, xyz);
	px = xyz[0];
	py = xyz[1];

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

	a = (p1x - px)*(p2y - py) - (p2x - px)*(p1y - py);
	b = (p2x - px)*(p3y - py) - (p3x - px)*(p2y - py);
	c = (p3x - px)*(p1y - py) - (p1x - px)*(p3y - py);	
	return (comp_sgn(sgn(a),sgn(b)) && comp_sgn(sgn(b),sgn(c)) && comp_sgn(sgn(a), sgn(c)));
}

struct cPoint
{
	pPoint pp;
	int id;	

	cPoint(pPoint a, int b) : pp(a), id(b) {};	

};

struct cEdge
{
	pair<cPoint, cPoint> edge;	
	int count;
	cEdge(pair<cPoint, cPoint> e) : edge(e) {count = 0;};		
};


struct cTriangle
{
	cPoint p0, p1, p2;
	

	cTriangle(cPoint a, cPoint b, cPoint c) : p0(a), p1(b), p2(c) {};
	
};

//metodo iterativo
double triangulate_cloud(list<pPoint> cloud_points, pFace backface){
	//pair of two points
	//ccPoint p0, p1, p2, point;	
	list< cEdge > edge_list;
	//new triangulation
	list< cTriangle > new_mesh;	
	list< pPoint > backup = cloud_points;
	list< cPoint > all_points;
	int sign1, sign2;
	double area1, area2, totalArea = 0, totalMass = 0, aux_area;

	//contructing cPoint list of pPoint
	int id = 1;
	for (std::list<pPoint>::iterator pit = cloud_points.begin(); pit != cloud_points.end(); pit++){
		all_points.push_back(cPoint(*pit, id));
		id++;
	}		

	cPoint p0 = all_points.back();
	all_points.pop_back();
	cPoint p1 = all_points.back();
	all_points.pop_back();
	cPoint p2 = all_points.back();
	all_points.pop_back();

	//inserindo edges a serem contabilizados
	edge_list.push_back(cEdge(make_pair(p0, p1)));
	edge_list.push_back(cEdge(make_pair(p1, p2)));
	edge_list.push_back(cEdge(make_pair(p2, p0)));

	//novo triangulo para a malha
	new_mesh.push_back(cTriangle(p0,p1,p2)); 

	//dado um edge encontrar sua face. jogar id do edge, localizar a face
	//printf("comeco triangulacao\n");
	//printf("%d %d %d\n", p0.id, p1.id, p2.id);
	//printf("%f, %f\n%f, %f\n%f, %f\n", P_x(p0), P_y(p0), P_x(p1), P_y(p1), P_x(p2), P_y(p2));
	area1 = signed_area(P_x(p0.pp), P_y(p0.pp), P_x(p1.pp), P_y(p1.pp), P_x(p2.pp), P_y(p2.pp));
	sign1 = sgn(area1);

	totalArea = abs(area1);	

	totalMass = calculate_elementMass(P_x(p0.pp), P_y(p0.pp), P_x(p1.pp), P_y(p1.pp), P_x(p2.pp), P_y(p2.pp), abs(area1));

	//lista de pontoss
	//printf("--------lista indice de pontos--------\n");
	//for (std::list<cPoint>::iterator pit = all_points.begin(); pit != all_points.end(); pit++)
	//	printf("%d\n",pit->id);
	//printf("fim da lista de indice de pontos--------\n");

	for (std::list<cEdge>::iterator eit = edge_list.begin(); eit != edge_list.end(); eit++){
		for (std::list<cPoint>::iterator pit = all_points.begin(); pit != all_points.end(); pit++){
			p0 = eit->edge.first;
			p1 = eit->edge.second;
			cPoint point = *pit;

			//printf("edge utilizado para comparacao\n");
			//printf("%lf %lf\n", P_x(p0), P_y(p0));
			//printf("%lf %lf\n", P_x(p1), P_y(p1));
			//printf("triangulo 2 utilizado para comparacao\n");
			//printf("%f, %f\n%f, %f\n%f, %f\n", P_x(p0), P_y(p0), P_x(p1), P_y(p1), P_x(point), P_y(point));
			area2 = signed_area(P_x(p0.pp), P_y(p0.pp), P_x(p1.pp), P_y(p1.pp), P_x(point.pp), P_y(point.pp));
			sign2 = sgn(area2);
			//printf("areas: %lf %lf\nsgns: %d %d\n", area1, area2, sign1, sign2);

			if(sign1 != sign2){
				p0 = eit->edge.second;
				p1 = eit->edge.first;

				new_mesh.push_back(cTriangle(p0, p1, point));
				//edge_list.push_back(cEdge(make_pair(p0, point)));
				//edge_list.push_back(cEdge(make_pair(point, p1)));
				edge_list.push_back(cEdge(make_pair(point, p0)));
				edge_list.push_back(cEdge(make_pair(p1, point)));
				//printf("-----------formou novo triangulo----------\n");
				//printf("%f, %f\n%f, %f\n%f, %f\n", P_x(p0), P_y(p0), P_x(p1), P_y(p1), P_x(point), P_y(point));
				//printf("%d %d %d\n", p0.id, p1.id, point.id);

				//calculating triangle mass
				aux_area =  abs(signed_area(P_x(p0.pp), P_y(p0.pp),  P_x(p1.pp), P_y(p1.pp), P_x(point.pp), P_y(point.pp)));
				totalArea += aux_area;
				totalMass += calculate_elementMass(P_x(p0.pp), P_y(p0.pp), P_x(p1.pp), P_y(p1.pp), P_x(point.pp), P_y(point.pp), aux_area);

				//printf("-----------end----------\n");
				pit = all_points.erase(pit);			
				break;									
				
			}		
			//caso encontre ponto, remove
			//caso nao encontra ponto, remove de qualquer forma
			//assim a lista nao vai ser infinita

			//encotnrou ponto, add no mesh
			//add mais dois edges a lista
		}

		eit = edge_list.erase(eit);
	}

	return totalMass;
}

double circle_func(double x, double y){

	return x + y;

}


double calculate_elementMass(pFace triangle){
	//triangle vertexes coordinates
	double xyz[3], p0x, p0y, p1x, p1y, p2x, p2y, area, centerx, centery;	

	V_coord(F_vertex(triangle, 0),xyz);
	p0x = xyz[0];
	p0y = xyz[1];
	V_coord(F_vertex(triangle, 1), xyz);
	p1x = xyz[0];
	p1y = xyz[1];
	V_coord(F_vertex(triangle, 2), xyz);	
	p2x = xyz[0];
	p2y = xyz[1];

	centerx = (p0x + p1x + p2x)/3;
	centery = (p0y + p1y + p2y)/3;
	area = abs(signed_area(p0x, p0y, p1x, p1y, p2x, p2y));

	return circle_func(centerx, centery)*area;

}


double calculate_elementMass(pFace triangle, double area){
	//triangle vertexes coordinates
	double xyz[3], p0x, p0y, p1x, p1y, p2x, p2y, centerx, centery;	

	V_coord(F_vertex(triangle, 0),xyz);
	p0x = xyz[0];
	p0y = xyz[1];
	V_coord(F_vertex(triangle, 1), xyz);
	p1x = xyz[0];
	p1y = xyz[1];
	V_coord(F_vertex(triangle, 2), xyz);	
	p2x = xyz[0];
	p2y = xyz[1];

	centerx = (p0x + p1x + p2x)/3;
	centery = (p0y + p1y + p2y)/3;
	
	return circle_func(centerx, centery)*area;

}


double calculate_elementMass(double p0x, double p0y, double p1x, double p1y, double p2x, double p2y, double area){
	//triangle vertexes coordinates
	double centerx, centery;	
	centerx = (p0x + p1x + p2x)/3;
	centery = (p0y + p1y + p2y)/3;
	
	return circle_func(centerx, centery)*area;

}
