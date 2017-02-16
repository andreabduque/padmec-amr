#include "interpolation.h"

//creates initial list of background elements
queue<pFace> create_initial_list(pFace new_meshFace, Octree* back_meshOctree, vector<int>* overlapped_IDelements){
	double xyz[3];
	queue<pFace> overlapped_elements;
	vector<int>::iterator it;

	for(int i = 0; i < 3; i++){
		pVertex vertex =  F_vertex(new_meshFace, i);

		V_coord(vertex,xyz);
		pFace background_face = (pFace)Octree_Search(xyz, back_meshOctree);
		if (!background_face){
			cout << "Finding element containing coordenate: " << xyz[0] << "\t" << xyz[1] << "\t" << xyz[2] << endl;
			throw Exception(__LINE__,__FILE__,"meshes dont overlap\n");
		}

		it = find(overlapped_IDelements->begin(), overlapped_IDelements->end(), EN_id(background_face));
		//add new face if it hasnt been visited
		if (it == overlapped_IDelements->end()){
			overlapped_elements.push(background_face);
			overlapped_IDelements->push_back(EN_id(background_face));
		}

	}
	return overlapped_elements;
}

//set entities ids
void set_entityID(pMesh mesh){
	FIter fit =  M_faceIter(mesh); //new mesh face iterator
	EIter eit =  M_edgeIter(mesh); //new mesh edge iterator
	unsigned int i = 1;

	while (pEntity edge = EIter_next(eit) ){
		EN_setID(edge, i++);
	}
	EIter_delete(eit);
	i = 1;
	while (pEntity face = FIter_next(fit) ){
		EN_setID(face, i++);
	}

	FIter_delete(fit);
}

double power(pPoint other, pPList  vertices_edge){
	pPoint P, Q;
	Vector* v;
	Vector* normal;

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

bool its_me(double x1, double y1, double x2, double y2){
	return (abs(x1 - x2) < EPSILON) && (abs(y1 - y2) < EPSILON);
}

list<pPoint> unique_points(list<pPoint> cloud_points){
	list<pPoint> aux_points;
	bool isThere;

	for (std::list<pPoint>::iterator it = cloud_points.begin(); it != cloud_points.end(); it++){
		isThere = false;
		for (std::list<pPoint>::iterator it2 = aux_points.begin(); it2 != aux_points.end(); it2++){
			if(its_me(P_x(*it), P_y(*it), P_x(*it2), P_y(*it2))){
					isThere = true;
			}
		}
		if(!isThere){
			aux_points.push_back(*it);
		}
	}
	return aux_points;
}

bool tolerance(double value){
	return abs(value - 0) < EPSILON;
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
	bool so_x, so_y;

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

				//intersection intervals in y
				aux_y[0] = max_point(min_point(P1, P2, 1), min_point(Q1, Q2, 1), 1);
				aux_y[1] = min_point(max_point(P1, P2, 1), max_point(Q1, Q2, 1), 1);

				so_x = false;
				so_y = false;

				if(abs(P_x(P1) - P_x(P2)) < EPSILON && abs(P_x(P1) - P_x(Q2)) < EPSILON && abs(P_x(P1) - P_x(Q1)) < EPSILON ){
					if (P_y(max_point(P1, P2, 1)) < P_y(min_point(Q1, Q2, 1))) {
				  		numPoints = 0;
				  		break;
					}
					else{
						A = aux_y[0];
						B = aux_y[1];
					}
					so_y = true;
				}

				if(abs(P_y(P1) - P_y(P2)) < EPSILON && abs(P_y(P1) - P_y(Q2)) < EPSILON && abs(P_y(P1) - P_y(Q1)) < EPSILON ){
					if (P_x(max_point(P1, P2, 0)) < P_x(min_point(Q1, Q2, 0))) {
					  numPoints = 0;
					  break;
					}
					else{
						A = aux_x[0];
						B = aux_x[1];
					}
					so_y = true;
				}

				if(!so_x && !so_y){
					if (P_x(max_point(P1, P2, 0)) < P_x(min_point(Q1, Q2, 0))) {
					  numPoints = 0;
					  break;
					}
					else{
						if (P_y(max_point(P1, P2, 1)) < P_y(min_point(Q1, Q2, 1))) {
					  		numPoints = 0;
					  		break;
						}
						else{
							A = aux_x[0];
							B = aux_y[1];
						}
					}

				}

				if(numPoints != 0){
					if(abs(P_x(A) - P_x(B)) < EPSILON && abs(P_y(A) - P_y(B)) < EPSILON){
						numPoints = 1;
					}
					else{
						numPoints = 2;
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
	return (p1x*p2y + p2x*p3y + p3x*p1y - p2x*p1y - p3x*p2y - p1x*p3y)*0.5;
}

double signed_area(pFace triangle){
	double p1x,  p1y,  p2x,  p2y,  p3x,  p3y, xyz[3];

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

	return (p1x*p2y + p2x*p3y + p3x*p1y - p2x*p1y - p3x*p2y - p1x*p3y)*0.5;
}

bool comp_sgn(int a, int b){
		return a == b || a == 0 || b == 0;
}

bool point_inside_triangle(pVertex vertex, pFace triangle) {
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

//somo useful stuff for triangulation
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


void show_face(pFace face){
	double xyz1[3], xyz2[3], xyz3[3];

	V_coord(F_vertex(face, 0), xyz1);
	V_coord(F_vertex(face, 1), xyz2);
	V_coord(F_vertex(face, 2), xyz3);

	printf("face of ID %d \n", EN_id(face));
	//x coords
	printf("x = [\n");
	printf("%lf\n", xyz1[0]);
	printf("%lf\n", xyz2[0]);
	printf("%lf\n", xyz3[0]);
	printf("]\n");
	//y coords
	printf("y = [\n");
	printf("%lf\n", xyz1[1]);
	printf("%lf\n", xyz2[1]);
	printf("%lf\n", xyz3[1]);
	printf("]\n");

}


//calculates centroid node solution value
void centroid_value(InterpolationDataStruct* pIData, pFace backface, double* field_val, double* grad_val){
		int vertex_id, dim = 3; //3 for triangle
		double scalar, val;

		printf("-----Calculando valor de sol para backface-------\n");
		show_face(backface);

		// interpolating for each scalar field: Sw and p
		for (int field=0; field < 2; field++){
			val = 0;
			printf("---Field %d-----\n", field);
				for(int i= 0; i < dim; i++){
					vertex_id = EN_id(F_vertex(backface, i));

					pIData->pGetDblFunctions[field](vertex_id - 1, scalar);		// get value from old mesh
					val += scalar*(1/3);
					printf("scalar %lf\n", scalar);
				}
				field_val[field] = val;
		}

		//interpolating gradients

}


//particular case of iterative delaunay.
//triangulates a convex cloud of points and calculates mass
double triangulate_cloud(InterpolationDataStruct* pIData, list<pPoint> cloud_points, pFace backface){
	//pair of two points
	list< cEdge > edge_list;
	//cloud of cPoints
	list< cPoint > all_points;
	int sign1, sign2, field = 0; //interpolate pressure and saturation and gradients
	double area1, area2, total_mass = 0, aux_area, field_val[2], grad_val[4];
	std::list<cEdge>::iterator aux_it;

	//calculating solution at centroid
	centroid_value(pIData, backface, field_val, grad_val);

	//constructing cPoint list of pPoint
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

	//adding edges that are not connected to a point yet
	edge_list.push_back(cEdge(make_pair(p0, p1)));
	edge_list.push_back(cEdge(make_pair(p1, p2)));
	edge_list.push_back(cEdge(make_pair(p2, p0)));

	//first random triangle from point cloud
	area1 = signed_area(P_x(p0.pp), P_y(p0.pp), P_x(p1.pp), P_y(p1.pp), P_x(p2.pp), P_y(p2.pp));
	sign1 = sgn(area1);

	total_mass = abs(area1);

	for (std::list<cEdge>::iterator eit = edge_list.begin(); eit != edge_list.end(); eit++){
		for (std::list<cPoint>::iterator pit = all_points.begin(); pit != all_points.end(); pit++){
			p0 = eit->edge.first;
			p1 = eit->edge.second;
			cPoint point = *pit;

			area2 = signed_area(P_x(p0.pp), P_y(p0.pp), P_x(p1.pp), P_y(p1.pp), P_x(point.pp), P_y(point.pp));
			sign2 = sgn(area2);

			//new triangle created
			if(sign1 != sign2){
				p0 = eit->edge.second;
				p1 = eit->edge.first;

				edge_list.push_back(cEdge(make_pair(point, p0)));
				edge_list.push_back(cEdge(make_pair(p1, point)));

				aux_area =  abs(signed_area(P_x(p0.pp), P_y(p0.pp),  P_x(p1.pp), P_y(p1.pp), P_x(point.pp), P_y(point.pp)));
				total_mass += aux_area;

				pit = all_points.erase(pit);
				break;
			}
		}
		eit = edge_list.erase(eit);
	}

	//cant forget this triangle
	if(edge_list.size() == 1 && all_points.size() == 1){

		p0 = edge_list.back().edge.first;
		p1 = edge_list.back().edge.second;
		edge_list.pop_back();

		cPoint point =  all_points.back();
		all_points.pop_back();

		aux_area =  abs(signed_area(P_x(p0.pp), P_y(p0.pp),  P_x(p1.pp), P_y(p1.pp), P_x(point.pp), P_y(point.pp)));
		total_mass += aux_area;

	}

	return total_mass*field_val[field];
}

double circle_func(double x, double y){

	return 1;

}


double gradx(){

	return 1;

}
double grady(){

	return 1;

}

double calculate_solution_center(pFace triangle){
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

	return circle_func(centerx, centery);

}

double calculate_element_mass(InterpolationDataStruct* pIData, pFace triangle){
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

double calculate_element_mass(InterpolationDataStruct* pIData, pFace triangle, double area){
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

double calculate_element_mass(double p0x, double p0y, double p1x, double p1y, double p2x, double p2y, double area){
	//triangle vertexes coordinates
	double centerx, centery;
	centerx = (p0x + p1x + p2x)/3;
	centery = (p0y + p1y + p2y)/3;

	return circle_func(centerx, centery)*area;

}

void debug_cloud(list<pPoint> cloud_points){
	//printf("calculating cloud for face of id %d and vertices \n", EN_id(new_mesh));
	printf("-------------------\n");
	printf("\ndebugging convex polygon\n");
	printf("number of points: %d\n", cloud_points.size());
	printf("x = \n");
	for (std::list<pPoint>::iterator it = cloud_points.begin() ; it != cloud_points.end(); it++){
		printf("%lf\n", P_x(*it));
	}
	printf("y = \n");

	for (std::list<pPoint>::iterator it = cloud_points.begin() ; it != cloud_points.end(); it++){
		printf("%lf\n", P_y(*it));
	}

	printf("-------------------\n");

}


double max_value(InterpolationDataStruct* pIData, pFace triangle){
	double values[3];
	int vertex_id, dim = 3, field = 0; //3 for triangle, 0 for pressure

	// printf("-----Calculando valor de sol para backface-------\n");
	// show_face(backface);

		for(int i= 0; i < dim; i++){
			vertex_id = EN_id(F_vertex(triangle, i));

			pIData->pGetDblFunctions[field](vertex_id - 1, values[i]);		// get value from old mesh
		//	printf("scalar %lf\n", scalar);
		}

	return max(values[0],max(values[1],values[2]));
}


double min_value(InterpolationDataStruct* pIData, pFace triangle){
	double values[3];
	int vertex_id, dim = 3, field = 0; //3 for triangle, 0 for pressure

	// printf("-----Calculando valor de sol para backface-------\n");
	// show_face(backface);

		for(int i= 0; i < dim; i++){
			vertex_id = EN_id(F_vertex(triangle, i));

			pIData->pGetDblFunctions[field](vertex_id - 1, values[i]);		// get value from old mesh
		//	printf("scalar %lf\n", scalar);
		}

	return min(values[0],min(values[1],values[2]));
}


double mesh_intersection(InterpolationDataStruct* pIData, pFace new_face, queue<pFace> overlapped_elements, vector<int> overlapped_IDelements, double* vals){
	//intersection polygon between 2 triangles
	list<pPoint> cloud_points;
	//intersection points between two edges
	vector<pPoint> aux_inter;

	//interpolated mass and real mass
	double interp_mass = 0, aux_mass = 0, aux_max, aux_min, max_val, min_val;
	bool debug = false, first = true;

	//auxiliar
	pPList faces_ofAvertex;
	pEntity edge1, edge2;
	pFace face, old_face;
	vector<int>::iterator it;

	if(debug){
		printf("New mesh ");
		show_face(new_face);
	}

	//compute intersection for all overlapped elements
	while(!overlapped_elements.empty()){

		old_face = overlapped_elements.front();
		overlapped_elements.pop();

		//calculating max and min value of vertices solutions from overlapped backmesh elements
		aux_max = max_value(pIData, old_face);
		aux_min = min_value(pIData, old_face);
		if(!first){
			max_val = max(aux_max, max_val);
			min_val = min(aux_min, min_val);
		}
		else{
			max_val = aux_max;
			min_val = aux_min;
		}

		first = false;

		if(debug){
			printf("Old face ");
			show_face(old_face);
		}

		//if vertice is inside new mesh, add vertex ball to list of overlapped elements
		for(int i = 0; i < 3; i++){
			if(point_inside_triangle(F_vertex(old_face, i), new_face)){
				faces_ofAvertex = V_faces((pVertex)F_vertex(old_face, i));

				for(int j = 0; j < PList_size(faces_ofAvertex); j++){
					face = (pFace)PList_item(faces_ofAvertex, j);
					//search if element of backmesh has already been visited
					it = find(overlapped_IDelements.begin(), overlapped_IDelements.end(), EN_id(face));
					//add new face if it hasnt
					if (it == overlapped_IDelements.end()){
						overlapped_elements.push(face);
						overlapped_IDelements.push_back(EN_id(face));
					}
				}
				cloud_points.push_back(V_point(F_vertex(old_face, i)));
			}
		}
		//marking visited element
		overlapped_IDelements.push_back(EN_id(old_face));

		//edge-edge intersections between element found in back mesh and element of new mesh
		for (int edge_face2 = 0; edge_face2 < 3; edge_face2++){
			for(int edge_face1 = 0; edge_face1 < 3; edge_face1++){

				//edge from new mesh
				edge2 = F_edge(new_face, edge_face2);
				//edge from old mesh
				edge1 = F_edge(old_face, edge_face1);

				//vector of intersection points can have zero, one or two intersection points
				aux_inter = edge_intersection(edge1, edge2);

				if (!aux_inter.empty()){
					cloud_points.insert( cloud_points.end(), aux_inter.begin(), aux_inter.end());

					//add the neighbouring triangle to the list
					for(unsigned int num = 0; num < E_numFaces(edge1); num++){
						face = E_face(edge1, num);
						if(EN_id(face) != EN_id(old_face)){
							face = E_face(edge1, num);
							break;
						}
					}
					it = find(overlapped_IDelements.begin(), overlapped_IDelements.end(), EN_id(face));
					//add new face if it hasnt been visited yet
					if (it == overlapped_IDelements.end()){
						overlapped_elements.push(face);
						overlapped_IDelements.push_back(EN_id(face));
					}
				}
			}
		}

		//add points that are inside of triangle of new mesh the cloud
		for(int i = 0; i < 3; i++){
			if(point_inside_triangle(F_vertex(new_face, i), old_face)){
				cloud_points.push_back(V_point(F_vertex(new_face, i)));
			}
		}
		//remove equal points
		cloud_points = unique_points(cloud_points);

		if(debug){
			debug_cloud(cloud_points);
		}
		if(cloud_points.size() >= 3){
			aux_mass = triangulate_cloud(pIData, cloud_points, old_face);
			interp_mass += aux_mass;

			if(debug){
				printf("polygon mass %lf\n", aux_mass);
			}
		}
		cloud_points.clear();
	}

	vals[0] = max_val;
	vals[1] = min_val;

	return interp_mass/abs(signed_area(new_face));
}



void get_center_coords(pFace triangle, double* center){
	double p0[3],  p1[3],  p2[3];

	V_coord(F_vertex(triangle, 0), p0);
	V_coord(F_vertex(triangle, 1), p1);
	V_coord(F_vertex(triangle, 2), p2);

	center[0] = (p0[0] + p1[0] + p2[0])/3;
	center[1] = (p0[1] + p1[1] + p2[1])/3;

}

vector< pair<int, double> > max_principle(InterpolationDataStruct* pIData, double interp_val, double max_val, double min_val, pFace new_face){
	double centerx, centery, uk_pi, center[2], phi, p[3];
	vector< pair<int, double> > solution;
	pVertex vertex;

	get_center_coords(new_face, center);
	centerx = center[0];
	centery = center[1];

	for(int i = 0; i < 3; i++){
		vertex = F_vertex(new_face, i);
		V_coord(vertex, p);

		uk_pi = interp_val + gradx()*(p[0] - centerx) + grady()*(p[1] - centery);

		//maximum principle correction test
		if(!((uk_pi >= min_val) && (uk_pi <= max_val))){
			//value needs to be corrected
			if(uk_pi - interp_val > EPSILON)
				phi = min(1.0, (max_val - interp_val)/(uk_pi - interp_val));
			else if (uk_pi - interp_val < EPSILON)
				phi = min(1.0, (min_val - interp_val)/(uk_pi - interp_val));
			else
				phi = 1;

			uk_pi = phi*(uk_pi - interp_val) + interp_val;

		}
		solution.push_back(make_pair(EN_id(vertex), uk_pi));
	}

	return solution;
}

void extrapolate_sol_to_vertices(SolutionMap sol_by_vertices, pMesh new_mesh){
	VIter vit = M_vertexIter(new_mesh);
	pair <SolutionMap::iterator, SolutionMap::iterator> ret;
	pair< pFace, double> aux;
	double num = 0, denom = 0, area = 0, aux_value, interp_value, real_value, p[3], norm;
	int count;

	while (pVertex vertex = VIter_next(vit) ){
   		 //ponderating corrected values by elements area
   		 ret = sol_by_vertices.equal_range(EN_id(vertex));
   		 count = 0;
   		 denom = 0;
   		 num = 0;
	   	 for (SolutionMap::iterator it = ret.first; it != ret.second; ++it){
	   	 	aux = it->second;
	   	 	aux_value = aux.second;

	   	 	area = abs(signed_area(aux.first));

	   	 	denom += area;
	   	 	num += area*aux_value;
	   	 	count++;
	   	 }

	   	 if(denom != 0){
	 		interp_value = num/denom;
	   	 }
	   	 else{
	   	 	num = 0;
	   	 }

	   	 V_coord(vertex, p);
	   	 real_value = circle_func(p[0], p[1]);
	   	 norm = abs(real_value*real_value - interp_value*interp_value)/(real_value*real_value);
	   	 printf("real value: %lf interpolated value: %lf\n\n\n", real_value, interp_value);
	   	 printf("norma para vertice de id %d eh %lf\n", EN_id(vertex), norm);

	}
}
