/*
 * interpolateBetweenMeshes.cpp
 *
 *  Created on: 11/02/2013
 *      Author: rogsoares
 */

#include "interpolation.h"

void calculate_GeometricCoefficients(InterpolationDataStruct* pIData, int dim){
	pEntity vertice_1,vertice_2,vertice_3;
	double ptn1[3],ptn2[3], ptn3[3], xyz[3];

	if(dim==2){
		VIter vit = M_vertexIter(pIData->m1);			// mesh to: m1
		while ( pEntity vertex = VIter_next(vit) ){
			 int row = EN_id(vertex);
			V_coord(vertex,xyz);

			// search entity in which the point above are
			pEntity element = (pEntity)Octree_Search(xyz,pIData->theOctree);
			if (!element){
				cout << "Finding element containing coordenate: " << xyz[0] << "\t" << xyz[1] << "\t" << xyz[2] << endl;
				throw Exception(__LINE__,__FILE__,"calculate_GeometricCoefficients: Entity not found! Exiting...\n");
			}

			//get coordinates
			vertice_1 = element->get(0,0);
			vertice_2 = element->get(0,1);
			vertice_3 = element->get(0,2);
			V_coord(vertice_1,ptn1);
			V_coord(vertice_2,ptn2);
			V_coord(vertice_3,ptn3);

			//calculate areas
			double a1 = F_area(xyz,ptn2,ptn3);
			double a2 = F_area(xyz,ptn1,ptn3);
			double a3 = F_area(xyz,ptn1,ptn2);
			double a = a1 + a2 + a3;
			a1 = (double)a1/a;
			a2 = (double)a2/a;
			a3 = (double)a3/a;

			pIData->pGeomCoeff(row,0) = a1;
			pIData->pGeomCoeff(row,1) = a2;
			pIData->pGeomCoeff(row,2) = a3;
		}
		VIter_delete(vit);
	}
	//	else {
	//		VIter vit = M_vertexIter(pIData->m2);
	//		while ( pEntity vertex = VIter_next(vit) ){
	//
	//			// get coordinate of a point from the mesh to be interpolated
	//			double xyz[3];
	//			V_coord(vertex,xyz);
	//			int row = EN_id(vertex);
	//			// search entity in which the point above are
	//
	//			pEntity element = (pEntity)Octree_Search(xyz,pIData->theOctree);
	//			if (!element){
	//				cout << "Entity not found! Exiting...\n";
	//				exit(1);
	//			}
	//
	//			//get coordinates
	//			vertice_1= element->get(0,0); V_coord(vertice_1,ptn1);
	//			vertice_2= element->get(0,1); V_coord(vertice_2,ptn2);
	//			vertice_3= element->get(0,2); V_coord(vertice_3,ptn3);
	//			vertice_4= element->get(0,3); V_coord(vertice_4,ptn4);
	//
	//			//calculate areas
	//			double v1 = R_Volume(xyz,ptn2,ptn3,ptn4);
	//			double v2 = R_Volume(xyz,ptn1,ptn3,ptn4);
	//			double v3 = R_Volume(xyz,ptn1,ptn2,ptn4);
	//			double v4 = R_Volume(xyz,ptn1,ptn2,ptn3);
	//			double v = v1 + v2 + v3 + v4;
	//			v1 = v1/v; v2= v2/v; v3= v3/v; v4= v4/v;
	//
	//			//calculate coefficients
	//			setGeomCoeff(row,0,v1);
	//			setGeomCoeff(row,1,v2);
	//			setGeomCoeff(row,2,v3);
	//			setGeomCoeff(row,3,v4);
	//		}
	//
	//		VIter_delete(vit);
	//	}
}

void calculate_LinearInterpolation(InterpolationDataStruct* pIData, int dim){
	double xyz[3], scalar, val, geom;
	int size = dim + 1;
	int field, rows[3];

	// Loop over new mesh's vertices. They will receive interpolated data.
	VIter vit = M_vertexIter(pIData->m1);
	while ( pEntity vertex = VIter_next(vit) ){
		int row = EN_id(vertex);
		V_coord(vertex,xyz);

		// search on old mesh the element that contains vertex of new mesh
		pEntity element = (pEntity)Octree_Search(xyz,pIData->theOctree);

		for (int i=0; i<size; i++){
			rows[i] = EN_id(element->get(0,i))-1;
		}

		for (int field=0; field<pIData->numFields; field++){		// for each scalar field: Sw and p
			val = 0;
			for(int i=0; i<size; i++){
				pIData->pGetDblFunctions[field](rows[i],scalar);		// get value from old mesh
				geom = pIData->pGeomCoeff(row,i);
				val += scalar*geom;
			}
			pIData->pSetDblFunctions[field](row-1,val);				// set value to new mesh
		}
	}
	VIter_delete(vit);
}

double calculate_QuadraticInterpolation(InterpolationDataStruct* pIData, int field){
	double GeomCoeff, residual, coord[3], linear, quadratic;
	pEntity v;
	double xyz[3];
	int dim = pIData->m2->getDim();
	int idx = 0;
	VIter vit = M_vertexIter(pIData->m1);	// m1: mesh to
	while ( pEntity vertex = VIter_next(vit) ){
		V_coord(vertex,xyz);
		int row = EN_id(vertex);
		//search entity in which the node is
		pEntity element = (pEntity)Octree_Search(xyz,pIData->theOctree);
		if (!element){
			cout << "WARNING:  Entity not found!\n";
		}

		//calculate residual
		residual = 0;
		for( int i = 0; i < dim+1; i++){
			v = element->get(0,i);			// v is a vertex of mesh from (m2)
			int id = EN_id(v);				// id is the vertex ID fo mesh from (m2)
			if (id > (M_numVertices(pIData->m2)+1)){
				throw Exception(__LINE__,__FILE__,"Index out of bound!");
			}
			V_coord(v,coord);
			GeomCoeff = pIData->pGeomCoeff(row,i);
			double dx = pIData->pGrad(id,0);
			double dy = pIData->pGrad(id,1);
			double dx2 = pIData->pGrad(id,2);
			double dy2 = pIData->pGrad(id,3);
			double dxdy = pIData->pGrad(id,4);

			double Res = (xyz[0]-coord[0])*dx + (xyz[1]-coord[1])*dy;
			Res += (xyz[0]-coord[0])*(xyz[0]-coord[0])*dx2/2;
			Res += (xyz[1]-coord[1])*(xyz[1]-coord[1])*dy2/2;
			Res += (xyz[0]-coord[0])*(xyz[1]-coord[1])*dxdy;
			residual += GeomCoeff*Res;
		}

		//calculate quadratic interpolation
		pIData->pGetDblFunctions[field](idx,linear);
		quadratic = linear + residual;
		pIData->pSetDblFunctions[field](idx,quadratic);
		idx++;
	}
	VIter_delete(vit);
	return 0;
}

//creates initial list of background elements
queue<pFace> create_initialList(pFace new_meshFace, Octree* back_meshOctree){
	double xyz[3];	
	queue<pFace> overlapped_elements;

	for(int i = 0; i < 3; i++){
		pVertex vertex =  F_vertex(new_meshFace, i);	

		V_coord(vertex,xyz);
		pFace background_face = (pFace)Octree_Search(xyz, back_meshOctree);
		if (!background_face){
			cout << "Finding element containing coordenate: " << xyz[0] << "\t" << xyz[1] << "\t" << xyz[2] << endl;		
			throw Exception(__LINE__,__FILE__,"meshes dont overlap\n");
		}
		overlapped_elements.push(background_face);
	}
	return overlapped_elements;
}

void debug_cloud(pFace new_mesh, vector<pPoint> cloud_points){
	double xyz[3];
	printf("calculating cloud for face of id %d and vertices \n", EN_id(new_mesh));

	for(int i = 0; i < 3; i++){
		pVertex vertex =  F_vertex(new_mesh, i);	

		V_coord(vertex,xyz);

		printf("(%lf, %lf) ", xyz[0], xyz[1]);
		
	}

	printf("\nCloud of points: \n");

	for(unsigned int i = 0; i < cloud_points.size(); i++){
			printf("(%lf, %lf)\n", P_x(cloud_points[i]), P_y(cloud_points[i]));
	} 

	printf("-------------------\n");

}

//set entities ids
void setEntityID(pMesh mesh){
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

	/*
	VIter vit = M_vertexIter(mesh) ;
	printf("----debugando id dos vertices-----\n");
	while (pEntity vertex = FIter_next(vit) ){	
		printf("vertex de id %d e coords (%lf,%lf) \n", EN_id(vertex), P_x(V_point(vertex)), P_y(V_point(vertex)));
	}

	VIter_delete(vit); */



}

//m2 is the mesh which sends data. m1 is the mesh to be interpolated.
void calculate_ConservativeInterpolation(InterpolationDataStruct* pIData, int dim){
	vector<pPoint> cloud_points;
	vector<pPoint> aux_inter;
	vector<int> cloud_IDpoints, overlapped_IDelements;
	vector< pair<int,int> > intersec_IDEdges;
	FIter fit2 =  M_faceIter(pIData->m1); //new mesh face iterator  
	queue<pFace> overlapped_elements;	
	pEntity face1, edge2, edge1;
	pPList faces_ofAvertex;
	pFace face;
	//iterators
	
	vector< pair<int,int> >::iterator it_pair, it_pair2;
	vector<int>::iterator it;

	freopen ("cloud_points.txt","w",stdout);


	setEntityID(pIData->m1);
	setEntityID(pIData->m2);

	

	//Loop over faces of new mesh to be interpolated
	while (pEntity face2 = FIter_next(fit2) ){	
		printf("id da face: %d\n",EN_id(face2));			

		//initial list of overlapped elements (elements from backmesh containing the vertices of the element of new mesh)
		overlapped_elements = create_initialList(face2, pIData->theOctree);
		
		
		//compute intersection for all overlapped elements
		while(!overlapped_elements.empty()){
			face1 = overlapped_elements.front();
			overlapped_elements.pop();
			//printf("id da face: %d\n",EN_id(face1));
			
			//marking visited element
			overlapped_IDelements.push_back(EN_id(face1));

			//edge-edge intersections between element found in back mesh and element of new mesh
			for (int edge_face2 = 0; edge_face2 < 3; edge_face2++){
				for(int edge_face1 = 0; edge_face1 < 3; edge_face1++){

					//edge from new mesh
					edge2 = F_edge(face2, edge_face2);
					//edge from old mesh
					edge1 = F_edge(face1, edge_face1);		


					//search if edge intersection has already been calculated
					it_pair = find(intersec_IDEdges.begin(), intersec_IDEdges.end(), make_pair(EN_id(edge2), EN_id(edge1)));
					it_pair2 = find(intersec_IDEdges.begin(), intersec_IDEdges.end(), make_pair(EN_id(edge1), EN_id(edge2)));


					//if it hasnt, calculate new intersection
					if(it_pair == intersec_IDEdges.end() && it_pair2 == intersec_IDEdges.end()){

						//printf("pair de edges: %d e %d \n",EN_id(edge2), EN_id(edge1));			

						//vector of intersection points can have zero, one or two intersection points					
						aux_inter = edge_intersection(edge1, edge2);

						for(unsigned int i = 0; i < aux_inter.size(); i++){					
							//if the point is a vertex, search for id
							if(EN_id(aux_inter[i]) != 0){
								//search if point has already been added to cloud							
								it = find(cloud_IDpoints.begin(), cloud_IDpoints.end(), EN_id(aux_inter[i]));
								//add new point if it hasnt
								if (it == cloud_IDpoints.end()){
									cloud_points.push_back(aux_inter[i]);
									//printf("inseri ponto a partir de intersecao(%lf, %lf)\n",P_x(aux_inter[i]), P_y(aux_inter[i]));
									cloud_IDpoints.push_back(EN_id(aux_inter[i]));
								}

								//point is a vertice, then add vertice ball to queue		
								faces_ofAvertex = V_faces((pVertex)aux_inter[i]);

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

							}
							else{

								//faces que dividem a aresta
								for(unsigned int num = 0; num < E_numFaces(edge1); num++){
									face = E_face(edge1, num);
									if(EN_id(face) != EN_id(face1)){
										face = E_face(edge1, num);
										break;
									}
								}

								 it = find(overlapped_IDelements.begin(), overlapped_IDelements.end(), EN_id(face));

								//add new face if it hasnt
								if (it == overlapped_IDelements.end()){
									overlapped_elements.push(face);	
									overlapped_IDelements.push_back(EN_id(face));								
								}	
								
								cloud_points.push_back(aux_inter[i]);	
								//printf("inseri vertice (%lf, %lf) de id %d\n",P_x(aux_inter[i]), P_y(aux_inter[i]), EN_id(aux_inter[i]));
			
							}
						}
						intersec_IDEdges.push_back(make_pair(EN_id(edge2), EN_id(edge1)));
					}		
				}
			}	
			//ver se o ponto ja foi adionado pra nuvem
			//acho que aqui vai precisar verificar os repetidos, pois nos casos degenerados pode haver
			//add points that are inside of both the triangles to the cloud
			for(int i = 0; i < 3; i++){
				if(point_insideTriangle(F_vertex(face2, i), face1)){
					//search if point has already been added to cloud							
					it = find(cloud_IDpoints.begin(), cloud_IDpoints.end(), EN_id(F_vertex(face2, i)));
					if (it == cloud_IDpoints.end()){
						cloud_points.push_back(V_point(F_vertex(face2, i)));
						//printf("inseri ponto a partir de intersecao(%lf, %lf)\n",P_x(aux_inter[i]), P_y(aux_inter[i]));
						cloud_IDpoints.push_back(EN_id(F_vertex(face2, i)));
					}

					//printf("ponto detrno do triang (%lf, %lf)\n", P_x(V_point(F_vertex(face2, i))), P_y(V_point(F_vertex(face2, i))));
				}
				if(point_insideTriangle(F_vertex(face1, i), face2)){

					it = find(cloud_IDpoints.begin(), cloud_IDpoints.end(), EN_id(F_vertex(face1, i)));
					if (it == cloud_IDpoints.end()){
						cloud_points.push_back(V_point(F_vertex(face1, i)));
						//printf("inseri ponto a partir de intersecao(%lf, %lf)\n",P_x(aux_inter[i]), P_y(aux_inter[i]));
						cloud_IDpoints.push_back(EN_id(F_vertex(face1, i)));
					}
					//printf("ponto detrno do triang (%lf, %lf)\n",P_x(V_point(F_vertex(face1, i))), P_y(V_point(F_vertex(face1, i))));
				}
			}
		}			
		//freopen ("/home/abd/cloud_points.txt","w", stdout);
		//printf("ID da face: %d \n", EN_id(face2));

		debug_cloud(face2, cloud_points);

		//triangulate




		cloud_points.clear();
		cloud_IDpoints.clear();
		intersec_IDEdges.clear();
		overlapped_IDelements.clear();		
	}
	FIter_delete(fit2);
	/*for(unsigned int i = 0; i < cloud_points.size(); i++){
			printf("(%lf, %lf)\n", P_x(cloud_points[i]), P_y(cloud_points[i]));
	} */
	printf("chegueiii");
	fclose (stdout);
	STOP();	
}

void calculate_DerivativesError(InterpolationDataStruct* pIData){
	double summ_dx = 0; double summ2_dx = 0; double MaxError_dx = 0;
	double summ_dy = 0; double summ2_dy = 0; double MaxError_dy = 0;
	double summ_dx2 = 0; double summ2_dx2 = 0; double MaxError_dx2 = 0;
	double summ_dy2 = 0; double summ2_dy2 = 0; double MaxError_dy2 = 0;
	double summ_dxdy = 0; double summ2_dxdy = 0; double MaxError_dxdy = 0;

	//loop on the mesh vertices
	VIter vit = M_vertexIter(pIData->m2);
	while (pEntity vertex = VIter_next(vit)){

		int id = EN_id(vertex);
		//get coordinates
		double xyz[3] = {0,0,0};
		V_coord(vertex, xyz);

		//ATTENTION: THE CORRECT GRADIENTS VALUES MUST BE PUT HERE BEFORE USE THIS FUNCTION
		double dx = 0;
		double dy = 2*xyz[1];
		double dx2 = 0;
		double dy2 = 2;
		double dxdy = 0;

		//calculate the global error (2-norm) related to the gradients
		double error_dx = pIData->pGrad(id,0) - dx;
		double error_dy = pIData->pGrad(id,1) - dy;
		double error_dx2 = pIData->pGrad(id,2) - dx2;
		double error_dy2 = pIData->pGrad(id,3) - dy2;
		double error_dxdy = pIData->pGrad(id,4) - dxdy;

		summ_dx = summ_dx + pow(error_dx,2);
		summ2_dx = summ2_dx + pow(dx,2);

		summ_dy = summ_dy + pow(error_dy,2);
		summ2_dy = summ2_dy + pow(dy,2);

		summ_dx2 = summ_dx2 + pow(error_dx2,2);
		summ2_dx2 = summ2_dx2 + pow(dx2,2);

		summ_dy2 = summ_dy2 + pow(error_dy2,2);
		summ2_dy2 = summ2_dy2 + pow(dy2,2);

		summ_dxdy = summ_dxdy + pow(error_dxdy,2);
		summ2_dxdy = summ2_dxdy + pow(dxdy,2);

		//calculate the maximmun node error
		error_dx = fabs(error_dx);
		if ( error_dx > MaxError_dx ) MaxError_dx = error_dx;

		error_dy = fabs(error_dy);
		if ( error_dy > MaxError_dy ) MaxError_dy = error_dy;

		error_dx2 = fabs(error_dx2);
		if ( error_dx2 > MaxError_dx2 ) MaxError_dx2 = error_dx2;

		error_dy2 = fabs(error_dy2);
		if ( error_dy2 > MaxError_dy2 ) MaxError_dy2 = error_dy2;

		error_dxdy = fabs(error_dxdy);
		if ( error_dxdy > MaxError_dxdy ) MaxError_dxdy = error_dxdy;
	}
}
