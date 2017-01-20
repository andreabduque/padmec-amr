
//caso queira usar delaunay

/*
void triangulate_cloud(list<pPoint> cloud_points){
	unsigned int num_points = cloud_points.size();
	del_point2d_t* vec_points = new del_point2d_t[num_points];
	del_point2d_t p0, p1, p2;

	int count = 0;
	for (std::list<pPoint>::iterator it = cloud_points.begin(); it != cloud_points.end(); it++){
		vec_points[count].x = P_x(*it);
		vec_points[count].y = P_y(*it);
		count++;
	}

	printf("-----------printando ai os triangulozinho do delaunay----------\n");
	
	delaunay2d_t* polygon = delaunay2d_from(vec_points, num_points);
	tri_delaunay2d_t* triangulation =tri_delaunay2d_from(polygon);

	printf("connectivity list %d\n", triangulation->num_triangles);
	
	printf("[");
	for(unsigned int i = 0; i < 3*triangulation->num_triangles; i = i+3){
		p0 = triangulation->points[triangulation->tris[i]];
		p1 = triangulation->points[triangulation->tris[i + 1]];
		p2 = triangulation->points[triangulation->tris[i + 2]];

		printf("%d %d %d;\n", triangulation->tris[i] + 1, triangulation->tris[i + 1] +1, triangulation->tris[i + 2]+1);		

	} printf("]");
	printf("coord x\n");
	for(unsigned int i = 0; i < triangulation->num_points; i++){
		p0 = triangulation->points[i];	
		
		printf("%lf\n", p0.x);
	}
	printf("coord y\n");
	for(unsigned int i = 0; i < triangulation->num_points; i++){
		p0 = triangulation->points[i];
		printf("%lf\n", p0.y);
	}

	printf("----------fim dos triangulozinho do delaunay----------\n");
	delaunay2d_release(polygon);
	tri_delaunay2d_release(triangulation);

}
*/
/*
list<pPoint> change_pointId(list<pPoint> cloud_points){
	unsigned int count = 1;
	pPoint A;
	list<pPoint> new_cloud;
	
	for (std::list<pPoint>::iterator it = cloud_points.begin(); it != cloud_points.end(); it++){
		A = P_new();
		P_setPos(A, P_x(*it), P_y(*it), 0);
		//printf("id da parada %d", EN_id(A));
		//count++;
		new_cloud.push_back(A);
	}

	return new_cloud;
} */
/*
void debug_intersect(list<pPoint> cloud_points, list< cTriangle > new_mesh){
	printf("printando coords  x da nuvem de pontos\n");
	for (std::list<pPoint>::iterator it = cloud_points.begin(); it != cloud_points.end(); it++){
		printf("%lf\n", P_x(*it));
	}
	printf("printando coords  y da nuvem de pontos\n");
	for (std::list<pPoint>::iterator it = cloud_points.begin(); it != cloud_points.end(); it++){
		printf("%lf\n", P_y(*it));
	}

	printf("printando lista de triangulacoes\n");
	for (std::list<cTriangle>::iterator it = new_mesh.begin(); it != new_mesh.end(); it++){
		printf("%d %d %d\n", it->p0, it->p1, it->p2);
	}
} */