#include <Edge.h>

Edge::Edge (Point _P, Point _Q) : Vector (_P, _Q){
	P = _P;
	Q = _Q;
	normal_to_edge = normal_unit_pos_vector();	
}
double Edge::power(Point other){	
	double proj = (P.x - other.x)*normal_to_edge->getX() + (P.y - other.y)*normal_to_edge->getY();
	return proj;
}

//pra debugar
double Edge::showMe(Point other){
	double power_value = power(other);	
	
	printf("\nPoint P: (%lf, %lf)\n", P.x, P.y);
	printf("Point Q: (%lf, %lf)\n", Q.x, Q.y);
	printf("Normal unit Vector: (%lf, %lf)\n", normal_to_edge->getX(), normal_to_edge->getY());
	printf("Power to respect to (%lf,%lf):  %lf\n", other.x, other.y, power_value);

	return power_value;

}

Point Edge::getP(){
	return P;
}

Point Edge::getQ(){
	return Q;
}