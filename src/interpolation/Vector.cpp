#include <interpolation.h>

// Vetor PQ
Vector::Vector (pPoint P, pPoint Q){
	x = P_x(Q) - P_x(P);
	y = P_y(Q) - P_y(P);
}

Vector::Vector (double _x, double _y){
	x = _x;
	y = _y;
}


Vector* Vector::normal_unit_pos_vector(){

 	Vector* normal = new Vector(y/module(), -x/module());
 	return normal;
}

double Vector::module(){
	double mod = sqrt(x*x + y*y);
	return mod;
}

double Vector::getX(){
	
	return x;
}

double Vector::getY(){
	
	return y;
}

