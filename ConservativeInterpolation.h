#ifndef __CONSERVATIVEINTERPOLATION_H_INCLUDED__
#define __CONSERVATIVEINTERPOLATION_H_INCLUDED__

//classe edge morre, permanece com a vector
//usa algum namespace pra nao dar conflito

vector<pPoint> edge_edge_Intersection(pEntity edge1, pEntity edge2);

//metodos de vector
double power(pPoint other);
Vector* normal_unit_pos_vector();
double module();




#endif