#include <Interpolation.h>
#include <Edge.h>

Interpolation::Interpolation(double m_mass){
		mass = m_mass;
}

vector<Point> Interpolation::edge_edge_Intersection(Edge *P, Edge *Q){
	bool P1, P2, Q1, Q2 = false;
	double Power_P1 = P->power(Q->getP());
	double Power_P2 = P->power(Q->getQ());
	double Power_Q1 = Q->power(P->getP());
	double Power_Q2 = Q->power(P->getQ());
	int numPoints = 1;
	vector<Point> intersect_points;
	//no maximo dois pontos de intersecao
	Point A, B;
	//intersect_points.push_back();
	//intersect_points.size(); 

	int* aux_x = new int[2];
	int* aux_y = new int[2];

	//number of powers that are zeros
	int count = 0;

	if(abs(Power_P1 - 0) < pow(10,-6)){
		P1= true;
		count++;
	}
	if(abs(Power_P2 - 0) < pow(10,-6)){
		P2 = true;
		count++;
	}
	if(abs(Power_Q1 - 0) < pow(10,-6)){
		Q1 = true;
		count++;
	}
	if(abs(Power_Q2 - 0) < pow(10,-6)){
		Q2 = true;
		count++;
	}

	switch (count){
		case 0:
				if(Power_P1*Power_P2 < 0 && Power_Q1*Power_Q2 < 0){
					A.x = P->getP().x +(Power_Q1/(Power_Q1-Power_Q2))*(P->getX());
					A.y = P->getP().y +(Power_Q1/(Power_Q1-Power_Q2))*(P->getY());
				}

				break;
		case 1:	
				if (P1 && Power_Q1*Power_Q2 < 0)
					A = Q->getP();
				else if (P2 && Power_Q1*Power_Q2 < 0)
					A = Q->getQ();
				else if (Q1 && Power_P1*Power_P2 < 0)
					A = P->getP();
				else if (Q2 &&  Power_P1*Power_P2 < 0) 
					A = P->getQ();
				else 
					numPoints = 0;

				break;
		case 2:	if((P1 & Q1) || (P1 & Q2))
					A = Q->getP();			
				else if ((P2 & Q1) || (P2 & Q2))
					A = Q->getQ();
				else
					numPoints = 0;
				
				break;
		case 4:
				//P->getP()
				//P->getQ()
				//Q->getP();
				//Q->getQ();

				//calculando intervalo de intersecao em x
				aux_x[0] = max(min(P->getP().x, P->getQ().x), min(Q->getP().x, Q->getQ().x));
				aux_x[1] = min(max(P->getP().x, P->getQ().x), max(Q->getP().x, Q->getQ().x));

				printf("\nolha nois aqui: %d %d\n", aux_x[0], aux_x[1]);

				if (aux_x[1] < aux_x[0]) {
				  numPoints = 0;
				  break;
				}
				else{
					//calculando intervalo de intersecao em y
					aux_y[0] = max(min(P->getP().y, P->getQ().y), min(Q->getP().y, Q->getQ().y));
					aux_y[1] = min(max(P->getP().y, P->getQ().y), max(Q->getP().y, Q->getQ().y));
					
					if (aux_y[1] < aux_y[0]) {
				  		numPoints = 0;
				  		break;
					}
					else{
						A.x = aux_x[0];
						A.y = aux_y[0];
						B.x = aux_x[1];
						B.y = aux_y[1];


						if(A.x == B.x && A.y == B.y)
							numPoints = 1;
						else
							numPoints = 2;
					}


				}
				//ve ai os casos
				//te vira
				break;
	}

	if(numPoints == 1){
		intersect_points.push_back(A);
	}
	else if (numPoints == 2){
		intersect_points.push_back(A);
		intersect_points.push_back(B);
	}

	return intersect_points;
}

