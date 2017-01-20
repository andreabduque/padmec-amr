#include <fstream>
#include <iostream>
#include <string>
#include "FileManagerFunctions.h"
#include <cstring>
#include <sstream>

#include <stdio.h>
#include <stdlib.h>
#include "Gmsh.h"
#include "GModel.h"
//#include "MElement.h"
#include "MVertexGMSH.h"
#include "meshGEdgeGMSH.h"
#include <sstream>
#include <iostream>
#include <set>
#include "SPoint3GMSH.h"


FileManagerFunctions::FileManagerFunctions (){

}

FileManagerFunctions::~FileManagerFunctions (){

}


void FileManagerFunctions::Get_Angular_Coef(set<pEdge> CommonEdges) { 
	set<pEdge>::iterator itped;

	for(itped=CommonEdges.begin(); itped!=CommonEdges.end(); itped++){
		double m; // Coeficiente angular da aresta

		pVertex Node1 = E_vertex(*itped,0);
		pVertex Node2 = E_vertex(*itped,1);

		double xyz1[3];
		double xyz2[3];
		V_coord(Node1, xyz1);
		V_coord(Node2, xyz2);

		double Denominador = xyz2[0]-xyz1[0];

		if (Denominador!=0){
			m = (xyz2[1]-xyz1[1])/(Denominador);
		}
		else{
			m = 1000;
		}


		EN_attachDataDbl (*itped,MD_lookupMeshDataId("Angular_Coef"), m);
	}
}

void FileManagerFunctions::CreateLineLoop(set<pEdge> CommonEdges, int c, string FileName, pMesh theMesh) { 

//void FileManagerFunctions::CreateLineLoop(set<pEdge> CommonEdges, int c, string FileName, pMesh theMesh, map<int, GVertexGMSH*> BoundaryPoints, GModel *a) { // "recebe" um set de Edges e cria uma matriz nx3 O set é o BEdges tantas vezes quantas forem necessárias

	// Preciso que esta funcão retorne as linhas de line loops para uma matriz externa, que servirá de parametro para outra funcao, "CreatePlaneSurface" e esta sim gravará no arquivo as plane surfaces.

	//map <int, GEdgeGMSH*> BoundaryLines;
	string pontobarra = "./";
	ostringstream os;
	os << pontobarra;
	os << FileName;

	std::string filename = os.str();

	ofstream Myfile(filename.c_str(), ios::app);



	set<pEdge>:: iterator tet;
	for(tet=CommonEdges.begin();tet!=CommonEdges.end(); tet++){
		pVertex Vx1 = E_vertex(*tet, 0);
		pVertex Vx2 = E_vertex(*tet, 1);
		Myfile<<"Line("<<EN_id(*tet)<<")"<<" "<<"="<<" "<<"{"<< EN_id(Vx1) <<","<<" "<< EN_id(Vx2) <<"};"<<endl;

		//GVertexGMSH* gv1 = BoundaryPoints.find(EN_id(Vx1))->second;
		//GVertexGMSH* gv2 = BoundaryPoints.find(EN_id(Vx2))->second;
		//BoundaryLines.insert ( std::pair< int, GEdgeGMSH* >(EN_id(*tet), a->addLine(*gv1, *gv2)) );
	}



	set<pEdge>::iterator itEds;
	pVertex Point1;
	pVertex Point2;

	int Dim=CommonEdges.size(); 

	//int Dimat3 = EdgesMap.size();
	int Dimat3 = CommonEdges.size();
	int matriz3[Dimat3][3];

	for (int m=1; m<=Dimat3; m++){
		matriz3[m][1] =0;
		matriz3[m][2] =0;
		matriz3[m][3] =0;
	}



	//cout << "Set CommonEdges CommonEdges.size() " << CommonEdges.size() << endl;
	for ( itEds=CommonEdges.begin() ; itEds != CommonEdges.end(); itEds++ ){

		pEdge Edgex = *itEds;

		Point1 = E_vertex(Edgex, 0);
		Point2 = E_vertex(Edgex, 1);

		int IDPoint1 = EN_id(Point1);
		int IDPoint2 = EN_id(Point2);
		//cout << "ID ARESTA " << EN_id(Edgex) << " ID PONTOS " << EN_id(Point1) << " " << EN_id(Point2) << endl;
	}

	int i=1;
	int matriz1[Dim][3];
	int matriz2[Dim][3];
	for (int m=1; m<=Dim; m++){
		matriz1[m][1] =0;
		matriz1[m][2] =0;
		matriz1[m][3] =0;
		matriz2[m][1] =0;
		matriz2[m][2] =0;
		matriz2[m][3] =0;
	}
	for ( itEds=CommonEdges.begin() ; itEds != CommonEdges.end(); itEds++ ){

		pEdge Edgex = *itEds;

		Point1 = E_vertex(Edgex, 0);
		Point2 = E_vertex(Edgex, 1);

		int IDPoint1 = EN_id(Point1);
		int IDPoint2 = EN_id(Point2);

		int k = EN_id(Edgex);
		matriz1 [i][1]=k;     // Criando a matriz Atribuindo ids às arestas usando o indice do conjunto map em vez de set
		matriz1 [i][2]=IDPoint1;
		matriz1 [i][3]=IDPoint2;
		i++;
		//cout << i <<" "<< IDPoint1 <<" "<< IDPoint2 << endl;

	}



	int b=1;
	for (int u=1; u<=3; u++){
		matriz2[1][u]=matriz1[1][u];
	}
	for(int t = 1; t<=3;t++){//zerando a 1a linha da matriz1
		matriz1[1][t]=0;
	}



	for(int f=1;f<Dim;f++){
		b=1;
		for(int h = 2; h<=3; h++){
			for (int k=2;k<=3;k++){
				if (b==1){
					for (int g=1;g<=Dim;g++){
						if (matriz2[f][h]==matriz1[g][k]){ //compara as colunas
							//copia a linha e zera a da primeira matriz
							for (int u=1; u<=3; u++){
								matriz2[f+1][u]=matriz1[g][u];

							}
							for(int r = 1; r<=3;r++){//zerando a linha correspondente da matriz1
								matriz1[g][r]=0;
							}
							b=0;
							break;
						}
					}
				}
			}
		}
	}


	for (int n=1; n<Dim;n++){
		if(matriz2[n][2]==matriz2[n+1][3]){
			matriz2[n+1][1]=-matriz2[n+1][1];
		}
		if(matriz2[n][3]==matriz2[n+1][3]){
			matriz2[n+1][1]=-matriz2[n+1][1];
		}
		if(n==Dim-1){
			if(matriz2[n+1][2]==matriz2[1][3]){
				matriz2[1][1]=-matriz2[1][1];
			}
			if(matriz2[n+1][3]==matriz2[1][3]){
				matriz2[1][1]=-matriz2[1][1];
			}
		}
	}

	//return matriz2;




	int LL = c;

	Myfile <<"Line Loop("<< LL <<") = {" << matriz2[1][1];

	for ( int k = 2 ; k<=Dim; k++ ){
		Myfile << ", " << matriz2[k][1]; // a matriz2 de cada vez guarda os ids e pontos das arestas já na ordem do line loop.
	}
	Myfile <<"};" << endl;

	// Agora criar matriz de duas colunas para armazenar PS e LL.

	int psreader = 0;

	set<pEdge>::iterator itCommonEdges;
	itCommonEdges = CommonEdges.begin();
	pEdge firstEdge = *itCommonEdges;

	EN_getDataInt (firstEdge, MD_lookupMeshDataId("PlaneSurface"), &psreader);

	//cout << "O line loop " << LL << " pertence ao PS " << psreader << endl;

	//Myfile << "Plane Surface(" << PS << ") = {" << PS << "};" << endl;

	Myfile.close();

}


	
void FileManagerFunctions::CreateLineLoop_3D(set<pEdge> CommonEdges, int c, string FileName, pMesh theMesh) { // "recebe" um set de Edges e cria uma matriz nx3 O set é o BEdges tantas vezes quantas forem necessárias

	// Preciso que esta funcão retorne as linhas de line loops para uma matriz externa, que servirá de parametro para outra funcao, "CreatePlaneSurface" e esta sim gravará no arquivo as plane surfaces.


	string pontobarra = "./";
	ostringstream os;
	os << pontobarra;
	os << FileName;

	std::string filename = os.str();

	ofstream Myfile(filename.c_str(), ios::app);



	set<pEdge>:: iterator tet;
	for(tet=CommonEdges.begin();tet!=CommonEdges.end(); tet++){
		pVertex Vx1 = E_vertex(*tet, 0);
		pVertex Vx2 = E_vertex(*tet, 1);
		//Myfile<<"Line("<<EN_id(*tet)<<")"<<" "<<"="<<" "<<"{"<< EN_id(Vx1) <<","<<" "<< EN_id(Vx2) <<"};"<<endl;
	}



	set<pEdge>::iterator itEds;
	pVertex Point1;
	pVertex Point2;

	int Dim=CommonEdges.size(); 

	//int Dimat3 = EdgesMap.size();
	int Dimat3 = CommonEdges.size();
	int matriz3[Dimat3][3];

	for (int m=1; m<=Dimat3; m++){
		matriz3[m][1] =0;
		matriz3[m][2] =0;
		matriz3[m][3] =0;
	}



	//cout << "Set CommonEdges CommonEdges.size() " << CommonEdges.size() << endl;
	for ( itEds=CommonEdges.begin() ; itEds != CommonEdges.end(); itEds++ ){

		pEdge Edgex = *itEds;

		Point1 = E_vertex(Edgex, 0);
		Point2 = E_vertex(Edgex, 1);

		int IDPoint1 = EN_id(Point1);
		int IDPoint2 = EN_id(Point2);
		//cout << "ID ARESTA " << EN_id(Edgex) << " ID PONTOS " << EN_id(Point1) << " " << EN_id(Point2) << endl;
	}

	int i=1;
	int matriz1[Dim][3];
	int matriz2[Dim][3];
	for (int m=1; m<=Dim; m++){
		matriz1[m][1] =0;
		matriz1[m][2] =0;
		matriz1[m][3] =0;
		matriz2[m][1] =0;
		matriz2[m][2] =0;
		matriz2[m][3] =0;
	}
	for ( itEds=CommonEdges.begin() ; itEds != CommonEdges.end(); itEds++ ){

		pEdge Edgex = *itEds;

		Point1 = E_vertex(Edgex, 0);
		Point2 = E_vertex(Edgex, 1);

		int IDPoint1 = EN_id(Point1);
		int IDPoint2 = EN_id(Point2);

		int k = EN_id(Edgex);
		matriz1 [i][1]=k;     // Criando a matriz Atribuindo ids às arestas usando o indice do conjunto map em vez de set
		matriz1 [i][2]=IDPoint1;
		matriz1 [i][3]=IDPoint2;
		i++;
		//cout << i <<" "<< IDPoint1 <<" "<< IDPoint2 << endl;

	}



	int b=1;
	for (int u=1; u<=3; u++){
		matriz2[1][u]=matriz1[1][u];
	}
	for(int t = 1; t<=3;t++){//zerando a 1a linha da matriz1
		matriz1[1][t]=0;
	}



	for(int f=1;f<Dim;f++){
		b=1;
		for(int h = 2; h<=3; h++){
			for (int k=2;k<=3;k++){
				if (b==1){
					for (int g=1;g<=Dim;g++){
						if (matriz2[f][h]==matriz1[g][k]){ //compara as colunas
							//copia a linha e zera a da primeira matriz
							for (int u=1; u<=3; u++){
								matriz2[f+1][u]=matriz1[g][u];

							}
							for(int r = 1; r<=3;r++){//zerando a linha correspondente da matriz1
								matriz1[g][r]=0;
							}
							b=0;
							break;
						}
					}
				}
			}
		}
	}


	for (int n=1; n<Dim;n++){
		if(matriz2[n][2]==matriz2[n+1][3]){
			matriz2[n+1][1]=-matriz2[n+1][1];
		}
		if(matriz2[n][3]==matriz2[n+1][3]){
			matriz2[n+1][1]=-matriz2[n+1][1];
		}
		if(n==Dim-1){
			if(matriz2[n+1][2]==matriz2[1][3]){
				matriz2[1][1]=-matriz2[1][1];
			}
			if(matriz2[n+1][3]==matriz2[1][3]){
				matriz2[1][1]=-matriz2[1][1];
			}
		}
	}

	//return matriz2;




	int LL = c;

	Myfile <<"Line Loop("<< LL <<") = {" << matriz2[1][1];

	for ( int k = 2 ; k<=Dim; k++ ){
		Myfile << ", " << matriz2[k][1]; // a matriz2 de cada vez guarda os ids e pontos das arestas já na ordem do line loop.
	}
	Myfile <<"};" << endl;

	// Agora criar matriz de duas colunas para armazenar PS e LL.

	int psreader = 0;

	set<pEdge>::iterator itCommonEdges;
	itCommonEdges = CommonEdges.begin();
	pEdge firstEdge = *itCommonEdges;

	EN_getDataInt (firstEdge, MD_lookupMeshDataId("PlaneSurface"), &psreader);

	//cout << "O line loop " << LL << " pertence ao PS " << psreader << endl;

	//Myfile << "Plane Surface(" << PS << ") = {" << PS << "};" << endl;

	Myfile.close();

}




set<pEdge> FileManagerFunctions::Separateset(set<pEdge> PSEdges, int PS){ // pega uma aresta e procura quais das outras arestas tem um ponto em comum com ela, separa UM set e deixa o resto das arestas no set original.

	set <pEdge> CommonEdges;

	CommonEdges.clear();

	//cout << "	separateset recebeu " << PSEdges.size() << " de PSEdges" << endl;

	pEdge AntEdge = *PSEdges.begin();  // Pegando uma aresta... Agora vai ver quais sao do mesmo lineloop, ou seja, quais tem pontos em comum...

	CommonEdges.insert(AntEdge);

	//cout << "Segmentation Fault 1" << endl;
	pVertex AntNode = E_vertex(AntEdge,1);

	pVertex firstVertex = AntNode;

	int c = 0;
	//cout << "Onde esta o laco 1" << endl;
	int g=1;

	int size = 0;

	while (g==1){

		//g=0;
		//cout << "Onde esta o laco 2" << endl;

		//	cout << "Segmentation Fault 2" << endl;
		pVertex VTX1 = E_vertex(AntEdge,0);
		//	cout << "Segmentation Fault 3" << endl;
		pVertex VTX2 = E_vertex(AntEdge,1);

		//	cout << "Aresta de vertices: " << EN_id(VTX1) <<" "<< EN_id(VTX2) << endl;

		if (VTX1!=AntNode){
			AntNode=VTX1;
		}
		else 	{
			AntNode=VTX2;
		}



		int edges = V_numEdges(AntNode);

		//	cout << "analisando o ponto " << EN_id(AntNode) << endl;
		for(int y=0;y<edges;y++){

			//		cout << "Segmentation Fault 4" << endl;
			pEdge edgex = V_edge(AntNode,y);

			int PlSfR=0;
			EN_getDataInt (edgex, MD_lookupMeshDataId("PlaneSurface"), &PlSfR);

			//		cout << "PlSfr = " << PlSfR << endl;


			if(PlSfR!=0){

				if (edgex==AntEdge){
					//				cout << "Essa era a aresta anterior" << endl;
				}

				if (edgex!=AntEdge){
					//				cout << "encontrada a outra aresta e inserida ao set, CommonEdges.size() = " << CommonEdges.size() << endl;
					if (PSEdges.count(edgex)==1){ // Só acrescenta se a aresta existir no PSEdges... - COLOQUEI ESSE IF MAS NAO ESTOU GOSTANDO, COLOQUEI A TITULO DE EMERGENCIA MESMO... PARA EVITAR QUE A FUNCAO PEGUE ARESTAS QUE NAO TÊM NADA A VER COM O LINE LOOP...
						y=edges;
						CommonEdges.insert(edgex);
						AntEdge=edgex;
						if(size==CommonEdges.size()){
							g=0;
						}
						size=CommonEdges.size();

					}
				}
			}
			//cout << "Onde esta o laco 5" << endl;
			if (y==edges) {
				//			cout << "terminou de analisar as arestas desse ponto, vamos ao proximo" << endl;
			}

		}


	}

	//cout << "	Separateset terminou de separar um set de arestas comuns com " << CommonEdges.size() <<" arestas" << endl;

	return CommonEdges;

}

string FileManagerFunctions::Create_Key_Point(pVertex Vt1){ 

	double coords1[3];

	V_coord(Vt1, coords1);

	string ponto = "."; 	//Criando a chave para inserir no novo map que servirá para atualizar BoundVertex

	ostringstream oss;
	oss << coords1[0];
	oss << ponto;
	oss << coords1[1];
	string XY = oss.str();

	return XY;
}



//set<pEdge> FileManagerFunctions::Merge_Edges(set<pEdge> CommonEdges, pMesh theMesh, map<int,pEdge> EdgesMap, int c, string FileName, GModel* a){

set<pEdge> FileManagerFunctions::Merge_Edges(set<pEdge> CommonEdges, pMesh theMesh, map<int,pEdge> EdgesMap, int c, string FileName){
		
	string pontobarra = "./";
	ostringstream os;
	os << pontobarra;
	os << FileName;

	std::string filename = os.str();

	ofstream Myfile(filename.c_str(), ios::app);


	set<pEdge> GeomBdryEdges;
	set<pEdge> newEdges;
	set<pEdge> SameEdge;
	set<pEdge> SameCoefEdges;

	set<pEdge>::iterator edgeiterator;
	set<pEdge>::iterator itEds;


	for(itEds=CommonEdges.begin(); itEds!=CommonEdges.end(); itEds++){
		if (E_numFaces(*itEds)==0){
			GeomBdryEdges.insert(*itEds);
		}
	}

	int y = 0;

	// iterando pelas arestas do contorno da geometria para separar o primeiro grupo de arestas com mesmo coeficiente angular

	while(GeomBdryEdges.size()>0){

		pEdge Edge1 = *GeomBdryEdges.begin();  //	PEGA UMA ARESTA

		double m1 = 0;

		EN_getDataDbl (Edge1, MD_lookupMeshDataId("Angular_Coef"), &m1); 
		SameCoefEdges.insert(Edge1);


		for(itEds=GeomBdryEdges.begin(); itEds!=GeomBdryEdges.end(); itEds++){  // E PEGA AS OUTRAS Q TIVEREM MESMO COEFICIENTE ANGULAR
			pEdge Edge2 = *itEds;

			double m2 = 0;
			EN_getDataDbl (Edge2, MD_lookupMeshDataId("Angular_Coef"), &m2);

			if (m1==m2){
				SameCoefEdges.insert(Edge2);
			}
		}

		while (SameCoefEdges.size()>0){

			// iniciando a funcao merge_edge
			pEdge Old_Edg1 = *SameCoefEdges.begin();
			pEdge Old_Edg2 = *SameCoefEdges.begin();

			pVertex Vert1 = E_vertex(Old_Edg1,0);
			pVertex Vert2 = E_vertex(Old_Edg2,1);

			SameEdge.insert(Old_Edg1);

			int v1 = 0;
			int v2 = 0;
			int V = 1;

			//agora está identificando as arestas que serao unidas	
			while (V==1){
				v1 =0;
				v2 =0;

				int Numb_Edges1 = V_numEdges(Vert1);
				int Numb_Edges2 = V_numEdges(Vert2);

				for(int i = 0 ; i < Numb_Edges1 ; i++){

					pEdge Edg1 = V_edge(Vert1, i);

					if (Edg1!=Old_Edg1 && SameCoefEdges.count(Edg1)==1){
						Vert1 = E_otherVertex(Edg1, Vert1);
						Old_Edg1 = Edg1;
						SameEdge.insert(Edg1);
						v1 = 1; // Preciso deletar da malha as arestas que serao mergidas...
						//cout << "v1 = " << v1 << endl;
					}
				}
				for(int i = 0 ; i < Numb_Edges2 ; i++){
					pEdge Edg2 = V_edge(Vert2, i);
					if (Edg2!=Old_Edg2 && SameCoefEdges.count(Edg2)==1){
						Vert2 = E_otherVertex(Edg2, Vert2);
						Old_Edg2 = Edg2;
						SameEdge.insert(Edg2);
						v2 = 1;
						//cout << "v2 = " << v2 << endl;
					}
				}
				//				cout << "v1 e v2 " << v1 << " " << v2 << endl;
				//				cout << "GeomBdryEdges.size() " << GeomBdryEdges.size() << endl;

				if (v1==0 && v2==0){
					V=0;
					// criando a nova aresta aqui dentro deste if e inseri-la no set newEdges e no set commonedges
					pGEntity ent = E_whatIn(Old_Edg2);
					pEdge edg = M_createE(theMesh, Vert1, Vert2, ent); 
					//newEdges.insert(edg);
					CommonEdges.insert(edg);
					int EMS = EdgesMap.size();
					//cout << "EdgesMap ta zerado?? EdgesMap.size(): " << EdgesMap.size() << endl;
					int chave = EMS+c;
					EdgesMap.insert (pair<int,pEdge>(chave,edg));

					Myfile<<"Line("<<chave<<")"<<" "<<"="<<" "<<"{"<< EN_id(Vert1) <<","<<" "<< EN_id(Vert2) <<"};"<<endl;

					//cout<< "Nova aresta criada" << endl;
				}
			}

			//terminou de identificar

			for(itEds=SameEdge.begin(); itEds!=SameEdge.end(); itEds++){
				GeomBdryEdges.erase(*itEds);
				SameCoefEdges.erase(*itEds);
				theMesh->DEL(*itEds); // Coloquei agora...
				CommonEdges.erase(*itEds);
			}
			SameEdge.clear();
		}
	}

	GeomBdryEdges.clear();
	SameCoefEdges.clear();
	return CommonEdges;


}













