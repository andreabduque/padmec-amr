/*
 * AdaptiveRemeshing.cpp
 *
 *  Created on: 26/04/2013
 *      Author: Saulo Menezes
 */
#include "AdaptiveRemeshing.h"
#include "ErrorEstimator.h"
#include "MRE.h"
#include "FileManager.h"
#include "Rebuilder2.h"
#include "Rebuilder3.h"
#include "Gmsh.h"
#include <stdio.h>
#include <stdlib.h>
#include "Gmsh.h"
#include <sstream>
#include <iostream>
#include <set>
#include "GModel.h"
#include "PView.h"
#include "PViewData.h"
#include "PluginManager.h"



AdaptiveRemeshing::AdaptiveRemeshing(int argc, char* argv[]) {
	GmshInitialize(argc, argv);
}

AdaptiveRemeshing::~AdaptiveRemeshing() {
	GmshFinalize();
}

void AdaptiveRemeshing::rodar(ErrorAnalysis  *pErrorAnalysis, pMesh & theMesh){
	#ifdef TRACKING_PROGRAM_STEPS
	cout << "TRACKING_PROGRAM_STEPS: AdaptiveRemeshing\tIN\n";
	#endif
	
	static ofstream FID;
	static bool key = true;
	static int counter = 0;
	
	if(key){
		FID.open("adaptation-monitor.txt");
		key = false;
	}
		
	if (!theMesh){
		throw Exception(__LINE__,__FILE__,"NULL Mesh!");
	}

	ErrorEstimator *E = new ErrorEstimator(theMesh);
	std::set<pEntity> nodesBGMesh;
	pErrorAnalysis->getRefUnrefElementsList(theMesh,E->elementList,nodesBGMesh);
	cout << "                                    : " << M_numFaces(theMesh) << endl;
	cout << "Number of elements to be (un)refined: " << E->elementList.size() << endl;
	system("rm View2D.pos");
	FID << "Adaptation: " << ++counter << endl;
	FID << "Before\n";
	FID << "Num. elements: " << M_numFaces(theMesh) << endl;
	FID << "Num. Vertices: " << M_numVertices(theMesh) << endl;
	FID << "Num. Elements to be removed: " << E->elementList.size() << endl;

	if (!E->elementList.size()){
		throw Exception(__LINE__,__FILE__,"No elements to be (un)refined");
	}
	else{
		theMesh->modifyState(3,2,1);
		theMesh->modifyState(3,1,1);
		theMesh->modifyState(3,0);

		theMesh->modifyState(2,1);
		theMesh->modifyState(2,3);
		theMesh->modifyState(2,0);

		theMesh->modifyState(1,3);
		theMesh->modifyState(1,2);
		theMesh->modifyState(1,0);

		theMesh->modifyState(0,2);
		theMesh->modifyState(0,1);
		theMesh->modifyState(0,3);

		E->LeitorDeLista(pErrorAnalysis,theMesh);
		
		int I=0;
		if (M_numRegions(theMesh)!=0){ // verificando se a malha é 3D ou 2D...
			I=3; // malha 3D
		}
		else{
			I=2; // malha 2D
		}

		cout << "Num tetras: " << M_numRegions(theMesh) << endl;
		cout << "Num vertices: " << M_numVertices(theMesh) << endl;
		cout << "Num faces: " << M_numFaces(theMesh) << endl;

		theMesh->modifyState(3,2,1);
		theMesh->modifyState(3,1,1);
		theMesh->modifyState(3,0);

		theMesh->modifyState(2,1);
		theMesh->modifyState(2,3);
		theMesh->modifyState(2,0);

		theMesh->modifyState(1,3);
		theMesh->modifyState(1,2);
		theMesh->modifyState(1,0);

		theMesh->modifyState(0,2);
		theMesh->modifyState(0,1);
		theMesh->modifyState(0,3);
		Remover *p = new Remover(theMesh);
		int GrtID = p->Iterador(E->ElementsIds, E->elementList,I, theMesh); // PAULO DIZ QUE O GMSH JA DA OS ELEMENTOS DE FRONTEIRA, VERIFICAR ISSO.
		
		if (I==2){
			/*
			 * recebe o set de faces a remover set<pFace> BoundaryFaces e a partir dele preenche os boundaryedges e boundarynodes,
			 * inicialmente recebe todos os faces (nao apenas os de contorno), mas depois das funcoes de remoção só sobram os boundary mesmo, daí o nome do set.
			 */
			p->BoundaryElements2D(theMesh,E->elementList);

			// Pegando o valor da maior aresta ANTES de rodar a funcao Boundaryelements2D novamente, quando ela rodar de novo, esse valor vai ser substituido, por isso tenho que guardá-lo agora.
			
			int numElementsToRemove = (int)E->elementList.size();			
			p->SaveBGMView1(E->elementList,numElementsToRemove); // Salvando o arquivo da View da malha de background.
			

			//   AGORA TENHO QUE CHAMAR AS FUNCOES IDENTIFICADORAS.
			int PSCounter = 0;  // Contador de plane surfaces
			p->removetetras(theMesh); // removendo os tetraedros indicados na lista
			p->removeinternalfaces_2D(theMesh,E->elementList); // removendo as faces indicadas na lista

			theMesh->modifyState(3,2,1);
			theMesh->modifyState(3,1,1);
			theMesh->modifyState(3,0);

			theMesh->modifyState(2,1);
			theMesh->modifyState(2,3);
			theMesh->modifyState(2,0);

			theMesh->modifyState(1,3);
			theMesh->modifyState(1,2);
			theMesh->modifyState(1,0);

			theMesh->modifyState(0,2);
			theMesh->modifyState(0,1);
			theMesh->modifyState(0,3);

			// Esta funcao remove elementos "estranhos" (dentes nos contornos, elementos que estejam causando uniao de lineloops, faces 		soltas no interior do buraco...)
			//p->removestgNodes(theMesh);
			p->removestgNodes(theMesh); // Ativando-a novamente, colocar uma autochamada dentro dela pra ficar mais elegante

			theMesh->modifyState(1,2,0);
			theMesh->modifyState(1,2,1);
			theMesh->modifyState(2,1,0);
			theMesh->modifyState(2,1,1);

			PSCounter=p->Identify_and_Remove_Edges_2D(theMesh);
			p->removeinternaledges(theMesh); // removendo as arestas das faces que foram removidas

			theMesh->modifyState(0,1,0);
			theMesh->modifyState(0,1);

			theMesh->modifyState(0,1,0);
			theMesh->modifyState(0,1,1);
			theMesh->modifyState(0,2,0);
			theMesh->modifyState(0,2,1);


			// Vou deixar as funcoes que removem elementos da fronteira desligadas por enquanto.


			p->SaveBGMView2(numElementsToRemove); // Salvando o arquivo da View da malha de background.
			p->removeinternalnodes(theMesh);

			float GE = p->GreatestEdge;  // ACHO QUE NÃO É PRECISO GUARDAR ESSE VALOR, BASTA FAZER O CÁLCULO DO TAMANHO DA ARESTA PARA CADA PAR DE NÓS
			
			// Agora salvando os arquivos:
			FileManager *fm;

			fm = new FileManager(theMesh, p->ReturnBoundaryEdges(), p->ReturnBoundaryNodes(),PSCounter);
			fm->SaveFileGeometry_2D(GE, theMesh,PSCounter);  // Neste momento o cl dos elementos é o comprimento da maior aresta  - GE
			fm->Clear_Containers();

			p->removeexternaledges(theMesh);
			p->removeexternalnodes(theMesh);

			theMesh->modifyState(3,2,1);
			theMesh->modifyState(3,1,1);
			theMesh->modifyState(3,0);

			theMesh->modifyState(2,1);
			theMesh->modifyState(2,3);
			theMesh->modifyState(2,0);

			theMesh->modifyState(1,3);
			theMesh->modifyState(1,2);
			theMesh->modifyState(1,0);

			theMesh->modifyState(0,2);
			theMesh->modifyState(0,1);
			theMesh->modifyState(0,3);

			// Agora a chamando as funcoes MakeMesh e RefineMesh para recriar a malha com a geometria recém salva em fm->SaveFileGeometry_2D();

			//Rebuilder2 *Q;
			//Q = new Rebuilder2();
			//Q->MakeMesh_2D("Geometria2D", 2);

			
			Rebuilder3 *Q;
			Q = new Rebuilder3();
			Q->MakeMesh_2D("Geometria2D");
			
			
			pMesh theMesh2;
			theMesh2 = MS_newMesh(0);
			
			cout << "Lendo 0.msh" << endl;
			
			M_load(theMesh2,"Geometria2D.msh");
			cout << "Lido o arquivo de malha gerado" << endl;
			cout << "Num edges: " << M_numEdges(theMesh2) << endl;


			theMesh2->modifyState(3,2,1);
			theMesh2->modifyState(3,1,1);
			theMesh2->modifyState(3,0);

			theMesh2->modifyState(2,1);
			theMesh2->modifyState(2,3);
			theMesh2->modifyState(2,0);

			theMesh2->modifyState(1,3);
			theMesh2->modifyState(1,2);
			theMesh2->modifyState(1,0);

			theMesh2->modifyState(0,2);
			theMesh2->modifyState(0,1);
			theMesh2->modifyState(0,3);

			theMesh->modifyState(3,2,1);
			theMesh->modifyState(3,1,1);
			theMesh->modifyState(3,0);

			theMesh->modifyState(2,1);
			theMesh->modifyState(2,3);
			theMesh->modifyState(2,0);

			theMesh->modifyState(1,3);
			theMesh->modifyState(1,2);
			theMesh->modifyState(1,0);

			theMesh->modifyState(0,2);
			theMesh->modifyState(0,1);
			theMesh->modifyState(0,3);

			double GID = (double)GrtID;
			fm->SaveFileMshRenumber(theMesh, "theMesh.msh",0); // Este terceiro parametro é o ponto inicial de contagem dos novos ids
			fm->SaveFileMshRenumber2(theMesh2, "theMesh2.msh",GID);

			Q->MakeMerge_2D("theMesh.msh", "theMesh2.msh","Final_01.msh");

// 			deleteMesh(theMesh);
// 			theMesh = MS_newMesh(0);
// 
// 			M_load(theMesh,"Final_01.msh");
// 
// 			theMesh->modifyState(3,2,1);
// 			theMesh->modifyState(3,1,1);
// 			theMesh->modifyState(3,0);
// 
// 			theMesh->modifyState(2,1);
// 			theMesh->modifyState(2,3);
// 			theMesh->modifyState(2,0);
// 
// 			theMesh->modifyState(1,3);
// 			theMesh->modifyState(1,2);
// 			theMesh->modifyState(1,0);
// 
// 			theMesh->modifyState(0,2);
// 			theMesh->modifyState(0,1);
// 			theMesh->modifyState(0,3);

			E->Clear_Containers(); // Classe Error Estimator
			fm->Clear_Containers(); // Classe filemanager
			p->Clear_Containers(); // Classe Remover
			
// 			FID << "After\n";
// 			FID << "Num. elements: " << M_numFaces(theMesh) << endl;
// 			FID << "Num. Vertices: " << M_numVertices(theMesh) << "\n\n";
			
			delete E; E = 0;
			delete p; p = 0;
			delete Q; Q = 0;
			delete fm; fm = 0;
		}
	}
	//system("rm View2D.pos Geometria2D Geometria2D.msh theMesh.msh theMesh2.msh");
	#ifdef TRACKING_PROGRAM_STEPS
	cout << "TRACKING_PROGRAM_STEPS: AdaptiveRemeshing\tOUT\n";
	#endif
}
