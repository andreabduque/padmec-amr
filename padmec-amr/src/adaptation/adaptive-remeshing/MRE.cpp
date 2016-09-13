#include "MRE.h"
#include <string>
#include <cstring>
#include <sstream>
#include <iostream>
#include <math.h>
#include <cmath>
#include <fstream>


// git add --all
// git commit -m " comentario da atualização "
// git push -u origin


Remover::Remover (pMesh theMesh){
	MtxTriangles = 0;
}

Remover::~Remover () {
}

void Remover::Clear_Containers(){
	PlaneSurfaceEdges.clear();
	GeometryBoundNodes.clear();
	GeomBoundNodes.clear();
	GeomBoundEdges.clear();
	MeshBoundEdges.clear();
	MeshBoundNodes.clear();
	AdaptNodes.clear();
	AdaptEdges.clear();
	BoundaryNodes.clear();
	BoundaryEdges.clear();
	RemovefromBoundEdges.clear(); 
	RemovefromBoundNodes.clear();
	RemovefromBoundFaces.clear();
	//BoundaryFaces.clear(); 
	GeomBoundFaces.clear();
	Tetras.clear();
	argi.clear();
	CommonVertex.clear();
	triangles.clear();
	VertexClsAux.clear();
	VertexCls.clear();
 	if (MtxTriangles){
 		delete[] MtxTriangles; MtxTriangles = 0;
 	}
 	GreatestEdge = 0;
}

void Remover::CriadordeLista3D(pMesh theMesh){ // funcao pra criar uma lista de elementos 3D  <----- Lembrar de retirar essa funcao depois

	int iddotetra = 40;

	cout << "Primeiro tetra: " << iddotetra << endl;


	////////// Atualizei ate aqui ///////////

	set <pRegion> Tetrasvizinhos;
	set <pFace> facesvizinhas;

	//set <pFace> facesremove;


	cout << "Primeiros vizinhos: " << endl;

	RIter rit = M_regionIter(theMesh);
	while (pRegion tity = RIter_next(rit)){

		if (EN_id(tity)==iddotetra){

			facesvizinhas.insert(R_face(tity, 0));
			facesvizinhas.insert(R_face(tity, 1));
			facesvizinhas.insert(R_face(tity, 2));
			facesvizinhas.insert(R_face(tity, 3));

			set<pFace>::iterator itfc;
			for ( itfc=facesvizinhas.begin() ; itfc != facesvizinhas.end(); itfc++ ){
				pFace FEC = *itfc;
				int numReg = F_numRegions(FEC);
				for(int t=0;t<numReg;t++){
					Tetrasvizinhos.insert(F_region(FEC, t));
				}
			}
		}
	}

	RIter_delete(rit);


	set<pRegion>::iterator itfcs;
	set<pFace>::iterator itfc;

	for (int v=2;v<=6;v++){
		cout << v << " vizinhos: " << endl;
		cout << "Tetrasvizinhos.size() " << Tetrasvizinhos.size() << endl;
		for ( itfcs=Tetrasvizinhos.begin() ; itfcs != Tetrasvizinhos.end(); itfcs++ ){
			if (v==9){
				cout << EN_id(*itfcs) << endl;
			}
			pRegion FEC = *itfcs;
			facesvizinhas.insert(R_face(FEC, 0));
			facesvizinhas.insert(R_face(FEC, 1));
			facesvizinhas.insert(R_face(FEC, 2));
			facesvizinhas.insert(R_face(FEC, 3));

		}

		cout << endl;
		cout << "facesvizinhas.size(): " << facesvizinhas.size() << endl;

		for ( itfc=facesvizinhas.begin() ; itfc != facesvizinhas.end(); itfc++ ){
			pFace FEC = *itfc;
			int numReg = F_numRegions(FEC);
			for(int t=0;t<numReg;t++){
				Tetrasvizinhos.insert(F_region(FEC, t));
			}
		}
		cout << endl;
	}

	cout << "Terminou funcao Criador de lista 3D" << endl;

	Tetrasvizinhos.set::~set();
	facesvizinhas.set::~set();
}

// funcao pra criar uma lista de elementos 2D  <----- Lembrar de retirar essa funcao depois
void Remover::CriadordeLista(pMesh theMesh){
	int iddaface = 754;
	cout << "Primeira face: " << iddaface << endl;
	set <pEdge> Edgesvizinhas;
	set <pFace> facesvizinhas;
	set <pFace> facesremove;
	cout << "Primeiros vizinhos: " << endl;

	FIter facit = M_faceIter(theMesh);
	while (pFace Fac = FIter_next(facit)){
		if (EN_id(Fac)==iddaface){

			Edgesvizinhas.insert(F_edge(Fac, 0));
			Edgesvizinhas.insert(F_edge(Fac, 1));
			Edgesvizinhas.insert(F_edge(Fac, 2));
			set<pEdge>::iterator itfc;
			for ( itfc=Edgesvizinhas.begin() ; itfc != Edgesvizinhas.end(); itfc++ ){
				pEdge FEC = *itfc;
				facesvizinhas.insert(E_face(FEC, 0));
				facesvizinhas.insert(E_face(FEC, 1));

				pFace face0 = E_face(FEC, 0);
				pFace face1 = E_face(FEC, 1);

				cout << EN_id(face0) << endl;
				cout << EN_id(face1) << endl;

			}

		}
	}

	FIter_delete(facit);

	set<pFace>::iterator itfcs;
	set<pEdge>::iterator itfc;

	for (int v=2;v<=8;v++){
		for ( itfcs=facesvizinhas.begin() ; itfcs != facesvizinhas.end(); itfcs++ ){
			pFace FEC = *itfcs;
			Edgesvizinhas.insert(F_edge(FEC, 0));
			Edgesvizinhas.insert(F_edge(FEC, 1));
			Edgesvizinhas.insert(F_edge(FEC, 2));
		}

		cout << endl;

		cout << v << " vizinhos: " << endl;

		for ( itfc=Edgesvizinhas.begin() ; itfc != Edgesvizinhas.end(); itfc++ ){
			pEdge FEC = *itfc;
			facesvizinhas.insert(E_face(FEC, 0));
			facesvizinhas.insert(E_face(FEC, 1));
		}
		if (v==3){
			for ( itfcs=facesvizinhas.begin() ; itfcs != facesvizinhas.end(); itfcs++ ){
				pFace RETG = *itfcs;
				facesremove.insert(*itfcs);
				cout << EN_id(RETG) << endl;

			}	
		}

		cout << endl;
	}
}



void Remover::SaveBGMView1(const std::list<pEntity>& elementList, int size){
	MtxTriangles = new eleMatriz[size];
	std::list<pEntity>::const_iterator itFaces;
	double xyz1[3];
	double xyz2[3];
	double xyz3[3];
	
	int i = 0;
	for(itFaces=elementList.begin(); (itFaces!=elementList.end()); itFaces++){
		pFace face = *itFaces;
// 		pVertex Point1 = F_vertex(*itFaces, 0);
// 		pVertex Point2 = F_vertex(*itFaces, 1);
// 		pVertex Point3 = F_vertex(*itFaces, 2);

		double height;
		EN_getDataDbl((pVertex)face->get(0,0),MD_lookupMeshDataId( "elem_height" ),&height);
		MtxTriangles[i].Cl1 = (float)height;
		EN_getDataDbl((pVertex)face->get(0,1),MD_lookupMeshDataId( "elem_height" ),&height);
		MtxTriangles[i].Cl2 = (float)height;
		EN_getDataDbl((pVertex)face->get(0,2),MD_lookupMeshDataId( "elem_height" ),&height);
		MtxTriangles[i].Cl3 = (float)height;

		V_coord((pVertex)face->get(0,0), xyz1);
		V_coord((pVertex)face->get(0,1), xyz2);
		V_coord((pVertex)face->get(0,2), xyz3);

		MtxTriangles[i].x1=xyz1[0];
		MtxTriangles[i].y1=xyz1[1];
		MtxTriangles[i].z1=xyz1[2];

		MtxTriangles[i].x2=xyz2[0];
		MtxTriangles[i].y2=xyz2[1];
		MtxTriangles[i].z2=xyz2[2];

		MtxTriangles[i].x3=xyz3[0];
		MtxTriangles[i].y3=xyz3[1];
		MtxTriangles[i].z3=xyz3[2];

		MtxTriangles[i].Pt1 = (pVertex)face->get(0,0);
		MtxTriangles[i].Pt2 = (pVertex)face->get(0,1);
		MtxTriangles[i].Pt3 = (pVertex)face->get(0,2);
		i++;
	}
	cout << "Concluido SaveBGMView1" << endl;
}

void Remover::SaveBGMView2(int size) {
	cout << "Iniciou SaveBGMView2" << endl;
	system("rm View2D.pos"); 
	cout << endl;
	cout << "Salvando View 2D - View2D.pos" << endl;
	ofstream Myfile("./View2D.pos", ios::app);
	if (!Myfile) {
		Myfile.open("./View2D.pos");
	}
	Myfile << endl;
	// inicio da rotina de gravacao da "View" que serve de background mesh
	Myfile << endl;
	Myfile << "View \"0\" {" << endl;  // o zero desta linha deve estar entre aspas...  <---- Alerta!!
	set<pVertex>::iterator itnodes;
	float GRTEDG = 0;
	float EdgeLength = 0;
	int NEDS = 0;

	for ( itnodes=BoundaryNodes.begin() ; itnodes != BoundaryNodes.end(); itnodes++ ){
		double xyz[3];
		V_coord(*itnodes, xyz);
		// Colcoar rotina para calcular para cada ponto o seu GE

		GRTEDG = 0;
		NEDS = V_numEdges(*itnodes);
		for (int h=0; h<NEDS; h++){
			pEdge Edge = V_edge(*itnodes,h);
			pVertex vertex1 = *itnodes;
			pVertex vertex2 = E_vertex(Edge,1); // tem uma funcao FMDB que pega o outro vertex
			
			if (vertex2 == vertex1){
				vertex2 = E_vertex(Edge,0);
			}
					
			double xyz1[3];
			double xyz2[3];

			V_coord(vertex1, xyz1);
			V_coord(vertex2, xyz2);
		
			EdgeLength = sqrt ((xyz1[0]-xyz2[0])*(xyz1[0]-xyz2[0]) + (xyz1[1]-xyz2[1])*(xyz1[1]-xyz2[1]) + (xyz1[2]-xyz2[2])*(xyz1[2]-xyz2[2]) );
			
			if (EdgeLength > GRTEDG){
				GRTEDG = EdgeLength;
			}
		}

		EN_attachDataDbl (*itnodes, MD_lookupMeshDataId("GRTEDG"), GRTEDG);
		Myfile<<"SP("<< xyz[0] <<","<<" "<<xyz[1]<<","<<" "<<xyz[2]<<")"<<"{"<<GRTEDG<<"};"<<endl;
		if (GreatestEdge<GRTEDG){
			GreatestEdge = GRTEDG;
		}
	}

	set<pEdge>::iterator itedges;
	for(itedges=BoundaryEdges.begin();itedges!=BoundaryEdges.end();itedges++){ // lendo o set e gravando no arquivo os dados.

		pVertex Point1 = E_vertex(*itedges, 0);
		pVertex Point2 = E_vertex(*itedges, 1);

		double xyz1[3];
		double xyz2[3];

		V_coord(Point1, xyz1);
		V_coord(Point2, xyz2);

		//SL(-0.5,0.5,0,0.8,0.5,0){0,0};
		Myfile << "SL(" << xyz1[0] << ", " << xyz1[1] << ", " << xyz1[2] << ", " << xyz2[0] << ", " << xyz2[1] << ", " << xyz2[2] << "){0,0};" << endl;

	}

	for(int c=0;c<size;c++){
		int g;
		int h;
		double Clx;
		h=V_numEdges(MtxTriangles[c].Pt1);  // Se o ponto tiver mais de duas arestas vizinhas entao nao é do contorno externo da geometria esta linha verifica isso
		g=BoundaryNodes.count(MtxTriangles[c].Pt1); // mas pode ser do contorno interno da geometria e esta verifica isso
		if (g==1 && h>2){
			
			EN_getDataDbl(MtxTriangles[c].Pt1,MD_lookupMeshDataId( "GRTEDG" ),&Clx);
			MtxTriangles[c].Cl1 = Clx;
		}
		h=V_numEdges(MtxTriangles[c].Pt2);
		g=BoundaryNodes.count(MtxTriangles[c].Pt2);
		if (g==1 && h>2){
			EN_getDataDbl(MtxTriangles[c].Pt2,MD_lookupMeshDataId( "GRTEDG" ),&Clx);
			MtxTriangles[c].Cl2 = Clx;
		}
		h=V_numEdges(MtxTriangles[c].Pt3);
		g=BoundaryNodes.count(MtxTriangles[c].Pt3);
		if (g==1 && h>2){
			EN_getDataDbl(MtxTriangles[c].Pt3,MD_lookupMeshDataId( "GRTEDG" ),&Clx);
			MtxTriangles[c].Cl3 = Clx;
		}

		Myfile << "ST(" << MtxTriangles[c].x1 << ", " << MtxTriangles[c].y1 << ", " << MtxTriangles[c].z1 << ", " 
		       << MtxTriangles[c].x2 << ", " << MtxTriangles[c].y2 << ", " << MtxTriangles[c].z2 << ", " 
			   << MtxTriangles[c].x3 << ", " << MtxTriangles[c].y3 << ", " << MtxTriangles[c].z3 << "){" 
			   << MtxTriangles[c].Cl1 << ", " << MtxTriangles[c].Cl2 << ", " << MtxTriangles[c].Cl3 << "};" << endl;
	}

	Myfile << "};" << endl;
	Myfile << endl;
	Myfile << "Background Mesh View[0];" << endl;


	Myfile.close(); // fim da gravação da View

	cout << "Concluido SaveBGMView2" << endl;
	cout << "Fim da gravação do arquivo View2D.pos" << endl;
	delete[] MtxTriangles; MtxTriangles = 0;
}


int Remover::Iterador(set <int> ListofElements, const list<pEntity>& BoundaryFaces, int I, pMesh theMesh){ // Esta funcao percorre toda a malha uma vez e levanta as informacoes necessárias

	int GreaterID=0;
	cout << "Entrou no iterador - ListofElements possui: " << ListofElements.size() << endl;

	//Considerando que seja uma malha 3D
// 	if(I==3){
// 		set<int>::iterator itElements;
// 		for ( itElements=ListofElements.begin() ; itElements != ListofElements.end(); itElements++ ){
// 
// 			RIter rit = M_regionIter(theMesh);     // O QUE???? UM ITERADOR DENTRO DE UM LOOOOOPPP????!!
// 
// 			while (pRegion reg = RIter_next(rit)){
// 
// 				if (EN_id(reg)==*itElements){
// 					Tetras.insert(reg);
// 					pFace F1 = R_face(reg , 0);
// 					pFace F2 = R_face(reg , 1);
// 					pFace F3 = R_face(reg , 2);
// 					pFace F4 = R_face(reg , 3);
// 				}
// 			}
// 
// 			RIter_delete(rit);
// 
// 			int FaceId = 1;
// 			FIter fit = M_faceIter(theMesh);
// 			while (pFace EntityFace = FIter_next(fit)){ //Atribuindo Ids às Faces
// 				EN_setID(EntityFace,FaceId);
// 				if (F_numRegions(EntityFace)==1){
// 					EN_attachDataInt (EntityFace, MD_lookupMeshDataId("BoundaryFace"), 1);
// 				}
// 
// 				FaceId++;
// 			}
// 			FIter_delete(fit);
// 
// 			EIter eiter = M_edgeIter(theMesh);
// 			int Edgeid=0;
// 			while (pEdge entedge = EIter_next(eiter)){
// 				Edgeid++;
// 				EN_setID(entedge, Edgeid);
// 
// 			}
// 
// 			EIter_delete(eiter);
// 
// 			// Atribuindo CL = 0 a todos os nodes
// 			VIter vit = M_vertexIter(theMesh);
// 			while (pVertex Vertex = VIter_next(vit)){
// 				EN_attachDataDbl (Vertex, MD_lookupMeshDataId("CL"), 0);
// 			}
// 			VIter_delete(vit);
// 
// 
// 		}
// 	}//SP(-0.5,-0.6,0){0.1};

	// Considerando que seja uma malha 2D
	if (I==2){
		//Malha sendo mapeada...
		//map <int, pFace> MeshMap;
// 		map<int, pFace>::iterator Mit;		
// 		FIter fit = M_faceIter(theMesh);
// 		while (pFace Fc = FIter_next(fit)){
// 			MeshMap.insert ( std::pair<int, pFace>(EN_id(Fc),Fc) );
// 		}
// 		FIter_delete(fit);
// 
// 		set<int>::iterator itElements;
// 		for ( itElements=ListofElements.begin() ; itElements != ListofElements.end(); itElements++ ){
// 			Mit=MeshMap.find(*itElements);
// 			BoundaryFaces.insert(Mit->second); // PODE MUITO BEM SER UM LIST
// 			triangles.insert(Mit->second);
// 		}

		int id = 0;		
		EIter eiter = M_edgeIter(theMesh);
		while (pEdge entedge = EIter_next(eiter)){
			id++;
			EN_setID(entedge, id);
			if(E_numFaces(entedge)==1){

				MeshBoundEdges.insert(entedge);
				pVertex VT1 = E_vertex(entedge , 0);
				if (EN_id(VT1)>GreaterID){
					GreaterID = EN_id(VT1);
				}
				pVertex VT2 = E_vertex(entedge , 1);
				if (EN_id(VT2)>GreaterID){
					GreaterID = EN_id(VT2);
				}
				MeshBoundNodes.insert(VT1);
				MeshBoundNodes.insert(VT2);
			}
		}

		EIter_delete(eiter);
		
		
		VIter vit = M_vertexIter(theMesh); // Atribuindo CL = 0 a todos os nodes
		while (pVertex Vertex = VIter_next(vit)){
			EN_attachDataDbl (Vertex, MD_lookupMeshDataId("CL"), 0);
		}
		VIter_delete(vit);
		
	}

	return GreaterID;
}


//          //          Funcao BoundaryElements         //           //

// void Remover::BoundaryElements3D(set <int> ListofElements, pMesh theMesh){
// 
// 
// 	cout << "          Iniciou BoundaryElements3D" << endl;
// 
// 	// preciso agora de uma funcao que a partir do id me devolva o objeto. Mas enquanto nao vem, eu vou usando o iterador mesmo...
// 
// 
// 	set<pRegion>::iterator itTetras;
// 	for ( itTetras=Tetras.begin() ; itTetras != Tetras.end(); itTetras++ ){
// 
// 		pRegion Tetraedro = *itTetras;
// 
// 		pFace Face1 = R_face(Tetraedro , 0);
// 		pFace Face2 = R_face(Tetraedro , 1);
// 		pFace Face3 = R_face(Tetraedro , 2);
// 		pFace Face4 = R_face(Tetraedro , 3);
// 
// 		BoundaryFaces.insert(Face1);
// 		BoundaryFaces.insert(Face2);
// 		BoundaryFaces.insert(Face3);
// 		BoundaryFaces.insert(Face4);
// 
// 
// 		if (F_numRegions(Face1)==1){
// 			GeomBoundFaces.insert(Face1);  // NAO PRECISA SER SET
// 		}
// 		if (F_numRegions(Face2)==1){
// 			GeomBoundFaces.insert(Face2);
// 		}
// 		if (F_numRegions(Face3)==1){
// 			GeomBoundFaces.insert(Face3);
// 		}
// 		if (F_numRegions(Face4)==1){
// 			GeomBoundFaces.insert(Face4);
// 		}
// 
// 
// 		pEdge Edge1 = Tetraedro->get(1, 0);
// 		pEdge Edge2 = Tetraedro->get(1, 1);
// 		pEdge Edge3 = Tetraedro->get(1, 2);
// 		pEdge Edge4 = Tetraedro->get(1, 3);
// 		pEdge Edge5 = Tetraedro->get(1, 4);
// 		pEdge Edge6 = Tetraedro->get(1, 5);
// 
// 		BoundaryEdges.insert(Edge1);
// 		BoundaryEdges.insert(Edge2);
// 		BoundaryEdges.insert(Edge3);
// 		BoundaryEdges.insert(Edge4);
// 		BoundaryEdges.insert(Edge5);
// 		BoundaryEdges.insert(Edge6);
// 
// 		pVertex Vertex1 = Tetraedro->get(0, 0);
// 		pVertex Vertex2 = Tetraedro->get(0, 1);
// 		pVertex Vertex3 = Tetraedro->get(0, 2);
// 		pVertex Vertex4 = Tetraedro->get(0, 3);
// 
// 		BoundaryNodes.insert(Vertex1);
// 		BoundaryNodes.insert(Vertex2);
// 		BoundaryNodes.insert(Vertex3);
// 		BoundaryNodes.insert(Vertex4);
// 
// 		if (E_numFaces(Edge1)==1){
// 			GeomBoundEdges.insert(Edge1);
// 			//double firstvalue = 1;
// 			//EN_attachDataInt (Edge1, MD_lookupMeshDataId("GeomBoundEdge"), firstvalue);
// 		};
// 		if (E_numFaces(Edge2)==1){
// 			GeomBoundEdges.insert(Edge2);
// 			//double firstvalue = 1;
// 			//EN_attachDataInt (Edge2, MD_lookupMeshDataId("GeomBoundEdge"), firstvalue);
// 		};
// 		if (E_numFaces(Edge3)==1){
// 			GeomBoundEdges.insert(Edge3);
// 			//double firstvalue = 1;
// 			//EN_attachDataInt (Edge3, MD_lookupMeshDataId("GeomBoundEdge"), firstvalue);
// 		};
// 		if (E_numFaces(Edge4)==1){
// 			GeomBoundEdges.insert(Edge4);
// 			//double firstvalue = 1;
// 			//EN_attachDataInt (Edge4, MD_lookupMeshDataId("GeomBoundEdge"), firstvalue);
// 		};
// 		if (E_numFaces(Edge5)==1){
// 			GeomBoundEdges.insert(Edge5);
// 			//double firstvalue = 1;
// 			//EN_attachDataInt (Edge4, MD_lookupMeshDataId("GeomBoundEdge"), firstvalue);
// 		};
// 		if (E_numFaces(Edge6)==1){
// 			GeomBoundEdges.insert(Edge6);
// 			//double firstvalue = 1;
// 			//EN_attachDataInt (Edge4, MD_lookupMeshDataId("GeomBoundEdge"), firstvalue);
// 		};
// 	}
// 
// }


void Remover::BoundaryElements2D(pMesh theMesh, const std::list<pEntity>& elementList){
	cout << "          Iniciou BoundaryElements2D" << endl;//ListofElements.size() << endl;
	std::list<pFace>::const_iterator itFaces;
	
	for ( itFaces=elementList.begin() ; itFaces != elementList.end(); itFaces++ ){
		pFace Face = *itFaces;

		pEdge Edge1 = Face->get(1, 0);
		pEdge Edge2 = Face->get(1, 1);
		pEdge Edge3 = Face->get(1, 2);

		BoundaryEdges.insert(Edge1);
		BoundaryEdges.insert(Edge2);
		BoundaryEdges.insert(Edge3);

		pVertex Vertex1 = Face->get(0, 0);
		pVertex Vertex2 = Face->get(0, 1);
		pVertex Vertex3 = Face->get(0, 2);

		BoundaryNodes.insert(Vertex1);
		BoundaryNodes.insert(Vertex2);
		BoundaryNodes.insert(Vertex3);

		if (E_numFaces(Edge1)==1){
			GeomBoundEdges.insert(Edge1);
		}

		if (E_numFaces(Edge2)==1){
			GeomBoundEdges.insert(Edge2);
		}

		if (E_numFaces(Edge3)==1){
			GeomBoundEdges.insert(Edge3);
		}

	}
}

// Nessa funcao as aprtes comentadas serviriam apenas para o caso 3D...
void Remover::removeinternalfaces_2D(pMesh theMesh, std::list<pEntity> &elementList){
	cout << "Iniciou removeinternalfaces" << endl;
	std::list<pFace>::iterator itfaces;
	for ( itfaces=elementList.begin() ; itfaces != elementList.end(); itfaces++ ){
		theMesh->DEL(*itfaces);
	}
}

void Remover::removestgNodes(pMesh theMesh){

	cout << "Iniciou RemoveStgNodes" << endl;

	theMesh->modifyState(2,1,0);
	theMesh->modifyState(2,1,1);

	theMesh->modifyState(1,2,0);
	theMesh->modifyState(1,2,1);

	int w=1;
	int conter = 0;

	list<pVertex> ControlNodes;
	set<pVertex> NewNodes;

	set<pFace>::iterator itfcs;
	set<pVertex>::iterator itnds;
	list<pVertex>::iterator litnds;

	for ( itnds=BoundaryNodes.begin() ; itnds != BoundaryNodes.end(); itnds++ ){
		ControlNodes.push_back(*itnds);
	}

	while (ControlNodes.size()>0){

		set<pFace> RemFACES;
		set<pEdge> VEDGES;
		conter++;

		for ( litnds=ControlNodes.begin() ; litnds != ControlNodes.end(); litnds++ ){
			pVertex Node = *litnds;

			if (MeshBoundNodes.count(Node)==0){ // Alterar isso para nao usar o set, usar attachdataint

				int numEdges = V_numEdges(Node);
				int numFaces = V_numFaces(Node); // O Modifystate entre nodes e faces está com problema...


				for (int i=0; i<numEdges; i++){
					pEdge VEDG = V_edge(Node,i);
					if(E_numFaces(VEDG)!=0){
						VEDGES.insert(VEDG);
					}

				}


				set<pEdge>::iterator itedgs;

				set<pFace> EFACES;

				for ( itedgs=VEDGES.begin() ; itedgs != VEDGES.end(); itedgs++ ){
					if (E_numFaces(*itedgs)>0){
						EFACES.insert(E_face(*itedgs,0));
						if(E_numFaces(*itedgs)==2){
							pEdge DG = E_face(*itedgs,1);
							EFACES.insert(DG);

						}
					}
				}

				// Porque o modifystate entre nos e faces esta com problema eu tenho que comparar o tamanho dos sets
				if ((VEDGES.size() > EFACES.size()+1)||((VEDGES.size() > EFACES.size())&& EFACES.size()==1) ){
					for(itfcs=EFACES.begin() ; itfcs != EFACES.end(); itfcs++ ){
						for(int y=0; y<=2; y++){

							pVertex V = F_vertex(*itfcs, y);
							BoundaryNodes.insert(V); // Inserir em outro set pra depois inserir nesse
							NewNodes.insert(V);

						}

						for(int y=0; y<=2; y++){

							pEdge E = F_edge(*itfcs, y);
							BoundaryEdges.insert(E);

						}

						RemFACES.insert(*itfcs);

					}
				}

				VEDGES.clear();
				EFACES.clear();

			}
			else{

				int numEdges = V_numEdges(Node);
				int numFaces = V_numFaces(Node);

				int Ctr=0;

				for (int i=0; i<numEdges; i++){
					pEdge VEDG = V_edge(Node,i);
					if(E_numFaces(VEDG)!=0){
						VEDGES.insert(VEDG);
					}
					if (MeshBoundEdges.count(VEDG)==1){
						Ctr = Ctr + E_numFaces(VEDG);
					}
				}
				if (Ctr==0){
					set<pEdge>::iterator itedgs;
					set<pFace>::iterator itfaces;
					set<pFace> EFACES;

					for ( itedgs=VEDGES.begin() ; itedgs != VEDGES.end(); itedgs++ ){
						if (E_numFaces(*itedgs)>0){
							EFACES.insert(E_face(*itedgs,0));
							if(E_numFaces(*itedgs)==2){
								pEdge DG = E_face(*itedgs,1);
								EFACES.insert(DG);
							}
						}
					}
					if(EFACES.size()>0){
						for(itfaces=EFACES.begin(); itfaces!=EFACES.end(); itfaces++){

							for(int y=0; y<=2; y++){

								pVertex V = F_vertex(*itfaces, y);
								BoundaryNodes.insert(V); // Inserir em outro set pra depois inserir nesse		
								NewNodes.insert(V);
							}

							for(int y=0; y<=2; y++){

								pEdge E = F_edge(*itfaces, y);
								BoundaryEdges.insert(E);

							}

							RemFACES.insert(*itfaces);

						}
					}
					EFACES.clear();
				}
				VEDGES.clear();
			}
		}

		ControlNodes.clear();

		for ( itnds=NewNodes.begin() ; itnds != NewNodes.end(); itnds++ ){
			ControlNodes.push_back(*itnds);
		}
		NewNodes.clear();


		for (itfcs=RemFACES.begin() ; itfcs != RemFACES.end(); itfcs++ ){
			theMesh->DEL(*itfcs);
		}
		RemFACES.clear();
		//w++;

		theMesh->modifyState(2,1,0);
		theMesh->modifyState(2,1,1);

		theMesh->modifyState(1,2,0);
		theMesh->modifyState(1,2,1);

	}

}

void Remover::removeinternaledges(pMesh theMesh){
	cout << "          Iniciou removeinternaledges" << endl;

	set<pEdge>::iterator itEdges;

	EIter edit = M_edgeIter(theMesh);
	while (pEdge entity = EIter_next(edit)){

		if (E_numFaces(entity)==0 && MeshBoundEdges.find(entity)==MeshBoundEdges.end()){	// NNAAOO PRECISA SER SET

			RemovefromBoundEdges.insert(entity);

		}
	}
	EIter_delete(edit);


	for ( itEdges=RemovefromBoundEdges.begin() ; itEdges != RemovefromBoundEdges.end(); itEdges++ ){
		BoundaryEdges.erase(*itEdges);
		theMesh->DEL(*itEdges);
	}

	RemovefromBoundEdges.clear();

}

void Remover::removeinternalnodes(pMesh theMesh){
	cout << "Iniciou removeinternalnodes" << endl;

	set <pVertex> ToRemoveVertex;
	VIter vit = M_vertexIter(theMesh);

	while (pVertex entity = VIter_next(vit)){
		if (V_numEdges(entity)==0){
			ToRemoveVertex.insert(entity);
		}
	}

	VIter_delete(vit);

	set<pVertex>::iterator itVts;
	for ( itVts=ToRemoveVertex.begin() ; itVts != ToRemoveVertex.end(); itVts++ ){
		theMesh->DEL(*itVts);
	}

}

set <pEdge> Remover::ReturnBoundaryEdges (){
	return BoundaryEdges;
}

set <pVertex> Remover::ReturnBoundaryNodes (){
	return BoundaryNodes;
}

// set <pFace> Remover::ReturnBoundaryFaces (){
// 	return BoundaryFaces;
// }

void Remover::removeexternaledges(pMesh theMesh){ // Essa funcao agora nao funciona por que o set GeomBoundEdges nao é mais preenchido

	cout << "Iniciou removeexternaledges" << endl;

	set<pEdge>::iterator itEdges;

	set<pEdge> EREMOVE;
	EIter eit = M_edgeIter(theMesh);
	while (pEdge entity = EIter_next(eit)){
		if (E_numFaces(entity)==0){
			EREMOVE.insert(entity);
		}
	}
	EIter_delete(eit);


	for ( itEdges=EREMOVE.begin() ; itEdges != EREMOVE.end(); itEdges++ ){
		theMesh->DEL(*itEdges);
	}
	EREMOVE.clear();

}

void Remover::removeexternalnodes(pMesh theMesh){
	cout << "Iniciou removeexternalnodes" << endl; 
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
	// Aqui é pra garantir que os pontos foram todos removidos
	set <pVertex> ToRemoveVertex;
	VIter vit = M_vertexIter(theMesh);
	while (pVertex entity = VIter_next(vit)){
		if (V_numEdges(entity)==0){
			ToRemoveVertex.insert(entity);
		}
	}
	VIter_delete(vit);
	set<pVertex>::iterator itVts;
	for ( itVts=ToRemoveVertex.begin() ; itVts != ToRemoveVertex.end(); itVts++ ){
		theMesh->DEL(*itVts);
	}
	cout << "terminou removeexternalnodes" << endl;
}

void Remover::removetetras(pMesh theMesh){

	cout << "Iniciou removetetras" << endl;

	set<pRegion>::iterator itTetras;
	for ( itTetras=Tetras.begin() ; itTetras != Tetras.end(); itTetras++ ){

		theMesh->DEL(*itTetras);

	}
}



void Remover::RemoveStgElements_2D(pMesh theMesh, float GE) { //detecta elementos estranhos numa malha 2D, os remove e atualiza set's de contorno

	set<pVertex>::iterator itVertex;
	set<pEdge>::iterator itedge;
	set<pEdge>::iterator itedge2;
	set<pFace>::iterator itface;

	set <pVertex> stgNodes;
	set <pEdge> stgEdges;
	set <pFace> stgFaces;

	stgNodes.clear();
	stgEdges.clear();
	stgFaces.clear();

	cout << "iniciou RemoveStgElements_2D" << endl;



	int h=1;
	while (h==1){
		h=0;  // Novo criterio: Todas as faces tem que ser vizinhas de ao menos 2 outras faces, se a face for vizinha de apenas 1 face, ou nenhuma outra face, então ela é uma face estranha e deve ser removida.
		for ( itVertex=BoundaryNodes.begin() ; itVertex != BoundaryNodes.end(); itVertex++ ){
			//cout << endl;
			//cout << "Analisando vertex " << EN_id(*itVertex) << " numero de faces deste vertice: " << V_numFaces(*itVertex) << endl;

			if(V_numFaces(*itVertex)>0){


				int NEVertex = V_numEdges(*itVertex);
				set<pEdge> NeighEdges;
				set<pFace> NeighFaces;

				for (int y = 0; y<NEVertex; y++){ // inserindo todas as arestas do no no set

					// para cada nó, verificar a qtd de faces e para cada face, verificar a quantidade de faces vizinhas dela...
					// se a qtd de faces vizinhas daquela face for <= 1, trata-se de face estranha.
					pEdge NeighEdge = V_edge(*itVertex,y);
					NeighEdges.insert(NeighEdge);
				}


				for(itedge=NeighEdges.begin();itedge!=NeighEdges.end();itedge++){ // inserindo todas as faces do no no set
					int Nfaces=E_numFaces(*itedge);
					if(Nfaces>0){
						pFace Nface = E_face(*itedge,0);
						NeighFaces.insert(Nface);
						if (Nfaces>1){
							NeighFaces.insert(E_face(*itedge,1));
						}
					}
				}


				NeighEdges.clear();

				//	cout<<"Terminou de inserir todas as faces vizinhas dos Nodes, agora analisado cada uma:"<<endl;

				int psr = 0;  // <--- psr é plane surface reader
				int stgpsr =0; // <--- plane surface reader dos elementos estranhos, guarda a informacao de qual ps estamos trabalha	ndo no momento

				for(itface=NeighFaces.begin();itface!=NeighFaces.end();itface++){ // analisando a qtd de vizinhos da face
					int FV=0;
					//		cout<< "Analisando face..." << endl;
					for (int i=0;i<=2;i++){
						pEdge Edg = F_edge(*itface,i);
						if (E_numFaces(Edg)==2){
							FV++;
						}
					}


					if (FV<2){ //todo este laco se aproveita... mas o criterio do if muda...
						h=1;

						if (V_numEdges(*itVertex)<=2){
							if(MeshBoundNodes.count(*itVertex)!=1){
								stgNodes.insert(*itVertex);
							}
						}



						stgFaces.insert(*itface);	


						//theMesh->DEL(stgFace); // removeu a face da malha, nao existe o set de faces para 2D


					}

				}

				NeighFaces.clear();
			}
		}		


		//cout << "todos os vertex analisados, procedendo remocao de elementos estranhos" << endl;

		// Agora, neste momento comecam a ser removidos os elementos da malha

		for (itface=stgFaces.begin() ;  itface != stgFaces.end(); itface++ ){
			for(int y=0; y<=2; y++){
				pVertex V = F_vertex(*itface, y);
				BoundaryNodes.insert(V);
			}
			for(int y=0; y<=2; y++){
				pEdge E = F_edge(*itface, y);
				BoundaryEdges.insert(E);
			}

			theMesh->DEL(*itface);

		}

		cout << endl;

		theMesh->modifyState(2,1,0);
		theMesh->modifyState(2,1,1);
		theMesh->modifyState(1,2,0);
		theMesh->modifyState(1,2,1);
		theMesh->modifyState(0,1,0);
		theMesh->modifyState(0,1,1);
		theMesh->modifyState(0,2,0);
		theMesh->modifyState(0,2,1);

		stgNodes.clear();
		stgEdges.clear();
		stgFaces.clear();

		if (h==1){
			//	cout << "Faces estranhas foram removidas, iniciando nova varredura" << endl;
			//	cout << endl;
		}
		if (h==0){
			//	cout << "Nao foram encontradas faces estranhas nesta varredura, finalizando limpeza" << endl;
		}
	}
	cout << "          finalizou RemoveStgElements_2D" << endl;
	cout << endl;


}

// void Remover::RemoveStgElements_3D(pMesh theMesh) {
// 	set<pVertex>::iterator itVertex;
// 	set<pEdge>::iterator itedge;
// 	set<pFace>::iterator itface;
// 	set<pRegion>::iterator ittetra;
// 
// 	set <pVertex> stgNodes;
// 	set <pEdge> stgEdges;
// 	set <pFace> stgFaces;
// 	set <pRegion> stgTetras;
// 
// 	int h=1;
// 	while (h==1){
// 		h=0;
// 		for ( itVertex=BoundaryNodes.begin() ; itVertex != BoundaryNodes.end(); itVertex++ ){ // este laco percorre os nós de contorno reconhecendo e capturando os elementos estranhos, depois de percorrido todo o set, os elementos sao apagados e o set é percorrido novamente para verificar se foi criado algum novo elemento estranho.
// 			if (V_numRegions(*itVertex)==1){
// 				h=1;
// 				pEdge Ed1 = V_edge(*itVertex, 0);
// 				pEdge Ed2 = V_edge(*itVertex, 1);
// 				pEdge Ed3 = V_edge(*itVertex, 2);
// 
// 				stgEdges.insert(Ed1); // capturando arestas "estranhas" a ser removidas
// 				stgEdges.insert(Ed2);
// 				stgEdges.insert(Ed3);
// 
// 				stgFaces.insert(E_face(Ed1, 0)); // capturando faces "estranhas" a ser removidas
// 				stgFaces.insert(E_face(Ed1, 1));
// 				stgFaces.insert(E_face(Ed2, 0));
// 				stgFaces.insert(E_face(Ed2, 1));
// 
// 				pRegion Tet = F_region(E_face(Ed2, 1), 0);
// 
// 				stgTetras.insert(Tet); // capturando tetra "estranho" a ser removido
// 
// 				for (int r = 0; r <= 3; r++){ // inserindo no boundary face a nova face
// 					pFace boundaryface = R_face(Tet , r);
// 					if (stgFaces.count(boundaryface) == 0){ // se a face nao for encontrada no stgfaces, ela é de contorno
// 						BoundaryFaces.insert(boundaryface);  // inserindo os novos boundary elements em seus set's
// 						BoundaryEdges.insert(F_edge(boundaryface, 0));
// 						BoundaryEdges.insert(F_edge(boundaryface, 1));
// 						BoundaryEdges.insert(F_edge(boundaryface, 2));
// 					}
// 				}
// 			}
// 
// 			for ( itVertex=stgNodes.begin() ; itVertex != stgNodes.end(); itVertex++ ){
// 				theMesh->DEL(*itVertex); //removendo nós estranhos da malha
// 				BoundaryNodes.erase(*itVertex);
// 			}
// 			for ( itedge=stgEdges.begin() ; itedge != stgEdges.end(); itedge++ ){
// 				theMesh->DEL(*itedge); //removendo arestas estranhas da malha
// 				BoundaryEdges.erase(*itedge);
// 			}
// 			for ( itface=stgFaces.begin() ; itface != stgFaces.end(); itface++ ){
// 				theMesh->DEL(*itface); //removendo faces da malha
// 				BoundaryFaces.erase(*itface);
// 			}
// 			for ( ittetra=stgTetras.begin() ; ittetra != stgTetras.end(); ittetra++ ){
// 				theMesh->DEL(*ittetra); //removendo tetraedros da malha
// 			}
// 			cout << endl;
// 			if (h==1){
// 				cout << "Tetras estranhos foram removidos, iniciando nova varredura" << endl;
// 			}
// 			if (h==0){
// 				cout << "Nao foram encontrados Tetras estranhos nesta varredura, finalizando limpeza" << endl;
// 			}
// 		}
// 
// 	}
// }
// 
// 
// int Remover::identifyvolumes_3D(pMesh Mesh){ // Versao 3D da funcao "identifysurfaces_2D"
// 
// 	FIter fit = M_faceIter(Mesh);  // Outro iterador que cedo ou tarde tem que sair daqui
// 	while (pFace face = FIter_next(fit)){
// 		double firstvalue = 0;
// 		EN_attachDataInt (face, MD_lookupMeshDataId("Volume"), firstvalue);
// 	}
// 	FIter_delete(fit);
// 
// 
// 
// 	int identify =1;
// 	set<pRegion> CpTetras; // uma copia de BoundaryTetras a ser preenchida e depois despreenchida...
// 	set<pRegion>::iterator itIdent;
// 
// 	for(itIdent=Tetras.begin(); itIdent!=Tetras.end(); itIdent++){ // Preenchendo a copia...
// 		CpTetras.insert(*itIdent);
// 	}
// 
// 	set <pRegion> StTetras;
// 	StTetras.insert(*CpTetras.begin());
// 
// 
// 	int size=1;
// 	int VLcounter = 0;
// 	while (size==1){
// 		int woutNInst = 1; // <-- Abreviacao de without new insertions...
// 		size=0;
// 
// 		///////////////////////// Pegando todos os elementos de UM volume ////////////////////////////
// 
// 		cout << "Pegando todos os elementos de UM volume" << endl;
// 
// 		for(itIdent=StTetras.begin(); itIdent!=StTetras.end(); itIdent++){ //comeca a percorrer o set StTetras que so tem inicialmente um elemento...
// 
// 			if (identify>=1){
// 				set <pFace> Faces;
// 				for (int r=0;r<=3;r++){// daqui...
// 					pFace face = R_face(*itIdent,r);
// 					Faces.insert(face); //inserindo os faces vizinho no set Faces...
// 				}
// 				set<pFace>::iterator itface;
// 				for( itface=Faces.begin();itface!=Faces.end();itface++){ // ate aqui, pegando vizinhos...
// 					set<pFace>::iterator itFC;
// 					for (itFC=CpTetras.begin(); itFC!=CpTetras.end(); itFC++){//verificando se o vizinho está na lista...
// 
// 
// 						pRegion Tetra1 = F_region(*itface,0);
// 						pRegion Tetra2 = F_region(*itface,1); // essa parte aqui da segmentation fault na fronteira
// 
// 
// 
// 						if (*itFC==Tetra1){ // se estiver, entra pro set surfaces na ordem...
// 							StTetras.insert(*itFC);
// 							woutNInst = 0;
// 						}
// 
// 						if (*itFC==Tetra2){
// 							StTetras.insert(*itFC);
// 							woutNInst = 0;
// 						}
// 					}
// 				}
// 				Faces.clear();
// 			}
// 			// identify...
// 		}
// 		////////////////// Terminou de pegar os elementos de UM volume ///////////////////////
// 
// 
// 
// 
// 
// 		///////////////// Retirando os elementos daquele plane surface da lista /////////////////////
// 		for(itIdent=StTetras.begin(); itIdent!=StTetras.end(); itIdent++){ //entra de novo no laco e remove de cpBdFaces
// 			if (CpTetras.count(*itIdent)==1){
// 				CpTetras.erase(*itIdent);
// 			}
// 		}
// 		//////////////////////////////// Terminou de retirar ////////////////////////////////////////
// 
// 
// 
// 		if (CpTetras.size()!=0){
// 			size=1;  // Se ainda existem elementos em CpBdFaces, entao o laco while deve rodar novamente
// 		}
// 
// 		identify++;
// 
// 
// 
// 
// 		if (size ==1){ // Se ainda tem elementos e nao houveram novas insercoes, entao foi separado uma casca de volume, anexando a informacao (qual o volume) aos elementos
// 			if (woutNInst == 1) {
// 				pRegion Tet = *CpTetras.begin();
// 				cout << "um Volume: " << endl;
// 
// 				VLcounter++;
// 
// 				//cout << "PScounter : " << PScounter << endl;
// 
// 				for(itIdent=StTetras.begin(); itIdent!=StTetras.end(); itIdent++){
// 					//cout << EN_id(*itIdent) << endl;
// 					//	Clcounter++; //colocar aqui a marcação do tag com o valor de PScounter...
// 					for (int f=0;f<=3;f++){
// 						pFace face = R_face(*itIdent , f); // pegando as faces vizinhas do tetra em questao...
// 
// 						double VL = VLcounter;// anexa a informacao de qual plane surface o elemento pertence
// 
// 						EN_attachDataInt (face, MD_lookupMeshDataId("Volume"), VL); // anexa a informacao de qual volume o elemento pertence
// 
// 
// 						// Atribuindo também às arestas e os vertex
// 
// 						for (int g=0; g<=2;g++){
// 
// 							pEdge edge = F_edge(face, g);
// 							pVertex vertex = F_vertex(face, g);
// 
// 							EN_attachDataInt (edge, MD_lookupMeshDataId("Volume"), VL);
// 							EN_attachDataInt (vertex, MD_lookupMeshDataId("Volume"), VL);
// 						}
// 
// 					}
// 
// 					int vlreader = 0;
// 					EN_getDataInt (*itIdent, MD_lookupMeshDataId("Volume"), &vlreader);
// 					//cout << "psreader: " << psreader << endl;
// 
// 				}
// 				StTetras.clear(); // limpa o set
// 				StTetras.insert(Tet); //para recomecar o ciclo de insersoes...
// 			}
// 		}
// 	}
// 	cout << "um Volume: " << VLcounter << endl;
// 	VLcounter++;
// 	//cout << "PScounter : " << PScounter << endl;
// 
// 
// 	/////////////////////
// 
// 
// 	// Agora, se o cpbdsurfaces esvaziou, anexa a informacao de qual é o plane surface e termina a funcao
// 	for(itIdent=StTetras.begin(); itIdent!=StTetras.end(); itIdent++){ // Apresenta os elementos que foram inseridos no laco...
// 		for (int f=0;f<=3;f++){
// 			pFace face = R_face(*itIdent , f); // pegando as faces vizinhas do tetra em questao...
// 
// 			double VL = VLcounter;// anexa a informacao de qual plane surface o elemento pertence
// 
// 			EN_attachDataInt (face, MD_lookupMeshDataId("Volume"), VL); // anexa a informacao de qual plane surface o elemento pertence
// 
// 
// 			// Atribuindo também às arestas e os vertex
// 
// 			for (int g=0; g<=2;g++){
// 
// 				pEdge edge = F_edge(face, g);
// 				pVertex vertex = F_vertex(face, g);
// 
// 				EN_attachDataInt (edge, MD_lookupMeshDataId("Volume"), VL);
// 				EN_attachDataInt (vertex, MD_lookupMeshDataId("Volume"), VL);
// 			}
// 
// 		}
// 
// 		int vlreader = 0;
// 		EN_getDataInt (*itIdent, MD_lookupMeshDataId("Volume"), &vlreader);
// 		//cout << "psreader: " << psreader << endl;
// 	}
// 	return VLcounter;
// 	cout << "Finalmente os volumes são " << VLcounter << endl;
// 
// 
// 
// }

int Remover::Identify_and_Remove_Edges_2D(pMesh theMesh){ // Nova versao de IdentifySurfaces_2D deve rodar antes de remover as arestas

	cout << "Iniciou Identify_and_Remove_Edges_2D" << endl;



	EIter eit = M_edgeIter(theMesh);
	while (pEdge edge = EIter_next(eit)){
		int firstvalue = 0; // atribuindo valor zero a todas as arestas no tag PlaneSurface.
		EN_attachDataInt (edge, MD_lookupMeshDataId("PlaneSurface"), firstvalue);

	}
	EIter_delete(eit);

	set<pEdge>:: iterator pite;

	set<pEdge> ControlEdges;
	for(pite=BoundaryEdges.begin(); pite!=BoundaryEdges.end(); pite++){
		ControlEdges.insert(*pite);
	}

	int PScounter = 1;
	while(ControlEdges.size()!=0){
		set <pEdge> PSFEdges;
		PSFEdges.insert(*ControlEdges.begin());

		int withNInst = 1;

		while (withNInst==1){  // Esse while separa as arestas pertencentes a um Plane Surface

			withNInst = 0; // <-- Abreviacao de with new insertions...

			for(pite=PSFEdges.begin(); pite!=PSFEdges.end(); pite++){

				pVertex EdVt1 = E_vertex(*pite, 0);
				pVertex EdVt2 = E_vertex(*pite, 1);

				int size = PSFEdges.size();

				for(int N=0; N<V_numEdges(EdVt1); N++){  // Pegando as arestas de um nó e inserindo no set
					pEdge E = V_edge(EdVt1, N);

					if(E_numFaces(E)<=1 && ControlEdges.count(E)==1){ // TA ESCAPANDO PELAS ARESTAS DO CONTORNO DA GEOMETRIA
						PSFEdges.insert(E);
						if (size != PSFEdges.size()){
							withNInst=1; // Se foi uma nova insercao seta -> withNInst=1;
						}
					}

				}

				for(int N=0; N<V_numEdges(EdVt2); N++){// Pegando as arestas de um nó e inserindo no set
					pEdge E = V_edge(EdVt2, N);
					if(E_numFaces(E)<=1 && ControlEdges.count(E)==1){
						PSFEdges.insert(E);
						if (size != PSFEdges.size()){
							withNInst=1; // Se foi uma nova insercao seta -> withNInst=1;
						}
					}
				}
			}
		} // Aqui fecha o while, todas as arestas de um plane surface foram separadas


		for(pite=PSFEdges.begin(); pite!=PSFEdges.end(); pite++){

			if(E_numFaces(*pite)==0 && MeshBoundEdges.count(*pite)==0){				

				theMesh->DEL(*pite);
				BoundaryEdges.erase(*pite);
				ControlEdges.erase(*pite);
			}
			if(MeshBoundEdges.count(*pite)==1 || E_numFaces(*pite)==1){
				//Esses aqui deixa e marca com o numero do PS
				EN_attachDataInt (*pite, MD_lookupMeshDataId("PlaneSurface"), PScounter);
				ControlEdges.erase(*pite);

			}
		}

		PSFEdges.clear();
		PScounter++;
	}

	PScounter = PScounter - 1;

	return PScounter;
}


// int Remover::identifysurfaces_2D(pMesh Mesh){
// 
// 	// Primeiro de tudo, percorrer todas as arestas da malha e setar o tag Plane Surface para zero.
// 	// ESTE É MAIS UM ITERADOR DESTES QUE EU TENHO QUE RETIRAR DEPOIS... NAO DA PRA FICAR PERCORRENDO TODA A MALHA DE VEZ EM QUANDO... PRECISO CRIAR UM FUNCAO QUE FACA ISSO UMA VEZ E CRIE UMA MATRIZ COM PONTEIROS, TAGS, IDS, E O QUE MAIS EU PRECISAR PARA USAR QUANDO TIVER QUE USAR.
// 
// 
// 	// Esta rotina é necessária para a funcao Removestgelements2D reconhecer qual aresta deve ser inserida no set BoundaryEdges
// 	// a que tiver o tag "PlaneSurface" setado com zero é a que nao pertencia ao planesurface e deve ser inserida entao...
// 
// 
// 	EIter eit = M_edgeIter(Mesh);
// 	while (pEdge edge = EIter_next(eit)){
// 		double firstvalue = 0; // atribuindo valor zero a todas as arestas no tag PlaneSurface.
// 		EN_attachDataInt (edge, MD_lookupMeshDataId("PlaneSurface"), firstvalue);
// 	}
// 	EIter_delete(eit);
// 
// 
// 	int identify =1;
// 	set<pFace> CpBdFaces; // uma copia de BoundaryFaces a ser preenchida e depois despreenchida...
// 	set<pFace>::iterator itIdent;
// 	for(itIdent=BoundaryFaces.begin(); itIdent!=BoundaryFaces.end(); itIdent++){ // Preenchendo a copia...
// 		CpBdFaces.insert(*itIdent);
// 	}
// 	cout << "Tamanho de cpbdFaces: " << CpBdFaces.size() << endl;
// 
// 
// 	set <pFace> surfaces;
// 	surfaces.insert(*CpBdFaces.begin());
// 
// 	cout << "Tamanho de CpBdFaces: " << CpBdFaces.size() << endl;
// 
// 	int size=1;
// 	int PScounter = 0;
// 	while (size==1){
// 		int woutNInst = 1; // <-- Abreviacao de without new insertions...
// 		size=0;
// 
// 		///////////////////////// Pegando todos os elementos de UM plane surface////////////////////////////
// 
// 		for(itIdent=surfaces.begin(); itIdent!=surfaces.end(); itIdent++){ //comeca a percorrer o set surfaces que so tem inicialmente um elemento...
// 			if (identify>=1){
// 				set <pEdge> Edges;
// 				for (int r=0;r<=2;r++){// daqui...
// 					pEdge edge = F_edge(*itIdent,r);
// 					Edges.insert(edge); //inserindo os edges vizinho no set Edges...
// 				}
// 
// 				set<pEdge>::iterator itedg;
// 				for( itedg=Edges.begin();itedg!=Edges.end();itedg++){ // ate aqui, pegando vizinhos...
// 					set<pEdge>::iterator itED;
// 					for (itED=CpBdFaces.begin(); itED!=CpBdFaces.end(); itED++){//verificando se o vizinho está na lista...
// 						pFace Face1 = E_face(*itedg,0);
// 						pFace Face2 = E_face(*itedg,1); // essa parte aqui da segmentation fault na fronteira
// 
// 
// 						if (*itED==Face1){ // se estiver, entra pro set surfaces na ordem...
// 							surfaces.insert(*itED);
// 							woutNInst = 0;
// 						}
// 
// 						if (*itED==Face2){
// 							surfaces.insert(*itED);
// 							woutNInst = 0;
// 						}
// 					}
// 				}
// 				Edges.clear();
// 			}
// 			// identify...
// 		}
// 		////////////////// Terminou de pegar os elementos de UM plane surface ///////////////////////
// 
// 
// 		///////////////// Retirando os elementos daquele plane surface da lista /////////////////////
// 		for(itIdent=surfaces.begin(); itIdent!=surfaces.end(); itIdent++){ //entra de novo no laco e remove de cpBdFaces
// 			if (CpBdFaces.count(*itIdent)==1){
// 				CpBdFaces.erase(*itIdent);
// 			}
// 		}
// 		//////////////////////////////// Terminou de retirar ////////////////////////////////////////
// 
// 
// 
// 		if (CpBdFaces.size()!=0){
// 			size=1;  // Se ainda existem elementos em CpBdFaces, entao o laco while deve rodar novamente
// 		}
// 
// 		identify++;
// 
// 
// 
// 		if (size ==1){ // Se ainda tem elementos e nao houveram novas insercoes, entao foi separado um plane surface, anexando a informacao (qual o plane surface) aos elementos
// 			if (woutNInst == 1) {
// 				pFace Fac = *CpBdFaces.begin();
// 				cout << "um Plane surface: " << endl;
// 				PScounter++;
// 				//cout << "PScounter : " << PScounter << endl;
// 
// 				for(itIdent=surfaces.begin(); itIdent!=surfaces.end(); itIdent++){
// 					//cout << EN_id(*itIdent) << endl;
// 					//	Clcounter++; //colocar aqui a marcação do tag com o valor de PScounter...
// 
// 					pEdge edge0 = F_edge(*itIdent , 0); // pegando as arestas vizinhas da face em questao...
// 					pEdge edge1 = F_edge(*itIdent , 1);
// 					pEdge edge2 = F_edge(*itIdent , 2);
// 
// 					double PS = PScounter;// anexa a informacao de qual plane surface o elemento pertence
// 
// 					EN_attachDataInt (edge0, MD_lookupMeshDataId("PlaneSurface"), PS); // anexa a informacao de qual 
// 					EN_attachDataInt (edge1, MD_lookupMeshDataId("PlaneSurface"), PS); // plane surface o
// 					EN_attachDataInt (edge2, MD_lookupMeshDataId("PlaneSurface"), PS); // elemento pertence
// 
// 					// Atribuindo também aos vertex
// 
// 					for (int g=0; g<=1;g++){
// 
// 						pVertex vertex0 = E_vertex(edge0,g);
// 						pVertex vertex1 = E_vertex(edge1,g);
// 						pVertex vertex2 = E_vertex(edge2,g);
// 
// 						EN_attachDataInt (vertex0, MD_lookupMeshDataId("PlaneSurface"), PS);
// 						EN_attachDataInt (vertex1, MD_lookupMeshDataId("PlaneSurface"), PS);
// 						EN_attachDataInt (vertex2, MD_lookupMeshDataId("PlaneSurface"), PS);
// 
// 					}
// 
// 					int psreader = 0;
// 					EN_getDataInt (edge0, MD_lookupMeshDataId("PlaneSurface"), &psreader);
// 					//cout << "psreader: " << psreader << endl;
// 
// 				}
// 				surfaces.clear(); // limpa o set
// 				surfaces.insert(Fac); //para recomecar o ciclo de insersoes...
// 			}
// 		}
// 	}
// 	cout << "um Plane surface: " << PScounter << endl;
// 	PScounter++;
// 	//cout << "PScounter : " << PScounter << endl;
// 
// 
// 
// 	// Agora, se o cpbdsurfaces esvaziou, anexa a informacao de qual é o plane surface e termina a funcao
// 	for(itIdent=surfaces.begin(); itIdent!=surfaces.end(); itIdent++){ // Apresenta os elementos que foram inseridos no laco...
// 		//	cout << EN_id(*itIdent) << endl;
// 
// 		pEdge edge0 = F_edge(*itIdent , 0); // pegando os elementos arestas vizinhos da face em questao...
// 		pEdge edge1 = F_edge(*itIdent , 1);
// 		pEdge edge2 = F_edge(*itIdent , 2);
// 
// 		double PS = PScounter;// anexa a informacao de qual plane surface o elemento pertence
// 		EN_attachDataInt (edge0, MD_lookupMeshDataId("PlaneSurface"), PS); // anexa a informacao de qual plane surface o
// 		EN_attachDataInt (edge1, MD_lookupMeshDataId("PlaneSurface"), PS); // elemento pertence
// 		EN_attachDataInt (edge2, MD_lookupMeshDataId("PlaneSurface"), PS);
// 
// 
// 		// Atribuindo também aos vertex
// 
// 		for (int g=0; g<=1;g++){
// 
// 			pVertex vertex0 = E_vertex(edge0,g);
// 			pVertex vertex1 = E_vertex(edge1,g);
// 			pVertex vertex2 = E_vertex(edge2,g);
// 
// 			EN_attachDataInt (vertex0, MD_lookupMeshDataId("PlaneSurface"), PS);
// 			EN_attachDataInt (vertex1, MD_lookupMeshDataId("PlaneSurface"), PS);
// 			EN_attachDataInt (vertex2, MD_lookupMeshDataId("PlaneSurface"), PS);
// 
// 		}
// 
// 
// 		int psreader = 0;
// 		EN_getDataInt (edge0, MD_lookupMeshDataId("PlaneSurface"), &psreader);
// 		//	cout << "psreader: " << psreader << endl;
// 	}
// 	return PScounter;
// 	cout << "Finalmente as Plane Surfaces são " << PScounter << endl;
// 
// }


// set <int> Remover::HoleNeighbours(pMesh theMesh2) { // Parece que está funcionando... Preenche uma nova lista de elementos a remover, no caso os vizinhos do buraco
// 
// 	cout << "Entrou em Holeneighbours" << endl;
// 	BoundaryEdges.clear();
// 
// 	EIter eiter = M_edgeIter(theMesh2);
// 	while (pEdge edge = EIter_next(eiter)){
// 		int N = E_numFaces(edge);
// 		if (N==1){
// 			BoundaryEdges.insert(edge);
// 		}
// 	}
// 
// 	EIter_delete(eiter);
// 
// 	cout << "BoundaryEdges.size(): " << BoundaryEdges.size() << endl;
// 
// 
// 	set <pFace> NewListofFaces;
// 	set <int> Listofids;
// 	set <pVertex> BoundNod;
// 
// 
// 	set<pVertex>::iterator iternodes;
// 	set<pEdge>::iterator iterBoundEdges;
// 
// 	for(iterBoundEdges=BoundaryEdges.begin(); iterBoundEdges!=BoundaryEdges.end(); iterBoundEdges++){
// 		pEdge Edge = *iterBoundEdges;
// 
// 		pVertex vert1 = E_vertex(Edge,0);
// 		pVertex vert2 = E_vertex(Edge,1);
// 
// 		BoundNod.insert(vert1);
// 		BoundNod.insert(vert2);
// 
// 	}
// 
// 	for(iternodes=BoundNod.begin(); iternodes!=BoundNod.end(); iternodes++){
// 		pVertex VT = *iternodes;
// 		int NumberEdges = V_numEdges(VT);
// 		for (int t=0;t<NumberEdges;t++){
// 			pEdge BEdge = V_edge(VT, t);
// 			int NumberFaces = E_numFaces(BEdge);
// 			for (int h=0;h<NumberFaces;h++){
// 				pFace BFace = E_face(BEdge, h);
// 				NewListofFaces.insert(BFace);
// 			}
// 		}
// 	}
// 
// 
// 	set<pFace>::iterator iterFaces;
// 	for ( iterFaces=NewListofFaces.begin() ; iterFaces != NewListofFaces.end(); iterFaces++ ){
// 		int id = EN_id(*iterFaces);
// 		Listofids.insert(id);	
// 	}
// 
// 	BoundaryEdges.clear();
// 
// 	return Listofids;
// 
// 
// 
// 
// }



void Remover::Merge_Edges(pMesh theMesh2) { // Copia de uma malha para a outra, as arestas que tem ponto em comum com a malha para onde vai, as arestas que vao ser copiadas sao SourceEdges e elas serao criadas na malha que contem os ReferenceNodes


	cout << "Iniciou Merge_Edges, BoundaryNodes.size() é " << BoundaryNodes.size() << endl;
	cout << "Iniciou Merge_Edges, AdaptEdges.size() é " << AdaptEdges.size() << endl;
	cout << "Iniciou Merge_Edges, BoundaryNodes.size() é " << BoundaryNodes.size() << endl;


	// Pego os BoundaryNodes de theMesh2 e mapeio (coloco em um map)
	set<pVertex>::iterator itRNodes;
	for ( itRNodes=BoundaryNodes.begin() ; itRNodes != BoundaryNodes.end(); itRNodes++ ){
		double xyz[3];
		V_coord(*itRNodes, xyz);
		string ponto = "."; 
		ostringstream oss;
		oss << xyz[0];
		oss << ponto;
		oss << xyz[1];

		string XY = oss.str();

		CommonVertex.insert ( pair<string, pVertex>(XY, *itRNodes) );
	}

	// quando criar a aresta em theMesh2, procuro quem é o ponto em theMesh2 (em CommonVertex) e o uso.

	// iniciando criação das arestas
	map<string,pVertex>::iterator itMapVert;
	set<pEdge>::iterator itSEdges;
	for ( itSEdges=AdaptEdges.begin() ; itSEdges != AdaptEdges.end(); itSEdges++ ){

		pGEntity ent = E_whatIn(*itSEdges);
		pVertex vtx1 = E_vertex(*itSEdges, 0);
		pVertex vtx2 = E_vertex(*itSEdges, 1);

		double xyz[3];

		// criando chave para o no 1
		V_coord(vtx1, xyz);
		string ponto = "."; 
		ostringstream oss;
		oss << xyz[0];
		oss << ponto;
		oss << xyz[1];
		string XY = oss.str();


		itMapVert = CommonVertex.find(XY);
		pVertex Node1 = (*itMapVert).second;

		// criando chave para o no 2
		V_coord(vtx2, xyz);
		ostringstream oss2;
		oss2 << xyz[0];
		oss2 << ponto;
		oss2 << xyz[1];
		XY = oss2.str();

		itMapVert = CommonVertex.find(XY);
		pVertex Node2 = (*itMapVert).second;

		// criando arestas na malha alvo

		pEdge edg = M_createE(theMesh2, Node1, Node2, ent);  // Esse theMesh2 aí ta suspeito...

		// Agora pegando o PS de um dos nós e atribuindo à aresta

		int psreader = 0;
		EN_getDataInt (Node1, MD_lookupMeshDataId("PlaneSurface"), &psreader);
		EN_attachDataInt (edg, MD_lookupMeshDataId("PlaneSurface"), psreader);

		cout << "Plane surface atribuido à aresta: " << psreader << endl;

		//Inserindo a aresta em BoundaryEdges
		BoundaryEdges.insert(edg);

	}

}


void Remover::ResetMRE() {  // Apaga o conteúdos de todos os sets da classe MRE para reinício e criação da malha de adaptação
	// Faz uma cópia de BoundaryNodes e BoundaryEdges em AdaptNodes e AdaptEdges para gerar o arquivo de geometria com eles.

	set<pVertex>::iterator itpVert;
	for ( itpVert=BoundaryNodes.begin() ; itpVert != BoundaryNodes.end(); itpVert++ ){
		AdaptNodes.insert(*itpVert);
	}

	set<pEdge>::iterator itpEdge;
	for ( itpEdge=BoundaryEdges.begin() ; itpEdge != BoundaryEdges.end(); itpEdge++ ){
		AdaptEdges.insert(*itpEdge);
	}


	GeometryBoundNodes.clear(); // acho que este está sem uso... verificar...
	GeomBoundNodes.clear();
	GeomBoundEdges.clear();
	MeshBoundEdges.clear();
	MeshBoundNodes.clear();
	BoundaryNodes.clear();
	BoundaryEdges.clear();
	RemovefromBoundEdges.clear(); // tive que colocar pra nao dar erro na contagem.
	RemovefromBoundNodes.clear();
	RemovefromBoundFaces.clear();
	//BoundaryFaces.clear(); // transformei o antigo set Faces em set BoundaryFaces
	GeomBoundFaces.clear(); // criado o GeomBoundFaces
	//	set<pFace> Facestoremove; //set do estimador de erro
	//	set<pFace> Facestoremove2; //set do estimador de erro
	//	set<pFace> Facestoremove3; //set do estimador de erro
	//	set<pFace> Facestoremove4; //set do estimador de erro
	Tetras.clear(); // Originalmente so este set é pRegion e o resto é pFace.
	argi.clear();


}











