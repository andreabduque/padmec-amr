#include "FileManager.h"
#include <fstream>
#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <iostream>
#include "FileManagerFunctions.h"


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


using std::ofstream;
using std::ifstream;
using std::ios;

using std::endl;
using std::cout;


FileManager::FileManager (pMesh theMesh, set <pFace> BoundaryFaces, set <pEdge> BoundaryEdges, set <pVertex> BoundaryNodes,int PSCounter){
	this->theMesh=theMesh;
	this->BoundaryFaces=BoundaryFaces;
	this->BoundaryEdges=BoundaryEdges;
	this->BoundaryNodes=BoundaryNodes;
	this->PSCounter=PSCounter;
}

FileManager::FileManager (pMesh theMesh, set <pEdge> BoundaryEdges, set <pVertex> BoundaryNodes, int PSCounter){
	this->theMesh=theMesh;
	this->BoundaryEdges=BoundaryEdges;
	this->BoundaryNodes=BoundaryNodes;
	this->PSCounter=PSCounter;
}

FileManager::FileManager (){
}

FileManager::~FileManager (){
	Clear_Containers();
}

void FileManager::FileReader(char *IDS){ // Esta funcao le o arquivo de entrada, e só... atribui os elementos a um set de inteiros chamado ListofElements... A funcao iterador vai usar esse set para encontrar pelos ids os elementos na malha e coloca-los em um set de faces ou tetras a depender do caso, outras funcoes irão removê-los da malha...
	ifstream infile;
	infile.open (IDS, ifstream::in);
	infile.seekg (0, ifstream::beg);
	string line;

	while (! infile.eof() ){
		getline (infile,line);
		if (!line.empty()){
			ListofElements.insert(atoi(line.c_str()));
		}
	}
	infile.close();
}

set<int> FileManager::ReturnListofElements(){
	return ListofElements;
}


//          //          Funcao SaveFileGeometry_2D         //           //
// Funcao para salvar arquivo de Geometria .geo


void FileManager::SaveFileGeometry_2D(float Cl, pMesh theMesh, int PSCounter){ 
	// Salva um ARQUIVO DE GEOMETRIA padrao do gmesh DO CONTORNO DO BURACO para gerar a malha com o cl escolhido   
	//General.Verbosity
	//Level of information printed during processing (0=no information)
	//Default value: 5
	//Saved in: General.OptionsFileName

	GmshSetOption("General", "Verbosity", 99.);
	GmshSetOption("General", "DetachedMenu", 0.);
	GmshSetOption("General", "Terminal", 1.);
	GmshSetOption("Mesh", "Algorithm", 5.);
	
	//Funções:
	//Recombine Surface {8};
	//Mesh.RecombineAll=1;
	//GmshSetOption("Mesh", "RecombineAll", 1.);//usar quadrangulos
	
	//const double lcMin = 1.e-3*0.01;
	const double lcMin = 1.e-2*0.1;
	const double lcMax = 10.0*0.01;
	GmshSetOption("Mesh","CharacteristicLengthMin",lcMin);
	//  GmshSetOption("Mesh","CharacteristicLengthMin",lcMax);
	GModel *a = new GModel();
	int NewEdgeID = 0;

	system("rm Geometria2D");
	system("rm Geometria2D.msh");

	cout << endl;
	cout << "Salvando geometria 2D - Geometria2D" << endl;
	ofstream Myfile("./Geometria2D", ios::app);

	if (!Myfile) {
		Myfile.open("./Geometria2D");
	}

	// ROTINA QUE CRIA A LISTA DE PONTOS INICIAL DO ARQUIVO - IDS E COORDENADAS DE TODOS OS PONTOS
	Myfile << "cl1 = " << Cl << ";" << endl;
	int i = 0;

	//É, não tá pegando porque tu fez um map de int pra GVertex.
	//Tem que ser:
	//map <int, GVertexGMSH*> BoundaryPoints;
	//E no insert:
	//BoundaryPoints.insert ( std::pair<int, GVertexGMSH>(EN_id(V), &arra[2] ) );
	//Perceba que agora to passando o endereço de arra[2], que é um vetor de GVertex.
	map <int, GVertexGMSH*> BoundaryPoints;
	BoundaryNodes.clear();
	set<pEdge>::iterator its;
	for ( its=BoundaryEdges.begin() ; its != BoundaryEdges.end(); its++ ){
		for(int y=0; y<2; y++){
			pVertex V = E_vertex(*its, y);
			double xyz[3];
			V_coord(V, xyz);
			BoundaryNodes.insert(V);

			BoundaryPoints.insert ( std::pair<int, GVertexGMSH*>(EN_id(V), a->addVertex(xyz[0], xyz[1], 0.0, Cl)) );
		}
	}
	
	set<pVertex>::iterator itnodes;
	for ( itnodes=BoundaryNodes.begin() ; itnodes != BoundaryNodes.end(); itnodes++ ){
		if (V_numEdges(*itnodes)>1){
			i=i+1;
			double xyz[3];
			V_coord(*itnodes, xyz);
			Myfile<<"Point("<<EN_id(*itnodes)<<")"<<" "<<"="<<" "<<"{"<<xyz[0]<<","<<" "<<xyz[1]<<","<<" "<<xyz[2]<<","<<" 				"<<Cl<<"};"<<endl;
		}
	}
	
	//map <int, GVertexGMSH>:: iterator itPoints;
	//for ( itPoints=BoundaryPoints.begin() ; itPoints != BoundaryPoints.end(); itPoints++ ){	
	//}

	// FIM DA ROTINA, PARTE INICIAL DO ARQUIVO GRAVADO


	//		//		//		//		//		//		//		//		//


	// ROTINA QUE SALVA NO ARQUIVO TODAS AS LINES E IDS DOS PONTOS QUE AS COMPOEM AINDA DENTRO DA FUNCAO SaveFileGeometry_2D
	
	set<pEdge> BEdges;
	set<pEdge>::iterator itEds;
	pVertex Point1;
	pVertex Point2;
	BEdges.clear();

	for (itEds=BoundaryEdges.begin(); itEds != BoundaryEdges.end(); itEds++){
		BEdges.insert(*itEds);
	}

	i=0;
	int QTPS=0;
	int FGH =0;

	set<pEdge>::iterator edgit;
	for(edgit=BEdges.begin();edgit!=BEdges.end();edgit++){ // lendo o map e gravando no arquivo os dados.

		EN_getDataInt (*edgit, MD_lookupMeshDataId("PlaneSurface"), &FGH);
		if (FGH>QTPS){
			QTPS = FGH;
		}

		pEdge Edgex = *edgit;
		Point1 = E_vertex(Edgex, 0);
		Point2 = E_vertex(Edgex, 1);

		int IDPoint1 = EN_id(Point1);
		int IDPoint2 = EN_id(Point2);

		if (NewEdgeID<EN_id(Edgex)){
			NewEdgeID=EN_id(Edgex);
		}
	}

	//FIM DA ROTINA QUE SALVA NO ARQUIVO A LISTA DE ARESTAS e cria a copia de BoundaryEdges e o map
	FileManagerFunctions *F = new FileManagerFunctions();
	int w=1;
	int c=BEdges.size();
	int PSurf;
	set<pEdge> ArraySet[PSCounter+1]; // Array de sets as arestas de cada plane surface preenchem um set do array
	for(itEds=BEdges.begin();itEds!=BEdges.end() ;itEds++){ // Percorre as BEdges a coloca cada aresta em um set específico correspondente a um Plane Surface.
		int PlaneSurfReader = 0;
		EN_getDataInt (*itEds, MD_lookupMeshDataId("PlaneSurface"), &PlaneSurfReader);
		if(PlaneSurfReader!=0){
			ArraySet[PlaneSurfReader].insert(*itEds);
		}
		else{
			cout<<"Tem aresta com PS zero"<< endl;
		}
	}

	int LL;
	int P=1;
	int PS=1;  // OLHA ATENCAO AQUI, NAO ESQUECA DE COLOCAR ISSO PRA UM DE NOVO
	while (P<=PSCounter){

		P++;
		w=1;
		set <int> lineloops;
		while (w==1){

			c++;
			set<pEdge> CommonEdges;

			set<int> nwEdges;

			CommonEdges=F->Separateset(ArraySet[PS], PS); // Separateset retornou o set commonEdges para ser usado por 	CreatePlaneSurface, retornou UM lineloop daquele plane surface dentro de CommonEdges


			// Iniciando a funcao de merge edges


			for(itEds=CommonEdges.begin();itEds!=CommonEdges.end() ;itEds++){ // Remove de PSEdges os elementos de CommonEdges
				ArraySet[PS].erase(*itEds);
			}

			if (ArraySet[PS].size()==0){ // quando ArraySet[PS] zerar, termina esta rotina e o flag w é setado para zero...
				w=0;
			}
			F->Get_Angular_Coef(CommonEdges);
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




				pEdge Edge1 = *GeomBdryEdges.begin();  //PEGA UMA ARESTA

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
							}
						}
						for(int i = 0 ; i < Numb_Edges2 ; i++){
							pEdge Edg2 = V_edge(Vert2, i);
							if (Edg2!=Old_Edg2 && SameCoefEdges.count(Edg2)==1){
								Vert2 = E_otherVertex(Edg2, Vert2);
								Old_Edg2 = Edg2;
								SameEdge.insert(Edg2);
								v2 = 1;
							}
						}

						if (v1==0 && v2==0){
							V=0;
							// criando a nova aresta aqui dentro deste if e inseri-la no set newEdges e no set commonedges
							NewEdgeID++;
							pGEntity ent = E_whatIn(Old_Edg2);
							pEdge edg = M_createE(theMesh, Vert1, Vert2, ent);
							EN_setID(edg, NewEdgeID);
							newEdges.insert(edg);
							CommonEdges.insert(edg);

							//		cout << "Criada mais uma aresta e acrescentada a CommonEdges" << endl;
							set<pEdge>:: iterator tytet;
							for(tytet=SameEdge.begin();tytet!=SameEdge.end(); tytet++){
								CommonEdges.erase(*tytet);
							}

							//cout<< "Nova aresta criada" << endl;
						}
					}

					//terminou de identificar

					for(itEds=SameEdge.begin(); itEds!=SameEdge.end(); itEds++){
						GeomBdryEdges.erase(*itEds);
						SameCoefEdges.erase(*itEds);
						//theMesh->DEL(*itEds); // Coloquei agora...
						//CommonEdges.erase(*itEds);
					}
					SameEdge.clear();

				}
			}

			GeomBdryEdges.clear();
			SameCoefEdges.clear();
			int dim = CommonEdges.size();
			//cout << "	Agora criando o lineloop, CommonEdges.size() " << CommonEdges.size() << endl; 
			LL=c;
			//F->CreateLineLoop(CommonEdges, c,"Geometria2D",theMesh, BoundaryPoints, a);  //matriz1 criada passa-la para a rotina de reordenacão para depois criar a curva. CreateLineLoop reordena estes edges 
			F->CreateLineLoop(CommonEdges, c,"Geometria2D",theMesh);  //matriz1 criada passa-la para a rotina de reordenacão para depois criar a curva. CreateLineLoop reordena estes edges 
			//como um lineloop - LL é um int ( o id do lineloop)

			set<pEdge>:: iterator nwedgs;
			for (nwedgs=SameEdge.begin();nwedgs!=SameEdge.end();nwedgs++){
				ArraySet[PS].erase(*nwedgs);
			}
			for (nwedgs=newEdges.begin(); nwedgs!=newEdges.end();nwedgs++){

				ArraySet[PS].erase(*nwedgs);

			}

			newEdges.clear();
			SameEdge.clear();

			//	cout << "Os line loops sao: " << LL << " ao total" << endl;

			lineloops.insert(LL);
			CommonEdges.clear();


		}
		Myfile << "Plane Surface(" << PS << ") = {"; // gravação no arquivo dos Plane surfaces

		set<int>::iterator itints;
		itints=lineloops.begin();
		Myfile << *itints;
		for (itints=lineloops.begin();itints!=lineloops.end();itints++){
			if (itints!=lineloops.begin()){
				Myfile << ", " << *itints;
			}
		}
		Myfile << "};" << endl;
		lineloops.clear();
		ArraySet[PS].clear(); // Acho que nao precisa mais disso
		PS++;
	}
	Myfile << "nan = 0.5;" << endl;
	Myfile << "Merge \'View2D.pos\';" << endl; // gravação no arquivo dos Plane surfaces
	Myfile.close(); // fim da gravação da geometria
	cout << "fim da gravacao da geometria da malha2 2D" << endl;
	cout << endl;
	delete F; F = 0;
	//delete a; a = 0;
}

// Funcoes para salvar os arquivos NNF em 2D e 3D

void FileManager::SaveFile1_2D(){ // As funcoes SaveFile criam os arquivos de saida, sao duas para o caso de 2D e duas para o caso de 3D ESSA AQUI SALVA O ARQUIVO ORIGINAL SEM OS ELEMENTOS QUE FORAM REMOVIDOS, A Savefile2_2D Salva em formato de malha apenas o contorno do buraco.
	cout << endl;
	cout << "Salvando malha 2D - Malha1_2D.dat" << endl;
	ofstream Myfile("./Malha1_2D.dat", ios::app);
	if (!Myfile) {
		Myfile.open("./Malha1_2D.dat");
	}

	Myfile << "%HEADER" << endl;
	Myfile << "File generated by UFPE" << endl;
	Myfile << endl;
	Myfile << "%NODE" << endl;
	Myfile << M_numVertices(theMesh) << endl;
	Myfile << endl;
	Myfile << "%NODE.COORD" << endl;
	Myfile << M_numVertices(theMesh) << endl;

	VIter vit = M_vertexIter(theMesh);
	while (pVertex Vert = VIter_next(vit)){

		int id = EN_id(Vert);

		double xyz[3];

		V_coord(Vert, xyz);

		Myfile << id << " " << xyz[0] << " " << xyz[1] << " " << xyz[2] << endl;  // saida de verificacao
	}

	VIter_delete(vit);

	Myfile << endl;
	Myfile << "%INTEGRATION.ORDER" << endl;
	Myfile << "1" << endl;
	Myfile << "1\t1\t1\t1\t1\t1\t1" << endl;
	Myfile << endl;
	Myfile << "%ELEMENT" << endl;
	Myfile << M_numFaces(theMesh) << endl;
	Myfile << endl;
	Myfile << "%ELEMENT.T3" << endl;
	Myfile << M_numFaces(theMesh) << endl;

	FIter fit = M_faceIter(theMesh);
	while (pFace FacE = FIter_next(fit)){

		int ID = EN_id(FacE);
		int Point1 = EN_id(F_vertex(FacE, 0));
		int Point2 = EN_id(F_vertex(FacE, 1));
		int Point3 = EN_id(F_vertex(FacE, 2));

		Myfile << ID << " " << Point1 << " " << Point2 << " " << Point3 << endl;
	}

	FIter_delete(fit);

	Myfile <<  endl;
	Myfile << "%END" << endl;

	Myfile.close(); // fim da gravação da malha1.
	cout << "fim da gravacao malha1_2D" << endl;
	cout << endl;
}


// funcao para salvar a malha2 essa aqui é pra salvar a MALHA DE CONTORNO Nesta funcao aparecem as arestas, corrigir o problema dos ids zerados

void FileManager::SaveFile2_2D(){  // Salva arquivo de contorno de malha NNF 2D

	cout << "Salvando malha 2D - Malha2_2D.dat" << endl;

	ofstream Myfile("./Malha2_2D.dat", ios::app);

	if (!Myfile) {
		Myfile.open("./Malha2_2D.dat");
	}

	Myfile << "%HEADER" << endl;
	Myfile << "File generated by UFPE" << endl;
	Myfile << endl;
	Myfile << "%NODE" << endl;
	Myfile << BoundaryNodes.size() << endl;
	Myfile << endl;
	Myfile << "%NODE.COORD" << endl;
	Myfile << BoundaryNodes.size() << endl;

	set<pVertex>::iterator itnodes;

	for ( itnodes=BoundaryNodes.begin() ; itnodes != BoundaryNodes.end(); itnodes++ ){
		int id = EN_id(*itnodes);
		double xyz[3];
		V_coord(*itnodes, xyz);
		Myfile << id << " " << xyz[0] << " " << xyz[1] << " " << xyz[2] << endl;
	}

	Myfile << endl;
	Myfile << "%INTEGRATION.ORDER" << endl;
	Myfile << "1" << endl;
	Myfile << "1\t1\t1\t1\t1\t1\t1" << endl;
	Myfile << endl;
	Myfile << "%ELEMENT" << endl;
	Myfile << BoundaryEdges.size() << endl;
	Myfile << endl;
	Myfile << "%EDGE" << endl;
	Myfile << BoundaryEdges.size() << endl;

	set<pEdge>::iterator itEdges;
	for ( itEdges=BoundaryEdges.begin() ; itEdges != BoundaryEdges.end(); itEdges++ ){

		int ID = EN_id(*itEdges);
		int Point1 = EN_id(E_vertex(*itEdges, 0));
		int Point2 = EN_id(E_vertex(*itEdges, 1));

		Myfile << ID << " " << Point1 << " " << Point2 << endl;

	}

	Myfile <<  endl;
	Myfile << "%END" << endl;
	cout << "fim da gravacao malha2_2D" << endl;
	Myfile.close();// fim da gravação.
}



void FileManager::SaveFile1_3D(){  // Salva a malha original sem os elementos removidos e em formato NNF 3D
	cout << "Salvando Malha 3D - Malha1_3D.dat" << endl;

	ofstream Myfile("./Malha1_3D.dat", ios::app);

	if (!Myfile) {
		Myfile.open("./Malha1_3D.dat");
	}


	Myfile << "%HEADER" << endl;
	Myfile << "File generated by UFPE" << endl;
	Myfile << endl;
	Myfile << "%NODE" << endl;
	Myfile << M_numVertices(theMesh) << endl;
	Myfile << endl;
	Myfile << "%NODE.COORD" << endl;
	Myfile << M_numVertices(theMesh) << endl;

	VIter vit = M_vertexIter(theMesh);
	while (pVertex Vert = VIter_next(vit)){

		int id = EN_id(Vert);
		double xyz[3];
		V_coord(Vert, xyz);
		Myfile << id << " " << xyz[0] << " " << xyz[1] << " " << xyz[2] << endl;  // saida de verificacao
	}

	VIter_delete(vit);

	Myfile << endl;
	Myfile << "%INTEGRATION.ORDER" << endl;
	Myfile << "1" << endl;
	Myfile << "1\t1\t1\t1\t1\t1\t1" << endl;
	Myfile << endl;
	Myfile << "%ELEMENT" << endl;
	Myfile << M_numRegions(theMesh) << endl; // Aqui comeca as alteracoes
	Myfile << endl;
	Myfile << "%ELEMENT.TETR4" << endl;
	Myfile << M_numRegions(theMesh) << endl;


	// Tambem nao estou gostando desse laco aqui, se eu conseguir uma funcao que me retorne os pontos associados a um tetraedro eu evito esses loops
	RIter rit = M_regionIter(theMesh);
	while (pRegion RegioN = RIter_next(rit)){
		pFace F1 = R_face(RegioN, 0);
		pFace F2 = R_face(RegioN, 1);
		set <pVertex> Points;
		int I=0;
		for (I=0; I<=2; I++){
			Points.insert(F_vertex(F1, I));
			Points.insert(F_vertex(F2, I));
		}
		set<pVertex>::iterator itPoints;
		Myfile << EN_id(RegioN);
		for (itPoints = Points.begin(); itPoints != Points.end(); itPoints++){
			Myfile << " " << EN_id(*itPoints);
		}
		Myfile << endl;
	}

	RIter_delete(rit);

	Myfile <<  endl;
	Myfile << "%END" << endl;

	cout << "fim da gravacao Malha1_3D.dat" << endl;
}


void FileManager::SaveFile2_3D(){  // Salva o arquivo de contorno de malha em NNF 3D
	cout << "Salvando malha 3D - Malha2_3D.dat" << endl;

	ofstream Myfile("./Malha2_3D.dat", ios::app);

	if (!Myfile) {
		Myfile.open("./Malha2_3D.dat");
	}

	Myfile << "%HEADER" << endl;
	Myfile << "File generated by UFPE" << endl;
	Myfile << endl;
	Myfile << "%NODE" << endl;
	Myfile << BoundaryNodes.size() << endl;
	Myfile << endl;
	Myfile << "%NODE.COORD" << endl;
	Myfile << BoundaryNodes.size() << endl;

	set<pVertex>::iterator itnodes;

	for ( itnodes=BoundaryNodes.begin() ; itnodes != BoundaryNodes.end(); itnodes++ ){
		int id = EN_id(*itnodes);
		double xyz[3];
		V_coord(*itnodes, xyz);
		Myfile << id << " " << xyz[0] << " " << xyz[1] << " " << xyz[2] << endl;
	}

	Myfile << endl;
	Myfile << "%INTEGRATION.ORDER" << endl;
	Myfile << "1" << endl;
	Myfile << "1\t1\t1\t1\t1\t1\t1" << endl;
	Myfile << endl;
	Myfile << "%ELEMENT" << endl;
	Myfile << BoundaryFaces.size() << endl;
	Myfile << endl;
	Myfile << "%ELEMENT.T3" << endl;
	Myfile << BoundaryFaces.size() << endl;

	set<pEdge>::iterator itFaces;
	for ( itFaces=BoundaryFaces.begin() ; itFaces != BoundaryFaces.end(); itFaces++ ){

		int ID = EN_id(*itFaces);
		int Point1 = EN_id(F_vertex(*itFaces, 0));
		int Point2 = EN_id(F_vertex(*itFaces, 1));
		int Point3 = EN_id(F_vertex(*itFaces, 2));

		Myfile << ID << " " << Point1 << " " << Point2 << " " << Point3 << endl;

	}

	Myfile <<  endl;
	Myfile << "%END" << endl;

	cout << "fim da gravacao Malha2_3D.dat" << endl;
}

void FileManager::Clear_Containers(){
	ListofElements.clear();
	BoundaryFaces.clear();
	BoundaryEdges.clear();
	BoundaryNodes.clear();
	BoundVertex.clear();
	Ident = 0;
	PSCounter = 0;
}

void FileManager::SaveFileMshRenumber(pMesh theMesh, string Filename, double newID){

	string pontobarra = "./";
	ostringstream os;
	os << pontobarra;
	os << Filename;

	std::string file_name = os.str();

	string rm = "rm ";
	ostringstream rms;
	rms << rm;
	rms << Filename;

	std::string Command_rm = rms.str();

	system(Command_rm.c_str());

	cout << "Salvando malha 2D - " << file_name << endl;

	ofstream Myfile(file_name.c_str(), ios::app);

	if (!Myfile) {
		Myfile.open(file_name.c_str());
	}

	int TotalVertex = M_numVertices(theMesh);

	Myfile << "$NOD" << endl;
	Myfile << TotalVertex << endl;


	// Criando os vertes no arquivo de malha
	//double newID = 0;
	VIter itVertex = M_vertexIter(theMesh);
	while (pVertex VT = VIter_next(itVertex)){
		newID++;
		EN_attachDataInt (VT, MD_lookupMeshDataId("RenumberID"), newID); // Atribuindo um Id renumerado, em ordem para todos os vertex
		double xyz[3];
		V_coord(VT, xyz);
		Myfile << newID << " " << xyz[0] << " " << xyz[1] << " " << xyz[2] << endl;
	}
	VIter_delete(itVertex);


	int TotalFaces = M_numFaces(theMesh);
	Myfile << "$ENDNOD" << endl;
	Myfile << "$ELM" << endl;
	Myfile << TotalFaces << endl;

	FIter itFace = M_faceIter(theMesh);
	while (pFace FC = FIter_next(itFace)){
		newID++;

		pVertex V1 = F_vertex(FC, 0);
		int ID1 = 0;
		EN_getDataInt (V1, MD_lookupMeshDataId("RenumberID"), &ID1);

		pVertex V2 = F_vertex(FC, 1);
		int ID2 = 0;
		EN_getDataInt (V2, MD_lookupMeshDataId("RenumberID"), &ID2);

		pVertex V3 = F_vertex(FC, 2);
		int ID3 = 0;
		EN_getDataInt (V3, MD_lookupMeshDataId("RenumberID"), &ID3);

		Myfile << newID << " 2 0 "<<newID<<" 3 " << ID1 << " " << ID2 << " " << ID3 << endl;

	}
	FIter_delete(itFace);

	Myfile << "$ENDELM" << endl;
}




void FileManager::SaveFileMshRenumber2(pMesh theMesh, string Filename, double newID){


	newID = newID+1;

	string pontobarra = "./";
	ostringstream os;
	os << pontobarra;
	os << Filename;

	std::string file_name = os.str();

	string rm = "rm ";
	ostringstream rms;
	rms << rm;
	rms << Filename;

	std::string Command_rm = rms.str();

	system(Command_rm.c_str());

	cout << "Salvando malha 2D - " << file_name << endl;

	ofstream Myfile(file_name.c_str(), ios::app);

	if (!Myfile) {
		Myfile.open(file_name.c_str());
	}

	int TotalVertex = M_numVertices(theMesh);

	Myfile << "$NOD" << endl;
	Myfile << TotalVertex << endl;


	// Criando os vertes no arquivo de malha
	//double newID = 0;
	VIter itVertex = M_vertexIter(theMesh);
	while (pVertex VT = VIter_next(itVertex)){
		newID++;
		EN_attachDataInt (VT, MD_lookupMeshDataId("RenumberID"), newID); // Atribuindo um Id renumerado, em ordem para todos os vertex
		double xyz[3];
		V_coord(VT, xyz);
		Myfile << newID << " " << xyz[0] << " " << xyz[1] << " " << xyz[2] << endl;
	}
	VIter_delete(itVertex);


	int TotalFaces = M_numFaces(theMesh);
	Myfile << "$ENDNOD" << endl;
	Myfile << "$ELM" << endl;
	Myfile << TotalFaces << endl;

	FIter itFace = M_faceIter(theMesh);
	while (pFace FC = FIter_next(itFace)){
		newID++;

		pVertex V1 = F_vertex(FC, 0);
		int ID1 = 0;
		EN_getDataInt (V1, MD_lookupMeshDataId("RenumberID"), &ID1);

		pVertex V2 = F_vertex(FC, 1);
		int ID2 = 0;
		EN_getDataInt (V2, MD_lookupMeshDataId("RenumberID"), &ID2);

		pVertex V3 = F_vertex(FC, 2);
		int ID3 = 0;
		EN_getDataInt (V3, MD_lookupMeshDataId("RenumberID"), &ID3);

		Myfile << newID << " 2 0 "<<newID<<" 3 " << ID1 << " " << ID2 << " " << ID3 << endl;

	}
	FIter_delete(itFace);

	Myfile << "$ENDELM" << endl;
}




void FileManager::Merge_msh(pMesh theMesh, pMesh theMesh2){

	FileManagerFunctions *H = new FileManagerFunctions();

	set <pVertex> QTDVIZINHOS;

	map <string, int> NewBoundVertex;

	int tinha = M_numFaces(theMesh);



	EIter edgiter = M_edgeIter(theMesh);
	while (pEdge entity = EIter_next(edgiter))
	{
		if (E_numFaces(entity)==1){
			pVertex Vertex1 = E_vertex(entity, 0);
			pVertex Vertex2 = E_vertex(entity, 1);

			string XY = H->Create_Key_Point(Vertex1);
			string XY2 = H->Create_Key_Point(Vertex2);

			BoundVertex.insert ( pair<string, pVertex>(XY, Vertex1) );
			BoundVertex.insert ( pair<string, pVertex>(XY2, Vertex2) );


		}
	}
	EIter_delete(edgiter);


	//////////////////////////////////////////////
	int newids = 100 * M_numVertices(theMesh);                      // Criar Rotina para definir newids de forma segura uma maneira é renumerar theMesh antes de mergê-la, assim fica int newids = M_numVertices(theMesh)
	cout << endl << " --------------------  Entrou em Merge_msh newids é " << newids << "   ----------------------" << endl << endl;					//////////////////////////////////////////////

	FIter fit = M_faceIter(theMesh2); // Notar que a malha percorrida é a theMesh2 <-- 2!

	while (pFace OrFc = FIter_next(fit)){

		pEdge Edg1 = F_edge(OrFc, 0);
		pEdge Edg2 = F_edge(OrFc, 1);
		pEdge Edg3 = F_edge(OrFc, 2);

		pVertex Vert1 = F_vertex(OrFc, 0);
		pVertex Vert2 = F_vertex(OrFc, 1);
		pVertex Vert3 = F_vertex(OrFc, 2); 

		EN_attachDataInt (Vert1, MD_lookupMeshDataId("BNode"), 0); //Colocando valor 0 em todos os vertices de theMesh2
		EN_attachDataInt (Vert2, MD_lookupMeshDataId("BNode"), 0);
		EN_attachDataInt (Vert3, MD_lookupMeshDataId("BNode"), 0);


		int NNeighEdg1 = E_numFaces(Edg1);

		//cout << "Qtd de vizinhos dessa aresta: " << NNeighEdg1 << endl;

		if (NNeighEdg1==1){
			pVertex Vt1 = E_vertex(Edg1, 0);
			pVertex Vt2 = E_vertex(Edg1, 1);
			EN_attachDataInt (Vt1, MD_lookupMeshDataId("BNode"), 1); //Substituindo o zero por 1 nas arestas de contorno
			EN_attachDataInt (Vt2, MD_lookupMeshDataId("BNode"), 1);


			string XY = H->Create_Key_Point(Vt1);

			string XY2 = H->Create_Key_Point(Vt2);

			NewBoundVertex.insert ( pair<string, int>(XY, EN_id(Vt1) + newids) ); // Esse map atualizará o map BoundVertex ao fim da função
			NewBoundVertex.insert ( pair<string, int>(XY2, EN_id(Vt2) + newids) );


		}

		int NNeighEdg2 = E_numFaces(Edg2);

		if (NNeighEdg2==1){
			pVertex Vt1 = E_vertex(Edg2, 0);
			pVertex Vt2 = E_vertex(Edg2, 1);
			EN_attachDataInt (Vt1, MD_lookupMeshDataId("BNode"), 1);//Substituindo o zero por 1 nas arestas de contorno
			EN_attachDataInt (Vt2, MD_lookupMeshDataId("BNode"), 1);


			//Criando a chave para inserir no novo map que servirá para atualizar BoundVertex

			string XY = H->Create_Key_Point(Vt1);

			string XY2 = H->Create_Key_Point(Vt2);

			NewBoundVertex.insert ( pair<string, int>(XY, EN_id(Vt1) + newids) ); // Esse map atualizará o map BoundVertex ao fim da função
			NewBoundVertex.insert ( pair<string, int>(XY2, EN_id(Vt2) + newids) );



		}
		int NNeighEdg3 = E_numFaces(Edg3);
		//cout << "Qtd de vizinhos dessa aresta: " << NNeighEdg3 << endl;
		if (NNeighEdg3==1){
			pVertex Vt1 = E_vertex(Edg3, 0);
			pVertex Vt2 = E_vertex(Edg3, 1);
			EN_attachDataInt (Vt1, MD_lookupMeshDataId("BNode"), 1);//Substituindo o zero por 1 nas arestas de contorno
			EN_attachDataInt (Vt2, MD_lookupMeshDataId("BNode"), 1);


			string XY = H->Create_Key_Point(Vt1);

			string XY2 = H->Create_Key_Point(Vt2);

			NewBoundVertex.insert ( pair<string, int>(XY, EN_id(Vt1) + newids) ); // Esse map atualizará o map BoundVertex ao fim da função
			NewBoundVertex.insert ( pair<string, int>(XY2, EN_id(Vt2) + newids) );


		}


		// Preparando a criação dos pontos

		pGEntity ent1 = V_whatIn(Vert1);
		pGEntity ent2 = V_whatIn(Vert2);
		pGEntity ent3 = V_whatIn(Vert3);

		double *xyz;
		xyz = (double*)malloc(3*sizeof(double));

		V_coord(Vert1, xyz);

		double *xyz2;
		xyz2 = (double*)malloc(3*sizeof(double));

		V_coord(Vert2, xyz);

		double *xyz3;
		xyz3 = (double*)malloc(3*sizeof(double));

		V_coord(Vert3, xyz);

		double *paramet;
		paramet = (double*)malloc(sizeof(double));
		paramet[0]=0;

		int BdNode = 0;

		map<string,pVertex>::iterator itBNd;

		int OriginalID;

		pVertex vetx;
		pVertex vetx2;
		pVertex vetx3;

		double coords[3];
		EN_getDataInt (Vert1, MD_lookupMeshDataId("BNode"), &BdNode);
		V_coord(Vert1, coords);


		if (BdNode==1) { // Se o vertex é do contorno da geometria, entao deve ser usado o id do ponto equivalente na malha original, procura-se no map BoundVertex, atraves das coordenadas que funcionam como key.


			string XY = H->Create_Key_Point(Vert1);

			itBNd=BoundVertex.find(XY); // Buscando no map o id original

			if (itBNd!=BoundVertex.end()) {
				vetx = (*itBNd).second;
			}
			else {
				BdNode=0;
			}
		}
		if (BdNode==0) { // Se o nó nao for de contorno, então é criado normalmente, com o newids mesmo...
			int IDS = EN_id(Vert1);
			int nids = IDS + newids;
			vetx = M_createVP2(theMesh, coords, paramet, nids, ent1);
		}

		EN_getDataInt (Vert2, MD_lookupMeshDataId("BNode"), &BdNode);
		V_coord(Vert2, coords);

		if (BdNode==1) { // Repetindo tudo para Vert2...

			string XY = H->Create_Key_Point(Vert2);

			itBNd=BoundVertex.find(XY); // Buscando no map o id original

			if (itBNd!=BoundVertex.end()) {
				vetx2 = (*itBNd).second;
			}
			else {
				BdNode=0;
			}

		}

		if (BdNode==0) { // Se o nó nao for de contorno, então é criado normalmente, com o newids mesmo...
			int IDS = EN_id(Vert2);
			int nids = IDS + newids;
			vetx2 = M_createVP2(theMesh, coords, paramet, nids, ent2);
		}

		EN_getDataInt (Vert3, MD_lookupMeshDataId("BNode"), &BdNode);
		V_coord(Vert3, coords);

		if (BdNode==1) { // Se o vertex é do contorno da geometria, entao deve ser usado o id do ponto equivalente na malha original, procura-se no map BoundVertex, atraves das coordenadas que funcionam como key.

			string XY = H->Create_Key_Point(Vert3);

			itBNd=BoundVertex.find(XY); // Buscando no map o id original


			if (itBNd!=BoundVertex.end()) {
				vetx3 = (*itBNd).second;
			}
			else {
				BdNode=0;
			}
		}
		if (BdNode==0) { // Se o nó nao for de contorno, então é criado normalmente, com o newids mesmo...
			int IDS = EN_id(Vert3);
			int nids = IDS + newids;
			vetx3 = M_createVP2(theMesh, coords, paramet, nids, ent3);
		}

		// Criando as arestas agora

		pEdge *edg;
		edg = (pEdge*)malloc(3*sizeof(pEdge));

		edg[0] = M_createE(theMesh, vetx, vetx2, ent1);
		edg[1] = M_createE(theMesh, vetx, vetx3, ent2);
		edg[2] = M_createE(theMesh, vetx2, vetx3, ent3);

		// Criando a face agora

		int *ide;
		ide = (int*)malloc(sizeof(int));
		ide[0]=2020;
		pFace fce = M_createF(theMesh, 3, edg, ide, ent1);

		// setando o id da face

		EN_setID(fce , 100000);

	}

	FIter_delete(fit);

	int TotalCount = 0;
	int NeighCount = 0;
	VIter vit = M_vertexIter(theMesh);
	while (pVertex ver = VIter_next(vit)){
		int ct = V_numEdges(ver);
		if (ct<=1){
			NeighCount++;
		}
		TotalCount++;
	}
	VIter_delete(vit);

	BoundVertex.clear();

	delete H; H = 0;
}



void FileManager::SaveFileGeometry_3D(float Cl){ // Salva um ARQUIVO DE GEOMETRIA padrao do gmesh DO CONTORNO DO BURACO para gerar a malha com o cl escolhido, SUBDIVIDINDO arestas

	// 1 - itera atraves de boundaryfaces para cada face cria um plane surface
	// 2 - Une todos os plane surfaces em um volume


	system("rm Geometria3D");
	system("rm Geometria3D.msh");

	cout << endl;
	cout << "Salvando geometria 3D - Geometria3D" << endl;
	ofstream Myfile("./Geometria3D", ios::app);

	if (!Myfile) {
		Myfile.open("./Geometria3D");
	}

	// ROTINA QUE CRIA A LISTA DE PONTOS INICIAL DO ARQUIVO - IDS E COORDENADAS DE TODOS OS PONTOS

	Myfile << "cl1 = " << Cl << ";" << endl;
	int i = 0;
	set<pVertex>::iterator itnodes;

	for ( itnodes=BoundaryNodes.begin() ; itnodes != BoundaryNodes.end(); itnodes++ ){

		i=i+1;
		double xyz[3];
		V_coord(*itnodes, xyz);
		Myfile<<"Point("<<EN_id(*itnodes)<<")"<<" "<<"="<<" "<<"{"<<xyz[0]<<","<<" "<<xyz[1]<<","<<" "<<xyz[2]<<","<<" 			"<<Cl<<"};"<<endl;
	}

	// FIM DA ROTINA, PARTE INICIAL DO ARQUIVO GRAVADO


	//		//		//		//		//		//		//		//		//




	// ROTINA QUE SALVA NO ARQUIVO TODAS AS LINES E IDS DOS PONTOS QUE AS COMPOEM AINDA DENTRO DA FUNCAO SaveFileGeometry_2D

	set<pEdge>::iterator itEds;
	pVertex Point1;
	pVertex Point2;


	i=0;
	map<int,pEdge> EdgesMap;  // criei um map com as arestas para atribuir ids às arestas... Uso para isso o indice do map
	map<int,pEdge>::iterator mapitedges;

	for(itEds=BoundaryEdges.begin(); itEds!=BoundaryEdges.end() ; itEds++ ){  //inserindo as arestas no map
		i=i+1;
		EdgesMap.insert ( pair<int,pEdge>(i,*itEds) );

	}

	for(mapitedges=EdgesMap.begin();mapitedges!=EdgesMap.end();mapitedges++){ // lendo o map e gravando no arquivo os dados.

		pEdge Edgex = (*mapitedges).second;
		Point1 = E_vertex(Edgex, 0);
		Point2 = E_vertex(Edgex, 1);

		int IDPoint1 = EN_id(Point1);
		int IDPoint2 = EN_id(Point2);

		Myfile<<"Line("<<EN_id(Edgex)<<")"<<" "<<"="<<" "<<"{"<<IDPoint1<<","<<" "<<IDPoint2<<"};"<<endl;

	}
	FileManagerFunctions *F = new FileManagerFunctions();
	int c=0;

	set<pFace>::iterator itFaces;


	for ( itFaces=BoundaryFaces.begin() ; itFaces != BoundaryFaces.end(); itFaces++ ){

		set<pEdge> PSEdges;

		PSEdges.insert(F_edge(*itFaces, 0));
		PSEdges.insert(F_edge(*itFaces, 1));
		PSEdges.insert(F_edge(*itFaces, 2));

		int PSurf;
		int LL;

		c++;

		int dim = PSEdges.size();

		int VFS=0;
		EN_getDataInt (*itFaces, MD_lookupMeshDataId("BoundaryFace"), &VFS);

		if ((F_numRegions(*itFaces)==0 && VFS == 1) || (F_numRegions(*itFaces)==1)){

			F->CreateLineLoop_3D(PSEdges, EN_id(*itFaces),"Geometria3D",theMesh); // troquei c por EN_id(*itFaces)

			LL=EN_id(*itFaces); // troquei c por EN_id(*itFaces)

			PSEdges.clear();

			Myfile << "Plane Surface(" << EN_id(*itFaces) << ") = { " << LL << "};" << endl; // gravação no arquivo dos Plane surfaces o planesurface é identificado pelo id da face, cada face é um plane surface.

		}


	}

	set <pFace> controlfaces;

	for ( itFaces=BoundaryFaces.begin() ; itFaces != BoundaryFaces.end(); itFaces++ ){
		controlfaces.insert(*itFaces);
	}

	int facescount = 0;
	int surfacescount = 0;
	while (controlfaces.size()>0) {

		set <pFace> commonVolumesurfaces;
		set <int> SurfacesCount;
		set <pFace> commonSurfacePlanes;

		facescount++;

		for ( itFaces=controlfaces.begin() ; itFaces != controlfaces.end(); itFaces++ ){
			int CSF = 0;
			EN_getDataInt (*itFaces, MD_lookupMeshDataId("Volume"), &CSF);
			if (CSF == facescount) {
				commonVolumesurfaces.insert(*itFaces);
			}
		}

		for ( itFaces=commonVolumesurfaces.begin() ; itFaces != commonVolumesurfaces.end(); itFaces++ ){
			controlfaces.erase(*itFaces);
		}

		int surfCount = facescount;



		while(commonVolumesurfaces.size()>0){ // Todas as faces sao do mesmo volume.

			surfacescount++;

			SurfacesCount.insert(surfacescount);

			// Colocando rotina que separa as plane surfaces
			commonSurfacePlanes.insert(*commonVolumesurfaces.begin());


			for ( itFaces=commonSurfacePlanes.begin() ; itFaces != commonSurfacePlanes.end(); itFaces++ ){

				pEdge E1 = F_edge(*itFaces , 0);
				pEdge E2 = F_edge(*itFaces , 1);
				pEdge E3 = F_edge(*itFaces , 2);

				int Q1 = E_numFaces(E1);
				int Q2 = E_numFaces(E2);
				int Q3 = E_numFaces(E3);

				for (int t=0; t<Q1;t++){ // Pegando as faces vizinhas da face através da aresta...
					pFace F1 = E_face(E1, t);
					int V = 0;
					EN_getDataInt (F1, MD_lookupMeshDataId("Volume"), &V);
					if (V==facescount){ // Se a face pertence ao volume em estudo, insere no set.
						commonSurfacePlanes.insert(F1);
					}
				}

				for (int t=0; t<Q2;t++){
					pFace F2 = E_face(E2, t);
					int V = 0;
					EN_getDataInt (F2, MD_lookupMeshDataId("Volume"), &V);
					if (V==facescount){
						commonSurfacePlanes.insert(F2);
					}
				}

				for (int t=0; t<Q3;t++){
					pFace F3 = E_face(E3, t);
					int V = 0;
					EN_getDataInt (F3, MD_lookupMeshDataId("Volume"), &V);
					if (V==facescount){
						commonSurfacePlanes.insert(F3);
					}
				}


			}

			// O volume é regido pelo facescount
			// as surface loops sao regidas pelo surfacescount

			Myfile << "Surface Loop(" << surfacescount << ") = { ";
			int VFS=0;

			for ( itFaces=commonSurfacePlanes.begin() ; itFaces != commonSurfacePlanes.end(); itFaces++ ){

				EN_getDataInt (*itFaces, MD_lookupMeshDataId("BoundaryFace"), &VFS);

				if ((F_numRegions(*itFaces)==0 && VFS == 1) || (F_numRegions(*itFaces)==1)){

					if(*itFaces!=*commonSurfacePlanes.begin()){
						Myfile << ", ";
					}

					Myfile << EN_id(*itFaces);
				}
			}
			Myfile << "};" << endl;

			for ( itFaces=commonSurfacePlanes.begin() ; itFaces != commonSurfacePlanes.end(); itFaces++ ){

				commonVolumesurfaces.erase(*itFaces);

			}

		}

		Myfile << "Volume(" << facescount << ") = { ";

		set<int>::iterator its;
		for ( its=SurfacesCount.begin() ; its != SurfacesCount.end(); its++ ){
			if(*its!=*SurfacesCount.begin()){
				Myfile << ", ";
			}
			Myfile << *its;
		}

		Myfile << "};" << endl;

		commonSurfacePlanes.clear();


	}


	Myfile.close(); // fim da gravação da geometria

	cout << "fim da gravacao da geometria da malha2 3D" << endl;
	cout << endl;

	delete F; F = 0;

}

