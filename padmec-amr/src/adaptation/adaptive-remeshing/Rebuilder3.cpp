#include "Rebuilder3.h"
#include <string>

#include "AMR.h"

#include "Gmsh.h"
#include "GModel.h"
#include "GRegionGMSH.h"
#include "MEdgeGMSH.h"
#include "MVertexGMSH.h"
#include "MFaceGMSH.h"

#include <cstring>
#include <sstream>



Rebuilder3::Rebuilder3 (){
}
Rebuilder3::~Rebuilder3 (){
}

void Rebuilder3::MakeMesh_2D(string FileName){ // Essa funcao pega a geometria que foi criada dos elementos de contorno e recria a malha

	string gmsh = "./mainSimple ";
	ostringstream os;
	os << gmsh;
	os << FileName;

	std::string Command = os.str();

	//ofstream Myfile(filename.c_str(), ios::app);

	cout << "Iniciando geração da malha a partir do arquivo de geometria - Geometria2D Cl escolhido: " << endl; 
    	//system("./gmsh -2 Geometria2D"); // ALTERAR ISSO PARA MUDAR O NOME DO ARQUIVO DE MALHA PARA MALHA2.MSH, ALGO ASSIM...

    	system(Command.c_str());

	cout << endl;
	cout << "Criado arquivo de malha Geometria2D.msh" << endl;
	cout << endl; 

}

/*void Rebuilder3::MakeMerge_2D(char * arq1,char * arq2,char * arq3){
	string gmsh = "./mainSimple ";
	ostringstream os;
	os << gmsh;
	os <<" " << arq1;
	os <<" " << arq2;
	os <<" " << arq3;

	std::string Command = os.str();

	//ofstream Myfile(filename.c_str(), ios::app);

	cout << "Iniciando geração da malha a partir do arquivo de geometria - Geometria2D Cl escolhido: " << endl; 
	//system("gmsh -2 Geometria2D"); // ALTERAR ISSO PARA MUDAR O NOME DO ARQUIVO DE MALHA PARA MALHA2.MSH, ALGO ASSIM...

	system(Command.c_str());
}*/

void Rebuilder3::MakeMerge_2D(char * Arquivo1,char * Arquivo2, char * ArquivoSaida){
	GmshSetOption("General", "Verbosity", 99.);
	GmshSetOption("Mesh", "Algorithm", 5.);
	GmshSetOption("General", "Terminal", 1.);
	GModel *m = new GModel();
	m->readMSH(Arquivo1);//M1.msh
	cout << "Argumentos 1 " << Arquivo1 << endl;
	cout << "Argumentos 2 " << Arquivo2 << endl;
	GmshMergeFile(Arquivo2);//M2.msh
	m->removeInvisibleElements();
	m->removeDuplicateMeshVertices(0.0000001);
	m->writeMSH(ArquivoSaida, 1.0); //test.msh
	delete m;
}


void Rebuilder3::RefineMesh_2D(){ // Refina o arquivo de malha Geometria2d.msh
	cout << "Refinando malha Geometria2D.msh" << endl;
	system("gmsh -refine Geometria2D.msh");
	cout << "concluido refinamento da malha Geometria2D.msh" << endl;
	cout << endl;
}


void Rebuilder3::MakeAdaptMesh_2D(){ // Essa funcao pega a geometria que foi criada dos elementos de contorno e recria a malha
	cout << "Iniciando geração da malha a partir do arquivo de geometria - Geometria2D" << endl;
	system("gmsh -2 AdaptGeometria2D"); // ALTERAR ISSO PARA MUDAR O NOME DO ARQUIVO DE MALHA PARA MALHA2.MSH, ALGO ASSIM...
	cout << endl;
	cout << "Criado arquivo de malha Geometria2D.msh" << endl;
	cout << endl;
}
