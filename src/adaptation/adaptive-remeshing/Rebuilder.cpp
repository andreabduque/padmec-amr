#include "Rebuilder.h"
#include <string>

//#include <fstream>
//#include <iostream>
//#include <string>

#include <cstring>
#include <sstream>



Rebuilder::Rebuilder (){
}

void Rebuilder::MakeMesh_2D(string FileName){ // Essa funcao pega a geometria que foi criada dos elementos de contorno e recria a malha

	string gmsh = "gmsh -2 ";
	//string gmsh = "./mainSimple ";
	ostringstream os;
	os << gmsh;
	os << FileName;

	std::string Command = os.str();

	//ofstream Myfile(filename.c_str(), ios::app);

	cout << "Iniciando geração da malha a partir do arquivo de geometria - Geometria2D Cl escolhido: " << endl; 
	//system("gmsh -2 Geometria2D"); // ALTERAR ISSO PARA MUDAR O NOME DO ARQUIVO DE MALHA PARA MALHA2.MSH, ALGO ASSIM...

	system(Command.c_str());

	MkMesh=1;
	cout << endl;
	cout << "Criado arquivo de malha Geometria2D.msh" << endl;
	cout << endl; 

	// O problema é que ele grava no formato 2.2 e nao no formato 1 que o M_Load lê, acho que vou me lascar fazendo um conversor pra isso.

	//Deve dar menos trabalho do que criar um modulo pra gravar em iges e depois ver como é que vou fazer pra abrir no open cascade pra gravar um formato de malha que eu nem sei ainda qual seria...

	// Tem a Library do netgen Nglib que vale a pena dar uma olhada...

}




void Rebuilder::RefineMesh_2D(){ // Refina o arquivo de malha Geometria2d.msh
	cout << "Refinando malha Geometria2D.msh" << endl;
	system("gmsh -refine Geometria2D.msh");
	cout << "concluido refinamento da malha Geometria2D.msh" << endl;
	cout << endl;
}


void Rebuilder::MakeAdaptMesh_2D(){ // Essa funcao pega a geometria que foi criada dos elementos de contorno e recria a malha
	cout << "Iniciando geração da malha a partir do arquivo de geometria - Geometria2D" << endl;
	system("gmsh -2 AdaptGeometria2D"); // ALTERAR ISSO PARA MUDAR O NOME DO ARQUIVO DE MALHA PARA MALHA2.MSH, ALGO ASSIM...
	MkMesh=1;
	cout << endl;
	cout << "Criado arquivo de malha Geometria2D.msh" << endl;
	cout << endl;
}

