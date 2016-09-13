#include "Rebuilder2.h"
#include <string>

Rebuilder2::Rebuilder2 (){
}

Rebuilder2::~Rebuilder2 () {
}

void Rebuilder2::MakeMesh_2D(char * ArquivoEntrada, int dimensao){ // Essa funcao pega a geometria que foi criada dos elementos de contorno e recria a malha

	const double lcMin = 0.0000001; //.e-2*0.1;
	GmshSetOption("Mesh","CharacteristicLengthMin",lcMin);
	GmshSetOption("General", "Terminal", 1.);
	GmshSetOption("Mesh", "Algorithm", 2.);

        char str[80];
        strcpy (str,ArquivoEntrada);
        strcat (str,".msh");
        cout << "Iniciando geração da malha a partir do arquivo de geometria - Geometria2D" << endl;
	GModel *m = new GModel();
        cout << "Criando arquivo de malha " << str << endl;
        m->readGEO(ArquivoEntrada);
        m->mesh(dimensao);//DIMENSAO DA MALHA
	m->writeMSH(str, 1.0); //Parâmetros: Nome do arquivo de saída, versão do arquivo
        cout << "Criacao de malha " << str << " completa" << endl;
        delete m;
	

}

void Rebuilder2::MakeMerge_2D(char * Arquivo1,char * Arquivo2, char * ArquivoSaida){

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
