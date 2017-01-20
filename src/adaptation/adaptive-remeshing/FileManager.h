//FileManager class

#include "includes.h"
#include <string>

class FileManager
{
public:
	// Contructor and Destructor

	pMesh theMesh;
	pMesh theMesh2;


	FileManager ();
	FileManager (pMesh theMesh, pMesh theMesh2);
	FileManager(pMesh theMesh, set<pFace> BoundaryFaces, set<pEdge> BoundaryEdges, set<pVertex> BoundaryNodes, int PSCounter);
	FileManager(pMesh theMesh, set<pEdge> BoundaryEdges, set<pVertex> BoundaryNodes, int PSCounter);
	~FileManager();

	void Clear_Containers(); // Reseta as variaveis da classe filemanager
	void SaveFile1_2D ();  // Essas funcoes gravam o arquivo em formato NNF
	void SaveFile2_2D ();
	void SaveFile1_3D ();
	void SaveFile2_3D ();
	void SaveFileMshRenumber (pMesh, string, double);
	void SaveFileMshRenumber2 (pMesh, string, double);
	void Merge_msh (pMesh, pMesh);  // Esta funcao une os arquivos de malha
	void SaveFileGeometry_2D(float Cl, pMesh theMesh, int PSCounter);  // Esta funcao Salva o arquivo de geometria necessario para o remeshing
	void SaveFileGeometry_3D(float Cl);  // Esta funcao Salva o arquivo de geometria necessario para o remeshing
	void SaveFile_AdaptGeometry_2D(string, float GE);  // Esta aqui salva a geometria da malha de adaptacao
	void FileReader (char *IDS);

	set <int> ReturnListofElements();
	int Ident;
	int PSCounter;


private:

	map <string,pVertex> BoundVertex; // mapa dos vertices da malha original necessário para merge_msh encontrar os BoundVertex e comparar com os pontos que ela está criando.
	set <int> ListofElements;
	set <pFace> BoundaryFaces;
	set <pEdge> BoundaryEdges;
	set <pVertex> BoundaryNodes;

};

