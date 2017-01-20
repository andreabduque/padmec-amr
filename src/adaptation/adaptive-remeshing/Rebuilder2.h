//Rebuilder class

#include <string>

#include "AMR.h"

#include "Gmsh.h"
#include "GModel.h"
#include "GRegionGMSH.h"
#include "MEdgeGMSH.h"
#include "MVertexGMSH.h"
#include "MFaceGMSH.h"

class Rebuilder2
{
public:
	// Contructor and Destructor
	pMesh theMesh;

	Rebuilder2 ();
	~Rebuilder2();

	void MakeMesh_2D (char * ArquivoEntrada,  int dimensao);
	void MakeMerge_2D(char * Arquivo1,char * Arquivo2, char * ArquivoSaida);

private:

};

