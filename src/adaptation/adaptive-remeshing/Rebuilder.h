//Rebuilder class

#include "includes.h"
#include <string>

class Rebuilder
{
public:
	// Contructor and Destructor

	int MkMesh;
	pMesh theMesh;

	Rebuilder ();
	~Rebuilder();

	void MakeMesh_2D (string FileName);
	void MakeAdaptMesh_2D ();
	//set <int> ReturnListofElements();
	//int Ident;
	void RefineMesh_2D();


	
	
private:

	//set <int> ListofElements;
	//set <pFace> BoundaryFaces;
	//set <pEdge> BoundaryEdges;
	//set <pVertex> BoundaryNodes;

};

