//Rebuilder class

#include "includes.h"
#include <string>

class Rebuilder3
{
public:
	// Contructor and Destructor

	pMesh theMesh;

	Rebuilder3 ();
	~Rebuilder3();

	void MakeMesh_2D (string FileName);
	void MakeAdaptMesh_2D ();
	//set <int> ReturnListofElements();
	//int Ident;
	void RefineMesh_2D();
	
	void MakeMerge_2D(char * , char * ,char * );

private:

	//set <int> ListofElements;
	//set <pFace> BoundaryFaces;
	//set <pEdge> BoundaryEdges;
	//set <pVertex> BoundaryNodes;

};

