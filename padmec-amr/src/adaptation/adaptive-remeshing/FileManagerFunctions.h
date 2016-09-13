//FileManager class

#include "includes.h"
#include <string>

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

class FileManagerFunctions
{
public:
	// Contructor and Destructor
	FileManagerFunctions ();
	~FileManagerFunctions();
	void Get_Angular_Coef(set<pEdge>);
	//void CreateLineLoop (set<pEdge>, int, string, pMesh, map<int, GVertexGMSH*>, GModel*  ); 
	void CreateLineLoop (set<pEdge>, int, string, pMesh);
	void CreateLineLoop_3D (set<pEdge>, int, string, pMesh); 
	set<pEdge> Separateset  (set<pEdge>, int);
	void CreatePlaneSurfMatrix (int,int);
	string Create_Key_Point(pVertex); 
	set<pEdge> Merge_Edges (set<pEdge>, pMesh, map<int,pEdge>, int, string);

private:

//	set <int> ListofElements;
//	set <pFace> BoundaryFaces;
//	set <pEdge> BoundaryEdges;
//	set <pVertex> BoundaryNodes;

};

