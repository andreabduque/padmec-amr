#include "ErrorAnalysis.h"
#include <string>


class ErrorEstimator
{
public:
	// Contructor and Destructor
	ErrorEstimator(pMesh theMesh);
	~ErrorEstimator();


	map<int,float> VCL;
	//set<pVertex> AllNodes;
	//set<pFace> ListedFaces;
	set<int> ElementsIds;
	//set<pEntity> elementSet;
	list<pEntity> elementList;

	void VertexList (pMesh theMesh);
	void LeitorDeLista(ErrorAnalysis *pErrorAnalysis, pMesh theMesh);
	void Clear_Containers();

private:
	char meshFilename[];

};

