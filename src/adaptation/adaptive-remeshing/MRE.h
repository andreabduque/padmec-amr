#include "includes.h"
#include <string>

class Remover
{
public:

	// Contructor and Destructor
	Remover(pMesh theMesh);
	~Remover();


	std::set<pEdge> PlaneSurfaceEdges;
	std::set<pVertex> GeometryBoundNodes; // acho que este está sem uso... verificar...
	std::set<pVertex> GeomBoundNodes;
	std::set<pEdge> GeomBoundEdges;
	std::set<pEdge> MeshBoundEdges;
	std::set<pEdge> MeshBoundNodes;
	std::set<pVertex> AdaptNodes;
	std::set<pEdge> AdaptEdges;
	std::set<pVertex> BoundaryNodes;
	std::set<pEdge> BoundaryEdges;
	std::set<pEdge> RemovefromBoundEdges; // tive que colocar pra nao dar erro na contagem.
	std::set<pVertex> RemovefromBoundNodes;
	std::set<pFace> RemovefromBoundFaces;
	//std::list<pFace> BoundaryFaces; // transformei o antigo std::set Faces em std::set BoundaryFaces
	std::set<pFace> GeomBoundFaces; // criado o GeomBoundFaces

	std::set<pRegion> Tetras; // Originalmente so este std::set é pRegion e o resto é pFace.
	std::vector <int> argi;
	std::map <string,pVertex> CommonVertex;
	float GreatestEdge;

	std::set <pFace> triangles;
	std::map <int,float> VertexClsAux; // Colocar esse como uma variavel local da funcao...
	std::map <int,float> VertexCls;

	typedef struct eleMatriz{
		int ID;
 		float x1,x2,x3,y1,y2,y3,z1,z2,z3;
 		float Cl1, Cl2, Cl3;
		pVertex Pt1, Pt2, Pt3;



	} EleMatriz;

	class operadorEleMatriz {
	public:
		bool operator()(const EleMatriz a, const EleMatriz b) {
			if(a.ID < b.ID) {
				return true;
			}
			return false;
		}
	};

	

	//void ClReader (char *CLS);
	
	void Clear_Containers();
	void SaveBGMView1(const std::list<pEntity>& elementList, int size); // Rotina para salvar view da malha de background
	void SaveBGMView2(int size); // Rotina para salvar view da malha de background
	void Merge_Edges(pMesh theMesh2);
	void ResetMRE();
	void removeinternalfaces_2D (pMesh theMesh, std::list<pEntity> &elementList);
	void removetetras (pMesh theMesh);
	int Iterador (std::set <int> ListofElements, const std::list<pEntity>& elementList ,int I, pMesh theMesh);
	void BoundaryElements3D (std::set <int> ListofElements, pMesh theMesh);
	void BoundaryElements2D (pMesh theMesh, const std::list<pEntity>& elementList);
	void removeinternaledges (pMesh theMesh);
	void removeinternalnodes (pMesh theMesh);
	void removeexternaledges(pMesh theMesh);
	void removeexternalnodes(pMesh theMesh);
	void RemoveStgElements_2D(pMesh theMesh, float GE);
	void RemoveStgElements_3D(pMesh theMesh);
	void removestgNodes(pMesh theMesh);
	
	void CriadordeLista(pMesh theMesh); // lembrar de retirar esse daqui depois
	void CriadordeLista3D(pMesh theMesh); // lembrar de retirar esse daqui depois
	//int identifysurfaces_2D(pMesh Mesh);
	int Identify_and_Remove_Edges_2D(pMesh theMesh);
	int identifyvolumes_3D(pMesh Mesh);

	std::set <EleMatriz, operadorEleMatriz> STriangles;
	eleMatriz * MtxTriangles; // As coordenadas e cls dos triangulos que formam a malha de background
	//std::set <int> HoleNeighbours(pMesh theMesh2);
	std::set <pFace> ReturnBoundaryEdges();
	std::set <pFace> ReturnBoundaryNodes();
	//std::set <pFace> ReturnBoundaryFaces();
	std::map <int, pFace> MeshMap;


private:
	///char meshFilename[];
	
	
};


