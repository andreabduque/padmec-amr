///*
// * MeshAdaptation.h
// *
// *  Created on: 16/02/2012
// *      Author: rogsoares
// */
//
//#ifndef MESHADAPTATION_H_
//#define MESHADAPTATION_H_
//
//#include "AMR.h"
//#include "Adaptation_SplitCallBack.h"
//
///*
// * Defining types:
// */
//typedef std::list<mFace*> SpecialElemList_2D;
//typedef std::list<mFace*>::iterator SpecialElemIterator_2D;
//
///*! brief For steady-state simulations. timeStep variable unnecessary
// */
//bool startMeshAdaptation(SimulatorParameters* pSimPar, PhysicPropData *pPPData, GeomData *pGCData,
//		                 AMR *pMAdapt, ErrorAnalysis *pErrorAnalysis, pMesh theMesh, int numFields);
//
///*! brief For all states (steady-state/transient) simulations. timeStep variable unnecessary
// */
//bool startMeshAdaptation(SimulatorParameters* pSimPar, PhysicPropData *pPPData, GeomData *pGCData,
//		                 AMR *pMAdapt, ErrorAnalysis *pErrorAnalysis, pMesh theMesh, double &timeStep, int numFields);
//
//
///*
// * Class H_Refinement is a base class where the right code to adapt 2D and 3D meshes will be used.
// */
//
//class H_Refinement : public AMR{
//
//public:
//
//	H_Refinement(){
//		//elem_id_special = MD_newMeshDataId("elem_id_special");
//	}
//
//	virtual ~H_Refinement(){
//		delete pCallback;
//	}
//
//	// it means that element was generated from a hang node
//	void setElementAsSpecial(pEntity elem, int i){
//		// unfortunately FMDB does note provide an attach function for booleans, which will save memory
//		EN_attachDataInt(elem,MD_lookupMeshDataId("elem_id_special"),i);
//	}
//
//	// says if element was generated from a hang node
//	static bool isElementSpecial(pEntity elem){
//		int flg = 0;
//		EN_getDataInt(elem,MD_lookupMeshDataId("elem_id_special"),&flg);
//		return flg;
//	}
//
//	/*! \brief: Eliminate special element's children (originated from hag-nodes)
//	 * \param pMesh FMDB mesh object
//	 */
//	virtual void removeSpecialElements(pMesh, ErrorAnalysis*)=0;
//
//	/*! \brief: Identify special elements (with hang-nodes) and store them on a list
//	 * \param theMesh FMDB mesh object
//	 */
//	virtual void createSpecialElementsList(pMesh theMesh, ErrorAnalysis *pErrorAnalysis)=0;
//
//	/*! \brief: eliminate hangnodes creating new elements
//	 * \param theMesh FMDB mesh object
//	 */
//	virtual void eliminateHangNodes(pMesh theMesh)=0;
//
//	/*! \brief: Guarantee that the refinement degree difference between two adjacent element is equal to one.
//	 *          The highest element refinement degree is limited to maxLevel to avoid excessive element subdivisions
//	 * \param pErrorAnalysis
//	 * \param theMesh FMDB mesh object
//	 */
//	virtual void unevenElements(ErrorAnalysis* pErrorAnalysis, pMesh theMesh, std::list<pEntity>& faceList, int max, int min)=0;
//
////	/*! \brief: Guarantee that the degree of freedom difference between two adjacent elements is 1.
////	 * \param pErrorAnalysis
////	 * \param theMesh FMDB mesh object
////	 */
////	virtual void unevenElements_part1(ErrorAnalysis* pErrorAnalysis, pMesh theMesh, std::list<pEntity>& faceList, int max, int min)=0;
////
////	/*! \brief: Correct A's element dor which is closed by two adjacent elements with dor+1. In this case, A element will be set to dor+1
////	 * \param pErrorAnalysis
////	 * \param theMesh FMDB mesh object
////	 */
////	virtual void unevenElements_part2(ErrorAnalysis* pErrorAnalysis, pMesh theMesh)=0;
//
//	/*
//	 * get maximum Degree of Refinement (DOR) difference between two adjacent elements
//	 */
//	virtual int checkMaxDORdifference(ErrorAnalysis* pErrorAnalysis, pMesh theMesh) = 0;
//
//	/*
//	 * perform mesh regularisation
//	 */
//	virtual void meshRegularization(ErrorAnalysis* pErrorAnalysis, pMesh theMesh, int min, std::list<pEntity>& elemList) = 0;
//	virtual bool leavesCorrectionMain(ErrorAnalysis*, pMesh, int, std::list<pEntity>&) = 0;
//	virtual void firstLeavesCorrection(ErrorAnalysis*, pMesh, pEntity, int, bool&) = 0;
//	virtual void secondLeavesCorrection(ErrorAnalysis*, pMesh, pEntity, int depth, bool&) = 0;
//
//	/*
//	 * After element mesh modification, recreate data structure
//	 */
//	virtual void refreshDataStructure(pMesh) = 0;
//
//	virtual int getNumElements(pMesh) = 0;
//
//	/*! \brief: Guarantee that all elements be flagged correctly (degree of refinement correction)
//	 * \param theMesh FMDB mesh object
//	 */
//	virtual void elementsFlagsCorrection(ErrorAnalysis*, pMesh) = 0;
//
//
//	/*
//	 * Common function for both 2-D and 3-D mesh adaptation
//	 */
//	void getChildren(pMesh, pEntity, std::vector<pFace>&);
//
//	/*! \brief: Perform element mesh subdivision and/or element leaves elimination
//	 * \param pErrorAnalysis
//	 * \param theMesh FMDB mesh object
//	 */
//	void adaptation(pMesh);
//
//	virtual void rodar(ErrorAnalysis *pErrorAnalysis, pMesh & theMesh)=0;
//
//protected:
//
//	Adaptation_SplitCallBack* pCallback;	// pointer passed to FMDB to create or delete mesh elements
//	SpecialElemList_2D SE_list_2D;			// list of special elements (containing hang-nodes)
//	SpecialElemIterator_2D SE_iter_2D;		// iterator for SE_list
//
//	/*
//	 * Do not adapt mesh without correct element flags from error analysis
//	 */
//	void setEFC(bool efc){ EFC = efc; }
//	bool EFC_isDone() const { return EFC; }
//
//private:
//
//	bool EFC;						// EFC - Element Flag Correction (true - OK, false - miss calling elementsFlagsCorrection()
//	//pMeshDataId elem_id_special;	// for elements on sigularity regions
//};
//
//#endif /* MESHADAPTATION_H_ */
