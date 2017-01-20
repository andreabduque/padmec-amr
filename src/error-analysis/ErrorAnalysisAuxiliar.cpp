/*
 * ErrorAnalysis_auxiliar.cpp
 *
 *  Created on: 29/02/2012
 *      Author: rogsoares
 */

#include "ErrorAnalysis.h"

void ErrorAnalysis::initialize(GeomData* pGCData, SimulatorParameters *pSimPar){

	// if initialize had already been called, it means pointers are still allocated from previous analysis.
	// Before proceed to a new analysis, delete them all.
	if (init){
		deletePointers();
	}

	int nnodes = pGCData->getNumNodes();
	int nelements = pGCData->getNumElements();

	// allocating memory
	alloc_DOUBLE_vector(__LINE__,__FILE__,pElemError,nelements);
	alloc_DOUBLE_vector(__LINE__,__FILE__,pWH_node,nnodes);
	alloc_BOOL_vector(__LINE__,__FILE__,isElementSingular,nelements);
	alloc_BOOL_vector(__LINE__,__FILE__,pElmToRemove,nelements);

	p_h_new = new Matrix<double>;
	p_h_new->allocateMemory(nelements,2);	// pressure and saturation

	// reseting array
	for (int i=0; i<nelements; i++){
		pElmToRemove[i] = false;
	}

	p_h_new->initialize(1e+10);

	// set to 0 the number of faces sharing the same vertex
	pGCData->calculateNumElemSharingVertex();

	// calculates Characteristic Dimension Length (CDL)
	pGCData->calculate_CDL();

	// weight CDL (h_old) on mesh vertices by the number of sharing elements
	pGCData->calculate_NodeAverage_CDL();

	// Use elements heights to define a minimum allowed element height for re-mesh.
	// This must be done for the very first time only.
	set_h_min(pGCData,pSimPar);

	// analysis error has been started
	init = true;
}

void ErrorAnalysis::resetVariables(int nelements){
	SGN = .0;
	SGN_sing = .0;
	averageError = .0;
	averageError_Singularity = .0;
	globalError = .0;
	globalError_Singularity = .0;
	for (int i=0; i<nelements; i++){
		pElemError[i] = .0;
	}
	setAllElementsAsNotSingular(nelements);
}

void ErrorAnalysis::deletePointers(){
	p_h_new->freeMemory();
	delete p_h_new; p_h_new = 0;
	dealloc_DOUBLE_vector(pElemError);
	dealloc_BOOL_vector(isElementSingular);
	dealloc_BOOL_vector(pElmToRemove);
	dealloc_DOUBLE_vector(pWH_node);
}

void ErrorAnalysis::setAllElementsAsNotSingular(int nelem){
	for (int i=0; i<nelem; i++){
		isElementSingular[i] = false;
	}
}

double ErrorAnalysis::calculate_ErrorSum(GeomData* pGCData, bool excludingSingularities=false){
	double error_sum = .0;
	double error;

	int nelements = pGCData->getNumElements();

	// loop over domains' faces
	if (!excludingSingularities){
		for (int i=0; i<nelements; i++){
			error = getElementError(i);
			error_sum -= error*error;
		}
	}
	else{
		for (int i=0; i<nelements; i++){
			if(!isSingular(i)){
				error = getElementError(i);
				error_sum -= error*error;
			}
		}
	}
	return error_sum;
}

void ErrorAnalysis::countElements(int nelem, bool seekforSingular=false){
	int numElem = 0;
	if (seekforSingular){
		for (int i=0; i<nelem; i++){
			if (isElementSingular[i]){
				numElem++;
			}
		}
		setNumElements_Singularity(numElem);
	}
	else{
		setNumElements(nelem);
	}
}

void ErrorAnalysis::set_h_min(GeomData* pGCData, SimulatorParameters *pSimPar){
	static bool h_min_set = true;
	if ( h_min_set ){
		int k = 0;
		h_min = 1e+10;
		for (int dom=0; dom<pGCData->getNumDomains(); dom++){
			int nfaces = pGCData->getNumElemPerDomain(dom);
			// loop over domains' faces
			for (int i=0; i<nfaces; i++){
				h_min = std::min(h_min, pGCData->getElem_CDL(k++) );
			}
		}
		h_min /= pSimPar->Remeshing_param3();
		h_min_set = false;
	}
}

/*************************************************************************************************************
 *************************************************************************************************************
 *************************************************************************************************************
 * Rogerio: 20/03/2015
 * Por questao de compatibilidade com os adaptadores (MAdlib e Gmsh) as funcoes abaixo:
 *
 * 	1)	update_h_new
 * 	2)  calculate_height_ratio
 * 	3)  getRefUnrefElementsList
 * 	4)  store_h_new
 *
 * continuarao usando o FMDB (argh!). Em um futuro proximo, o FMDB devera ser expurgado daqui.
 *************************************************************************************************************
 *************************************************************************************************************
 *************************************************************************************************************
 */
//void ErrorAnalysis::update_h_new(pMesh theMesh){
//	int i = 0;
//	pEntity node;
//	VIter vit = M_vertexIter(theMesh);
//	while ( (node = VIter_next(vit)) ){
//		//cout << EN_id(node)<< "\t" << pStore_h_new[i] << "\th_min = " << get_h_min() << "\tsingular? " << isSingular(node) << endl;
//		EN_attachDataDbl(node,MD_lookupMeshDataId( "elem_height" ),pStore_h_new[i]);
//		if (pStore_h_new[i]<.0){
//			throw Exception(__LINE__,__FILE__,"Nengative h_new!!!");
//		}
//		if (pStore_h_new[i] < get_h_min()){
//			char msg[256]; sprintf(msg,"h_new too small: %f\th_min = %f",pStore_h_new[i],get_h_min());
//			throw Exception(__LINE__,__FILE__,msg);
//		}
//		if (pStore_h_new[i] > 1e+9){
//			char msg[256]; sprintf(msg,"h_new too big: %f\th_min = %f",pStore_h_new[i],get_h_min());
//			throw Exception(__LINE__,__FILE__,msg);
//		}
//		i++;
//	}
//	VIter_delete(vit);
//}

//void ErrorAnalysis::calculate_height_ratio(pMesh theMesh){
//	double h_new_node, h_old_node;
//	pEntity node;
//	VIter vit = M_vertexIter(theMesh);
//	while ( node = VIter_next(vit) ){
//		EN_getDataDbl(node,MD_lookupMeshDataId( "h_old_node" ),&h_old_node);
//		EN_getDataDbl(node,MD_lookupMeshDataId( "elem_height" ),&h_new_node);
//		EN_attachDataDbl(node,MD_lookupMeshDataId( "height_ratio" ),(double)(h_new_node/h_old_node));
//	}
//	VIter_delete(vit);
//}

// New version of calculate_height_ratio
// ================================================================================
//void ErrorAnalysis::calculate_height_ratio(GeomData* pGCData){
//	for(i=0; i<nnodes; i++){
//
//	}
//	double h_new_node, h_old_node;
//	pEntity node;
//	VIter vit = M_vertexIter(theMesh);
//	while ( node = VIter_next(vit) ){
//		EN_getDataDbl(node,MD_lookupMeshDataId( "h_old_node" ),&h_old_node);
//		EN_getDataDbl(node,MD_lookupMeshDataId( "elem_height" ),&h_new_node);
//		EN_attachDataDbl(node,MD_lookupMeshDataId( "height_ratio" ),(double)(h_new_node/h_old_node));
//	}
//	VIter_delete(vit);
//}

//void ErrorAnalysis::getRefUnrefElementsList(pMesh theMesh, std::list<pEntity> &elementList, std::set<pEntity>& nodesBGMesh){
//	#ifdef TRACKING_PROGRAM_STEPS
//	cout << "TRACKING_PROGRAM_STEPS: getRefUnrefElementsList\tIN\n";
//	#endif
//
//	// avoid messing things
//	if (elementList.size()){
//		elementList.clear();
//		nodesBGMesh.clear();
//	}
//
//	pEntity elem;
//	std::list<double> heightList;
//	double param1;// = _pSimPar->Remeshing_param1();
//	double param2;// = _pSimPar->Remeshing_param2();
//	printf("param1: %.5f\tparam2: %.5f\n",param1,param2);
//	double h_old, h_new;
//	if (!theMesh){
//		throw Exception(__LINE__,__FILE__,"NULL Mesh!");
//	}
//
//	int dim = theMesh->getDim();
//	if (dim==2){
//// 		double hmin = get_h_min();
//// 		avgError = getAverageError();
//		FIter fit = M_faceIter (theMesh);
//		while ( (elem = FIter_next(fit)) ){
//// 			error = getElementError(elem);
//// 			h_old = getElement_CDL(elem);
//// 			h_new = h_old;
//			EN_getDataDbl(elem,MD_lookupMeshDataId( "h_new" ),&h_new);
//// 			if (fabs(error) > 1e-8){
//// 				h_new = h_old*(avgError/error);
//// 			}
//			EN_attachDataInt(elem,MD_lookupMeshDataId( "elem_to_remove" ),0);
//			if ( h_new < param1*h_old){
//				elementList.push_back(elem);
//				//heightList.push_back(h_new);
//				EN_attachDataInt(elem,MD_lookupMeshDataId( "elem_to_remove" ),1);
//			}
//			if (h_new > param2*h_old){
//				elementList.push_back(elem);
//				//heightList.push_back(h_new);
//				EN_attachDataInt(elem,MD_lookupMeshDataId( "elem_to_remove" ),1);
//			}
//		}
//		FIter_delete(fit);
//	}
//	else if (dim==3){
//		throw Exception(__LINE__,__FILE__,"Under construction!");
//	}
//	else{
//		throw Exception(__LINE__,__FILE__,"Mesh dimension unknown.");
//	}
//
//	heightList.clear();
//	#ifdef TRACKING_PROGRAM_STEPS
//	cout << "TRACKING_PROGRAM_STEPS: getRefUnrefElementsList\tOUT\n";
//	#endif
//}

//void ErrorAnalysis::store_h_new(pMesh theMesh){
//	int i = 0;
//	double height;
//	pEntity node;
//	VIter vit = M_vertexIter(theMesh);
//	while ( (node = VIter_next(vit)) ){
////		if (!isSingular(node)){
//			EN_getDataDbl(node,MD_lookupMeshDataId( "elem_height" ),&height);
//			pStore_h_new[i] = std::min(pStore_h_new[i],height);
////		}
//		i++;
//	}
//	VIter_delete(vit);
//}
