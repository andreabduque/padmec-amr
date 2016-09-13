///*
// * Adaptation_SplitCallBack.h
// *
// *  Created on: 04/06/2012
// *      Author: rogsoares
// */
//
//#ifndef ADAPTATION_SPLITCALLBACK_H_
//#define ADAPTATION_SPLITCALLBACK_H_
//
//#include "ErrorAnalysis.h"
//
//class Adaptation_SplitCallBack : public mSplitCallbacks{
//public:
//	Adaptation_SplitCallBack(mMesh*){
//	}
//
//	~Adaptation_SplitCallBack(){};
//
//	int operator()(mEntity *e){
//		//int ret;
//		// If a leave is not flagged to be refined, maybe its parent be.
//		// If the parent is also not flagged to refine, keep going down looking for a flagged parent up to entity's root.
////		pEntity e_tmp = e;
////		if (refine){														// REFINEMENT
////			// Don't look to an element with children to refine it
////			if (theMesh->getRefinementDepth(e)) {
////				return 0;
////			}
////			while (1){
////				int dor = pEA->getLevelOfRefinement(e_tmp);
////				int parent_depth = theMesh->getRefinementDepth(e_tmp);
////				if ( dor-parent_depth>0 ){
////				//	cout << "returning 1!!!\n";
////					return 1;
////				}
////				else{
////					if (e_tmp == e->root()){
////						return 0;
////					}
////					e_tmp = e_tmp->parent();
////				}
////			}
////		}
////		else{																// UNREFINEMENT
////			if (theMesh->getRefinementDepth(e)==1){
////				int IDs[3];
////				getTriVerticesIDs(e,IDs);
////			//	cout << "IDs: " << IDs[0] << " " << IDs[1] << " " << IDs[2] << endl;
////				list<mEntity*> leaves;
////				e->getLeaves(leaves);
////				int numLeaves = (int)leaves.size();
////				// if element has two leaves, it's a special element. It's not handled here
////				if (numLeaves==4){
////					/*
////					 * leaf with highest level will decide if parent element must be unrefined or keep as it is.
////					 * initialize highestLevel variable to take the highest level to refine/unrefine
////					 */
////					int highestLevel = -1000;
////					for(std::list<mEntity*>::iterator lit=leaves.begin(); lit!=leaves.end();lit++){
////						highestLevel = std::max(highestLevel,pEA->getLevelOfRefinement( *lit ));
////					}
////					if (highestLevel>=0){
////						return 0;
////					}
////					else{
////						//pEA->setLevelOfRefinement((pEntity)e,highestLevel+1);
////						return -1;
////					}
////				}
////			}
////		}
//	}
//
//	void setUnrefine(){
//		refine = false;
//	}
//
//	void setRefine(){
//		refine = true;
//	}
//
//	virtual void splitCallback(mEntity*){};
//	virtual void unsplitCallback(mEntity*){};
//
//private:
//	pMesh theMesh;
//	ErrorAnalysis* pEA;
//	bool refine;
//};
//
//
//
//
//#endif /* ADAPTATION_SPLITCALLBACK_H_ */
