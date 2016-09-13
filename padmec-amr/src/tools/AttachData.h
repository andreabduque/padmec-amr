/*
 * AttachData.h
 *
 *  Created on: 25/12/2008
 *      Author: rogerio
 */

#ifndef ATTACHDATA_H_
#define ATTACHDATA_H_

#include "auxiliar.h"

namespace PRS{

	/**
	 * A set of data can be associated to mesh entities like physical properties
	 * or geometric parameters. User should derive a class from AttacData to set/get
	 * data to vertices, edges, faces or tetrahedrals. User should define a pointer
	 * to a structure that holds all data needed. The job to attach and to get
	 * this pointer to some mesh entity is performed by AttacData class. Before
	 * Invoke getAttachedData_pointer or setAttachedData_pointer, setMeshDataId
	 * MUST be called!
	 * */
	class AttachData
	{
	public:
		AttachData(){}
		~AttachData(){}

	protected:

		// give a name to the set of data to be attached to mesh entity
		void setMeshDataId(string mdId){
			dataAttached_id = MD_lookupMeshDataId( mdId.c_str() );
		}

		template <class T>
		T* getAttachedData_pointer(pEntity ent){
			void *p; p=0;
			EN_getDataPtr(ent,dataAttached_id,&p);
			return (p) ? (T*)p : new T;
		}

		template <class T>
		void setAttachedData_pointer(pEntity ent, T* ptr){
			void *p; p=0;
			p = (T*)ptr;
			EN_attachDataPtr(ent,dataAttached_id,p);
		}

		template <class T>
		void deleteAttachedData_pointer(pEntity ent, T* ptr){
			EN_deleteData(ent,dataAttached_id);
		}

	private:

		pMeshDataId dataAttached_id;
	};

}


#endif /* ATTACHDATA_H_ */
