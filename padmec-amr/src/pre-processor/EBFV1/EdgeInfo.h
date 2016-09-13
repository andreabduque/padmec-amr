/*
 * EdgeInfo.h
 *
 *  Created on: Apr 21, 2011
 *      Author: rogsoares
 *
 *      Manages informations about Cij coefficients of all edges only with remote
 *      copies. It's was built to make easier Cij unification over partition
 *      boundaries.
 */

#ifndef EDGEINFO_H_
#define EDGEINFO_H_

#include "includes.h"


class IdPair{
public:
	IdPair(int id0,int id1):_id0(id0),_id1(id1){}
	~IdPair(){}
	int getId0() const { return _id0; }
	int getId1() const { return _id1; }

private:
	int _id0;
	int _id1;
};


// local class for unify Cij's
class EdgeInfo{
public:

	enum xyzLists {xlist, ylist, zlist};
	typedef list<IdPair>::iterator IDIter;
	typedef list<double>::iterator CijIter;
	typedef map<int, int >::const_iterator MCI; // MCI - map const iterator

	EdgeInfo() {}

	~EdgeInfo() {}

	IDIter IdPair_beginIterator() { return listEdgeIds.begin(); }
	IDIter IdPair_endIterator() { return listEdgeIds.end(); }

	void printEdgeIds(){
		CijIter Cijx = xList.begin();
		CijIter Cijy = yList.begin();
		CijIter Cijz = zList.begin();
		for (IDIter iter = listEdgeIds.begin(); iter != listEdgeIds.end(); iter++,Cijx++,Cijy++,Cijz++)
			printf("[%d]aresta: %d %d :: %f %f %f\n",P_pid(),iter->getId0()+1, iter->getId1()+1,*Cijx,*Cijy,*Cijz);
	}

	CijIter Cij_beginIterator(xyzLists whichone){
		whichList = whichone;
		return Cij_beginIterator();
	}

	CijIter Cij_beginIterator(){
		switch (whichList){
		case xlist:  return xList.begin(); break;
		case ylist:  return yList.begin(); break;
		case zlist:  return zList.begin(); break;
		}
	}

//	CijIter Cij_endIterator(xyzLists whichone){
//		whichList = whichone;
//		return Cij_endIterator();
//	}
//
//	CijIter Cij_endIterator(){
//		switch (whichList){
//		case xlist:  xList.end(); break;
//		case ylist:  yList.end(); break;
//		case zlist:  zList.end(); break;
//		}
//		return
//	}

	void printListVals(){
		CijIter Cijx = xList.begin();
		CijIter Cijy = yList.begin();
		CijIter Cijz = zList.begin();
		for (; Cijx != xList.end(); Cijx++,Cijy++,Cijz++)  printf("[%d] Cij: %f %f %f\n",P_pid(),*Cijx,*Cijy,*Cijz);
	}
	int numEdges() const {  return listEdgeIds.size();  }

	void setWhichList(xyzLists wl) { whichList = wl; }

	void setEdgeIds(IdPair ids){
		listEdgeIds.push_back(ids);
	}

	const list<IdPair>& getRemoteEdgeIds() { return  listEdgeIds; }

	const list<double>& getRemoteEdgeCijList()  { return getRemoteEdgeCijList(whichList); }

	void setRemoteEdgeCij(const double &val){
		setRemoteEdgeCij(whichList,val);
	}

	void setRemoteEdgeCij(xyzLists whichone, const double &val){
		switch (whichone){
		case xlist:  xList.push_back(val); break;
		case ylist:  yList.push_back(val); break;
		case zlist:  zList.push_back(val); break;
		}
	}

	int sizeRemoteEdgeCij(xyzLists whichone){
		switch (whichone){
		case xlist:  return xList.size(); break;
		case ylist:  return yList.size(); break;
		case zlist:  return zList.size(); break;
		}
	}

	const list<double>& getRemoteEdgeCijList(xyzLists whichone){
		switch (whichone){
		case xlist:  return xList; break;
		case ylist:  return yList; break;
		case zlist:  return zList; break;
		}
	}

	void updateRemoteEdgeCijList(CijIter iter, double val){
		switch (whichList){
		case xlist:  xList.insert(iter, val); break;
		case ylist:  yList.insert(iter, val); break;
		case zlist:  zList.insert(iter, val); break;
		}
	}

	// after everything is done, clean the room
	void deleteEdgeIds() { listEdgeIds.clear(); }

	// after update each coordenate of Cij delete it
	void deleteCij_list(){
		deleteCij_list(whichList);
	}

	void deleteCij_list(xyzLists whichone){
		switch (whichone){
		case xlist:  xList.clear(); break;
		case ylist:  yList.clear(); break;
		case zlist:  zList.clear(); break;
		}
	}

	//mapRows() must be called BEFORE use numRowstoImport() and rowstoImport()
	void mapRows(){
		set<int> sTemp;
		for ( IDIter iter = listEdgeIds.begin(); iter != listEdgeIds.end(); iter++)
			sTemp.insert( iter->getId0() );

		int i = 0;
		for( set<int>::iterator sIter = sTemp.begin(); sIter != sTemp.end(); i++, sIter++)
			mRows.insert(make_pair(*sIter, i) );
		sTemp.clear();
	}

	void showMapping(){
		MCI miter = mRows.begin();
		for (; miter !=mRows.end(); miter++)
			printf("[%d] ido=%d -> %d\n",P_pid(),miter->first, miter->second);
	}

	int numRowstoImport() const { return mRows.size(); }

	void rowstoImport(int *rows) const{
		MCI mapIter = mRows.begin();
		for ( int i=0; mapIter != mRows.end(); mapIter++, i++)
			rows[i] = mapIter->first;
	}

	// which row must be read from parallel matrix
	int findMappedRow(const int &row){
		MCI mapIter = mRows.find(row);
		return mapIter->second;
	}

	void deleteMapRows() { mRows.clear(); }

private:

	xyzLists whichList;               // enum to specify which list of coordenate must be returned
	list<IdPair> listEdgeIds;     // list for edges' ids
	// list to store Cij coordenates
	list<double> xList;
	list<double> yList;
	list<double> zList;

	map<int,int> mRows;   // mapping to chose which rows must be import from parallel matrix

};


#endif /* EDGEINFO_H_ */
