/*
 * Matrix.h
 *
 *  Created on: May 26, 2011
 *      Author: rogsoares
 */

#include "auxiliar.h"

#ifndef MATRIX_H_
#define MATRIX_H_

/*
 * This template class provides a simple way to create a matrix of any type and to access its values efficiently.
 * Values are stored in a array of pointers to arrays and not like objects, which will lead to a big overhead for
 * problems where data are continuously read like those found in numeric simulations. Mesh information like nodes
 * id and coordinates MUST NOT be stored like objects.
 */

template <class T>
class Matrix {

public:

	Matrix (){
		_rows = _cols = 0;
	}

	Matrix (int rows, int cols):_rows(rows),_cols(cols){
		mat = new T*[_rows];
		for (int i=0; i<_rows; i++){
			mat[i] = new T[_cols];
		}
	}
	
	Matrix (int rows, int cols,T val):_rows(rows),_cols(cols){
		mat = new T*[_rows];
		for (int i=0; i<_rows; i++){
			mat[i] = new T[_cols];
		}
		initialize(val);
	}

	~Matrix (){
		if (mat){
			for (int i=0; i<_rows; i++){
				delete mat[i];
			}
			delete [] mat; mat = 0;
		}
	}

	void getsize(int &m, int&n){
		m = _rows;
		n = _cols;
	}

	void allocateMemory(int rows, int cols){
		_rows = rows;
		_cols = cols;
		allocateMemory();
	}

	void allocateMemory(int rows){ // we are supposing matrix has one column: matrix is  a vector
		_rows = rows;
		_cols = 1;
		allocateMemory();
	}

	void allocateMemory(){
		mat = new T*[_rows];
		for (int i=0; i<_rows; i++){
			mat[i] = new T[_cols];
		}
	}

	void freeMemory(){
		for (int i=0; i<_rows; i++){
			delete[] mat[i];
		}
		delete[] mat;
		mat = 0;
	}

	void initialize(T val){
		if (!_rows || !_cols)
			cout << "WARNING: Cannot initialize matrix with null size.\n";
		for (int i=0; i<_rows; i++){
			for (int j=0; j<_cols; j++){
				mat[i][j] = val;
			}
		}
	}

	void printfMatrix() const{
		for (int i=0; i<_rows; i++){
			for (int j=0; j<_cols; j++){
				cout << mat[i][j] << "  ";
			}
		cout << endl;
		}
				
	}

	void print(const char* filename) const{
		ofstream fid;
		fid.open(filename);
		fid << setprecision(6) << fixed;
		for (int i=0; i<_rows; i++){
			for (int j=0; j<_cols; j++){
				fid << mat[i][j] << "  ";
			}
			fid << endl;
		}
		fid.close();
	}

	T& operator()( int i, int j ) {
		if ((i<0 && i>=_rows) || (j<0 && j>=_cols))
			throw Exception(__LINE__,__FILE__,"Attempt of getting value out of bound\n");
		else
			return mat[i][j];
	}

	T& operator()( int i) {
		if ((i<0 && i>=_rows)){
			throw Exception(__LINE__,__FILE__,"Attempt of getting value out of bound\n");
		}
		else{
			return mat[i][0];
		}
	}

	void setValue(int i, int j, const T &val){
		if ((i<0 && i>=_rows) || (j<0 && j>=_cols)){
			throw Exception(__LINE__,__FILE__,"Attempt of getting value out of bound\n");
		}
		else{
			mat[i][j] = val;
		}
	}

	void setValue(int i, const T &val){
		mat[i][0] = val;
	}

	void getRow(int row,const T **vec){
		*vec = mat[row]; 
	}

	const T* getrowconst(int row){
		return mat[row];
	}

	T* getrow(int row){
		return mat[row];
	}

	const T getValue(int i, int j) const{
		if (i>=_rows || j>=_cols){
			cerr << "ERROR: Index out of bounds!\n";
			exit(1);
		}
		return mat[i][j];
	}

	T getValue(int i){
		return mat[i][0];
	}

	void setRow(int row,const T *vec){ 
		for (int i=0; i<_cols; i++){
			mat[row][i] = vec[i];
		}
	} 

	void printNumRowsCols() const{
		cout << "rows: " << _rows <<"\tcols: " << _cols << endl;
	}


private:

	int _rows;
	int _cols;
	T **mat;

};

#endif /* MATRIX_H_ */
