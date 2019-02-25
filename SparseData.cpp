#include <iostream>
#include <vector>
using namespace std;

/*
 * In case the matrices don't have the same #rows and #cols.
 */
template <class DT>
class ExceptionAdd {};

/*
 * In case the #cols of the first matrix doesn't equal the #rows of the second.
 */
template <class DT>
class ExceptionMultiply {};

/*
 * In case the matrices don't have the same commonValue.
 */
template <class DT>
class ExceptionCV {};

/*
 * SparseMatrix
 *
 *	Created on: February 23, 2019
 *		Author: David Gilstad
 */

 /*
  * This class is a row in sparse matrix form. A row is made up of the row index,
  * column index, and value of an uncommon element in the original matrix.
  */
template <class DT>
class SparseRow {
protected:
	int row, col; // row and column index
	DT value; // the element's value
public:
	SparseRow();
	SparseRow(int r, int c, DT val);
	virtual ~SparseRow();
	int getRow();
	int getCol();
	DT getValue();
	void add(int val);
	void display();
	void set(int r, int c, int val);
};

/*
 * This class is a sparse matrix. A sparse matrix has a majority of values that are
 * equal, and is made up of the uncommon values, scattered throughout the matrix. A
 * sparse matrix is a more simple and efficient way to represent and operate on
 * the original matrix.
 */
template <class DT>
class SparseMatrix {
friend ostream& operator<< (ostream& s, const SparseMatrix& M);
protected:
	int noRows, noCols;	// #rows and #cols in the original matrix
	int	noNonSparseValues; // # of uncommon values in the matrix
	DT commonValue; // the common value in the original matrix
	vector<SparseRow<DT> >* myMatrix; // sparse matrix containing uncommon elements
public:
	SparseMatrix();
	SparseMatrix(int n, int m, int cv);
	virtual ~SparseMatrix();
	int getRow(int i);
	int getCol(int i);
	DT getVal(int i);
	int getCols();
	int getRows();
	DT valueOf(int r, int c);
	bool add(int r, int c, DT val);
	vector<SparseRow<DT> >* getMatrix();
	void setNoNSV(int size);
	void setSparseRow(int pos, int r, int c, DT& v);
	void displayMatrix();
	void readInMatrix();
	SparseMatrix* operator*(SparseMatrix& M);
	SparseMatrix* operator+(SparseMatrix& M);
	SparseMatrix* operator!();
};

/*
 * Displays the sparse rows in the sparse matrix in order.
 */
template <class DT>
ostream& operator<< (ostream& s, SparseMatrix<DT>& M) {
	for (int i = 0; i < M.getMatrix()->size(); ++i)
		s << M.getRow(i) << ", " << M.getCol(i) << ", " << M.getVal(i) << endl;
	return s;
}

/*
 * Sparse matrix multiplication: Multiplies this matrix with M. Assumes #cols is
 * equal to M's #rows, and common value of both are 0.
 *
 * Return the product of the multiplication.
 */
template <class DT>
SparseMatrix<DT>* SparseMatrix<DT>::operator*(SparseMatrix<DT>& M) {
	if (noCols != M.getRows()) throw ExceptionMultiply<int>();
	if (commonValue != M.commonValue) throw ExceptionCV<int>();

	SparseMatrix<DT>* P = new SparseMatrix<DT>(noRows, M.getCols(), 0);
	int v;
	for (int i = 0; i < noNonSparseValues; i++)
		// go through each element and try to find multplication match in M
		for (int mCol = 0; mCol < M.noCols; mCol++) {
			v = getVal(i) * M.valueOf(getCol(i), mCol);
			// if match in M exists, put the product in product matrix P
			if (!P->add(getRow(i), mCol, v) && v != 0) {
				SparseRow<DT>* p = new SparseRow<DT>(getRow(i), mCol, v);
				P->getMatrix()->push_back(*p);
			}
		}
	P->setNoNSV(P->getMatrix()->size());
	return P;
}

 /*
  * Sparse matrix addition: Add the two matrices. #rows, #cols, and common
  * value of both must all be equal.
  *
  * Returns the sum of the matrices.
  */
template <class DT>
SparseMatrix<DT>* SparseMatrix<DT>::operator+(SparseMatrix<DT>& M) {
	if (noCols != M.getCols() && noRows != M.getRows())	throw ExceptionAdd<int>();
	if (commonValue != M.commonValue) throw ExceptionCV<int>();

	SparseMatrix<DT>* S = new SparseMatrix<DT>(noRows, noCols, commonValue);
	// assume all nonsparse elements don't have a matching nonsparse value in M
	for (int i = 0; i < noNonSparseValues; i++) {
		SparseRow<DT>* s = new SparseRow<DT>(getRow(i), getCol(i), getVal(i) + commonValue);
		S->getMatrix()->push_back(*s);
	}
	for (int i = 0; i < M.noNonSparseValues; i++) {
		// if M has matching nonsparse value, subtract common value added above
		if (!S->add(M.getRow(i), M.getCol(i), M.getVal(i) - commonValue)) {
			SparseRow<DT>* s = new SparseRow<DT>(M.getRow(i), M.getCol(i), M.getVal(i) + commonValue);
			S->getMatrix()->push_back(*s);
		}
		else {
			cout << "\n Added " << M.getVal(i) - commonValue << " to element ";
			S->getMatrix()->at(i).display();
		}
	}
	S->setNoNSV(S->getMatrix()->size());
	return S;
}

 /*
  * Returns the transpose matrix. A transpose has inversed row sizes, column sizes,
  * and indexes, a.k.a. flipping itself on the diagnol that's from [0,0] to [n,n].
  */
template <class DT>
SparseMatrix<DT>* SparseMatrix<DT>::operator!() {
	SparseMatrix<DT>* T = new SparseMatrix<DT>(noCols, noRows, commonValue);
	for (int i = 0; i < noNonSparseValues; i++) {
		SparseRow<DT>* t = new SparseRow<DT>(getCol(i), getRow(i), getVal(i));
		T->getMatrix()->push_back(*t);
	}
	setNoNSV(noNonSparseValues);
	return T;
}

/*
 * Default sparse row constructor. Sets row and col to -1, value to 0.
 */
template <class DT>
SparseRow<DT>::SparseRow() {
	row = -1, col = -1, value = 0;
}

/*
 * Non-default constructor. Sets row, col and val to the respective arguments.
 */
template <class DT>
SparseRow<DT>::SparseRow(int r, int c, DT val) {
	row = r, col = c, value = val;
}

/*
 * Sparse row destructor.
 */
template <class DT>
SparseRow<DT>::~SparseRow() {}

/*
 * Default sparse matrix constructor. Sets noRows, noCols to -1, all else to 0.
 */
template <class DT>
SparseMatrix<DT>::SparseMatrix() {
	noRows = -1, noCols = -1, commonValue = 0, noNonSparseValues = 0;
	myMatrix = new vector<SparseRow<DT> >();
}

/*
 * Non-default sparse matrix constructor. Sets variables to the given parameters 
 * (noRows = n, noCols = m), and sets noNonSparseValues to 0.
 */
template <class DT>
SparseMatrix<DT>::SparseMatrix(int n, int m, int cv) {
	noRows = n, noCols = m, commonValue = cv, noNonSparseValues = 0;
	myMatrix = new vector<SparseRow<DT> >();
}

/*
 * Sparse matrix destructor.
 */
template <class DT>
SparseMatrix<DT>::~SparseMatrix() {
	delete myMatrix;
}

/************** SPARSE ROW METHODS ***************/
/*
 * Returns the row index of the original element.
 */
template <class DT>
int SparseRow<DT>::getRow() {
	return row;
}

/*
 * Returns the column index of the original element.
 */
template <class DT>
int SparseRow<DT>::getCol() {
	return col;
}

/*
 * Returns the value of the element.
 */
template <class DT>
DT SparseRow<DT>::getValue() {
	return value;
}

/*
 * Adds the given value to the elements value.
 */
template <class DT>
void SparseRow<DT>::add(int val) {
	value += val;
}

/*
 * Prints out the sparse row in the form: row#, col#, value.
 */
template <class DT>
void SparseRow<DT>::display() {
	cout << row << ", " << col << ", " << value << endl;
}

/*
 * Sets row, col, and value to r, c, and val respectively.
 */
template <class DT>
void SparseRow<DT>::set(int r, int c, int val) {
	row = r, col = c, value = val;
}

/*************** SPARSE MATRIX METHODS ***************/
/*
 * Return the number of columns in the matrix.
 */
template <class DT>
int SparseMatrix<DT>::getCols() {
	return noCols;
}

/*
 * Return the number of rows in the matrix.
 */
template <class DT>
int SparseMatrix<DT>::getRows() {
	return noRows;
}

/*
 * Returns the column of the element at index i of myMatrix.
 */
template <class DT>
int SparseMatrix<DT>::getCol(int i) {
	return myMatrix->at(i).getCol();
}

/*
 * Returns the row of the element at index i of myMatrix.
 */
template <class DT>
int SparseMatrix<DT>::getRow(int i) {
	return myMatrix->at(i).getRow();
}

/*
 * Returns the value of the element at index i of myMatrix.
 */
template <class DT>
DT SparseMatrix<DT>::getVal(int i) {
	return myMatrix->at(i).getValue();
}

/*
 * Returns the value at element [r,c], or the common value if [r,c] not found.
 */
template <class DT>
DT SparseMatrix<DT>::valueOf(int r, int c) {
	for (int i = 0; i < myMatrix->size(); i++)
		if (getRow(i) == r && getCol(i) == c)
			return getVal(i);
	return commonValue;
}

/*
 * Goes through the matrix from 0-size and tries to add val to an element [r,c].
 *
 * Returns true if added succesfully, false if [r,c] isn't found.
 */
template <class DT>
bool SparseMatrix<DT>::add(int r, int c, DT val) {
	for (int i = 0; i < (int)(myMatrix->size()); i++)
		if (getRow(i) == r && getCol(i) == c) {
			myMatrix->at(i).add(val);
			return true;
		}
	return false;
}

/*
 * Returns myMatrix
 */
template <class DT>
vector<SparseRow<DT> >* SparseMatrix<DT>::getMatrix() {
	return myMatrix;
}

/*
 * Sets the noNonSpareValues to the given argument size.
 */
template <class DT>
void SparseMatrix<DT>::setNoNSV(int size) {
	noNonSparseValues = size;
}

/*
 * Sets the sparse row at index pos in myMatrix to the given r, c, and v.
 */
template <class DT>
void SparseMatrix<DT>::setSparseRow(int pos, int r, int c, DT& v) {
	myMatrix->at(pos) = SparseRow<DT>(r, c, v);
}

/*
 * Recreates and prints the original matrix, with a tab between each element, and a
 * new line between each row.
 */
template <class DT>
void SparseMatrix<DT>::displayMatrix() {
	for (int row = 0; row < noRows; row++) {
		for (int col = 0; col < noCols; col++) {
			cout << (*this).valueOf(row, col) << "\t";
		}
		cout << endl;
	}
}

/*
 * Takes in user input to create the sparse matrix.
 */
template <class DT>
void SparseMatrix<DT>::readInMatrix() {
	int val;
	for (int row = 0; row < noRows; ++row)
		for (int col = 0; col < noCols; ++col) {
			cin >> val;
			if (val != commonValue)
				myMatrix->push_back(SparseRow<DT>(row, col, val));
		}
	noNonSparseValues = myMatrix->size();
}

/*************** MAIN PROGRAM ***************
 *
 * Reads in two sparse matrices from the user, performs operations on the matrices
 * and displays them in.
 */
int main() {
	int n, m, cv;
	SparseMatrix<int>* temp;

	// read in first matrix
	cout << "First Matrix:" << endl;
	cin >> n >> m >> cv;
	SparseMatrix<int>* firstOne = new SparseMatrix<int>(n, m, cv);
	(*firstOne).readInMatrix();

	cout << "First one in sparse matrix format" << endl;
	cout << (*firstOne);

	cout << "First one in normal matrix format" << endl;
	(*firstOne).displayMatrix();

	// read in second matrix
	cout << "Second Matrix" << endl;
	cin >> n >> m >> cv;
	SparseMatrix<int>* secondOne = new SparseMatrix<int>(n, m, cv);
	(*secondOne).readInMatrix();

	cout << "Second one in sparse matrix format" << endl;
	cout << (*secondOne);

	cout << "Second one in normal matrix format" << endl;
	(*secondOne).displayMatrix();

	cout << "After Transpose first one in normal format" << endl;
	temp = !(*firstOne);
	(*temp).displayMatrix();

	cout << "After Transpose second one in normal format" << endl;
	temp = !(*secondOne);
	(*temp).displayMatrix();

	try {
		cout << "Multiplication of matrices in sparse matrix form:" << endl;
		temp = (*secondOne)*(*firstOne);
		cout << (*temp);
	}
	catch (ExceptionAdd<int> e1) {
		cout << "Error: Matrices have invalid size for addition." << endl;
	}
	catch (ExceptionCV<int> e2) {
		cout << "Error: Matrices must have same common value" << endl;
	}

	try {
		cout << "Addition of matrices in sparse matrix form:" << endl;
		temp = (*secondOne) + (*firstOne);
		cout << (*temp);
	}
	catch (ExceptionMultiply<int> e1) {
		cout << "Error: Matrices have invalid size for multiplication." << endl;
	}
	catch (ExceptionCV<int> e2) {
		cout << "Error: Matrices must have same common value" << endl;
	}

	delete firstOne;
	delete secondOne;
	delete temp;

	system("pause");
	return 0;
}
