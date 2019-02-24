#include <iostream>
#include <vector>
using namespace std;

/*
 * SparseMatrix
 *
 *	Created on: January 23, 2019
 *		Author: David Gilstad
 */

 /*
  * This class is a row in sparse matrix form. A row is made up of the row index,
  * column index, and value of an uncommon element in the original matrix.
  */
class SparseRow {
protected:
	int row, col, value; // row and column index, and the element's value
public:
	SparseRow();
	SparseRow(int r, int c, int val);
	int getRow();
	int getCol();
	int getValue();
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
class SparseMatrix {
protected:
	int noRows, noCols;	// #rows and #cols in the original matrix
	int	noNonSparseValues; // # of uncommon values in the matrix
	int	commonValue; // the common value in the original matrix
	vector<SparseRow>* myMatrix; // sparse matrix containing uncommon elements
public:
	SparseMatrix();
	SparseMatrix(int n, int m, int cv, int noNSV);
	int getRow(int i);
	int getCol(int i);
	int getVal(int i);
	int getCols();
	int getRows();
	int valueOf(int r, int c);
	bool add(int r, int c, int val);
	vector<SparseRow>* getMatrix();
	SparseMatrix* Transpose();
	SparseMatrix* Multiply(SparseMatrix& M);
	SparseMatrix* Add(SparseMatrix& M);
	void sort();
	void display();
	void displayMatrix();
	void readInMatrix();
};

/*
 * Default sparse row constructor. Sets row and col to -1, value to 0.
 */
SparseRow::SparseRow() {
	row = -1, col = -1, value = 0;
}

SparseRow::SparseRow(int r, int c, int val) {
	row = r, col = c, value = val;
}

/*
 * Default sparse matrix constructor. Sets noRows, noCols to -1, all else to 0.
 */
SparseMatrix::SparseMatrix() {
	noRows = -1, noCols = -1, commonValue = 0, noNonSparseValues = 0;
	myMatrix = new vector<SparseRow>();
}

/*
 * Non-default sparse matrix constructor. Sets all variables to the given
 * parameters (noRows = n, noCols = m).
 */
SparseMatrix::SparseMatrix(int n, int m, int cv, int noNSV) {
	noRows = n, noCols = m, commonValue = cv, noNonSparseValues = noNSV;
	myMatrix = new vector<SparseRow>();
}

/************** SPARSE ROW METHODS ***************/
/*
 * Returns the row index of the original element.
 */
int SparseRow::getRow() {
	return row;
}

/*
 * Returns the column index of the original element.
 */
int SparseRow::getCol() {
	return col;
}

/*
 * Returns the value of the element.
 */
int SparseRow::getValue() {
	return value;
}

/*
 * Adds the given value to the elements value.
 */
void SparseRow::add(int val) {
	value += val;
}

/*
 * Prints out the sparse row in the form: row#, col#, value.
 */
void SparseRow::display() {
	cout << row << ", " << col << ", " << value << endl;
}

/*
 * Sets row, col, and value to r, c, and val respectively.
 */
void SparseRow::set(int r, int c, int val) {
	row = r, col = c, value = val;
}

/*************** SPARSE MATRIX METHODS ***************/
/*
 * Return the number of columns in the matrix.
 */
int SparseMatrix::getCols() {
	return noCols;
}

/*
 * Return the number of rows in the matrix.
 */
int SparseMatrix::getRows() {
	return noRows;
}

int SparseMatrix::getCol(int i) {
	return myMatrix->at(i).getCol();
}

int SparseMatrix::getRow(int i) {
	return myMatrix->at(i).getRow();
}

int SparseMatrix::getVal(int i) {
	return myMatrix->at(i).getValue();
}

/*
 * Returns the value at element [r,c], or the common value if [r,c] not found.
 */
int SparseMatrix::valueOf(int r, int c) {
	for (int i = 0; i < noNonSparseValues; i++)
		if (getRow(i) == r && getCol(i) == c)
			return getVal(i);
	return commonValue;
}

/*
 * Goes through the matrix from 0-size and tries to add val to an element [r,c].
 *
 * Returns true if added succesfully, false it the element isn't found.
 */
bool SparseMatrix::add(int r, int c, int val) {
	for (int i = 0; i < (int)(myMatrix->size()); i++)
		if (getRow(i) == r && getCol(i) == c) {
			myMatrix->at(i).add(val);
			return true;
		}
	return false;
}

vector<SparseRow>* SparseMatrix::getMatrix() {
	return myMatrix;
}

/*
 * Returns the transpose matrix. A transpose had inversed row sizes, column sizes,
 * and indexes, a.k.a. flipping itself on the diagnol that from [0,0] to [n,n].
 */
SparseMatrix* SparseMatrix::Transpose() {
	SparseMatrix* T = new SparseMatrix(noCols, noRows, commonValue, noNonSparseValues);
	for (int i = 0; i < noNonSparseValues; i++)
		T->getMatrix()->push_back(SparseRow(getCol(i), getRow(i), getVal(i)));
	return T;
}

/*
 * Sparse matrix multiplication: Multiplies this matrix with M. Assumes #cols is
 * equal to M's #rows, and common value of both are 0.
 *
 * Return the product of the multiplication.
 */
SparseMatrix* SparseMatrix::Multiply(SparseMatrix& M) {
	if (noCols != M.getRows()) {
		cout << "Error: Matrices have invalid size for multiplication." << endl;
		return new SparseMatrix();
	} // #cols of the first must be equal to #rows of the second matrix

	SparseMatrix* P = new SparseMatrix(noRows, M.getCols(), 0, noRows*M.getCols());
	int v;

	for (int i = 0; i < noNonSparseValues; i++)
		// go through each element and try to find multplication match in M
		for (int mCol = 0; mCol < M.noCols; mCol++) {
			v = getVal(i) * M.valueOf(getCol(i), mCol);
			// if match in M exists, put the product in product matrix P
			if (!(*P).add(getRow(i), mCol, v) && v != 0)
				(*P).getMatrix()->push_back(SparseRow(getRow(i), mCol, v));
		}
	(*P).noNonSparseValues = P->getMatrix()->size();
	return P;
}

/*
 * Sparse matrix addition: Add the two matrices. Assumes #rows, #cols, and common
 * value of both are equal.
 *
 * Returns the sum of the matrices.
 */
SparseMatrix* SparseMatrix::Add(SparseMatrix& M) {
	if (noRows != M.getRows() || noCols != M.getCols()) {
		cout << "Error: Matrices have invalid size for addition" << endl;
		return new SparseMatrix();
	} // matrices must have equal sizes

	SparseMatrix* S = new SparseMatrix(noRows, noCols, commonValue, noRows*noCols);

	// assume all nonspares elements don't have a matching nonsparse value in M
	for (int i = 0; i < noNonSparseValues; i++)
		(*S).getMatrix()->push_back(SparseRow(getRow(i), getCol(i), getVal(i) + commonValue));
	for (int i = 0; i < M.noNonSparseValues; i++) {
		// if M has matching nonsparse value, subtract common value added above
		if (!(*S).add(M.getRow(i), M.getCol(i), M.getVal(i) - commonValue))
			(*S).getMatrix()->push_back(SparseRow(M.getRow(i), M.getCol(i), M.getVal(i) + commonValue));
	}
	(*S).noNonSparseValues = S->getMatrix()->size();
	return S;
}

/*
 * Sorts the sparse matrix lowest to highest using row and then column. Uses
 * insertion sorting.
 */
/*
void SparseMatrix::sort() {
	SparseRow temp;
	for (int i = 0; i < noNonSparseValues; i++)
		for (int j = i; j > 0; j--) {
			int r1 = myMatrix[j].getRow(), c1 = myMatrix[j].getCol(),
				r2 = myMatrix[j - 1].getRow(), c2 = myMatrix[j - 1].getCol();
			if (r1 < r2 || r1 == r2 && c1 < c2) {
				temp = myMatrix[j];
				myMatrix[j] = myMatrix[j - 1];
				myMatrix[j - 1] = temp;
			}
			else j = 0;
		}
}
*/

/*
 * Displays the sparse rows in the sparse matrix in order.
 *
 * @see SparseRow::display()
 */
void SparseMatrix::display() {
	for (int i = 0; i < noNonSparseValues; ++i)
		myMatrix->at(i).display();
}

/*
 * Recreates and prints the original matrix, with a tab between each element, and a
 * new line between each row.
 */
void SparseMatrix::displayMatrix() {
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
void SparseMatrix::readInMatrix() {
	int val;
	for (int row = 0; row < noRows; ++row)
		for (int col = 0; col < noCols; ++col) {
			cin >> val;
			if (val != commonValue)
				myMatrix->push_back(SparseRow(row, col, val));
		}
}

/*************** MAIN PROGRAM ***************
 *
 * Reads in two sparse matrices from the user, performs operations on the matrices
 * and displays them in.
 */
int main() {
	int n, m, cv, noNSV;
	SparseMatrix* temp;

	cin >> n >> m >> cv >> noNSV;
	SparseMatrix* firstOne = new SparseMatrix(n, m, cv, noNSV);
	(*firstOne).readInMatrix();
	// read in first and second matrix
	cin >> n >> m >> cv >> noNSV;
	SparseMatrix* secondOne = new SparseMatrix(n, m, cv, noNSV);
	(*secondOne).readInMatrix();

	cout << "First one in sparse matrix format" << endl;
	(*firstOne).display();

	cout << "First one in normal matrix format" << endl;
	(*firstOne).displayMatrix();

	cout << "Second one in sparse matrix format" << endl;
	(*secondOne).display();

	cout << "Second one in normal matrix format" << endl;
	(*secondOne).displayMatrix();

	cout << "After Transpose first one in normal format" << endl;
	temp = (*firstOne).Transpose();
	(*temp).displayMatrix();

	cout << "After Transpose second one in normal format" << endl;
	temp = (*secondOne).Transpose();
	(*temp).displayMatrix();

	cout << "Multiplication of matrices in sparse matrix form:" << endl;
	temp = (*secondOne).Multiply(*firstOne);
//	(*temp).sort();
	(*temp).display();

	cout << "Addition of matrices in sparse matrix form:" << endl;
	temp = (*secondOne).Add(*firstOne);
//	(*temp).sort();
	(*temp).display();

	system("pause");
	return 0;
}
