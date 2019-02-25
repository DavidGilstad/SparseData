#include <iostream>
#include <vector>
using namespace std;

class ExceptionAdd {};

class ExceptionMultiply {};

class ExceptionCV {};

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
	virtual ~SparseRow();
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
friend ostream& operator<< (ostream& s, const SparseMatrix& M);
protected:
	int noRows, noCols;	// #rows and #cols in the original matrix
	int	noNonSparseValues; // # of uncommon values in the matrix
	int	commonValue; // the common value in the original matrix
	vector<SparseRow>* myMatrix; // sparse matrix containing uncommon elements
public:
	SparseMatrix();
	SparseMatrix(int n, int m, int cv);
	virtual ~SparseMatrix();
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
//	void sort();
	void setNoNSV(int size);
	void setSparseRow(int pos, int r, int c, int v);
	void display();
	void displayMatrix();
	void readInMatrix();

	SparseMatrix* operator*(SparseMatrix& M);
	SparseMatrix* operator+(SparseMatrix& M);
	SparseMatrix* operator!();
};

ostream& operator<< (ostream& s, SparseMatrix& M) {
	for (int i = 0; i < M.getMatrix()->size(); ++i)
		s << M.getRow(i) << ", " << M.getCol(i) << ", " << M.getVal(i) << endl;
	return s;
}

SparseMatrix* SparseMatrix::operator*(SparseMatrix& M) {
	try {
		if (noCols != M.getRows())
			throw ExceptionMultiply();
		if (commonValue != M.commonValue)
			throw ExceptionCV();
	}
	catch (ExceptionMultiply e1) {
		cout << "Error: Matrices have invalid size for multiplication." << endl;
		return new SparseMatrix();
	}
	catch (ExceptionCV e2) {
		cout << "Error: Matrices must have same common value" << endl;
		return new SparseMatrix();
	}

	SparseMatrix* P = new SparseMatrix(noRows, M.getCols(), 0);
	int v;

	for (int i = 0; i < noNonSparseValues; i++)
		// go through each element and try to find multplication match in M
		for (int mCol = 0; mCol < M.noCols; mCol++) {
			v = getVal(i) * M.valueOf(getCol(i), mCol);
			// if match in M exists, put the product in product matrix P
			if (!P->add(getRow(i), mCol, v) && v != 0)
				P->getMatrix()->push_back(SparseRow(getRow(i), mCol, v));
		}
	P->setNoNSV(P->getMatrix()->size());
	return P;
}

SparseMatrix* SparseMatrix::operator+(SparseMatrix& M) {
	try {
		if (noCols != M.getCols() && noRows != M.getRows())
			throw ExceptionAdd();
		if (commonValue != M.commonValue)
			throw ExceptionCV();
	}
	catch (ExceptionAdd e1) {
		cout << "Error: Matrices have invalid size for addition." << endl;
		return new SparseMatrix();
	}
	catch (ExceptionCV e2) {
		cout << "Error: Matrices must have same common value" << endl;
		return new SparseMatrix();
	}

	SparseMatrix* S = new SparseMatrix(noRows, noCols, commonValue);

	// assume all nonsparse elements don't have a matching nonsparse value in M
	for (int i = 0; i < noNonSparseValues; i++)
		S->getMatrix()->push_back(SparseRow(getRow(i), getCol(i), getVal(i) + commonValue));
	for (int i = 0; i < M.noNonSparseValues; i++) {
		// if M has matching nonsparse value, subtract common value added above
		if (!S->add(M.getRow(i), M.getCol(i), M.getVal(i) - commonValue))
			S->getMatrix()->push_back(SparseRow(M.getRow(i), M.getCol(i), M.getVal(i) + commonValue));
	}
	S->setNoNSV(S->getMatrix()->size());
	return S;
}

SparseMatrix* SparseMatrix::operator!() {
	SparseMatrix* T = new SparseMatrix(noCols, noRows, commonValue);
	for (int i = 0; i < noNonSparseValues; i++) {
		T->getMatrix()->push_back(SparseRow(getCol(i), getRow(i), getVal(i)));
	}
	setNoNSV(noNonSparseValues);
	return T;
}

/*
 * Default sparse row constructor. Sets row and col to -1, value to 0.
 */
SparseRow::SparseRow() {
	row = -1, col = -1, value = 0;
}

SparseRow::SparseRow(int r, int c, int val) {
	row = r, col = c, value = val;
}

SparseRow::~SparseRow() {}

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
SparseMatrix::SparseMatrix(int n, int m, int cv) {
	noRows = n, noCols = m, commonValue = cv, noNonSparseValues = 0;
	myMatrix = new vector<SparseRow>();
}

SparseMatrix::~SparseMatrix() {
	delete myMatrix;
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
	for (int i = 0; i < myMatrix->size(); i++)
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
	SparseMatrix* T = new SparseMatrix(noCols, noRows, commonValue);
	for (int i = 0; i < noNonSparseValues; i++) {
		T->getMatrix()->push_back(SparseRow(getCol(i), getRow(i), getVal(i)));
	}
	setNoNSV(noNonSparseValues);
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

	SparseMatrix* P = new SparseMatrix(noRows, M.getCols(), 0);
	int v;

	for (int i = 0; i < noNonSparseValues; i++)
		// go through each element and try to find multplication match in M
		for (int mCol = 0; mCol < M.noCols; mCol++) {
			v = getVal(i) * M.valueOf(getCol(i), mCol);
			// if match in M exists, put the product in product matrix P
			if (!P->add(getRow(i), mCol, v) && v != 0)
				P->getMatrix()->push_back(SparseRow(getRow(i), mCol, v));
		}
	P->setNoNSV(P->getMatrix()->size());
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

	SparseMatrix* S = new SparseMatrix(noRows, noCols, commonValue);

	// assume all nonsparse elements don't have a matching nonsparse value in M
	for (int i = 0; i < noNonSparseValues; i++)
		S->getMatrix()->push_back(SparseRow(getRow(i), getCol(i), getVal(i) + commonValue));
	for (int i = 0; i < M.noNonSparseValues; i++) {
		// if M has matching nonsparse value, subtract common value added above
		if (!S->add(M.getRow(i), M.getCol(i), M.getVal(i) - commonValue))
			S->getMatrix()->push_back(SparseRow(M.getRow(i), M.getCol(i), M.getVal(i) + commonValue));
	}
	S->setNoNSV(S->getMatrix()->size());
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

void SparseMatrix::setNoNSV(int size) {
	noNonSparseValues = size;
}

void SparseMatrix::setSparseRow(int pos, int r, int c, int v) {
	myMatrix->at(pos) = SparseRow(r, c, v);
}

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
	noNonSparseValues = myMatrix->size();
}

/*************** MAIN PROGRAM ***************
 *
 * Reads in two sparse matrices from the user, performs operations on the matrices
 * and displays them in.
 */
int main() {
	int n, m, cv;
	SparseMatrix* temp;

	cin >> n >> m >> cv;
	SparseMatrix* firstOne = new SparseMatrix(n, m, cv);
	(*firstOne).readInMatrix();
	// read in first and second matrix
	cin >> n >> m >> cv;
	SparseMatrix* secondOne = new SparseMatrix(n, m, cv);
	(*secondOne).readInMatrix();

	cout << "First one in sparse matrix format" << endl;
	cout << (*firstOne);

	cout << "First one in normal matrix format" << endl;
	(*firstOne).displayMatrix();

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

	cout << "Multiplication of matrices in sparse matrix form:" << endl;
	temp = (*secondOne)*(*firstOne);
	cout << (*temp);

	cout << "Addition of matrices in sparse matrix form:" << endl;
	temp = (*secondOne)+(*firstOne);
	cout << (*temp);

	system("pause");
	return 0;
}

// do we need display method
// what should I do about method that needs commonValue as DT?