#include <iostream>
#include <cstring>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <sys/time.h>

using namespace std;

vector<int> make_vector(const int i, const int j, const int v) {
	vector<int> result;
	result.clear();
	result.push_back(i);
	result.push_back(j);
	result.push_back(v);
	return result;
}

double mean(const vector<int> &arr) {
	double result = 0.0;
	for (int i = 0; i < arr.size(); ++i) {
		result += arr[i];
	}
	return result / arr.size();
}

double standard_deviation(const vector<int> &arr) {
	const double mu = mean(arr);
	double result = 0.0;
	for (int i = 0; i < arr.size(); ++i) {
		result += (arr[i] - mu) * (arr[i] - mu);
	}
	return sqrt(result / arr.size());
}

// Get the millisecond
void time_ms(long int &ms) {
	struct timeval tp;
	gettimeofday(&tp, NULL);
	ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
}

// Difference in milliseconds
long int difftime_ms(long int &start, long int &end) {
	return end - start;
}

// +--------------+
// | Class Matrix |
// +--------------+

class Matrix {
public:
	int nRows;
	int nColumns;
	int size;

	vector< pair<int, int> > row_by_row_access() {
		vector< pair<int, int> > result;
		result.clear();
		for (int row = 0; row < nRows; ++row) {
			for (int column = 0; column < nColumns; ++column) {
				result.push_back(make_pair(row, column));
			}
		}
		return result;
	}

	vector< pair<int, int> > column_by_column_access() {
		vector< pair<int, int> > result;
		result.clear();
		for (int column = 0; column < nColumns; ++column) {
			for (int row = 0; row < nRows; ++row) {
				result.push_back(make_pair(row, column));
			}
		}
		return result;
	}

	vector< pair<int, int> > random_access(const int &N) {
		vector< pair<int, int> > result;
		result.clear();
		for (int i = 0; i < N; ++i) {
			const int r = rand() % nRows;
			const int c = rand() % nColumns;
			result.push_back(make_pair(r, c));
		}
		return result;
	}
};

// +------------------------------------------+
// | Class Matrix1 - Single Dimensional Array |
// +------------------------------------------+

template<typename TYPE>
class Matrix1 : public Matrix {
public:
	Matrix1(const int &nRows, const int &nColumns) {
		this -> nRows = nRows;
		this -> nColumns = nColumns;
		size = this -> nRows * this -> nColumns;

		arr = new TYPE [size];
	}

	TYPE at(const int i, const int j) {
		return arr[i * nRows + j];
	}

	~Matrix1() {
		delete[] arr;
	}

	// Single dimensional array
	TYPE *arr;

	TYPE sum(const vector< pair<int, int> > &access) {
		TYPE result = 0;
		for (int row = 0; row < nRows; ++row) {
			for (int column = 0; column < nColumns; ++column) {
				result += *this(row, column);
			}
		}
		return result;
	}
};

// +--------------------------------+
// | Class Matrix2 - Array Of Array |
// +--------------------------------+

template<typename TYPE>
class Matrix2 : public Matrix {
public:
	Matrix2(const int &nRows, const int &nColumns) {
		this -> nRows = nRows;
		this -> nColumns = nColumns;
		size = this -> nRows * this -> nColumns;

		arr = new TYPE* [this -> nRows];
		for (int row = 0; row < this -> nRows; ++row) {
			arr[row] = new TYPE [this -> nColumns];
		}
	}

	TYPE at(const int i, const int j) {
		return arr[i][j];
	}

	~Matrix2() {
		for (int row = 0; row < nRows; ++row) {
			delete[] arr[row];
		}
		delete[] arr;
	}

	// Array of array
	TYPE **arr;
};

// +------------------------+
// | Sum Operation - Matrix |
// +------------------------+

template<typename TYPE>
TYPE sum_op(const Matrix1<TYPE> &mat1, const vector< pair<int, int> > &access) {
	TYPE result = 0;
	for (int i = 0; i < access.size(); ++i) {
		result += mat1.at(access[i].first, access[i].second);
	}
	return result;
}

template<typename TYPE>
TYPE sum_op(const Matrix2<TYPE> &mat2, const vector< pair<int, int> > &access) {
	TYPE result = 0;
	for (int i = 0; i < access.size(); ++i) {
		result += mat2.at(access[i].first, access[i].second);
	}
	return result;
}

// +--------------+
// | Class Tensor |
// +--------------+

class Tensor {
public:
	int nRows;
	int nColumns;
	int nDepth;
	int size;

	vector< vector<int> > sequential_access() {
		vector< vector<int> > result;
		result.clear();
		for (int row = 0; row < nRows; ++row) {
			for (int column = 0; column < nColumns; ++column) {
				for (int depth = 0; depth < nDepth; ++depth) {
					result.push_back(make_vector(row, column, depth));
				}
			}
		}
		return result;
	}

	vector< vector<int> > random_access(const int &N) {
		vector< vector<int> > result;
		result.clear();
		for (int i = 0; i < N; ++i) {
			const int row = rand() % nRows;
			const int column = rand() % nColumns;
			const int depth = rand() % nDepth;
			result.push_back(make_vector(row, column, depth));
		}
		return result;
	}
};

// +------------------------------------------+
// | Class Tensor1 - Single Dimensional Array |
// +------------------------------------------+

template<typename TYPE>
class Tensor1 : public Tensor {
public:
	Tensor1(const int &nRows, const int &nColumns, const int &nDepth) {
		this -> nRows = nRows;
		this -> nColumns = nColumns;
		this -> nDepth = nDepth;
		stride = this -> nColumns * this -> nDepth;
		size = this -> nRows * this -> nColumns * this -> nDepth;

		arr = new TYPE [size];
	}

	TYPE at(const int i, const int j, const int v) {
		return arr[i * stride + j * nDepth + v];
	}

	TYPE& operator()(const int i, const int j, const int v) {
		return arr[i * stride + j * nDepth + v];
	}

	~Tensor1() {
		delete[] arr;
	}

	// Single dimensional array
	TYPE *arr;

	int stride;
};

// +-----------------------------------------+
// | Class Tensor2 - Array Of Array Of Array |
// +-----------------------------------------+

template<typename TYPE>
class Tensor2 : public Tensor {
public:
	Tensor2(const int &nRows, const int &nColumns, const int &nDepth) {
		this -> nRows = nRows;
		this -> nColumns = nColumns;
		this -> nDepth = nDepth;
		size = this -> nRows * this -> nColumns * this -> nDepth;

		arr = new TYPE** [this -> nRows];
		for (int row = 0; row < this -> nRows; ++row) {
			arr[row] = new TYPE* [this -> nColumns];
			for (int column = 0; column < this -> nDepth; ++column) {
				arr[row][column] = new TYPE [this -> nDepth];
			}
		}
	}

	TYPE at(const int i, const int j, const int v) {
		return arr[i][j][v];
	}

	TYPE& operator()(const int i, const int j, const int v) {
		return arr[i][j][v];
	}

	~Tensor2() {
		for (int row = 0; row < nRows; ++row) {
			for (int column = 0; column < nColumns; ++column) {
				delete[] arr[row][column];
			}
			delete[] arr[row];
		}
		delete[] arr;
	}

	// Array of array of array
	TYPE ***arr;
};

// +------------------------+
// | Sum Operation - Tensor |
// +------------------------+

template<typename TYPE>
TYPE sum_op(const Tensor1<TYPE> &ten1, const vector< vector<int> > &access) {
	TYPE result = 0;
	for (int i = 0; i < access.size(); ++i) {
		result += ten1.at(access[i][0], access[i][1], access[i][2]);
	}
	return result;
}

template<typename TYPE>
TYPE sum_op(const Tensor2<TYPE> &ten2, const vector< vector<int> > &access) {
	TYPE result = 0;
	for (int i = 0; i < access.size(); ++i) {
		result += ten2.at(access[i][0], access[i][1], access[i][2]);
	}
	return result;
}

// +---------------------+
// | Matrix Benchmarking |
// +---------------------+

template<typename TYPE>
void matrix_benchmarking(const int &nRows, const int &nColumns, const int &nTimes) {
	Matrix1<TYPE> mat1(nRows, nColumns);
	Matrix2<TYPE> mat2(nRows, nColumns);

	vector<int> rr1;
	vector<int> rr2;
	vector<int> cc1;
	vector<int> cc2;
	vector<int> random1;
	vector<int> random2;

	for (int t = 0; t < nTimes; ++t) {
		vector< pair<int, int> > row_by_row_access = mat1.row_by_row_access();
		vector< pair<int, int> > column_by_column_access = mat1.column_by_column_access();
		vector< pair<int, int> > random_access = mat1.random_access(nRows * nColumns);

		long int start, end;

		time_ms(start);
		sum_op<TYPE>(mat1, row_by_row_access);
		time_ms(end);
		rr1.push_back(difftime_ms(start, end));

		time_ms(start);
		sum_op<TYPE>(mat2, row_by_row_access);
		time_ms(end);
		rr2.push_back(difftime_ms(start, end));

		time_ms(start);
		sum_op<TYPE>(mat1, column_by_column_access);
		time_ms(end);
		cc1.push_back(difftime_ms(start, end));

		time_ms(start);
		sum_op<TYPE>(mat2, column_by_column_access);
		time_ms(end);
		cc2.push_back(difftime_ms(start, end));

		time_ms(start);
		sum_op<TYPE>(mat1, random_access);
		time_ms(end);
		random1.push_back(difftime_ms(start, end));

		time_ms(start);
		sum_op<TYPE>(mat2, random_access);
		time_ms(end);
		random2.push_back(difftime_ms(start, end));
	}

	cout << "*** Row by Row ***" << endl;
	cout << "Single dimensional array: mu = " << mean(rr1) << ", standard_deviation = " << standard_deviation(rr1) << endl;
	cout << "Array of array: mu = " << mean(rr2) << ", standard_deviation = " << standard_deviation(rr2) << endl << endl;

	cout << "*** Column by Column ***" << endl;
	cout << "Single dimensional array: mu = " << mean(cc1) << ", standard_deviation = " << standard_deviation(cc1) << endl;
	cout << "Array of array: mu = " << mean(cc2) << ", standard_deviation = " << standard_deviation(cc2) << endl << endl;

	cout << "*** Random ***" << endl;
	cout << "Single dimensional array: mu = " << mean(random1) << ", standard_deviation = " << standard_deviation(random1) << endl;
	cout << "Array of array: mu = " << mean(random2) << ", standard_deviation = " << standard_deviation(random2) << endl << endl;
}

// +---------------------+
// | Tensor Benchmarking |
// +---------------------+

template<typename TYPE>
void tensor_benchmarking(const int &nRows, const int &nColumns, const int &nDepth, const int &nTimes) {
	Tensor1<TYPE> ten1(nRows, nColumns, nDepth);
	Tensor2<TYPE> ten2(nRows, nColumns, nDepth);

	vector<int> seq1;
	vector<int> seq2;
	vector<int> random1;
	vector<int> random2;

	for (int t = 0; t < nTimes; ++t) {
		vector< vector<int> > sequential_access = ten1.sequential_access();
		vector< vector<int> > random_access = ten1.random_access(nRows * nColumns * nDepth);

		long int start, end;

		time_ms(start);
		sum_op<TYPE>(ten1, sequential_access);
		time_ms(end);
		seq1.push_back(difftime_ms(start, end));

		time_ms(start);
		sum_op<TYPE>(ten2, sequential_access);
		time_ms(end);
		seq2.push_back(difftime_ms(start, end));

		time_ms(start);
		sum_op<TYPE>(ten1, random_access);
		time_ms(end);
		random1.push_back(difftime_ms(start, end));

		time_ms(start);
		sum_op<TYPE>(ten2, random_access);
		time_ms(end);
		random2.push_back(difftime_ms(start, end));
	}

	cout << "*** Sequential access ***" << endl;
	cout << "Single dimensional array: mu = " << mean(seq1) << ", standard_deviation = " << standard_deviation(seq1) << endl;
	cout << "Array of array of array: mu = " << mean(seq2) << ", standard_deviation = " << standard_deviation(seq2) << endl << endl;

	cout << "*** Random access ***" << endl;
	cout << "Single dimensional array: mu = " << mean(random1) << ", standard_deviation = " << standard_deviation(random1) << endl;
	cout << "Array of array of array: mu = " << mean(random2) << ", standard_deviation = " << standard_deviation(random2) << endl << endl;
}

// +--------------+
// | Main Program |
// +--------------+

int main(int argc, char **argv) {
	cout << "--- Small size 100x100 - Integer ------------------------------" << endl;
	matrix_benchmarking<int>(100, 100, 10);

	cout << "--- Medium size 1000x1000 - Integer ------------------------------" << endl;
	matrix_benchmarking<int>(1000, 1000, 10);

	cout << "--- Large size 10000x10000 - Integer ------------------------------" << endl;
	matrix_benchmarking<int>(10000, 10000, 10);

	cout << "--- Small size 100x100 - Float ------------------------------" << endl;
	matrix_benchmarking<float>(100, 100, 10);

	cout << "--- Medium size 1000x1000 - Float ------------------------------" << endl;
	matrix_benchmarking<float>(1000, 1000, 10);

	cout << "--- Large size 10000x10000 - Float ------------------------------" << endl;
	matrix_benchmarking<float>(10000, 10000, 10);

	cout << "--- Small size 100x100 - Double ------------------------------" << endl;
	matrix_benchmarking<double>(100, 100, 10);

	cout << "--- Medium size 1000x1000 - Double ------------------------------" << endl;
	matrix_benchmarking<double>(1000, 1000, 10);

	cout << "--- Large size 10000x10000 - Double ------------------------------" << endl;
	matrix_benchmarking<double>(10000, 10000, 10);

	cout << "--- Tensor 100x100x100 - Integer -----------------------------------" << endl;
	tensor_benchmarking<int>(100, 100, 100, 10);

	cout << "--- Tensor 100x100x100 - Float -----------------------------------" << endl;
	tensor_benchmarking<float>(100, 100, 100, 10);

	cout << "--- Tensor 100x100x100 - Double -----------------------------------" << endl;
	tensor_benchmarking<double>(100, 100, 100, 10);

	return 0;
}