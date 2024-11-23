#include <cmath>
#include <complex>
#include <iostream>
#include <vector>
using namespace std;

#define matrix_type complex<double>

class Matrix {
   public:
    int rows, cols, size = rows * cols;
    vector<matrix_type> data;

    // Constructor
    Matrix(int r = 0, int c = 0, const vector<matrix_type>& d = {})
        : rows(r), cols(c), data(d) {
        if (data.empty()) {
            data.resize(rows * cols, 0);
        }
    }

    // get matrix element
    matrix_type getElement(int r, int c) const {
        if (r < 0 || r >= rows || c < 0 || c >= cols) {
            cerr << "Invalid indices for matrix element access" << endl;
            return NAN;  // Return NaN for invalid access
        }
        return data[r * cols + c];
    }

    // print matrix
    void print() const {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                cout << getElement(i, j) << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

    // matrix multiply
    Matrix operator*(const Matrix& other) const {
        if (cols != other.rows) {
            cerr << "Matrix dimensions don't match for multiplication" << endl;
            return Matrix(0, 0, {});
        }
        vector<matrix_type> result(rows * other.cols, 0);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < other.cols; ++j) {
                for (int k = 0; k < cols; ++k) {
                    result[i * other.cols + j] +=
                        data[i * cols + k] * other.data[k * other.cols + j];
                }
            }
        }
        return Matrix(rows, other.cols, result);
    }

    Matrix& operator*=(const Matrix& other) {
        if (cols != other.rows) {
            cerr << "Matrix dimensions don't match for multiplication" << endl;
            return *this;
        }
        vector<matrix_type> result(rows * other.cols, 0);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < other.cols; ++j) {
                for (int k = 0; k < cols; ++k) {
                    result[i * other.cols + j] +=
                        data[i * cols + k] * other.data[k * other.cols + j];
                }
            }
        }
        data = result;
        cols = other.cols;
        return *this;
    }

    // matrix multiply(scalar)
    Matrix operator*(matrix_type scalar) const {
        vector<matrix_type> result(data);
        for (auto& value : result) {
            value *= scalar;
        }
        return Matrix(rows, cols, result);
    }

    friend Matrix operator*(matrix_type scalar, const Matrix& matrix) {
        return matrix * scalar;
    }

    Matrix& operator*=(complex<double> scalar) {
        *this = *this * scalar;
        return *this;
    }

    // matrix division(scalar)
    Matrix operator/(matrix_type scalar) const {
        vector<matrix_type> result(data);
        for (auto& value : result) {
            value /= scalar;
        }
        return Matrix(rows, cols, result);
    }

    Matrix& operator/=(complex<double> scalar) {
        *this = *this / scalar;
        return *this;
    }

    // matrix addition
    Matrix operator+(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            cerr << "Matrix dimensions don't match for addition" << endl;
            return Matrix(0, 0, {});
        }
        vector<matrix_type> result(rows * cols);
        for (int i = 0; i < rows * cols; ++i) {
            result[i] = data[i] + other.data[i];
        }
        return Matrix(rows, cols, result);
    }

    Matrix& operator+=(const Matrix& other) {
        *this = *this + other;
        return *this;
    }

    // matrix subtraction
    Matrix operator-(const Matrix& other) const {
        if (rows != other.rows || cols != other.cols) {
            cerr << "Matrix dimensions don't match for subtraction" << endl;
            return Matrix(0, 0, {});
        }
        vector<matrix_type> result(rows * cols);
        for (int i = 0; i < rows * cols; ++i) {
            result[i] = data[i] - other.data[i];
        }
        return Matrix(rows, cols, result);
    }

    Matrix& operator-=(const Matrix& other) {
        *this = *this - other;
        return *this;
    }

    // transpose
    Matrix T() const {
        vector<matrix_type> result(cols * rows);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                result[j * rows + i] = data[i * cols + j];
            }
        }
        return Matrix(cols, rows, result);
    }

    // matrix shape
    void shape() const { cout << rows << "x" << cols << endl; }
};

// dagger
Matrix operator!(Matrix a) {
    a = a.T();
    vector<matrix_type> result(a.rows * a.cols);
    for (int i = 0; i < a.cols * a.rows; i++) {
        result[i] = conj(a.data[i]);
    }
    return Matrix(a.rows, a.cols, result);
}

// trace
complex<double> Tr(Matrix a) {
    complex<double> sum = 0.0;
    if (a.rows != a.cols) {
        cerr << "Matrix dimensions don't match for Trace" << endl;
        return 0;
    } else {
        for (int i = 0; i < a.rows; i++) {
            sum += a.data[i * a.cols + i];
        }
        return sum;
    }
}

// 2x2 determinant
complex<double> det(Matrix a) {
    if (a.rows != a.cols || a.rows != 2) {
        cerr << "can't calculate determinant, only for 2x2" << endl;
        return 0;
    } else {
        return a.data[0] * a.data[3] - a.data[1] * a.data[2];
    }
}

/*
2x2 pauli matrix
input = 0,1,2,3
correp. unit,x,y,z
 */
Matrix pauli(int mu) {
    matrix_type i(0, 1);
    switch (mu) {
        case 0:
            return Matrix(2, 2, {1, 0, 0, 1});
        case 1:
            return Matrix(2, 2, {0, 1, 1, 0});
        case 2:
            return Matrix(2, 2, {0, -i, i, 0});
        case 3:
            return Matrix(2, 2, {1, 0, 0, -1});
        default:
            cerr << "no such matrix" << endl;
            return Matrix(0, 0);
    }
}

// zero matrix
Matrix zeros(int row, int col = 0) {
    if (col == 0) {
        col = row;
    }
    return Matrix(row, col, {});
}

// one matrix
Matrix ones(int row, int col = 0) {
    if (col == 0) {
        col = row;
    }
    return Matrix(row, col, vector<matrix_type>(row * col, 1));
}

// identity mateix
Matrix I(int row) {
    vector<matrix_type> result(row * row, 0);
    for (int i = 0; i < row * row; i += (row + 1)) {
        result[i] = 1;
    }
    return Matrix(row, row, result);
}