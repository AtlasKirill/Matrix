//
// Created by kirill on 20.11.16.
//

#ifndef CLASS_MATRIX_MATRIX_H
#define CLASS_MATRIX_MATRIX_H

#include <iostream>
#include <cstdio>
#include <cassert>
#include <utility>
#include <cmath>
#include <vector>
#include <cstring>

using namespace std;
namespace Arrays {
    template<typename T=double>
    class Matrix final {
    private:
        T **ptr;
        size_t lines;
        size_t columns;
    public:
        size_t get_lines() const { return lines; };

        size_t get_columns() const { return columns; };

        void get_and_put_ptr(T *);

        /*void get_and_put_column(double* buffer,int index){            //to get culumn number of index
            for (int i = 0; i <lines ; ++i)                             // it is unused
                buffer[i]=ptr[i][index];
        }*/

        Matrix(size_t, size_t, T *array);

        Matrix(size_t, size_t, vector<T> &array);

        Matrix(size_t, size_t);

        Matrix();    //default constructor
        Matrix(const Matrix<T> &rhs);    //copy constructor
        Matrix(Matrix<T> &&rhs);   //move-constructor

        ~Matrix();

        bool operator==(const Matrix<T> &rhs);

        bool operator!=(const Matrix<T> &rhs);

        Matrix<T> &operator+=(const Matrix<T> &rhs);

        Matrix<T> &operator-=(const Matrix<T> &rhs);

        Matrix<T> &operator*=(double a);

        Matrix<T> &operator=(const Matrix<T> &rhs);

        Matrix<T> &operator=(Matrix<T> &&rhs); //move-assign

        template<class X>
        friend ostream &operator<<(ostream &os, const Matrix<X> &rhs);

        Matrix<T> &transpose();

        Matrix<T> &inverse();

        Matrix<T> &mult_eq(const Matrix<T> &rhs);

        void cleanup();            //is used for deleting class field "ptr" in different methods
        void matrix_build(T *);      //is used in constructors
        void print() const;

        void kroneker(const Matrix<T> &b);

        double determinant();

        double Hilbert_generate(size_t n);   //makes the matrix of the Hilbert matrix and count its determinant

    };

    template<typename T>
    Matrix<T>::Matrix() {
        ptr = nullptr;
        lines = 0;
        columns = 0;
    }     //default ctor

    template<typename T>
    Matrix<T>::Matrix(size_t lines, size_t columns, T *array) : lines(lines), columns(columns) {
        assert(array);
        matrix_build(array);
    }

    template<typename T>
    Matrix<T>::Matrix(size_t lines, size_t columns, vector<T> &vector) : lines(lines), columns(columns) {
        assert(vector.size() != 0);
        double array[vector.size()];
        for (int i = 0; i < vector.size(); ++i)
            array[i] = vector[i];
        matrix_build(array);
    }

    template<typename T>
    Matrix<T>::Matrix(size_t lines, size_t columns) : lines(lines), columns(columns) {
        ptr = new T *[lines];
        for (size_t i = 0; i < lines; ++i) {
            ptr[i] = new T[columns];
        }

        for (size_t i = 0; i < lines; ++i)
            for (size_t j = 0; j < columns; ++j)
                ptr[i][j] = 0;
    }

    template<typename T>
    Matrix<T>::Matrix(const Matrix<T> &rhs) : lines(rhs.lines), columns(rhs.columns) {
        if (rhs.ptr) {
            ptr = new T *[lines];
            for (size_t i = 0; i < lines; ++i) {
                ptr[i] = new T[columns];
            }

            for (size_t i = 0; i < lines; ++i)
                for (size_t j = 0; j < columns; ++j)
                    ptr[i][j] = rhs.ptr[i][j];
        } else
            ptr = rhs.ptr;
    }

    template<typename T>
    Matrix<T>::Matrix(Matrix<T> &&rhs):lines(0), columns(0), ptr(nullptr) {
        lines = rhs.lines;
        columns = rhs.columns;
        ptr = rhs.ptr;
        rhs.lines = rhs.columns = 0;
        rhs.ptr = nullptr;
    }

    template<typename T>
    Matrix<T>::~Matrix() {
        cleanup();
    }

//------------------------------------------------------------ Operators
    template<typename T>
    bool Matrix<T>::operator==(const Matrix<T> &rhs) {
        if ((lines == rhs.lines) && (columns == rhs.columns)) {
            for (size_t i = 0; i < lines; ++i)
                for (size_t j = 0; j < columns; ++j)
                    if (ptr[i][j] != rhs.ptr[i][j])
                        return false;
            return true;
        }
        return false;
    }

    template<typename T>
    bool Matrix<T>::operator!=(const Matrix<T> &rhs) {
        return !(*this == rhs);
    }

/////////////////////////////////////////////
    template<typename T>
    Matrix<T> &Matrix<T>::operator=(Matrix<T> &&rhs) {     //move-assign
        if (this != &rhs) {
            this->cleanup();
            ptr = rhs.ptr;
            lines = rhs.lines;
            columns = rhs.columns;
            rhs.lines = rhs.columns = 0;
            rhs.ptr = nullptr;
        }
        return *this;
    }

    template<typename T>
    Matrix<T> &Matrix<T>::operator=(const Matrix<T> &rhs) {
        if ((this == &rhs) || (*this == rhs))
            return *this;
        Matrix<T> tmp(rhs);
        *this = std::move(tmp);
        return *this;
    }

/////////////////////////////////////////////
    template<typename T>
    Matrix<T> &Matrix<T>::operator+=(const Matrix<T> &rhs) {
        assert((lines == rhs.lines) && (columns == rhs.columns));
        for (size_t i = 0; i < lines; ++i)
            for (size_t j = 0; j < columns; ++j)
                ptr[i][j] = ptr[i][j] + rhs.ptr[i][j];
        return *this;
    }

    template<typename T>
    Matrix<T> &Matrix<T>::operator-=(const Matrix<T> &rhs) {
        assert((lines == rhs.lines) && (columns == rhs.columns));
        for (size_t i = 0; i < lines; ++i)
            for (size_t j = 0; j < columns; ++j)
                ptr[i][j] = ptr[i][j] - rhs.ptr[i][j];
        return *this;
    }

    template<typename T>
    Matrix<T> &Matrix<T>::operator*=(double a) {
        for (size_t i = 0; i < lines; ++i)
            for (size_t j = 0; j < columns; ++j)
                ptr[i][j] *= static_cast<T>(a);
        return *this;
    }

    template<typename T>
    Matrix<T> operator+(const Matrix<T> &lhs, const Matrix<T> &rhs) {
        assert((lhs.get_lines() == rhs.get_lines()) && (lhs.get_columns() == rhs.get_columns()));
        Matrix<T> tmp(std::move(lhs));
        tmp += rhs;
        return tmp;
    }

    template<typename T>
    Matrix<T> operator-(const Matrix<T> &lhs, const Matrix<T> &rhs) {
        assert((lhs.get_lines() == rhs.get_lines()) && (lhs.get_columns() == rhs.get_columns()));
        Matrix<T> tmp(std::move(lhs));
        tmp -= rhs;
        return tmp;
    }

    template<typename T>
    Matrix<T> operator*(double a, const Matrix<T> &rhs) {
        Matrix<T> tmp(std::move(rhs));
        tmp *= static_cast<T>(a);
        return tmp;
    }

    template<typename T>
    Matrix<T> operator*(const Matrix<T> &rhs, double a) {
        Matrix<T> tmp(std::move(rhs));
        tmp *= static_cast<T>(a);
        return tmp;
    }

    template<typename T>
    ostream &operator<<(ostream &os, const Matrix<T> &rhs) {
        assert((rhs.columns != 0) && (rhs.lines != 0));
        for (size_t i = 0; i < rhs.lines; ++i) {
            for (size_t j = 0; j < rhs.columns; ++j) {
                os.width(6);
                os << rhs.ptr[i][j] << " ";
            }
            os << endl;
        }
        os << endl;
        return os;
    }

//--------------------------------------------  methods

    template<typename T>
    void Matrix<T>::get_and_put_ptr(T *bufer) {
        assert(bufer);
        T *m_ptr = new T[lines * columns]{};
        int k = 0;
        for (int i = 0; i < lines; ++i)
            for (int j = 0; j < columns; ++j) {
                m_ptr[k] = ptr[i][j];
                k++;
            }
        for (int l = 0; l < lines * columns; ++l)
            bufer[l] = m_ptr[l];
        delete[] m_ptr;
    }

    template<typename T>
    void Matrix<T>::cleanup() {
        if (ptr) {
            for (size_t i = 0; i < lines; ++i)
                delete[] ptr[i];
        };
        delete[] ptr;
    }

    template<typename T>
    void Matrix<T>::matrix_build(T *array) {
        ptr = new T *[lines];
        for (size_t i = 0; i < lines; ++i) {
            ptr[i] = new T[columns];
        }
        size_t k = 0;
        for (size_t i = 0; i < lines; ++i)
            for (size_t j = 0; j < columns; ++j) {
                ptr[i][j] = array[k];
                k++;
            }
    }

    template<typename T>
    Matrix<T> &Matrix<T>::transpose() {
        Matrix<T> tmp(columns, lines);
        for (unsigned i = 0; i < columns; i++)
            for (unsigned j = 0; j < lines; j++) {
                tmp.ptr[i][j] = ptr[j][i];
            }
        *this = move(tmp);
        return *this;
    }

    template<typename T>
    Matrix<T> &Matrix<T>::inverse() {
        assert(lines == columns);
        assert(determinant() != 0);

// create the same matrix
        Matrix<T> H(lines, lines);
        for (unsigned i = 0; i < lines; i++)
            for (unsigned j = 0; j < lines; j++)
                H.ptr[i][i] = 1.0;

// direct course of Gauss method
        double temp;
        for (unsigned k = 0; k < lines; k++) {
            temp = ptr[k][k];

            for (unsigned j = 0; j < lines; j++) {
                ptr[k][j] /= temp;
                H.ptr[k][j] /= temp;
            }

            for (unsigned i = k + 1; i < lines; i++) {
                temp = ptr[i][k];

                for (unsigned j = 0; j < lines; j++) {
                    ptr[i][j] -= ptr[k][j] * temp;
                    H.ptr[i][j] -= H.ptr[k][j] * temp;
                }
            }
        }
// inverse course of Gauss method
        for (unsigned k = lines - 1; k > 0; k--)
            for (int i = k - 1; i >= 0; i--) {
                temp = ptr[i][k];
                for (unsigned j = 0; j < lines; j++) {
                    ptr[i][j] -= ptr[k][j] * temp;
                    H.ptr[i][j] -= H.ptr[k][j] * temp;
                }
            }
        *this = move(H);
        return *this;
    }

    template<typename T>
    double Matrix<T>::determinant() {
        assert(lines == columns);
        Matrix<T> tmp(*this);
        size_t n = lines;
        double det = 1;
        for (size_t k = 0; k < n; k++)    // Forward Elimination
        {
            for (size_t j = k + 1; j < n; ++j) {
                size_t m = k + 1;
                while (tmp.ptr[k][k] == 0) {     //swapping
                    if (m == n)
                        break;
                    if (tmp.ptr[m][k] != 0) {
                        std::swap(tmp.ptr[k], tmp.ptr[m]);
                        det *= -1;
                    }
                    m++;
                }
                double t1 = tmp.ptr[k][k];
                double t2 = tmp.ptr[j][k];
                for (size_t i = k; i < n; ++i) {
                    if (t1 != 0)
                        tmp.ptr[j][i] = (t1 * tmp.ptr[j][i] - t2 * tmp.ptr[k][i]) / t1;
                    else
                        tmp.ptr[j][i] = t1 * tmp.ptr[j][i] - t2 * tmp.ptr[k][i];
                }
            }
        }
        for (size_t i = 0; i < n; ++i)   //counting determinant
            det *= tmp.ptr[i][i];
        if (det == -0)
            det = 0;
        return det;
    }

    template<typename T>
    Matrix<T> &Matrix<T>::mult_eq(const Matrix<T> &rhs) {
        assert(columns == rhs.lines);
        Matrix<T> product(lines, rhs.columns);
        for (size_t i = 0; i < lines; ++i)
            for (size_t j = 0; j < rhs.columns; ++j)
                for (size_t inner = 0; inner < columns; inner++)
                    product.ptr[i][j] += ptr[i][inner] * rhs.ptr[inner][j];
        *this = std::move(product);
        return *this;
    }

    template<typename T>
    Matrix<T> mult(const Matrix<T> &lhs, const Matrix<T> &rhs) {
        assert(lhs.get_columns() == rhs.get_lines());
        Matrix<T> tmp(std::move(lhs));
        tmp.mult_eq(rhs);
        return tmp;
    }

    template<typename T>
    void Matrix<T>::print() const {
        assert((columns != 0) && (lines != 0));
        for (size_t i = 0; i < lines; ++i) {
            for (size_t j = 0; j < columns; ++j) {
                cout.width(6);
                cout << ptr[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

    template<typename T>
    void Matrix<T>::kroneker(const Matrix<T> &b) {

        Matrix<T> **block = new Matrix<T> *[lines];                         //"block" contain blocks of matrices
        for (size_t k = 0; k < lines; ++k) {
            block[k] = new Matrix<T>[columns];
        }
        for (size_t i = 0; i < lines; ++i)
            for (size_t j = 0; j < columns; ++j)
                block[i][j] = b * ptr[i][j];

        Matrix<T> temp(lines * b.lines, columns * b.columns);
        *this = temp;

        for (size_t i = 0; i < lines; ++i)
            for (size_t j = 0; j < columns; ++j)
                ptr[i][j] = block[i * b.lines / lines][j * b.columns / columns].ptr[i % b.lines][j % b.columns];


        for (size_t i = 0; i < block[0][0].lines; i++) {
            delete[] block[i];
        }
        delete[] block;

    }

    auto factorial = [](size_t n) {
        unsigned long long r;
        for (r = 1; n > 1; r *= (n--));
        return r;
    };

    template<typename T>
    double Matrix<T>::Hilbert_generate(size_t n) {
        Matrix<T> tmp(n, n);
        *this = std::move(tmp);
        for (size_t i = 0; i < n; ++i)
            for (size_t j = 0; j < n; ++j)
                ptr[i][j] = 1.0 / (i + j + 1.0);
        double det = 1.0;
        for (size_t i = 2; i < n; ++i)
            det *= factorial(i);
        det = pow(det, 3);
        for (size_t i = n; i < 2 * n; ++i)
            det /= static_cast<double>(factorial(i));
        return det;
    }
}
#endif //CLASS_MATRIX_MATRIX_H