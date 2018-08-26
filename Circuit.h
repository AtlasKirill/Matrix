//
// Created by kirill on 09.02.17.
//

#ifndef CLASS_MATRIX_CIRCUIT_H
#define CLASS_MATRIX_CIRCUIT_H

#include "Matrix.h"
#include <fstream>
#include <string>
#include <cctype>
#include <cstdlib>

using std::vector;
using std::cout;

using namespace Arrays;
namespace Current{
    class Circuit final{
    private:
        Matrix<> R;Matrix<> Y;Matrix<> A;Matrix<> A_FULL;Matrix<> E;Matrix<> J;Matrix<> U0;Matrix<> I0;
        size_t nodes, edge;
    public:
        Circuit(vector<double> &a, vector<double> &a_full, vector<double> &y, vector<double> &r, vector<double> &e, vector<double> &j,size_t & nodes, size_t & edge);
        Circuit(const Circuit&lhs);

        void calculate_voltage();
        void calculate_currents();
        void print(){cout<<"Potentials are:"<<endl; cout<<U0; cout<<"Currents are:"<<endl; cout<<I0;}

    };
}

#endif //CLASS_MATRIX_CIRCUIT_H
