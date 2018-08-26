//
// Created by kirill on 09.02.17.
//
#include "Circuit.h"

using namespace std;

namespace Current{
    Circuit::Circuit(const Circuit &lhs) { A=lhs.A; R=lhs.R; Y=lhs.Y; E=lhs.E; J=lhs.J; U0=lhs.U0; I0=lhs.I0;}

    Circuit::Circuit(vector<double> &a, vector<double> &a_full, vector<double> &y, vector<double> &r, vector<double> &e, vector<double> &j,size_t & nodes, size_t & edge): nodes(nodes), edge(edge) {
        //assert(a); assert(a_full); assert(y); assert(r); assert(e); assert(j);
        A_FULL={nodes,edge,a_full};
        A={(nodes-1),edge,a};      //here matrix of connections (A) aren't full
        Y={edge,edge,y};
        R={edge,1,r}; E={edge,1,e}; J={edge,1,j}; U0={(nodes-1),1}; I0={edge,1};             //here matrix of potentials (U0) aren't full
    }
    void Circuit::calculate_voltage() {
        /*using formula     A*Y*At*U0=-A(J+Y*E)    */
        Matrix<>A_copy=A;
        Matrix<> At=A_copy.transpose(); Matrix<> a1=mult(A,Y); a1=mult(a1,At); a1=a1.inverse();
        Matrix<> a2=mult(Y,E); a2+=J; a2=mult(A,a2); a2=(-1)*a2;
        U0=mult(a1,a2);
        double * u0=new double [nodes]{};
        U0.get_and_put_ptr(u0+1);    //here pass u0+1 because we want to write from second cell, the first cell is the first node with zero potential
        U0={nodes,1,u0};     //here matrix of potentials is full
        delete[]u0;
    }
    void Circuit::calculate_currents(){
        Matrix<> A_F=A_FULL;    //A_F is needed not to change matrix A_FULL
        A_F.transpose();
        Matrix<> A_fi=mult(A_F,U0);

        A_fi+=E;
        size_t lines=A_fi.get_lines();
        double A_fi_ptr [lines]{};    //here will be stored the result of sum of potentials and voltages of voltage sources
        double Y_ptr [Y.get_lines()*Y.get_columns()]{};     //here will be stored the result of division of potentials with voltages by resistance
        A_fi.get_and_put_ptr(A_fi_ptr);
        Y.get_and_put_ptr(Y_ptr);
        double A_fi_ptr2 [lines*lines]{};
        double Y_ptr2 [Y.get_lines()]{};

        for (int i = 0; i <lines; ++i)
            A_fi_ptr2[i*lines+i]=A_fi_ptr[i];
        for (int i = 0; i <Y.get_lines(); ++i)
            Y_ptr2[i]=Y_ptr[i*Y.get_lines()+i];

        A_fi={lines,lines,A_fi_ptr2};
        Y={lines,1,Y_ptr2};                      //????

        A_fi=mult(A_fi,Y);
        A_fi+=J;
        I0=A_fi;
    }
}

