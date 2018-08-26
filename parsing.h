//
// Created by kirill on 15.03.17.
//

#ifndef CLASS_MATRIX_PARSING_H
#define CLASS_MATRIX_PARSING_H

#include <vector>
#include <dirent.h>
#include "Circuit.h"


using namespace Arrays;
using namespace Current;
using namespace std;


namespace Parse{
    class Parsing {
    private:
        string str_;
        int index;

    public:
        Parsing();
        void skip_spaces();
        bool parse_digit(size_t& DIGIT);
        bool parse_number(size_t & NUM);
        bool parse_dashes();
        bool parse_comma();
        bool parse_real(double&R);
        bool parse_scolon();
        bool parse_string(size_t &Q1, size_t &Q2, double &R, double &J, double &E);
        Circuit * parsing(const string);
    };
    void fill_conduct_array(size_t edge, vector<double> &r, vector<double> &y);                          //this...
    void fill_connection_array(size_t edge, vector <double> &a, vector<size_t> &q1, vector<size_t> &q2);              //...this...
    void fill_full_connection_array(size_t edge, vector <double> &a, vector<size_t> &q1, vector<size_t> &q2);         //...and this methods are not necessary in the class Parsing
}

#endif //CLASS_MATRIX_PARSING_H
