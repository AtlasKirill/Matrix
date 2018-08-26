//
// Created by kirill on 08.02.17.
//
#include "parsing.h"

using namespace Current;
using namespace Arrays;
using namespace std;

namespace Parse {
    Parsing::Parsing(){index=0;}

    void Parsing::skip_spaces() {
        while (isspace(str_[index]))
            index++;
    }

    bool Parsing::parse_digit(size_t &DIGIT) {
        const char *sample[] = {"0", "1", "2", "3", "4", "5", "6", "7", "8", "9"};
        for (int i = 0; i < 10; ++i)
            if (str_[index] == *sample[i]) {
                DIGIT = static_cast<size_t >(atoi(sample[i]));
                index++;
                return true;
            }
        return false;
    }

    bool Parsing::parse_number(size_t &NUM) {
        size_t DIGIT = 0;
        skip_spaces();
        if (!parse_digit(DIGIT))
            return false;
        NUM = DIGIT;
        while (parse_digit(DIGIT)) {
            NUM = NUM * 10 + DIGIT;
            index++;
        }
        return true;
    }

    bool Parsing::parse_dashes() {
        skip_spaces();
        if ((str_[index] == '-') && (str_[index + 1] == '-')) {
            index += 2;
            return true;
        }
        return false;
    }

    bool Parsing::parse_comma() {
        skip_spaces();
        if (str_[index] == ',') {
            index++;
            return true;
        }
        return false;
    }

    bool Parsing::parse_real(double &R) {
        skip_spaces();
        size_t digit = 0;

        int signum = 1;
        if (str_[index] == '-') {
            signum = -1;
            index++;
        }
        if (!parse_digit(digit))
            return false;
        R = digit;
        while (parse_digit(digit))
            R = R * 10 + digit;
        int fraction = -1;
        if (str_[index] == '.') {
            index++;
            while (parse_digit(digit)) {
                R += digit * pow(10, fraction);
                fraction--;
            }
            R *= signum;
            return true;
        }
        return false;
    }

    bool Parsing::parse_scolon() {
        skip_spaces();
        if (str_[index] == ';') {
            index++;
            return true;
        }
        return false;

    }

    bool Parsing::parse_string(size_t &Q1, size_t &Q2, double &R, double &J, double &E) {
        index = 0;
        if (!parse_number(Q1))
            return false;
        if (!parse_dashes())
            return false;
        if (!parse_number(Q2))
            return false;
        if (!parse_comma())
            return false;
        if (!parse_real(R))
            return false;
        assert(R != 0);
        if (!parse_scolon())
            return false;
        if ((str_.find('V') == -1) & (str_.find('A') != -1)) {
            if (!parse_real(J))
                return false;
            E = 0;
            skip_spaces();
            if (str_[index] == 'A')
                return true;
        }
        if ((str_.find('V') != -1) & (str_.find('A') == -1)) {
            if (!parse_real(E))
                return false;
            J = 0;
            skip_spaces();
            if (str_[index] == 'V')
                return true;
        }
        return (str_.find('V') == -1) & (str_.find('A') == -1);    //it is equivalent to // if ((str_.find('V') == -1) & (str_.find('A') == -1))
                                                                                         //     return true;
                                                                                         // return false;
    }

    void fill_conduct_array(size_t edge, vector<double> &r, vector<double> &y) {
        for (size_t i = 0; i < edge; i++)
            for (size_t j = 0; j < edge; j++)
                if (i == j)
                    y[j + i * edge] = 1 / r[i];
                else
                    y[j + i * edge] = 0;
    }

    void fill_connection_array(size_t edge, vector <double> &a, vector<size_t> &q1, vector<size_t> &q2) {

        for (int i = 0; i < edge; i++) {
            if (q1[i] != 1)
                a[i + q1[i] * edge - 2 * edge] = 1;
            if (q2[i] != 1)
                a[i + q2[i] * edge - 2 * edge] = -1;
        }
    }

    void fill_full_connection_array(size_t edge, vector <double> &a, vector<size_t> &q1, vector<size_t> &q2) {
        for (int i = 0; i < edge; i++) {
            a[i + q1[i] * edge - edge] = 1;
            a[i + q2[i] * edge - edge] = -1;
        }
    }


    Circuit * Parsing::parsing(const string File_path) {

        ifstream in(File_path);
        assert(in.is_open());
        size_t edge = 0;    //count a number of edges
        while (!in.eof())
            while (getline(in, str_))
                ++edge;

        //input.seekg(0, ios::beg);   //it doesn't work correct
        in.close();

        vector<size_t> q1(edge); // first colomn of nodes
        vector<size_t> q2(edge); // second colomn of nodes
        vector<double> r(edge); // colomn of resistances
        vector<double> y(edge*edge); // array of conductivities
        vector<double> j(edge); // colomn of currents
        vector<double> e(edge); // colomn of voltage

        ifstream input(File_path);
        assert(input.is_open());

        for (int k = 0; k < edge; ++k) {
            getline(input, str_);
            parse_string(q1[k], q2[k], r[k], j[k], e[k]);
        }

        input.close();
        //filling of conductivities array
        fill_conduct_array(edge, r, y);

        //search of a number of nodes
        size_t nodes = 0;
        for (size_t i = 0; i < edge; ++i) {
            if (nodes < q1[i])
                nodes = q1[i];
            if (nodes < q2[i])
                nodes = q2[i];
        }
        //filling of connections array
        vector <double> a (edge * (nodes - 1)); // array of connections       мы учитываем, что узлы нумеруются с единицы
        fill_connection_array(edge, a, q1, q2);
        vector <double> a_full (edge * nodes); // full array of connections
        fill_full_connection_array(edge, a_full, q1, q2);

        Circuit *T = new Circuit(a, a_full, y, r, e, j, nodes, edge);
        return T;
    }
}