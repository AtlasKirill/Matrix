#include <sstream>
#include "parsing.h"

//using namespace Arrays;
//using namespace Current;
using namespace Parse;
using namespace std;


int main(int argc, char *argv[]) {
    int c=6;//atoi(argv[2]);
    if(c>6){
        cout<<"The input number of circuits is greater than number of examples"<<endl;
        return 1;
    }
    Parsing A;
    Circuit*D= nullptr;
    for (int i = 1; i <=c; ++i) {
        ostringstream file;
        file<<"/home/kirill/ClionProjects/class_matrix/input/input"<<i<<".txt";//string(argv[1])<<"/input"<<i<<".txt";
        cout<<"CIRCUIT NUMBER "<<i<<endl;
        D =A.parsing(file.str());
        D->calculate_voltage();
        D->calculate_currents();
        D->print();
        cout<<endl;
    }
    delete D;
//    vector<double> b = {1, 2, 3, 4, 5, 6};
//    Matrix<> A(2, 3);
//    Matrix<> B;
//    B = A;


    return 0;
}
