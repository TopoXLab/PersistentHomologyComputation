%module persis_homo_optimal

%include "std_string.i"
%include "std_vector.i"

using namespace std;

%template(cpp_nested1_IntVector) vector<int>;
%template(cpp_nested2_InvVector) vector<vector<int>>;
%template(cpp_nested3_IntVector) vector<vector<vector<int>>>;
%template(cpp_nested4_IntVector) vector<vector<vector<vector<int>>>>;

%template(cpp_nested1_DoubleVector) vector<double>;
%template(cpp_nested2_DoubleVector) vector<vector<double>>;
%template(cpp_nested3_DoubleVector) vector<vector<vector<double>>>;

%{
#include "PersistenceComputer.h"
%}

//double-check that this is indeed %include !!!
%include "PersistenceComputer.h"