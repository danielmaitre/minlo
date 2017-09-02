%module MINLOReader

%include "stl.i"
%include "std_string.i"

%include "std_vector.i"

%include "carrays.i"
%array_functions(float, floatArray);


%{
#include "pdf.h"
#include "MinloReader.h"
%}

namespace std {
   %template(vector_str) vector<string>;
   %template(vector_d) vector<double>;
   %template(vector_i) vector<int>;
};

%include <typemaps.i>



%include "ntuplereader/NtupleInfo.h"
%include "ntuplereader/EventReaderBase.h"
%include "pdf.h"
%include "MinloInfo.h"

%apply double & OUTPUT { double& q0 };
%apply double & OUTPUT { double &scaleForNLO };

%include "MinloReader.h"


	
