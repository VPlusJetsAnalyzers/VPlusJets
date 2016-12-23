#include "vector"
#ifdef __CINT__ 
#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs;
#pragma link C++ class vector<double*>+;
#pragma link C++ class vector<double*>::*;
#ifdef G__VECTOR_HAS_CLASS_ITERATOR
#pragma link C++ operators vector<double*>::iterator;
#pragma link C++ operators vector<double*>::const_iterator;
#pragma link C++ operators vector<double*>::reverse_iterator;
#endif
#endif
