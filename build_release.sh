#!/bin/tcsh
setenv PREFIX "/usr/local"
setenv LIBDIR "$PREFIX/lib"
setenv LIB64DIR "$PREFIX/lib64"
setenv CC "gcc"
setenv CXX "g++"
setenv CPPFLAGS "-I$PREFIX/include"
setenv LDFLAGS "-L$LIBDIR -L$LIB64DIR"
cmake -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX -DCMAKE_BUILD_TYPE=Release 
