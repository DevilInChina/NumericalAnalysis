gcc eigensolve.c -o eigen -lm -lpthread -fopenmp
./eigen $1 $2
rm eigen
