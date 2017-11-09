export OMP_NUM_THREADS=4
g++ -std=c++11 -fopenmp -O3 main.cpp -o omp_test