//
// Created by devilinchina on 10/5/19.
//
#include "baseOpt.h"
const int maxs=1024*1024;

int main() {
    freopen("/home/devilinchina/CLionProjects/number/data.txt", "r", stdin);
    int n;
    scanf("%d", &n);

    MALLOC(x, type, n * n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            scanf("%f", x + i * n + j);
        }
    }
    MALLOC(ans, type, n);
    int maxss = 5000;
    qrqeigensolver(ans, x, n, &maxss, 0.0001);

    //descend_power(ans ,x, n, &maxss, 0.0001);

    showMtx(ans, 1, n);
//    showMtx(vec,1,n);
}
/// gcc eigensolve.c -o main -lm -fopenmp -lpthread
