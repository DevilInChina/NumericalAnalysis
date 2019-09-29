#include "baseOpt.h"
const int maxn = 1024*1024;
int main(int argc,char **argv) {

    init();
    int cnts = 0;
    int n, m;
    type *A, *b;
    type *x;
    loadA_b(argv[2], &A, &b, &n);
    x = malloc(sizeof(type) * n);
    struct timeval tick, tickend;
    for (int i = 0; i < n; ++i)x[i] = 0.0;
    clock_t tim = 0;
    int cnt = 0;
    int maxcnt = MAX_ITER;
    if (strcmp(argv[1], "doo") == 0) {
        tim = clock();
        gettimeofday(&tick, NULL);
        lu_doolittle(A, x, b, n);
    } else if (strcmp(argv[1], "gau") == 0) {
        tim = clock();
        gettimeofday(&tick, NULL);
        gauss(A, x, b, n);

    } else if (strcmp(argv[1], "cro") == 0) {
        tim = clock();
        gettimeofday(&tick, NULL);
        lu_crout(A, x, b, n);
    } else if (strcmp(argv[1], "cho") == 0) {
        tim = clock();
        gettimeofday(&tick, NULL);
        cholesky(A, x, b, n);
    } else if (strcmp(argv[1], "jac") == 0) {
        tim = clock();
        gettimeofday(&tick, NULL);
        jacobi(A, x, b, n, &cnt, &maxcnt, THRESHOLD);
    } else if (strcmp(argv[1], "gs") == 0) {
        tim = clock();
        gettimeofday(&tick, NULL);
        gs(A, x, b, n, &cnt, &maxcnt, THRESHOLD);
    } else if (strcmp(argv[1], "sor") == 0) {
        tim = clock();
        gettimeofday(&tick, NULL);
#ifdef MULTHREAD
        type w = 1;
        type dep = 0.1;
        MatInfo k;
        k.A = A;
        k.b = b;
        k.x = NULL;
        k.n = n;
        k.w = w;
        pthread_t tid = pthread_self();
        for(int i = 0 ; i < 11 ; ++i) {
            pthread_attr_t s;
            double* ret;
            int err = pthread_create(&tid, ret, (void *) (Thread), &k);
            k.w+=dep;
            printf("%lf\n",ret);
            if(err){
                printf("err\n");
            }
        }
#else
        sor(A, x, b, n, &cnt, &maxcnt, THRESHOLD, 1.5);
#endif
    } else if (strcmp(argv[1], "cg") == 0) {
        tim = clock();
        gettimeofday(&tick, NULL);
        cg(A, x, b, n, &cnt, &maxcnt, THRESHOLD);
    }

    tim = clock() - tim;
    gettimeofday(&tickend, NULL);
    double s = dealRes(A, x, b, cnt, maxcnt, n);
    ///printf("%lf\n",s);
    if (isnormal(s)) {
        if (s < RESPRE) {
            printf("%s %s ", argv[1], argv[2]);
            printf("TrueTime:%.2lfms CpuTime:%.2fms %.7f %d\n",
                   (tickend.tv_usec - tick.tv_usec) / 1000.0 + 1000.0 * (tickend.tv_sec - tick.tv_sec),
                   1000.0 * tim / CLOCKS_PER_SEC, s, cnt);
        }
    }
    free(A);
    free(b);
    free(x);
}