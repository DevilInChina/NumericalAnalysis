#include "baseOpt.h"


int main(int argc,char **argv){
    int cnts = 0;
    int n, m;

    type *A,*b;
    type *x;
    loadA_b(argv[2],&A,&b,&n);
    x= malloc(sizeof(type)*n);

    for(int i = 0 ; i < n ; ++i)x[i] = 0.0;
    clock_t tim =0;

    int cnt = 0;
    int maxcnt = MAX_ITER;
    if (strcmp(argv[1], "doo") == 0) {
        tim = clock();
        lu_doolittle(A,x,b,n);
    } else if (strcmp(argv[1], "gau") == 0) {
        tim = clock();
        gauss(A,x,b,n);

    } else if (strcmp(argv[1], "cro") == 0) {
        tim = clock();
        lu_crout(A,x,b,n);
    } else if(strcmp(argv[1],"cho")==0){
        tim = clock();
        cholesky(A,x,b,n);
    } else if(strcmp(argv[1],"jac")==0){
        tim = clock();
        jacobi(A, x, b, n, &cnt, &maxcnt, THRESHOLD);
    }else if(strcmp(argv[1],"gs")==0){
        tim=clock();
        gs(A, x, b, n, &cnt, &maxcnt, THRESHOLD);
    }else if(strcmp(argv[1],"sor")==0){
        tim = clock();
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
    }else if(strcmp(argv[1],"cg")==0){
        tim = clock();
        cg(A, x, b, n, &cnt, &maxcnt, THRESHOLD);
    }

    tim = clock() - tim;
    double s = dealRes(A,x,b,cnt,maxcnt,n);
    ///printf("%lf\n",s);
    if(isnormal(s)) {
        if(s<RESPRE) {
            printf("%s %s ", argv[1], argv[2]);
            printf("%.5f %.9f\n", 1000.0 * tim / CLOCKS_PER_SEC,s);
        }

    }
    free(A);
    free(b);
    free(x);
}