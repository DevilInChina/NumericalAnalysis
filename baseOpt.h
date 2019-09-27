//
// Created by devilinchina on 9/12/19.
//
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <sys/time.h>
#include <unistd.h>
#define CORE_NUMBERS sysconf(_SC_NPROCESSORS_ONLN)
#ifndef NUMBER_BASEOPT_H
#define NUMBER_BASEOPT_H

#define type float
#define tempType int
#define true 1
#define false 0

#define USEPARALLEL
///#define DOT_PARALLEL

#define EPS 1e-5

#define THRESHOLD 1e-8
#define RESPRE 1e-4
#define N 10005
type temp[N];
#define CORE_NUMBER 1
#define MAX_ITER 1000

void init(){
   // CORE_NUMBER = CORE_NUMBERS;
}

//// one
void swap(type *a,type*b){
    *(tempType*)a^=*(tempType*)b^=*(tempType*)a^=*(tempType*)b;
}
void swapIJ(type *A,int i,int j,int beg,int n){
    type *a1 = A+n*j+beg;
    type *a2 = A+n*i+beg;
    n-=beg;
#ifdef USEPARALLEL
    #pragma omp parallel for
    for(int ii =0  ; ii < n ; ++ii){
        swap(a1+ii,a2+ii);
    }
#else
    while (a2!=en){
        swap(a1,a2);
        ++a1,++a2;
    }
#endif
}
void addKi_to_j(type *A ,type k, int i,int j,int beg,int n){
    type *a1 = A+n*j+beg;
    type *a2 = A+n*i+beg;
    type *en = A+n*i+n;
    n-=beg;
#ifdef USEPARALLEL
#pragma omp parallel for
    for(int ii =0  ; ii < n ; ++ii){
        a1[ii]+=k*a2[ii];
    }
#else
    while (a2!=en){
        *a1+=k*(*a2);
        ++a1,++a2;
    }
#endif
}
int isZero(const type Zero){
    return fabs(Zero)<EPS;
}
type temps[5000];


void ll_dotprod(type *result,type *a,type *b ,int len){
    *result = 0.0;
    type *ba = a;
    type *bb = b;
    type *ea = a+len;
    while (ba!=ea){
        *result+=(*ba)*(*bb);
        ++ba;
        ++bb;
    }
}


void parallel_dotprod(type *result,type*a,type*b,int len){
    type Tmps[CORE_NUMBER] ;//= malloc(sizeof(type)*(CORE_NUMBER));
    // memset(Tmps,0, sizeof(type)*(CORE_NUMBER));
    int tot = len/CORE_NUMBER;
#pragma omp parallel for
    for(int i = 0 ; i < CORE_NUMBER ; ++i){
        ll_dotprod(Tmps+i,a+i*tot,b+i*tot,tot);
    }
    *result = 0.0;
    ll_dotprod(result,a+CORE_NUMBER*tot,b+CORE_NUMBER*tot,len%CORE_NUMBER);
    for(int i = 0 ; i < CORE_NUMBER ; ++i){
        *result+=Tmps[i];
    }
}

void dotprod(type *result,type*a,type*b,int len){
#ifdef DOT_PARALLEL
    parallel_dotprod(result,a,b,len);
#else
    ll_dotprod(result,a,b,len);
#endif
}

type sumLrkUki( type*l, type*u,int r, int i,int n) {
    type re = 0.0;
    type *bel = l+r*n;
    type *beu = u+i*n;
    dotprod(&re,bel,beu,r);
    return re;
}

type sumLikUkr( type*l, type*u,int r, int i,int n){
    type re = 0.0;
    type *bel = l+i*n;
    type *beu = u+r*n;
    dotprod(&re,bel,beu,r);
    return re;
}

void calLow(const type*A,type*x,const type *b,int n){///no parallel
    for(int i = 0 ; i < n ; ++i){
        type k = 0;
        for(int j = 0 ; j < i ; ++j){
            k+=x[j]*A[i*n+j];
        }
        x[i] = (b[i]-k)/A[i*n+i];
    }
}

void calUp(const type *A,type*x,const type *b,int n){///no parallel
    for(int i = n-1 ; i >= 0 ; --i){
        type k = 0;
        for(int j = i ; j < n ; ++j){
            k+=x[j]*A[i+j*n];
        }
        x[i] = (b[i]-k)/A[i*n+i];
    }
}
type calSum(type *l,int i,int k,int n){
    type ret = 0.0;
    type *bel = l+i*n;
    type *enl = l+i*n+k;
    type *belk= l+k*n;
    while (bel!=enl){
        ret+=(*bel)*(*belk);
        ++bel;
        ++belk;
    }
    return ret;
}
type* loadMtx(const char * filePath,int *n,int *m){
    FILE *fp = fopen(filePath,"r");
    fscanf(fp,"%d %d",n,m);
    // printf("%d %d\n",*n,*m);
    long long k = 1ll*(*n)*(*m);
    type *mem = (type*)malloc(sizeof(type)*k);
    if(mem==NULL)return mem;
    type *en = mem+k;
    type *st = mem;
    while (st!=en){
        fscanf(fp,"%f",st);
        ++st;
    }
    return mem;
}

void loadA_b(char*filePath,type **A,type **b,int *n) {///alloc a memory
    FILE *fp = fopen(filePath,"r");
    fscanf(fp,"%d",n);
    int tot = *n*(*n);
    *A = malloc(sizeof(type)*tot);
    *b = malloc(sizeof(type)*(*n));
    for(int i = 0 ; i < tot ; ++i){
        fscanf(fp,"%f",*A+i);
    }
    for(int i = 0 ; i < *n ; ++i){
        fscanf(fp,"%f",*b+i);
    }
}

void showMtx(type *a,int n,int m){
    for(int i = 0 ; i < n ; ++i){
        for(int j = 0 ; j < m ; ++j){
            printf("%5.2f",*(a+i*m+j));
        }
        printf("\n");
    }
}
//// two



void ll_matvec(type *y,type *A,type *x,int m,int n){
    type *bA = A;
    type *eA = A+m*n;
    type *bx = x;
    type *ex = x+n;
    type *by = y;
    for(;bA!=eA;){
        if(bx==ex){
            bx = x;
            ++by;
        }
        *by+=(*bx)*(*bA);
        ++bx;
        ++bA;
    }
}

void matvec(type *y,type *A,type *x,int m,int n){
#pragma omp parallel for
    for(int i = 0 ; i < m ; ++i){
        dotprod(y+i,A+i*n,x,n);
    }
}
void residual(type *A,type*x,type *b,type *y,int n){
    for(int i = 0 ; i < n ; ++i)y[i] = 0.0;
    matvec(y,A,x,n,n);
    for(int i = 0 ; i < n ; ++i)y[i] = b[i] - y[i];
}

type getF(type *a,type *A,int n){////using temp storage
    type ret = 0;
    memset(temp,0, sizeof(type)*n);
    matvec(temp,A,a,n,n);
    dotprod(&ret,a,temp,n);
    return ret;
}

double dealRes(type *A,type *x,type *b,int cnt,int maxiter,int n){
    /*if(maxiter==-1){
        printf("true,with deal %d\n",cnt);
    }else if(cnt==-1){
        printf("error:");
    }*/
    double s = 0;
    type *y = malloc(sizeof(type)*n);
    residual(A,x,b,y,n);

    for(int i = 0 ; i < n ; ++i){
        s+=y[i]*y[i];
    }
    free(y);
    return sqrt(s);
}

void gauss(type *A,type *x,type *b,int n){
    int i,j;
    int *stai = (int*)malloc(sizeof(int)*n);
    int *staj = (int*)malloc(sizeof(int)*n);
    int top = 0;
    for(i = 0,j = 0 ; i < n && j < n ; ++i,++j){
        int flag = true;
        do {
            if (isZero(*(A + n * i + j))) {///zero
                for(int ti = i+1 ; ti < n ; ++ti){
                    if(!isZero(*(A+ti*n+j))){
                        swapIJ(A,i,ti,j,n);
                        swapIJ(b,i,ti,0,1);
                        flag = false;
                        break;
                    }
                }
            }else{
                flag = false;
            }
            if(flag){
                ++j;
            }
            if(j>=n)break;
        }while (flag);
        if(isZero(*(A + n * i + j)) && !isZero(b[i])){
            printf("Error-0-1\n");/// no solve
            free(stai);
            free(staj);
            return;
        }
        if(j==n)break;
        stai[top] = i;
        staj[top] = j;///store the loc
        ++top;
        type cur = *(A+n*i+j);
        for(int ti = i+1 ; ti < n ; ++ti){
            if(!isZero(*(A+n*ti+j))){
                type ks = -1.0**(A+n*ti+j)/cur;
                addKi_to_j(A,ks,i,ti,j,n);

                addKi_to_j(b,ks,i,ti,0,1);
            }
        }
    }
    int multi = false;

    if(i==j){/// has the only result
        --j,--i;
        while (i>=0 && j >=0){
            type k = 0;
            for(int ks = j+1 ; ks < n ; ++ks){
                k+=x[ks]*A[i*n+ks];
            }
            x[i] = (b[i]-k)/A[i*n+j];
            --i,--j;
        }
    }else {///multi result
        multi = true;
        while (top--){
            i = stai[top];
            --j;
            while (j!=staj[top]) {
                x[j] = 1;
                --j;
            }
            type k = 0;
            for(int ks = staj[top]+1 ; ks < n ; ++ks){
                k+=x[ks]*A[i*n+ks];
            }
            x[j] = (b[i] - k)/A[i*n+staj[top]];
        }
    }
    free(stai);
    free(staj);
}

void lu_crout(type *A,type *x,type *b,int n){
    type *l = (type *)malloc(sizeof(type)*n*n);
    type *u = (type *)malloc(sizeof(type)*n*n);
    for(int i = 0 ; i < n ; ++i){
        l[i*n] = A[i*n];
        u[i*n] = A[i]/(*l);
    }

    for(int r = 1; r < n ; ++r){
        for(int i = 1 ; i <= r; ++i){
            l[r*n+i] = A[r*n+i] - sumLrkUki(l,u,r, i,n);
            if(i==r) u[r*n+r]=1;
            else u[i+r*n] = (A[i*n+r] - sumLikUkr(l,u,r, i,n)) / l[i*n+i];
        }
    }
    type *y = (type*)malloc(sizeof(type)*n);
    calLow(l,y,b,n);
    calUp(u,x,y,n);
    free(y);
    free(l);
    free(u);
}

void lu_doolittle(type *A,type *x,type *b,int n){
    type *l = (type *)malloc(sizeof(type)*n*n);
    type *u = (type *)malloc(sizeof(type)*n*n);
    for(int i = 0 ; i < n ; ++i){
        u[i*n] = A[i];
        l[i*n] = A[i*n]/(*u);
    }

    for(int r = 1; r < n ; ++r){
        for(int i = r ; i < n; ++i){
            u[r+i*n] = A[r*n+i] - sumLrkUki(l,u,r, i,n);
            if(i==r) l[r*n+r]=1;
            else l[i*n+r] = (A[i*n+r] - sumLikUkr(l,u,r, i,n)) / u[r*n+r];
        }
    }
    type *y = (type*)malloc(sizeof(type)*n);
    calLow(l,y,b,n);
    calUp(u,x,y,n);
    free(y);
    free(l);
    free(u);
}

void cholesky(type *A,type *x,type *b,int n){
    type *l = (type *)malloc(sizeof(type)*n*n);
    for(int k = 0 ; k < n ; ++k){
        type sum = calSum(l,k,k,n);
        sum = A[k*n+k]-sum;
        if(sum>0){
            l[k*n+k] = sqrt(1.0*sum);
        }else l[k*n+k] = 0;
        //l[k*n+k] = sum>0?(sqrt(1.0*sum)):0;
        for(int i = k+1 ; i < n ; ++i){
            sum = calSum(l,i,k,n);
            l[i*n+k]=(A[i*n+k]-sum)/l[k*n+k];
        }
    }
    for(int i = 0 ; i < n ; ++i){
        for(int j = i+1 ; j < n ; ++j){
            l[j*n+i] = l[i*n+j];
        }
    }
    type *y = (type*)malloc(sizeof(type)*n);
    calLow(l,y,b,n);
    calUp(l,x,y,n);
    free(l);
    free(y);
}

void jacobi(type *A, type *x, type *b, int n, int *iter, int *maxiter, type threshold){
    int canExit = true;
    if(*iter>*(maxiter))
        return;
    *iter = *iter+1;
    type dot;
    type curS;
    for(int i = 0 ; i < n ; ++i){
        dot = 0.0;
        curS = A[i*n+i];
        dotprod(&dot,A+i*n,x,n);
        // dotprod(&dot,A+i*n,x,n);
        temp[i] = b[i] - dot+curS*x[i];
        //   type xx = x[i];
        if(!isZero(curS))x[i]/=curS;
        else temp[i] = x[i];
        if(fabs(temp[i] - x[i]) > threshold){
            canExit = false;
        }
    }

    memcpy(x,temp, sizeof(type)*n);

    if(canExit) {
        *maxiter = -1;//// marks for calculate a result
        return;
    }
    else{
        jacobi(A, x, b, n, iter, maxiter, threshold);
    }
}

void gs(type *A, type *x, type *b, int n, int *iter, int *maxiter, type threshold){
    int canExit = true;
    if(*iter>*(maxiter)) {
        *iter = -1;
        return;
    }
    *iter = *iter+1;
    type dot;
    type curS;
    for(int i = 0 ; i < n ; ++i){
        dot = 0.0;
        curS = A[i*n+i];
        dotprod(&dot,A+i*n,x,n);
        // dotprod(&dot,A+i*n,x,n);
        type c = b[i] - dot;
        //   type xx = x[i];
        if(!isZero(curS)) {
            x[i] = x[i] + 1.0*c/curS;
        }
        else x[i] = 1.0;
        //    printf("%.5f ",x[i]);
        if(fabs(temp[i] - x[i]) > threshold){
            canExit = false;
        }
    }
    //  printf("\n");
    if(canExit) {
        *maxiter = -1;//// marks for calculate a result
        return;
    }
    else{
        gs(A, x, b, n, iter, maxiter, threshold);
    }
}

void sor(type *A, type *x, type *b, int n, int *iter, int *maxiter, type threshold, type w){
    int canExit = true;
    if(*iter>*(maxiter)) {
        *iter = -1;
        return;
    }
    *iter = *iter+1;
    type dot;
    type curS;
    for(int i = 0 ; i < n ; ++i){
        dot = 0.0;
        curS = A[i*n+i];
        dotprod(&dot,A+i*n,x,n);
        // dotprod(&dot,A+i*n,x,n);
        type c = b[i] - dot;
        //   type xx = x[i];
        if(!isZero(curS)) {
            x[i] = x[i] + w*c/curS;
        }
        else x[i] = 1.0;
        //    printf("%.5f ",x[i]);
        if(fabs(temp[i] - x[i]) > threshold){
            canExit = false;
        }
    }
    //  printf("\n");
    if(canExit) {
        *maxiter = -1;//// marks for calculate a result
        return;
    }
    else{
        sor(A, x, b, n, iter, maxiter, threshold, w);
    }
}

void cg(type *A, type *x, type *b, int n, int *iter, int *maxiter, type threshold) {
    type *p = malloc(sizeof(type) * n);
    residual(A, x, b, p, n);
    type *r = malloc(sizeof(type) * n);
    memcpy(r, p, sizeof(type) * n);
    for (int i = 0; i < n; ++i)x[i] = 0.0;
    type rdot, rdot_1 = 1.0;
    for (*iter = 0; *iter < *maxiter; ++*iter) {
        dotprod(&rdot, r, r, n);
        if (sqrt(rdot) < RESPRE) {
            break;
        }
        type beta;
        if (!*iter);
        else {
            beta = rdot / rdot_1;
            for (int i = 0; i < n; ++i) {
                p[i] = r[i] + beta * p[i];
            }
        }
        type alpha = getF(p, A, n);
        alpha = rdot / alpha;

        //if(isnormal(alpha)){
        for (int i = 0; i < n; ++i) {
            x[i] += p[i] * alpha;
            r[i] -= temp[i] * alpha;
        }
        rdot_1 = rdot;
    }
}


type bisection(type (*func)(type),type left,type right,int* maxiter,type threshold){
    type ret = (right+left)/2;
    int cnt = 0;
    while (fabs(right-left)>threshold){
        if(cnt>*maxiter){
            *maxiter = -1;///no solve
            break;
        }
        ret = (right+left)/2;
        if(func(right)*func(ret)<=0){
            left = ret;
        }else{
            right = ret;
        }
        ++cnt;
    }
    *maxiter = cnt;
    return ret;
}

type fixedpoint(type (*func)(type),type init,int *maxiter,type threshold){
    type p =init;
    int cnt = 0;
    do{
        if(cnt>*maxiter){
            *maxiter = -1;
            break;
        }
        ++cnt;
        init = p;
        p = func(p)+p;
    }while (fabs(p-init)>threshold);
    *maxiter = cnt;
    return p;
}

type newtonraphson(type (*func)(type),type (*funcderivative)(type),type init,int *maxiter,type threshold){
    type x = init;
    int cnt = 0;
    do{
        if(cnt>*maxiter){
            *maxiter = -1;
            break;
        }
        ++cnt;
        init = x;
        x = x-func(x)/funcderivative(x);
    }while (fabs(x-init)>threshold);
    *maxiter = cnt;
    return x;
}

type secant(type (*func)(type),type initx0,type initx1,int *maxiter,type threshold){
    int cnt = 0;
    type k,y;
    do{
        if(cnt>*maxiter){
            *maxiter = -1;
            break;
        }
        ++cnt;
        y = func(initx0);
        k = (y-func(initx1))/(initx0-initx1);
        initx1 = initx0;
        initx0 = initx0 - y/k;
    }while (fabs(initx0-initx1)>threshold);
    *maxiter = cnt;
    return (initx0+initx1)/2;
}

#endif //NUMBER_BASEOPT_H
