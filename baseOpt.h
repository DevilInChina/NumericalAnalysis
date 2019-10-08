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
#define CORE_NUMBER 4
#define MAX_ITER 1000
#define MALLOC(name,T,SIZE)\
typeof(T)* name = (typeof(T)*)malloc(sizeof(T)*(SIZE))
void init(){
    int c;
    typeof(c) x = 2;

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


void l2_dotprod(type *result,type *a,type *b ,int len){
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
void ll_dotprod(type *result,type *a,type *b ,int len){
    *result = 0.0;
    type *ba = a;
    type *bb = b;
    type *ea = a+len/CORE_NUMBER*CORE_NUMBER;
    type *tea = a+len;
    type res[CORE_NUMBER] = {0};
    for(;ba!=ea; ba+=CORE_NUMBER,bb+=CORE_NUMBER){
        res[0]+=ba[0]*bb[0];
        res[1]+=ba[1]*bb[1];
        res[2]+=ba[2]*bb[2];
        res[3]+=ba[3]*bb[3];
        /*res[4]+=ba[4]*bb[4];
        res[5]+=ba[5]*bb[5];
        res[6]+=ba[6]*bb[6];
        res[7]+=ba[7]*bb[7];*/
    }
    while (ba!=tea){
        *result+=*ba*(*bb);
        ++ba;
        ++bb;
    }
    for(int i = 0 ;i < CORE_NUMBER ; ++i)*result+=res[i];
}

void trans_to_T_matrix(type*A,int n,int m){
    MALLOC(tempA,type,m*n);

    for(int i = 0 ; i < n ; ++i){
        for(int j = 0 ; j < m ;++j){
            tempA[j*n+i] = A[i*m+j];
        }
    }
    memcpy(A,tempA, sizeof(type)*n*m);
    free(tempA);
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
    printf("[");
    for(int i = 0 ; i < n ; ++i){
        printf("[");
        for(int j = 0 ; j < m ; ++j){
            printf("%7.2f",*(a+i*m+j));
            if(j!=m-1)printf(",");
        }
        printf("]");
        if(i!=n-1)printf(",\n");
    }
    printf("]\n");
}
//// two

void axpby(type a,const type *x,type b,type*y,int n){
#pragma omp parallel for
    for(int i = 0 ; i < n ; ++i){
        y[i]=a*x[i]+b*y[i];
    }
}

void eable(type *A,int n){
    type res;
    dotprod(&res,A,A,n);
    res = sqrt(res);
    for(int i = 0 ; i < n ; ++i){
        A[i]/=res;
    }
}

void one_able(type *A,int n){
    type f=A[0];
    for(int i = 0 ; i < n ; ++i){
        A[i]/=f;
    }
}

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

void matmat(type *C,type *A,type *B,int m,int k,int n){
    trans_to_T_matrix(B,k,n);
    memset(C,0, sizeof(type)*m*n);
#pragma omp parallel for
    for(int i = 0 ; i < m ; ++i){
        for(int j = 0 ; j < n ; ++j){
            for(int kk = 0 ; kk < k ; ++kk){
                C[i*n+j]+=A[i*k+kk]*B[j*k+kk];
            }
        }
    }
    trans_to_T_matrix(B,n,k);
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

void polynomial_interopolation(int n,const type *x,type *y,type *a,type threshold){
    MALLOC(A,type,n*n);
    type *Ae = A+n*n;
    int ii = 0;
    for(type*col = A ; col!=Ae ;){
        type k = 1;
        for(int j = 0 ; j < n ; ++j,++col){
            *col = k;
            k*=x[ii];
        }
        ++ii;
    }
    lu_crout(A,a,y,n);
    free(A);
}

type lagrange_interopolation(int n,const type *x,type *y,type a,type threshold){
    type A = 1.0;
    for(int i = 0 ; i < n ; ++i) {
        A *= (a - x[i]);
        if(fabsf(a-x[i])<threshold){
            return y[i];
        }
    }
    type ret = 0.0;
    for(int i = 0 ; i < n ; ++i){
        type B = y[i];
        for(int j = 0 ; j < n ; ++j){
            if(i==j)continue;
            B/=(x[i]-x[j]);
        }
        ret+=B*A/(a-x[i]);
    }
    return ret;
}

type newtown_interopolation(int n,const type *x,type *y,type a,type threshold){
    MALLOC(dif,type,n*(n-1)/2);
    int las = n;
    type *beg = dif;
    type ret = 0;
    for(int i = 1 ; i < n ; ++i){///n-1
        beg[i-1] = (y[i] - y[i-1])/(x[i]-x[i-1]);
    }
    type B = a-x[0];
    ret += y[0]+beg[0]*B;
    int cnt = 1;
    for(int i = 2 ; i < n ; ++i){///n-2
        --las;
        beg+=las;
        B*=(a-x[i-1]);
        for(int j = 0 ; j < n - i ; ++j){
            beg[j] = (*(beg+j-las+1)-*(beg-las+j))/(x[j+i]-x[j]);
        }
        ret+=(*beg)*B;
    }
    free(dif);
    return ret;
}

void deal_origin(type *a,type x,int n,int mark){
    a[0] = mark;
    for(int i = 1 ; i < n ; ++i){
        a[i] = a[i-1]*x;
    }
}

void deal_dif1(type *a,type x, int n,int mark){
    type k = 1;
    int fac = 1;
    a[0] = 0;
    for(int i = 1; i < n ; ++i){
        a[i] = k*fac*mark;
        k*=x;
        ++fac;
    }
}

void deal_dif2(type *a,type x,int n,int mark){
    type k = 1;
    int fac = 2;
    a[0] = a[1] =0;
    for(int i = 2 ; i < n ; ++i){
        a[i] = k*fac*mark;
        k*=x;
        fac*=(i+1);
        fac/=(i-1);
    }
}

void cubic_spline_interpolation(int n,const type *x,type *y,type*a,type threshold){
    MALLOC(A,type,16*(n-1)*(n-1));
    MALLOC(Y,type,4*(n-1));
    memset(A,0, sizeof(type)*16*(n-1)*(n-1));
    memset(Y,0, sizeof(type)*4*(n-1));
    int len = 4*(n-1);
    type *curL = A;
    for(int i = 0 ; i < n-1 ; ++i){
        deal_origin(curL+i*4,x[i],4,1);
        Y[i] = y[i];
        curL+=len;
    }
    Y[n-1] = y[n-1];
    deal_origin(curL+(n-2)*4,x[n-1],4,1);
    curL+=len;
    for(int i = 1 ; i < n-1 ; ++i){/// 3*n-6 equ
        deal_origin(curL+i*4-4,x[i],4,1);
        deal_origin(curL+i*4,x[i],4,-1);
        curL+=len;
        deal_dif1(curL+i*4-4,x[i],4,1);
        deal_dif1(curL+i*4,x[i],4,-1);
        curL+=len;
        deal_dif2(curL+i*4-4,x[i],4,1);
        deal_dif2(curL+i*4,x[i],4,-1);
        curL+=len;
    }
    deal_dif1(curL,x[0],4,1);///1
    curL+=len;
    deal_dif2(curL+4*(n-2),x[n-1],4,1);////1

   // showMtx(A,len,len);
  //  showMtx(Y,len,1);
    /////1+1+n+3*n-6 = 4*n-4

    int cnt=0,maxcnt = 1000;
    gauss(A,a,Y,len);
    //showMtx(A,len,len);
}

void power(type *eigenvalues,type*eigen_vec,type *A,int n,int *maxiter,type threshold) {
    MALLOC(eigenvector, type, n);
    MALLOC(eigenvector_new, type, n);
    for (int i = 0; i < n; ++i)eigenvector[i] = 1;
    float evalue = eigenvector[0];
    float evalue_new;
    int iter = 0;
    float error = 0;
    long long sz = sizeof(type) * n;
    do {
        matvec(eigenvector_new,A, eigenvector,n,n);
        evalue_new = eigenvector_new[0];
        for (int i = 0; i < n; ++i)eigenvector_new[i] /= evalue_new;
        memcpy(eigenvector, eigenvector_new, sz);
        error = fabsf((evalue_new - evalue) / evalue_new);
        evalue = evalue_new;
        ++iter;
    } while (iter < *maxiter && error > threshold);
    *eigenvalues = evalue;
    memcpy(eigen_vec,eigenvector_new, sizeof(type)*n);
    one_able(eigen_vec,n);
    free(eigenvector);
    free(eigenvector_new);
}

void descend_power(type *ans,type *A,int n,int *maxiter,type threshold){
    MALLOC(As,type,n*n);
    memcpy(As,A, sizeof(type)*n*n);
    MALLOC(vecA,type,n);
    MALLOC(ev,type,n);
    MALLOC(va,type,n*n);
    for(int i = 0 ; i < n ; ++i){
        int maxs = *maxiter;
        power(ans+i,ev,As,n-i,&maxs,threshold);
        memcpy(vecA,As, sizeof(type)*(n-i));
        int tot = n-i-1;
        matmat(va,ev,vecA,tot+1,1,tot+1);
        axpby(-1.0,va,1.0,As,(n-i)*(n-i));
        for(int j = 0 ; j < tot ; ++j){
            for(int k = 0 ; k < tot ; ++k){
                As[j*(tot)+k] = As[(j+1)*(tot+1)+k+1];
            }
        }
    //    showMtx(As,n,n);
    }
    free(vecA);
    free(As);
    free(ev);
    free(va);
}

#define show(x,ch)\
printf("%.5lf%c",x,ch);fflush(stdout)
void schmidt_orthogonalization(type *A,int n,int m){
    int flag = n<m;
    trans_to_T_matrix(A,n,m);
    MALLOC(beta_vec,(A),m);
    MALLOC(alph_vec,(A),m);
    MALLOC(dot_beta,type,m);
    MALLOC(dot_alpha,type,m);
    for(int i = 0 ; i < m ; ++i){
        beta_vec[i] = A+i*n;
        alph_vec[i] = A+i*n;
    }
    for(int i = 1 ; i < m ; ++i){////m
        dotprod(dot_beta+i-1,beta_vec[i-1],beta_vec[i-1],n);////cal dot beta[i-1]
#pragma omp parallel for
        for(int j = 0 ; j < i ; ++j){/////m*m
            dotprod(dot_alpha+j,alph_vec[i],beta_vec[j],n);/////m*m*n
            //// cal dot a[i]*b[j] to alp[j]
        }
        for(int j = 0 ; j < i ; ++j){////non-parallel
            axpby(-dot_alpha[j] / dot_beta[j],beta_vec[j],1,alph_vec[i],n);

        }
     //  exit(0);
    }
    for(int i = 0 ; i < m ;++i){
        eable(alph_vec[i],n);
    }
    free(dot_beta);
    free(dot_alpha);
    free(beta_vec);
    free(alph_vec);
    trans_to_T_matrix(A,m,n);
}


void qr(type*Q,type*R,type *A,int n,int m){
    memcpy(Q,A, sizeof(type)*n*m);
    schmidt_orthogonalization(Q,n,m);
    trans_to_T_matrix(Q,n,m);
    matmat(R,Q,A,m,n,m);
   // showMtx(R,m,m);
    trans_to_T_matrix(Q,m,n);
}

#define Sp(x) ((x)*(x))

void qrqeigensolver(type*ans,type*A,int n,int *maxiter,type threshold){
    MALLOC(newA1,type,n*n);
    MALLOC(newA2,type,n*n);
    MALLOC(Q,type,n*n);
    MALLOC(R,type,n*n);
    memcpy(newA1,A, sizeof(type)*n*n);
    type error = 0;

    int iterator = 0;
    do{
        error = 0;
        if(iterator>*maxiter){
            *maxiter = -1;
            break;
        }
        qr(Q,R,newA1,n,n);
        matmat(newA2,R,Q,n,n,n);
        for(int i = 0 ; i < n ; ++i){
            error+=Sp(newA2[i*n+i]-newA1[i*n+i]);
        }
        memcpy(newA1,newA2, sizeof(type)*n*n);
        ++iterator;
        error = sqrt(error);
        printf("%lf\n",error);
    }while (error>threshold);

    if(*maxiter!=-1){
        *maxiter = iterator;
    }
    for(int i = 0 ; i < n ; ++i){
        ans[i] = newA1[i*n+i];
    }

    free(Q);
    free(R);
    free(newA1);
    free(newA2);

}

#endif //NUMBER_BASEOPT_H
