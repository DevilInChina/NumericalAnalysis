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
#define CORE_NUMBER 4
#define MAX_ITER 1000
double mem = 0;
int CNT = 0;

#define MALLOC(name,T,SIZE)\
typeof(T)* name = (typeof(T)*)malloc(sizeof(T)*(SIZE))
#define MALLCP(name,from,T,SIZE)\
MALLOC(name,T,SIZE);\
memcpy(name,from,sizeof(T)*(SIZE))
#define CALLOC(name,T,SIZE)\
typeof(T)* name = (typeof(T)*)calloc(sizeof(T),SIZE)

void init(){
    int c;
    typeof(c) x = 2;

   // CORE_NUMBER = CORE_NUMBERS;
}

//// one
void swap(type *a,type*b){
    (a!=b?(*(tempType*)a^=*(tempType*)b^=*(tempType*)a^=*(tempType*)b):*a);
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
#define trans_Next(i,m,n)\
((i)%(n)*(m))+((i)/(n))

#define trans_Prev(i,m,n)\
((i)%(m)*(n))+((i)/(m))
void transMoveData(type *mtx,int i,int m,int n){
    type s = mtx[i];
    int cur = i;
    int pre = trans_Prev(i,m,n);
    while (pre!=i){
        mtx[cur]=mtx[pre];
        cur = pre;
        pre = trans_Prev(cur,m,n);
    }
    mtx[cur] = s;
}
void trans_to_T_matrix(type*A,int m,int n){
  //  MALLOC(tempA,type,m*n);
    int k = m*n;
    for(int i = 0 ; i < k ; ++i){
        int next = trans_Next(i,m,n);
        while (next>i){
            next = trans_Next(next,m,n);
        }
        if(next==i){
            transMoveData(A,i,m,n);
        }
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

void setXtoA(type *X,type A,int len){
#pragma omp parallel for
    for(int i = 0 ; i < len ; ++i){
        X[i] = A;
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
            printf("%7.9f",*(a+i*m+j));
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
    if(!isZero(res))
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

type getF(type *a,type *A,int n,type *temp){////using temp storage
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

void matmat_transB(type *C,type *A,type *B,int m,int k,int n){
    memset(C,0, sizeof(type)*m*n);
#pragma omp parallel for
    for(int i = 0 ; i < m ; ++i){
        for(int j = 0 ; j < n ; ++j){
            for(int kk = 0 ; kk < k ; ++kk){
                C[i*n+j]+=A[i*k+kk]*B[j*k+kk];
            }
        }
    }
}
void matmat_transA(type *C,type *A,type *B,int m,int k,int n){
    trans_to_T_matrix(A,m,k);
    matmat(C,A,B,k,m,n);
    trans_to_T_matrix(A,k,m);
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
    MALLOC(temp,type,n);
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
    free(temp);
    memcpy(x,temp, sizeof(type)*n);

    if(canExit) {
        *maxiter = -1;//// marks for calculate a result
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
        type xx = x[i];
        if(!isZero(curS)) {
            x[i] = x[i] + 1.0*c/curS;
        }
        else x[i] = 1.0;
        //    printf("%.5f ",x[i]);
        if(fabsf(xx - x[i]) > threshold){
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
        type xx = x[i];
        type c = b[i] - dot;
        if(!isZero(curS)) {
            x[i] = x[i] + w*c/curS;
        }
        else x[i] = 1.0;

        if(fabsf(xx - x[i]) > threshold){
            canExit = false;
        }
    }
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
    MALLOC(temp,type,n);
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
        type alpha = getF(p, A, n,temp);
        alpha = rdot / alpha;

        for (int i = 0; i < n; ++i) {
            x[i] += p[i] * alpha;
            r[i] -= temp[i] * alpha;
        }
        rdot_1 = rdot;
    }
    free(temp);
    free(r);
    free(p);
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
    ///get The biggest eigenvalue and it's vector
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

void geteigenVec(type *eigenValues,type *eigenVec,type *A,int n){
    MALLOC(teM,type,n*n);
    MALLOC(Ans,type,n*n);
    memset(Ans,0, sizeof(type)*n*n);
    for(int i = 0 ; i < n ; ++i){
        Ans[i*n+i] = eigenValues[i];
    }
    printf("vectors:\n");
    showMtx(eigenVec,n,n);
    printf("value:\n");
    matmat(teM,Ans,eigenVec,n,n,n);
    showMtx(teM,n,n);
    printf("vector\n");
    /*
    for(int i = 0 ; i < n ; ++i){
        matvec(teM+i*n,A,eigenVec+i*n,n,n);
    }
     */
    matmat_transB(teM,eigenVec,A,n,n,n);
   // trans_to_T_matrix(teM,n,n);
    showMtx(teM,n,n);
}

void descend_power(type *ans,type *eigenVec,type *A,int n,int *maxiter,type threshold){
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
    }
    free(vecA);
    free(As);
    free(ev);
    free(va);
}

#define show(x,ch)\
printf("%.5lf%c",x,ch);fflush(stdout)
void schmidt_orthogonalization(type *A,int m,int n){
  ///  int flag = n<m;
    trans_to_T_matrix(A,m,n);
    MALLOC(beta_vec,(A),n);
    MALLOC(alph_vec,(A),n);
    MALLOC(dot_beta,type,n);
    MALLOC(dot_alpha,type,n);
    for(int i = 0 ; i < n ; ++i){
        beta_vec[i] = A+i*m;
        alph_vec[i] = A+i*m;
    }
    type temp;
    for(int i = 1 ; i < n ; ++i){////m
        dotprod(dot_beta+i-1,beta_vec[i-1],beta_vec[i-1],m);////cal dot beta[i-1]
#pragma omp parallel for
        for(int j = 0 ; j < i ; ++j){/////m*m
            dotprod(dot_alpha+j,alph_vec[i],beta_vec[j],m);/////m*m*n
            //// cal dot a[i]*b[j] to alp[j]
        }
        for(int j = 0 ; j < i ; ++j) {////non-parallel
            temp = dot_beta[j];
            if (isZero(temp))
                axpby(0, beta_vec[j], 1, alph_vec[i], m);
            else {
                if(!isZero(dot_beta[j]))
                    axpby(-dot_alpha[j] / dot_beta[j], beta_vec[j], 1, alph_vec[i], m);
            }
        }
     //  exit(0);
    }
    for(int i = 0 ; i < n ;++i){
        eable(alph_vec[i],m);
    }
    free(dot_beta);
    free(dot_alpha);
    free(beta_vec);
    free(alph_vec);
    trans_to_T_matrix(A,n,m);
}


void qr(type*Q,type*R,type *A,int m,int n){
    /// can means calculator vector
    memcpy(Q,A, sizeof(type)*n*m);
    schmidt_orthogonalization(Q,m,n);
  //  printf("ins:");
   // showMtx(Q,m,n);
    trans_to_T_matrix(Q,m,n);
    matmat(R,Q,A,n,m,n);
 //   trans_to_T_matrix(Q,m,n);
}

#define Sp(x) ((x)*(x))
void toOne(type *A,int m,int n){
    for(int i = 0 ; i < m ; ++i){
        type sum;
        dotprod(&sum,A+i*n,A+i*n,n);
        if(A[i*n]<0)sum=-sum;
        for(int j = 0 ; j < n ; ++j){
            A[i*n+j]/=sum;
        }
    }
}
void toSne(type *A,int m,int n){
    for(int i = 0 ; i < m ; ++i){
        type k = A[i*n];
        for(int j = 0 ; j < n ; ++j){
            A[i*n+j]/=k;
        }
    }
}
void qrqeigensolver(type*ans,type *eigenVector,type*A,int n,int *maxiter,type threshold){
    MALLOC(newA1,type,n*n);
    MALLOC(newA2,type,n*n);
    MALLOC(Q,type,n*n);
    MALLOC(R,type,n*n);
    memcpy(newA1,A, sizeof(type)*n*n);
    MALLOC(teM,type,n*n);
    CALLOC(Qp,type,n*n);

    for(int i = 0 ; i < n ; ++i){
        Qp[i*n+i] = 1;
    }
    type error = 0;

    int iterator = 0;
    do{
        error = 0;
        if(iterator>*maxiter){
            *maxiter = -1;
            break;
        }

        qr(Q,R,newA1,n,n);///get Qt and R
        matmat(eigenVector,Q,Qp,n,n,n);////eigenVec = Q(n)t * Q(n-1)t * ... * Q(1)t
        memcpy(Qp,eigenVector, sizeof(type)*n*n);
        matmat_transB(newA2,R,Q,n,n,n);

        for(int i = 0 ; i < n ; ++i){
            error+=Sp(newA2[i*n+i]-newA1[i*n+i]);
        }
        memcpy(newA1,newA2, sizeof(type)*n*n);
        ++iterator;
        error = sqrt(error);
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
    free(teM);
    free(Qp);
}

void mid_integral(type *result,type (*func)(type),type beg,type end,int scale){
    MALLOC(scales,type,scale);
    type dif = (end-beg)/scale;
    setXtoA(scales,dif,scale);
    MALLOC(fx,type,scale);
#pragma omp parallel for
    for(int i = 0 ; i < scale ; ++i){
        fx[i] = func(beg+(dif*(2*i+1)/2));
    }

    dotprod(result,fx,scales,scale);
    free(fx);
    free(scales);
}

void trapezoid_integral(type *result,type (*func)(type),type beg,type end,int scale){
    type ans = 0;
    type dif = (end-beg)/scale;
    for(int i = 1; i < scale ; ++i){
        ans+=func(dif*i+beg);
    }
    ans+=(func(beg)+func(end))/2;
    *result = ans*dif;
}

void simpson13_integral(type *result,type (*func)(type),type beg,type end,int scale){
    type dif = (end-beg)/scale;
    type ans = 0;
    for(int i = 1 ; i < scale ; ++i){
        switch (i&1){
            case 1:
                ans+=4*func(beg+i*dif);
                break;
            default:
                ans+=2*func(beg+i*dif);
                break;
        }
    }
    ans+=func(end)+func(beg);
    ans*=dif/3;
    *result = ans;
}


void simpson38_integral(type *result,type (*func)(type),type beg,type end,int scale){
    type dif = (end-beg)/scale;
    type ans = 0;
    for(int i = 1 ; i < scale ; ++i){
        switch (i%3){
            case 0:
                ans+=2*func(beg+i*dif);
                break;
            default:
                ans+=3*func(beg+i*dif);
                break;
        }
    }
    ans+=func(end)+func(beg);
    ans*=3.0/8*dif;
    *result = ans;
}

void simpson_integral(type *result,type (*func)(type),type beg,type end,int scale) {
    type a = 0, b = 0;
    if (scale > 1) {
        if (scale % 3 != 0) {
            int res = scale % 3;
            type mid ;
            if (res & 1) {
                res += 3;
            }
            mid = end - (end - beg) * (scale - res) / scale;

            simpson13_integral(&a, func, beg, mid, res);
            if (scale - res)
                simpson38_integral(&b, func, mid, end, scale - res);
            *result = a + b;
        } else {
            simpson38_integral(result, func, beg, end, scale);
        }
    } else simpson13_integral(result, func, beg, end, scale);
}

void romberg_integral(type *result,type (*func)(type),type beg,type end,int scale){
    MALLOC(results,type,scale*scale);
    setXtoA(results,0,scale*scale);
    type dif = (end-beg)/2;
    int nesg = 2;
    for(int i = 0 ; i < scale ; ++i){
        type Temp = 0;
        for (int j = 1 ; j < nesg ; ++j){
            Temp+=2*func(beg+j*dif);
        }
        Temp+=func(beg)+func(beg+nesg*dif);
        results[i*scale+0] = Temp*dif/2;
        //Temp/=2;
        nesg*=2;
        dif/=2;
    }
    for(int j = 1; j < scale ; ++j){
        for(int i = 0 ; i < scale - j ; ++i){
            results[i*scale+j] = (4.0*results[(i+1)*scale+j-1]-1.0*results[i*scale+j-1])/3.0;
        }
    }
    *result = results[scale-1];
    free(results);
}

typedef struct csr{
    int *ia;
    int *ja;
    type *va;
}csr;
void csrmm(const int *ia,const int *ja, const type *va,
        const int *ib, const int *jb, const type *vb,
        csr*ans,
        int m, int k, int n) {
    MALLOC(anscol,type,n);
    ans->ia = (int*)malloc((m+1)* sizeof(int));
    MALLOC(ctemp,ans->ja,m);
    MALLOC(cval,ans->va,m);
    int nnz = 0;
    ans->ia[0] = 0;
    int *ic = ans->ia;

    for(int i = 0 ; i < m ; ++i){
        setXtoA(anscol,0,n);
        int mnnz = ia[i+1];
        for(int j = ia[i] ; j < mnnz ; ++j){
            int ks = ja[j];
            for(int pos = ib[ks] ; pos != ib[ks+1] ; ++pos){
                anscol[jb[pos]]+=va[j]*vb[pos];
            }
        }
        for(int j = 0;  j < n ; ++j){
            if(!isZero(anscol[j])){
                ++nnz;
            }
        }
        ic[i+1]=nnz;
        ctemp[i] = malloc((ic[i+1]-ic[i])* sizeof(int));
        cval[i] = malloc((ic[i+1]-ic[i])* sizeof(type));
        int cnt = 0;
        for(int j = 0;  j < n ; ++j){
            if(!isZero(anscol[j])){
                ctemp[i][cnt] = j;
                cval[i][cnt] = anscol[j];
                ++cnt;
            }
        }
    }
    int cnt =0;
    ans->ja = (int*)malloc(sizeof(int)*nnz);
    int *jc = ans->ja;
    ans->va = (type*)malloc(sizeof(int)*nnz);

    type*vc=ans->va;
    for(int i = 0 ; i < m ; ++i){
        memcpy(jc+ic[i],ctemp[i], (ic[i+1]-ic[i])*sizeof(int));
        memcpy(vc+ic[i],cval[i],(ic[i+1]-ic[i])* sizeof(type));
        free(ctemp[i]);
        free(cval[i]);
    }

    free(ctemp);
    free(cval);
}


void tridiagonalization(type *A,type *T,type *P,int n){
    MALLOC(p1,type,n);
    MALLOC(p2,type,n);
    MALLOC(w,type,n);
    MALLOC(wp,type,n);
    memset(p1,0,sizeof(n));
    p1[0]=1.0;
    matvec(wp,A,p1,n,n);
    type alpha;
    dotprod(&alpha,wp,p1,n);
    T[0]=alpha;
    for(int i=0;i<n;i++){
        w[i]=wp[i]-alpha*p1[i];
    }
    for(int i=0;i<n;i++){
        P[i*n]=p1[i];
    }
    for(int j=1;j<n;j++){
        type beta;
        dotprod(&beta,w,w,n);
        beta = sqrt(beta);
        T[(j-1)*n+j]=beta;
        T[j*n+j-1]=beta;
        for(int i=0;i<n;i++){
            if(!isZero(beta))
            p2[i]=w[i]/beta;
            else p2[i] = w[i];
        }
        matvec(wp,A,p2,n,n);

        dotprod(&alpha,wp,p2,n);
        T[j*n+j]=alpha;
        for(int i=0;i<n;i++){
            w[i]=wp[i]-alpha*p2[i]-beta*p1[i];
        }
        for(int i=0;i<n;i++){
            P[i*n+j]=p2[i];
        }
        memcpy(p1,p2,sizeof(type)*n);
    }
    free(wp);
    free(w);
    free(p2);
    free(p1);
}

typedef struct Image{
    type *F;
    type *T;
    int m,n;
}Image;
void toDiag(type*A,int n){
    for(int i = 0 ; i < n ; ++i){
        swap(A+i,A+i*n+i);
    }
}
void matmatMatTrans(type *P,type*T,int n,int m){
    MALLOC(Temp,type,n*m);
    matmat(Temp,P,T,n,m,m);
    MALLOC(Ret,type,n*n);
    matmat_transB(Ret,Temp,P,n,m,n);
    showMtx(Ret,n,n);
}
void lanczos(type *image,type*U,type*sigma,int n) {

    MALLOC(P,type,n*n);
    MALLOC(T,type,n*n);

    tridiagonalization(image, T, P, n);

    int maxs = 20;
   // showMtx(T,n,n);
    qrqeigensolver(sigma, U, T, n, &maxs, 0.000000001);

    matmat_transB(T,P,U,n,n,n);

    memcpy(U,T, sizeof(type)*n*n);

    free(T);
    free(P);
}
void drill(type*T,const type *A,int n,int m){///n>m
    for(int i = 0 ; i < n ; ++i){
        for(int j = 0 ; j < m ; ++j){
            T[i*m+j] = A[i*n+j];
        }
    }
}
void matMdiag(type*A,type *diag,int n,int m){////n*m m
    for(int i = 0 ; i < n ; ++i){
        for(int j = 0 ; j < m ; ++j){
            A[i*m+j]*=diag[j];
        }
    }
}

void lanczosOne(type*image,type*compress,type compression_rate,int n){
    CALLOC(sigma,type,n*n);

    lanczos(image,compress,sigma,n);

    type sum = 0;
    for(int i = 0 ; i < n ; ++i){
        sum+=fabsf(sigma[i]);
    }
    int m = 1;
    type cur = fabsf(*sigma);
    for(; m < n ; ++m){
        if(cur/sum>compression_rate)break;
        cur+=fabsf(sigma[m]);
    }

    MALLOC(Temp,type,n*m);
    drill(Temp,compress,n,m);///n m
    matMdiag(Temp,sigma,n,m);///
    MALLOC(TTS,type,n*n);
    matmat_transB(TTS,Temp,compress,n,m,n);
    memcpy(compress,TTS, sizeof(type)*n*n);
    free(TTS);
    free(Temp);
    free(sigma);
}

void svd(type *img,int n){
    MALLOC(U,type,n*n);
    MALLOC(V,type,n*n);
    matmat_transB(V,img,img,n,n,n);
    trans_to_T_matrix(img,n,n);
    matmat_transB(U,img,img,n,n,n);
    trans_to_T_matrix(img,n,n);
    MALLOC(ans,type,n*n);
    MALLOC(Ue,type,n*n);
    MALLOC(Ve,type,n*n);
    int maxs = 20;
    qrqeigensolver(ans,Ue,U,n,&maxs,THRESHOLD);


    qrqeigensolver(ans,Ve,V,n,&maxs,THRESHOLD);

    for(int i = 0 ; i < n ; ++i){
        ans[i] = sqrt(ans[i]);
    }
    toDiag(ans,n);
    MALLOC(kk,type,n*n);

   // trans_to_T_matrix(Ue,n,n);
    matmat(kk,Ue,ans,n,n,n);
    matmat(ans,kk,Ve,n,n,n);
    showMtx(img,n,n);
    showMtx(ans,n,n);
}
/*
 *
A= U SIGME V
 */

#endif //NUMBER_BASEOPT_H
