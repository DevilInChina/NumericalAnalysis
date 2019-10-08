//
// Created by devilinchina on 9/29/19.
//

#include "baseOpt.h"
int main(int argc,char *argv[]){
    int n;
    freopen(argv[2],"r",stdin);
    freopen(argv[3],"w",stdout);
    scanf("%d",&n);
    MALLOC(a,type,n);
    MALLOC(x,type,n);
    MALLOC(y,type,n);
    for(int i = 0 ; i < n ; ++i){
        scanf("%f %f",x+i,y+i);
    }
    int showed = 0;
    if(strcmp(argv[1],"poly")==0) {
        polynomial_interopolation(n, x, y, a, THRESHOLD);
    }else if(strcmp(argv[1],"lag")==0){
        for(type xx = -1.0 ; xx<1.0 ; xx+=0.01){
            printf("%.5f\n",lagrange_interopolation(n,x,y,xx,THRESHOLD));
        }
        showed = 1;
    }else if(strcmp(argv[1],"newtown")==0){
        for(type xx = -1.0 ; xx<1.0 ; xx+=0.01){
            printf("%.5f\n",newtown_interopolation(n,x,y,xx,THRESHOLD));
        }
        showed = 1;
    }else if(strcmp(argv[1],"cubic")==0){
        free(a);
        a = (type*)malloc(sizeof(type)*4*(n-1));
        cubic_spline_interpolation(n,x,y,a,THRESHOLD);
        for(int i = 0 ; i < n-1 ; ++i){
            for(int j = 0 ; j < 4 ; ++j){
                if(j!=3)
                printf("%.5f ",a[i*4+j]);
                else printf("%.5f",a[i*4+j]);
            }
            printf("\n");
        }
        showed = 1;
    }
    if(!showed){
        for (int i = 0; i < n; ++i) {
            printf("%.5f\n", a[i]);
        }
    }
    free(a);
    free(x);
    free(y);
}