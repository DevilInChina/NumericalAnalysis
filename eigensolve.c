//
// Created by devilinchina on 10/5/19.
//
#include "baseOpt.h"
const int maxs=1024*1024;

int main(int argc,char *argv[]) {
    int n,m;
    type*x=loadMtx(argv[2],&n,&m);

    MALLOC(ans, type, n);
    int maxss = 5000;
    if(strcmp(argv[1],"qr")==0) {
        qrqeigensolver(ans, x, n, &maxss, 0.0001);
    }else if(strcmp(argv[1],"desc")==0){
        descend_power(ans ,x, n, &maxss, 0.0001);
    }else{
        printf("error\n");
        return 0;
    }

    showMtx(ans, 1, n);
//    showMtx(vec,1,n);
}
///