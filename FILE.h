#ifndef FILE_H_
#define FILE_H_

#include<stdio.h>

//n:要素数
//a:n個の配列
void int_fscanf(int n,int *a, FILE* fp)
{
    int tmp;
    for(int i=0;i<n;i++){
        fscanf(fp,"%d",&tmp);
        a[i] = tmp;
    }
}

//n:要素数
//a:n個の配列
void double_fscanf(int n,double *a, FILE* fp)
{
    double tmp;
    for(int i=0;i<n;i++){
        fscanf(fp,"%lf",&tmp);
        a[i] = tmp;
    }
}

#endif //FILE_H_
