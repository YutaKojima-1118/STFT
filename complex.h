#ifndef COMPLEX_H_
#define COMPLEX_H_

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

typedef struct complex{
    double re;
    double im;
}complex;

complex e_powered(double r)
{
    complex cmp;
    cmp.re = cos(r);
    cmp.im = sin(r);
    return cmp;
}

complex sum(complex c1, complex c2)
{
    complex answer;
    answer.re = c1.re + c2.re;
    answer.im = c1.im + c2.im;
    return answer;
}

complex c_c_multi(complex c1, complex c2)
{
    complex answer;
    answer.re = c1.re*c2.re - c1.im*c2.im;
    answer.im = c1.im*c2.re + c1.re*c2.im;
    return answer;
}

complex r_c_multi(double r, complex c)
{
    c.im *= (double)r;
    c.re *= (double)r;
    return c;
}

complex divide(complex c, int n)
{
    c.im /= (double)n;
    c.re /= (double)n;
    return c;
}

//絶対値
double power(complex c){
    return sqrt(c.im*c.im + c.re*c.re);
}

void print_complex(complex c)
{
    if(c.re >= 0.05 || c.re < -0.05 /*c.re!=0*/){
        if(c.im >= -0.05 && c.im < 0.05){
            printf("%.1f",c.re);
        }
        else if(c.im > 0){
            printf("%.1f + j%.1f",c.re,c.im);
        }
        else{
            printf("%.1f - j%.1f",c.re,-c.im);
        }
    }
    else{
        if(c.im >= -0.05 && c.im < 0.05){
            printf("0");
        }
        else if(c.im > 0){
            printf("j%.1f",c.im);
        }
        else{
            printf("- j%.1f",-c.im);
        }
    }
}

void init(int n, complex *c)
{
    for(int i=0;i<n;i++){
        c[i].re=0;
        c[i].im=0;
    }
}

complex*** gen_triple_pointer_complex(int a,int b,int c)
{
    complex ***tmp = (complex ***)malloc(sizeof(complex**)*a);
    for(int i=0;i<a;i++){
        tmp[i] = (complex**)malloc(sizeof(complex*)*b);
        for(int j=0;j<b;j++){
            tmp[i][j] = (complex*)malloc(sizeof(complex)*c);
            init(c,tmp[i][j]);
        }
    }
    return tmp;
}

complex** gen_double_pointer_complex(int col, int row)
{
    complex **tmp = (complex **)malloc(sizeof(complex*)*col);
    for(int i=0;i<col;i++){
        tmp[i] = (complex*)malloc(sizeof(complex)*row);
        init(row,tmp[i]);
    }
    return tmp;
}

complex* gen_pointer_complex(int n)
{
    complex *tmp = (complex *)malloc(sizeof(complex)*n);
    init(n,tmp);
    return tmp;
}

//n:要素数
//a:n個の配列
void real_fscanf(int n,complex *c, FILE* fp)
{
    double tmp;
    for(int i=0;i<n;i++){
        fscanf(fp,"%lf",&tmp);
        c[i].re = tmp;
    }
}

#endif //COMPLEX_H_
