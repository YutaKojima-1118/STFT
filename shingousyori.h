#ifndef SHINGOUSYORI_H_
#define SHINGOUSYORI_H_

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"complex.h"
#define MY_PI (atan(1.0)*4)

//DCT

int DCT(int n, double *input, double *output)
{
    for(int k=0;k<n;k++){
        output[k] = 0;
        for(int i=0;i<n;i++){
            output[k] += input[i] * cos(((2.0*i+1)*k*MY_PI)/2.0/n);
        }
        //printf("%.1f ",output[k]);
    }
    //printf("\n");
}

int iDCT(int n, double *input, double *output)
{
    double sum=0;
    for(int k=0;k<n;k++){
        sum += 1.0 / 2.0 * input[0];
        for(int i=1;i<n;i++){
            double tmp;
            tmp = input[i]*cos(((2.0*(double)k+1.0)*(double)i*MY_PI)/(2*n));
            //printf("%f ",tmp);
            sum += tmp;


        }
        output[k] = sum * 2/n;
        //printf("\n");
        sum = 0;
    }
/*
    for(int i=0;i<n;i++){
        printf("%.1f ",output[i]);
    }
    printf("\n");
*/
}

void DCT_2d(int col, int row, double **input, double **output)
{
    double **tmp1 = (double**)malloc(sizeof(double*)*col); //tmp1[COL][ROW]
    for(int i=0;i<col;i++){
        tmp1[i] = (double*)malloc(sizeof(double)*row);
        for(int j=0;j<row;j++){
            tmp1[i][j]=0;
        }
    }

    double **tmp2 = (double**)malloc(sizeof(double*)*row); //tmp2[ROW][COL]
    double **tmp3 = (double**)malloc(sizeof(double*)*row);
    for(int i=0;i<row;i++){
        tmp2[i] = (double*)malloc(sizeof(double)*col);
        tmp3[i] = (double*)malloc(sizeof(double)*col);
        for(int j=0;j<row;j++){
            tmp2[i][j]=0;
            tmp3[i][j]=0;
        }
    }

    for(int i=0;i<col;i++){ //横向きにDCT
        DCT(row,input[i],tmp1[i]);
    }

    for(int i=0;i<col;i++){
        for(int j=0;j<row;j++){
            tmp2[j][i] = tmp1[i][j];
        }
    }

    for(int i=0;i<row;i++){ //縦向きDCT
        DCT(col,tmp2[i],tmp3[i]);
    }

    for(int i=0;i<col;i++){
        for(int j=0;j<row;j++){
            output[i][j] = tmp3[j][i];
        }
    }

    for(int i=0;i<row;i++){
        free(tmp1[i]);
        free(tmp2[i]);
        free(tmp3[i]);
    }
    free(tmp1);
    free(tmp2);
    free(tmp3);
}

void iDCT_2d(int col, int row, double **input, double **output)
{
    double **tmp1 = (double**)malloc(sizeof(double*)*col); //tmp1[COL][ROW]
    for(int i=0;i<col;i++){
        tmp1[i] = (double*)malloc(sizeof(double)*row);
        for(int j=0;j<row;j++){
            tmp1[i][j]=0;
        }
    }

    double **tmp2 = (double**)malloc(sizeof(double*)*row); //tmp2[ROW][COL]
    double **tmp3 = (double**)malloc(sizeof(double*)*row);
    for(int i=0;i<row;i++){
        tmp2[i] = (double*)malloc(sizeof(double)*col);
        tmp3[i] = (double*)malloc(sizeof(double)*col);
        for(int j=0;j<row;j++){
            tmp2[i][j]=0;
            tmp3[i][j]=0;
        }
    }

    for(int i=0;i<col;i++){ //横向きにDCT
        iDCT(row,input[i],tmp1[i]);
    }

    for(int i=0;i<col;i++){
        for(int j=0;j<row;j++){
            tmp2[j][i] = tmp1[i][j];
        }
    }

    for(int i=0;i<row;i++){ //縦向きDCT
        iDCT(col,tmp2[i],tmp3[i]);
    }

    for(int i=0;i<col;i++){
        for(int j=0;j<row;j++){
            output[i][j] = tmp3[j][i];
        }
    }

    for(int i=0;i<row;i++){
        free(tmp1[i]);
        free(tmp2[i]);
        free(tmp3[i]);
    }
    free(tmp1);
    free(tmp2);
    free(tmp3);
}

//DFT

void DFT(int n,  complex *input, complex *output)
{
    complex s,tmp;
    s.re=0; s.im=0; tmp.re=0; tmp.im=0;
    for(int k=0;k<n;k++){
        for(int i=0;i<n;i++){
            tmp = c_c_multi(input[i],e_powered((-2.0*MY_PI*(double)i*(double)k)/(double)n));
            s = sum(tmp,s);
        }
        //printf(" F[%d] = ",k);
        output[k] = s;
        //print_complex(output[k]);
        //printf("\n");
        s.re=0; s.im=0;
    }
}

void iDFT(int n,  complex *input, complex *output)
{
    complex s,tmp;
    s.re=0; s.im=0; tmp.re=0; tmp.im=0;
    for(int k=0;k<n;k++){
        for(int i=0;i<n;i++){
            tmp = c_c_multi(e_powered((2.0*MY_PI*(double)i*(double)k)/(double)n),input[i]);
            s = sum(tmp,s);
        }
        s = divide(s,n);
        output[k] = s;
        //printf(" f[%d] = %.1f\n",k,output[k].re);
        s.re=0; s.im=0;
    }
}

void DFT_2d(int col, int row, complex **input, complex **output)
{
    complex **tmp1 = (complex**)malloc(sizeof(complex*)*col); //tmp1[COL][ROW]
    for(int i=0;i<col;i++){
        tmp1[i] = (complex*)malloc(sizeof(complex)*row);
        init(row,tmp1[i]);
    }

    complex **tmp2 = (complex**)malloc(sizeof(complex*)*row); //tmp2[ROW][COL]
    complex **tmp3 = (complex**)malloc(sizeof(complex*)*row);
    for(int i=0;i<row;i++){
        tmp2[i] = (complex*)malloc(sizeof(complex)*col);
        init(row,tmp2[i]);
        tmp3[i] = (complex*)malloc(sizeof(complex)*col);
        init(row,tmp3[i]);
    }

    for(int i=0;i<col;i++){ //横向きにDCT
        DFT(row,input[i],tmp1[i]);
    }

    for(int i=0;i<col;i++){
        for(int j=0;j<row;j++){
            tmp2[j][i] = tmp1[i][j];
        }
    }

    for(int i=0;i<row;i++){ //縦向きDCT
        DFT(col,tmp2[i],tmp3[i]);
    }

    for(int i=0;i<col;i++){
        for(int j=0;j<row;j++){
            output[i][j] = tmp3[j][i];
        }
    }

    for(int i=0;i<row;i++){
        free(tmp1[i]);
        free(tmp2[i]);
        free(tmp3[i]);
    }
    free(tmp1);
    free(tmp2);
    free(tmp3);
}

void iDFT_2d(int col, int row, complex **input, complex **output)
{
    complex **tmp1 = (complex**)malloc(sizeof(complex*)*col); //tmp1[COL][ROW]
    for(int i=0;i<col;i++){
        tmp1[i] = (complex*)malloc(sizeof(complex)*row);
        init(row,tmp1[i]);
    }

    complex **tmp2 = (complex**)malloc(sizeof(complex*)*row); //tmp2[ROW][COL]
    complex **tmp3 = (complex**)malloc(sizeof(complex*)*row);
    for(int i=0;i<row;i++){
        tmp2[i] = (complex*)malloc(sizeof(complex)*col);
        init(row,tmp2[i]);
        tmp3[i] = (complex*)malloc(sizeof(complex)*col);
        init(row,tmp3[i]);
    }

    for(int i=0;i<col;i++){ //横向きにDCT
        iDFT(row,input[i],tmp1[i]);
    }

    for(int i=0;i<col;i++){
        for(int j=0;j<row;j++){
            tmp2[j][i] = tmp1[i][j];
        }
    }

    for(int i=0;i<row;i++){ //縦向きDCT
        iDFT(col,tmp2[i],tmp3[i]);
    }

    for(int i=0;i<col;i++){
        for(int j=0;j<row;j++){
            output[i][j] = tmp3[j][i];
        }
    }

    for(int i=0;i<row;i++){
        free(tmp1[i]);
        free(tmp2[i]);
        free(tmp3[i]);
    }
    free(tmp1);
    free(tmp2);
    free(tmp3);
}

//FFT

void make_W(int I,int n, complex *w)
{
    if(I==0){
        for(int i=0;i<n;i++){
            w[i] = e_powered(-2.0*(double)i*MY_PI/(double)n);
        }
    }
    else{
        for(int i=0;i<n;i++){
            w[i] = e_powered(2.0*(double)i*MY_PI/(double)n);
        }
    }
}

int calc_m(int n)
{
    int m=0;
    while(n>0){
        n/=2;
        m++;
    }
    return m-1;
}

void sort(int n,complex *c)
{
    complex *c2 = (complex *)malloc(sizeof(complex)*n);
    init(n,c2);
    for(int i=0;i<n;i++){
        c2[i] = c[i];
    }
    int m = calc_m(n);
    for(int i=0;i<n;i++){
        int index=0;
        int tmp=i;
        for(int j=0;j<m;j++){
            if(tmp == 0) break;
            if(tmp%2==1){
                index += pow(2,m-j-1);
            }
            tmp = tmp >> 1;
        }
        c[index] = c2[i];
    }
    free(c2);
}

void block(int I,int n, complex *input, complex *output)
{
    complex *w = (complex*)malloc(sizeof(complex)*n);
    init(n,w);
    make_W(I,n,w);
    init(n,output);
    int cmd=0;
    for(int i=0;i<n;i++){
        if(i<n/2){
            output[i] = sum(output[i],input[i]);
            output[i+n/2] = sum(output[i+n/2],input[i]);
        }
        else{
            output[i-n/2] = sum(output[i-n/2],c_c_multi(input[i],w[(i-n/2)%n]));
            output[i] = sum(output[i],c_c_multi(input[i],w[i%n]));
        }
    }
}

void step(int I,int n, complex *input, complex *output)
{
    int m=calc_m(n);
    complex *in = (complex *)malloc(sizeof(complex)*n);
    init(n,in);
    complex *out = (complex *)malloc(sizeof(complex)*n);
    init(n,out);
    complex *tmp=NULL;
    for(int i=0;i<n;i++){
        in[i] = input[i];
    }
    for(int i=1;i<=m;i++){  //ステップは1,2,3...
        for(int j=0;j<pow(2,m-i);j++){
            block(I,pow(2,i),(in+j*(int)pow(2,i)),(out+j*(int)pow(2,i)));
        }
        tmp = in;
        in = out;
        out = tmp;
    }
    for(int i=0;i<n;i++){
        output[i] = in[i];
    }

    free(in);
    free(out); //free(tmp)は，out=tmpのためしない
}

void FFT(int n,complex* input,complex* output)
{
    sort(n,input);
    step(0,n,input,output);
}

void iFFT(int n,complex* input,complex* output)
{
    sort(n,input);
    step(1,n,input,output);
    for(int i=0;i<n;i++){
        output[i] = divide(output[i],n);
    }
}

void FFT_2d(int col, int row, complex **input, complex **output)
{
    complex **tmp1 = (complex**)malloc(sizeof(complex*)*col); //tmp1[COL][ROW]
    for(int i=0;i<col;i++){
        tmp1[i] = (complex*)malloc(sizeof(complex)*row);
        init(row,tmp1[i]);
    }

    complex **tmp2 = (complex**)malloc(sizeof(complex*)*row); //tmp2[ROW][COL]
    complex **tmp3 = (complex**)malloc(sizeof(complex*)*row);
    for(int i=0;i<row;i++){
        tmp2[i] = (complex*)malloc(sizeof(complex)*col);
        tmp3[i] = (complex*)malloc(sizeof(complex)*col);
        init(row,tmp2[i]);
        init(row,tmp3[i]);
    }

    for(int i=0;i<col;i++){ //横向きにDCT
        FFT(row,input[i],tmp1[i]);
    }

    for(int i=0;i<col;i++){
        for(int j=0;j<row;j++){
            tmp2[j][i] = tmp1[i][j];
        }
    }

    for(int i=0;i<row;i++){ //縦向きDCT
        FFT(col,tmp2[i],tmp3[i]);
    }

    for(int i=0;i<col;i++){
        for(int j=0;j<row;j++){
            output[i][j] = tmp3[j][i];
        }
    }

    for(int i=0;i<row;i++){
        free(tmp1[i]);
        free(tmp2[i]);
        free(tmp3[i]);
    }
    free(tmp1);
    free(tmp2);
    free(tmp3);
}

void iFFT_2d(int col, int row, complex **input, complex **output)
{
    complex **tmp1 = (complex**)malloc(sizeof(complex*)*col); //tmp1[COL][ROW]
    for(int i=0;i<col;i++){
        tmp1[i] = (complex*)malloc(sizeof(complex)*row);
    }

    complex **tmp2 = (complex**)malloc(sizeof(complex*)*row); //tmp2[ROW][COL]
    complex **tmp3 = (complex**)malloc(sizeof(complex*)*row);
    for(int i=0;i<row;i++){
        tmp2[i] = (complex*)malloc(sizeof(complex)*col);
        tmp3[i] = (complex*)malloc(sizeof(complex)*col);
    }

    for(int i=0;i<col;i++){ //横向きにDCT
        iFFT(row,input[i],tmp1[i]);
    }

    for(int i=0;i<col;i++){
        for(int j=0;j<row;j++){
            tmp2[j][i] = tmp1[i][j];
        }
    }

    for(int i=0;i<row;i++){ //縦向きDCT
        iFFT(col,tmp2[i],tmp3[i]);
    }

    for(int i=0;i<col;i++){
        for(int j=0;j<row;j++){
            output[i][j] = tmp3[j][i];
        }
    }

    for(int i=0;i<row;i++){
        free(tmp1[i]);
        free(tmp2[i]);
        free(tmp3[i]);
    }
    free(tmp1);
    free(tmp2);
    free(tmp3);
}

//その他処理

//最大周波数の計算
//fsはサンプリング周波数
double calc_frequency(int n, double fs, complex *input)
{
    double max=0.0;
    int index=0;
    for(int i=0;i<n/2+1;i++){
        double tmp = power(input[i]);
        //printf("%d:%.1f\n",i,tmp);
        if(tmp > max){
            max = tmp;
            index = i;
            //printf("%d %f\n",index,tmp);
        }
    }
    //printf("index=%d, power=%f ",index, max);
    //printf(" %f\n",(double)index/n*fs);
    //printf("\nindex1 = %d\n",index);
    return (double)index/n*fs;
}

//2番目に高い周波数の計算
//fsはサンプリング周波数
//要チェック
double calc_second_frequency(int n, double fs, complex *input)
{
    double max1=0.0, max2=0.0;
    int index1=0, index2=0;
    for(int i=0;i<n/2-1;i++){
        double tmp_m = power(input[i]);
        int tmp_i = 0;
        //printf("%d:%.1f\n",i,tmp);
        if(tmp_m > max2){
            max2 = tmp_m;
            index2 = i;
            //printf("%d %f\n",index,tmp);
        }
        if(max1 < max2){
            tmp_m = max2;
            max2 = max1;
            max1 = tmp_m;
            tmp_i = index2;
            index2 = index1;
            index1 = tmp_i;
        }
    }
    //printf("index=%d, power=%f ",index, max);
    //printf(" %f\n",(double)index/n*fs);
    printf("\nindex2 = %d\n",index2);
    return (double)index2/n*fs;
}

//最大周波数のインデックス計算
int calc_max_index(int n, complex *input)
{
    double max=0.0;
    int index=0;
    for(int i=0;i<n/2+1;i++){
        double tmp = power(input[i]);
        if(tmp > max){
            max = tmp;
            index = i;
        }
    }
    return index;
}

//2番目に大きい周波数のインデックス計算
//要チェック
int calc_second_index(int n, complex *input)
{
    double max1=0.0,max2=0.0;
    int index1=0, index2=0;
    for(int i=0;i<n/2+1;i++){
        double tmp_m = power(input[i]);
        int tmp_i = 0;
        if(tmp_m > max2){
            max2 = tmp_m;
            index2 = i;
        }
        if(max1 < max2){
            tmp_m = max2;
            max2 = max1;
            max1 = tmp_m;
            tmp_i = index2;
            index2 = index1;
            index1 = tmp_i;
        }
    }
    return index2;
}

//オーバーラップしないSTFT
//n:全体の要素数
//d:幅
//outputは，n/d×dの二重ポインタ型の二次元配列
void STFT(int n,int d,complex *input,complex **output)
{
    complex **tmp = (complex **)malloc(sizeof(complex*)*n/d);
    for(int i=0;i<n/d;i++){
        tmp[i] = (complex*)malloc(sizeof(complex)*d);
        init(d,tmp[i]);
    }
    for(int i=0;i<n/d;i++){
        for(int j=0;j<d;j++){
            tmp[i][j] = input[i*d+j];
        }
    }
    for(int i=0;i<n/d;i++){
        FFT(d,tmp[i],output[i]);
    }

    for(int i=0;i<n/d;i++){
        free(tmp[i]);
    }
    free(tmp);
}

//50%でオーバーラップするSTFT
//n:全体の要素数
//d:幅
//outputは，n/d*2 × dの二重ポインタ型の二次元配列
//outputは，n/d × 2 × dの三重ポインタ型の三次元配列
void STFT_over(int n,int d,complex *input,complex ***output)
{
    complex ***tmp = gen_triple_pointer_complex(n/d,2,d);

//inputからtmpに代入
    for(int i=0;i<n/d;i++){
        for(int j=0;j<2;j++){
            for(int k=0;k<d;k++){
                if(i*d+k+j*d/2 < n){
                    tmp[i][j][k] = input[i*d+k+j*d/2];
                }
                else{
                    tmp[i][j][k].re = 0; tmp[i][j][k].im = 0;
                    //printf("caution! ---in STFT_over, input[%d]\n",i*d+k+j*d/2);
                }
            }
        }
    }
//tmpをoutputにFFT
    for(int i=0;i<n/d;i++){
        for(int j=0;j<2;j++){
            FFT(d,tmp[i][j],output[i][j]);
        }
    }

    for(int i=0;i<n/d;i++){
        for(int j=0;j<2;j++){
            free(tmp[i][j]);
        }
        free(tmp[i]);
    }
    free(tmp);
    //printf("\n");
}

//cutする部分を作る
void STFT_cut(int n,int d,complex *input, complex *output)
{
    complex **stft = gen_double_pointer_complex(n/d,d);
    STFT(n,d,input,stft);
    for(int i=0;i<n/d;i++){
        double sum=0;
        int x=400000;
        for(int j=0;j<d;j++){
            if(stft[i][j].re > x || stft[i][j].im > x){
                stft[i][j].re *= 100;
                stft[i][j].im += 100;
            }
            else{
                stft[i][j].re /= 100;
                stft[i][j].im /= 100;
            }
        }
    }
/*
    for(int i=0;i<n/d;i++){
        int index = calc_max_index(d,stft[i]);
        printf("index = %d, %f\n",index,(double)index/d*44100);
        for(int j=0;j<d;j++){
            if(j != index){
                stft[i][j].re=0; stft[i][j].im=0;
            }
            //printf("%d,%d: %f %f\n",i,j,stft[i][j].re,stft[i][j].im);
        }
    }
*/
    for(int i=0;i<n/d;i++){
        iFFT(d,stft[i],(output+i*d));
    }
}

//DCTに対するhi_cut(low_pass)
//pはカットオフ周波数
void hi_cut_DCT(double p, double* c, int n, double fs)
{
    int index = p*(n/(fs/2));
    for(int i=0;i<n;i++){
        if(i > index){
            c[i]=0.0;
        }
    }
}

//DCTに対するlow_cut(hi_pass)
//pはカットオフ周波数
void low_cut_DCT(double p, double* c, int n, double fs)
{
    int index = p*(n/(fs/2));
    for(int i=0;i<n;i++){
        if(i < index){
            c[i]=0.0;
        }
    }
}

//DCTに対するband_cut
void band_cut_DCT(double high,double low,double* c, int n, double fs)
{
    low_cut_DCT(low,c,n,fs);
    hi_cut_DCT(high,c,n,fs);
}

//FFTに対するhi_cut(low_pass)
//pはカットオフ周波数
void hi_cut_FFT(double p, complex* input, complex * output, int n, double fs)
{
    int index = p*(n/fs);
    for(int i=0;i<n;i++){
        if(i>index && i<n-index){
            output[i].re=0.0; output[i].im=0.0;
            //c[n-i-1].re=0.0; c[n-i-1].im=0.0;
        }
        else{
            output[i] = input[i];
        }
    }
}

//FFTに対するlow_cut(hi_pass)
//pはカットオフ周波数
void low_cut_FFT(double p, complex* input, complex* output, int n, double fs)
{
    int index = p*(n/fs);
    //printf("index=%d\n",index);
    for(int i=0;i<n;i++){
        if(i<index || i>n-index){
            output[i].re=0.0; output[i].im=0.0;
            //c[n-i-1].re=0.0; c[n-i-1].im=0.0;
        }
        else{
            output[i] = input[i];
        }
    }
}

//FFTに対するband_cut
//high,lowは周波数
void band_cut_FFT(double high,double low,complex* input, complex* output, int n, double fs)
{
    low_cut_FFT(low, input, output, n, fs);
    hi_cut_FFT(high, input, output, n, fs);
}

/*
//スペクトルサブトラクション
void subtraction(int n, complex *input, complex *output)
{
    //band_cut(1000,100,f,N,FS);
    int cnt = 5;
    int X = 1024;
    complex *f_silent = gen_pointer_complex(X);
    complex *F_silent = gen_pointer_complex(X);
    complex *F = gen_pointer_complex(X);

    for(int j=0;j<cnt;j++){
        for(int i=0;i<X;i++){
            f_silent[i].re += input[i+j*X].re/cnt;
        }
    }
    FFT(X,f_silent,F_silent);

    for(int i=0;i<n/X;i++){
        FFT(X,(input+i*X),(F+i*X));
        for(int j=0;j<X;j++){
            F[i*X+j].re -= F_silent[j].re;
            F[i*X+j].im -= F_silent[j].im;
            if(F[i*X+j].re < 0) F[i*X+j].re=0;
            if(F[i*X+j].im < 0) F[i*X+j].im=0;
        }
        iFFT(X,(F+i*X),(output+i*X));
    }
}

*/

/*
//移動平均フィルタ
void move_ave_filter(int n, complex *input, complex *output)
{
    int FILTER_SIZE = 64;
    for(int i=0;i<n;i++){
        int start = i-FILTER_SIZE/2;
        int end = i+FILTER_SIZE/2;
        if(start < 0) start=0;
        if(end > n-1) end=n-1;
        complex sum;
        sum.re=0.0; sum.im=0.0;

        for(int j=0;j<=end;j++){
            sum.re += input[j].re;
            sum.im += input[j].im;
        }
        output[i].re = sum.re / (end-start+1);
        output[i].im = sum.im / (end-start+1);
    }
}

*/

#endif //SHINGOUSYORI_H_
