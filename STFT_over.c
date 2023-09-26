#include"shingousyori.h"

#define N 4096
#define FS 100
#define T 40.96
#define D 256

int main()
{
    FILE *fp=fopen("FFT_a.txt","r");
    if(fp==NULL){
        printf("オープン失敗\n");
        return 0;
    }
    complex *f = gen_pointer_complex(N);
    complex ***F = gen_triple_pointer_complex(N/D,2,D);
    real_fscanf(N,f,fp);
    STFT_over(N,D,f,F);

    for(int i=0;i<N/D;i++){
        for(int j=0;j<2;j++){
            double frequency = calc_frequency(D,FS,F[i][j]);
            printf("%.2fHz\n",frequency);
        }
    }

    fclose(fp);
}
