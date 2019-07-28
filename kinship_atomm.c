#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "nrutil.c"


#define MAXLEN 2048
#define NTIMES 2
#define MIN 10e-5

char genotypefile[MAXLEN]="genotype.txt",sizefile[MAXLEN]="size.txt",kinshipfile[MAXLEN]="kinship.txt",freqfile[MAXLEN]="freqfile.txt";

FILE *genotypef,*sizef,*memberf,*kinshipf,*freqf;

int n1;
int **G;
double **center_G;
int nSNP_total;



int main(int argc, char *argv[])
{
    int i,t,k;
    int gfile=0,afile=0,ofile;
    if (argc>1){
        for (i=1;i<argc && argv[i][0]=='-';i++){
            switch(argv[i][1])
            {
                case 'g':
                    strncpy(genotypefile,argv[++i],MAXLEN);
                    printf("genoype input: %s\n",genotypefile);
                    gfile=1;
                    break;
                case 's':
                    strncpy(sizefile,argv[++i],MAXLEN);
                    printf("size input: %s\n",sizefile);
                    gfile=1;
                    break;
                case 'k':
                    strncpy(kinshipfile,argv[++i],MAXLEN);
                    printf("user specified output: %s\n",kinshipfile);
                    ofile=1;
                    break;
                case 'f':
                    strncpy(freqfile,argv[++i],MAXLEN);
                    printf("user specified output: %s\n",freqfile);
                    ofile=1;
                    break;
                default:
                    printf ("Unknown option");
                    exit(1);
            }    
        }
    }

    
    
    if((genotypef=fopen(genotypefile,"r"))==NULL){
        printf("Cannot open genotype file.\n");
        exit(1);
    }
    
    if((sizef=fopen(sizefile,"r"))==NULL){
        printf("Cannot open size file.\n");
        exit(1);
    }
    
    if((kinshipf=fopen(kinshipfile,"w"))==NULL){
        printf("Cannot open kinship file.\n");
        exit(1);
    }
    
    if((freqf=fopen(freqfile,"w"))==NULL){
        printf("Cannot open freq file.\n");
        exit(1);
    }                                                                       
                                  
    
    double thresh=0.01;
    
    fscanf(sizef,"%d %d %lf",&n1, &nSNP_total,&thresh);
    
    G=imatrix(1,n1,1,nSNP_total);
    center_G=dmatrix(1,n1,1,nSNP_total);
    
    int j=0;
    int p=0;

    int *nonmissing=ivector(1,n1);
    double **totalfreq=dmatrix(1,nSNP_total,1,2);
    
    int chr, pos;
    double value;
    int nonmissingcount;
    char SNP[n1][MAXLEN];
    
    double **cov=dmatrix(1,n1,1,n1);
    double **nonmissing_pair=dmatrix(1,n1,1,n1);
    
    char char1[MAXLEN];
    
    double x1,x2,x3;
    int nallele;
    

    for(i=1;i<=n1+1;i++){
        fscanf(genotypef,"%s ",char1);
    }
    
    char firstallele[MAXLEN];
    
    int snp=0;
    int nclass;

    
    fprintf(freqf,"index\t allel1_freq\t allel2_freq\t deletion_freq\n");
    
      printf("Start computing genetic relatedness matrix (GRM) using both mutation and deletion polymorphisms...\nVariants with minor allele frequency (MAF) less than %.3lf are removed from the estimation...\n", thresh);
    
    
    while(fscanf(genotypef,"%d",&chr)==1){
        
        nallele=0;
        nclass=1;
        
        p++;
        totalfreq[p][1]=0;
        totalfreq[p][2]=0;
        nonmissingcount=0;
        
        for(i=1;i<=n1;i++){
            nonmissingcount++;
            fscanf(genotypef,"%s",char1);
            strcpy(SNP[i],char1);
            
            
            if ((strcmp(SNP[i],"N")!=0)&(strcmp(SNP[i],"-")!=0)&(nallele==0)){
                strcpy(firstallele,SNP[i]);
                nallele=1;
                G[i][p]=1;//type 1, mutation allele
                totalfreq[p][1]++;
            }
            
            
            else if(strcmp(SNP[i],"-")==0){
                G[i][p]=2;//type 2 (optional; deletion)
                nclass=2;
                totalfreq[p][2]++;
            }
            
            else if (strcmp(SNP[i],"N")==0){
                nonmissingcount--;
                G[i][p]=-9;
            }
            
            
            else if ((strcmp(SNP[i],firstallele)==0)){
                G[i][p]=1;//type 1
                totalfreq[p][1]++;
            }
            
            else if ((strcmp(SNP[i],firstallele)!=0)){
                G[i][p]=3;//type 3, mutation allele
            }
        }
        
        totalfreq[p][1]=(totalfreq[p][1])/(nonmissingcount);
        totalfreq[p][2]=(totalfreq[p][2])/(nonmissingcount);
        
        
        fprintf(freqf,"%d\t%.3lf\t%.3lf\t%.3lf\n",chr,totalfreq[p][1],1-totalfreq[p][1]-totalfreq[p][2],totalfreq[p][2]);
        
        x1=fmax(totalfreq[p][1],(totalfreq[p][1]==0));
        x2=fmax(totalfreq[p][2],(totalfreq[p][2]==0));
        x3=fmax(1-totalfreq[p][1]-totalfreq[p][2],((1-totalfreq[p][1]-totalfreq[p][2])<MIN));
      
       
        //consider genetic variant with mutation polymorphism; deletion polymorphism is optional
        if((totalfreq[p][1]!=0)*((totalfreq[p][1]+totalfreq[p][2])!=1)*(fmin(fmin(x1,x2),x3)!=1)*(fmin(fmin(x1,x2),x3)>=thresh)){
            snp++;
            
            for(i=1;i<=n1;i++){
                for(j=1;j<=n1;j++){
                    
                    if((G[i][p]!=(-9))&(G[j][p]!=(-9))){
                        nonmissing_pair[i][j]++;
                        
                        if((G[i][p]==G[j][p])&(G[i][p]==1))
                            cov[i][j]+=(1-totalfreq[p][1])/totalfreq[p][1];
                        
                        
                        
                        else if ((G[i][p]==G[j][p])&(G[i][p]==2))
                            
                            cov[i][j]+=(1-totalfreq[p][2])/totalfreq[p][2];
                        
                        
                        else if ((G[i][p]==G[j][p])&(G[i][p]==3))
                            
                            
                            cov[i][j]+=(totalfreq[p][1]+totalfreq[p][2])/(1-totalfreq[p][1]-totalfreq[p][2]);
                        
                        else
                            
                            cov[i][j]--;                       
                        }
                    }
            }
        }
    }
    
    
    
    for(i=1;i<=n1;i++){
        for(j=1;j<=n1;j++){
            cov[i][j]=cov[i][j]/nonmissing_pair[i][j];
            fprintf(kinshipf,"%lf\t",cov[i][j]);
            
        }
        fprintf(kinshipf,"\n");
    }
    
    printf("GRM estimation is done.\n");
    
    printf("Total number of genetic variants considered in the GRM estimation: %d\n",snp);
    
}
