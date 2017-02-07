#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>



int quant[8][8];
double qScale;

// This function returns the value of Cu*Cv according to the value of u and v
double coeff (int u, int v) {
 double Cu =1.0, Cv =1.0;
double pdt;    
    if (u==0)
      Cu = ((double)1.0)/(sqrt(2));
    if( v==0)
      Cv= ((double)1.0)/(sqrt(2));
    pdt = (Cu*Cv);
    return pdt;

}

// This function is used to print the dct output in zig zag reorder
void printZZ( FILE *oFp, int val, int itr)
{
        if (itr%8 ==0 && itr !=0)
            fprintf(oFp,"\n");    
	fprintf(oFp,"%5d",val);

}

// This function performs DCT, Quantization and calls the function to print the output in zig-zag reorder
void encode (FILE *oFp, int array[8][8]){
	//qt is the matrix that holds the final integer value of a pixel after DCT and quantization
	int u,v,i,j , qt[8][8];
	// matrix F represents the result of DCT.
	double F[8][8], Sum[8][8];
	for (v=0;v<8;v++){
	   for(u=0;u<8;u++){
	   	    Sum[u][v]=0.0;
	   	    for (i=0;i<8;i++){
			    for(j=0;j<8;j++){
	   	            Sum[u][v] = Sum[u][v] + (array[i][j] * cos((double)((((2*i) +1)* u* M_PI)/16))* cos((double)((((2*j) +1)* v* M_PI)/16)));
				}
			}
	   	    F[u][v]= ((float)(coeff(u,v) * Sum[u][v]))/4;
	   	    if(quant[u][v] !=0)
	            qt[u][v] = (int)round(F[u][v]/(qScale*quant[u][v]));
	        else
	            qt[u][v]= (int)round(F[u][v]);
	        // set the range of qt to be [0,255]
	        if(qt[u][v]>128)
	           qt[u][v]=128;
	        if (qt[u][v] < -127)
	           qt[u][v]= -127;
	        qt[u][v] = qt[u][v]+127;
	   }
    }
  // This part of the function zig-zag reorders the matrix qt for output to file
   
    int itr =0;
	u=0; 
    v=0;
    printZZ(oFp, qt[u][v], itr);
    itr++;
    while ( u<8 && v<8){
        v++;
        if(u<8 && v<8){
    		printZZ(oFp, qt[u][v], itr);
    		itr++;
    	}
    	else if (u<7 && v>0){
	        v--;
	        u++ ;
	        printZZ(oFp, qt[u][v], itr);
	    	itr++;
		}
		else break;
    	while (v>0 && u<7 && v<8){
    		v--;
        	u++;
			printZZ(oFp, qt[u][v], itr);
    		itr++;
		}
        u++;
        if(u<8 && v<8){
          	printZZ(oFp, qt[u][v], itr);
   			itr++;
		}
	    else if (u>0 && v<7) {
            u--;
	    	v++ ;
			printZZ(oFp, qt[u][v], itr);
    		itr++;
		}
		while(u>0 && u<8 && v<7){
			u--;
			v++;
			printZZ(oFp, qt[u][v], itr);
    		itr++;
		}
	}
    fprintf(oFp,"\n");
}







int main(int argc, char ** argv){

    FILE *iImg, *qFile,  *oFp;
    iImg = fopen(argv[1], "r");  // input image file
    qFile = fopen(argv[2],"r"); // Quantization file
    qScale = atof(argv[3]);  // The value of quantization factor
    

    if (iImg == NULL) {
        fprintf(stderr, "Cannot open input file %s",argv[1]);
        exit(1);
    }
     
    if (qFile == NULL) {
        fprintf(stderr, "Cannot open input file %s",argv[2]);
        exit(1);
    }

    int count, i,j;
    // read in the quantization matrix
    for (i =0; i< 8; i++){
      	for (j=0; j<8; j++){
        	fscanf(qFile, "%d", &quant[i][j]);
        }
    }
	fclose(qFile);

    int xsize, ysize, maxVal;
    char line[40];
    // read in the header of the input pgm file.
    fscanf(iImg, "%s", line);
    fscanf(iImg, "%d" "%d", &xsize, &ysize);
    fscanf(iImg, "%d", &maxVal);
    // y corresponds to the row (vertical variation) and x corresponds to column
    int input[ysize][xsize];  
    int c;
    // read the characters in the pgm file, convert them to corresponding ASCII code 
	// and save it in the matrix for input image pixels
    c = fgetc(iImg);
    for (i=0; i<ysize; i++){
        for (j =0; j<xsize;j++){
            c = fgetc(iImg);
            unsigned char ch = (unsigned char) c;
            input[i][j] = ch;
        }
    }
    fclose(iImg);

    // open the output file and write its header
    oFp = fopen(argv[4],"w");
    if(oFp ==NULL) {
          fprintf(stderr, "Cannot open output file %s", argv[3]);
          exit(1);
      }
    char *str = "MYDCT\n" ;
    fprintf(oFp, "%s%d %d\n%lf\n", str,xsize,ysize,qScale);    

    // a 16X16 macroblock 'a' is split into 4 8X8 macroblocks: a1, a2,a3,a4
    int a[16][16], a1[8][8], a2[8][8], a3[8][8], a4[8][8], k,l;
        for (k=0;k<ysize;k=k+16){
        	for (l=0;l<xsize; l=l+16){
    		for (i=0; i<16; i++){	
                for (j=0; j<16; j++){
                    a[i][j]= input[k+i][l+j];
		        }
    	   }
    	   for (i=0; i<8; i++){	   
                for (j=0; j<8; j++){
			        a1[i][j] = a[i][j];
			        a3[i][j] = a[i+8][j];
			        a2[i][j] = a[i][j+8];
			        a4[i][j] = a[i+8][j+8]; 
               }
           }
           // encode each of the four 8X8 blocks with the help of encode function
    	    fprintf(oFp,"%d %d\n", l,k);
            encode(oFp,a1);
    	    fprintf(oFp,"%d %d\n", l+8, k);	
    	    encode(oFp,a2);
    	    fprintf(oFp,"%d %d\n", l,k+8);	
    	    encode(oFp,a3);
    	    fprintf(oFp,"%d %d\n", l+8,k+8);	
    	    encode(oFp,a4);
    	        		
		}
  }
    fclose(oFp);
return 0;
}

