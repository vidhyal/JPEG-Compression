#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


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

int quant[8][8];
double qScale;
int main (int argc, char ** argv){
	FILE * iImg, *qFile, *oFp;
	iImg = fopen(argv[1], "r");       // open input image file read only mode
       qFile = fopen(argv[2],"r");    // open quantization file read only mode
    oFp = fopen(argv[3],"w");         // open output file write only mode
    if (iImg == NULL) {
          fprintf(stderr, "Cannot open input file %s",argv[1]);
          exit(1);
      }
    if (qFile == NULL) {
          fprintf(stderr, "Cannot open input file %s",argv[2]);
          exit(1);
    }
    if (oFp == NULL) {
          fprintf(stderr, "Cannot open input file %s",argv[2]);
          exit(1);
    }
    
    int count, i,j;
	// read and save quantization matrix in quant
    for (i =0; i< 8; i++){
      for (j=0; j<8; j++){
        fscanf(qFile, "%d", &quant[i][j]);
      }
    }
    fclose(qFile);
    
    int xsize, ysize, maxVal=255;
    char line[40];
    // read header of input dct file
    fscanf(iImg, "%s", line);
    fscanf(iImg, "%d" "%d", &xsize, &ysize);
    fscanf(iImg, "%lf", &qScale);
    // write the header to teh output pgm file
    fprintf(oFp,"%s\n%d %d\n%d\n", "P5", xsize, ysize, maxVal);
    int qt[8][8],fun[8][8];
    int output[ysize][xsize];
    double F[8][8], Sum[8][8];
    int c=0, m,n,u,v;
    
    //loop for all input and variable c is the loop variable
    while(c<(xsize*ysize/256*4)) {
    	//read the starting points (pixel) of the 8X* blocks
    	fscanf(iImg, "%d %d", &n, &m);
    	
    	// read and unzig-zag the input matrix into qt 
    	u=0,v=0;
    	fscanf(iImg, "%d", &qt[u][v]);
    	while ( u<8 && v<8){
        	v++;
        	if(u<8 && v<8){
    			fscanf(iImg,"%d", &qt[u][v]);
       		}
       		else if (u<7 && v>0){
            	v--;
            	u++ ;
	            fscanf(iImg,"%d", &qt[u][v]);
        	}
			else break;
    	    while (v>0 && u<7 && v<8){
    			v--;
        		u++;
	    		fscanf(iImg,"%d", &qt[u][v]);
			}
    	    u++;
            if(u<8 && v<8){
            	fscanf(iImg,"%d", &qt[u][v]);
    		}
	    	else if (u>0 && v<7) {
            	u--;
	    		v++ ;
				fscanf(iImg,"%d", &qt[u][v]);
			}
			while(u>0 && u<8 && v<7){
				u--;
				v++;
				fscanf(iImg,"%d", &qt[u][v]);
			}
		}
		c++;
        // multiplying the coefficients with the quantization matrix values and quantization factor 
		for (v =0; v<8; v++){
			for (u=0; u<8; u++){
				qt[v][u] = qt[v][u] -127;
	        	F[v][u] = (float)qt[v][u] * qScale* quant[v][u];
			}
		}
		
		// code for IDCT
		for (i=0; i<8; i++){
       		for (j =0; j<8;j++){
         		Sum[i][j]= 0;
         		for (v=0; v<8; v++){ 
                    for(u=0; u<8;u++){
                        if(F[v][u] != 0){
	                        Sum[i][j] = Sum[i][j]+ (F[v][u] *coeff(u,v)* cos((double)((((2*i) +1)* v* M_PI)/16))* cos((double)((((2*j) +1)* u* M_PI)/16)));// summations inside the IDCT formula
						}
					}
        	 	}	
				// Result of IDCT  
	    	    fun[i][j] =  (int)round(Sum[i][j]/4);
			// set the range of function value to be [-127,128]
                if (fun[i][j]>255)
                    fun[i][j]=255;
                if (fun[i][j]<0)
                    fun[i][j] =0;
                output[i+m][j+n]= fun[i][j];
            }
        }    	
    }
	// writing out the 8 bit per pixel values of the entire image to the pgm file
    for(i=0; i<ysize;i++){
		for (j=0; j<xsize; j++){
    		char ch = output[i][j];
			fputc(ch, oFp);
		}	
	}
	
	fclose(iImg);
	fclose(oFp);
}
