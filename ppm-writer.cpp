#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <malloc.h>  
#include <memory.h>

typedef struct {
  unsigned char r;
  unsigned char g;
  unsigned char b;
} RGBDATA;


void WritePPM(int XDIM, int YDIM, RGBDATA *dataset, FILE* fp)
{
  int i,j;
  //printf("we are in writeppm function!");

  fprintf(fp, "P6\n%d %d\n%d\n", XDIM, YDIM, 255);
  for (j=0; j<YDIM; j++)
    for (i=0; i<XDIM; i++) {
      fputc(dataset[j*XDIM+i].r, fp);
      fputc(dataset[j*XDIM+i].g, fp);
      fputc(dataset[j*XDIM+i].b, fp);
    }

  fclose(fp);
  
}
