#include <stdio.h>
#include <complex.h>
#include "mpi.h"
#include <stdlib.h>
#include <string.h>
#define P 2

unsigned char cal_pixel(double complex d, double b,unsigned char N)
{
  complex double z = 0;
  unsigned char count = 0;

  while (cabs(z) < b && count < N)
  {
    z = z * z + d;
    count++;
  }
  return count;
}

int main(int argc, char **argv)
{
  int w= 2048;
  int h= 2048;
  int tag = 100;
  int rank,x,y,i,rc;
  int num_proc=2;
  int count = 0;
  int row = 0;
  double b=2;
  int dx = 2 * b / (w - 1);
  int dy = 2 * b / (h - 1);
  unsigned char N = 255;
  MPI_Status status;
  rc = MPI_Init(&argc, &argv);
  rc = MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
  rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int image[1024][2048] = {0};
  /*
  image = (int**)malloc(w*sizeof(int *));
  for (i=0; i<w; i++){
    image[i]=(int *)malloc(h/num_proc*sizeof(int));
  }
  */
  if (rank == 0){//master
    for(i=1;i<num_proc;i++){
      rc = MPI_Send(&image, w*(h/num_proc), MPI_INT, i, 0, MPI_COMM_WORLD);
      printf("%d\n",i);
    }
  }

  else{ //slave 
    rc = MPI_Recv(&image, w*(h/num_proc), MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    printf("%p",image);
    /*
    for (x = 0; x < w - 1; x++){
      float dreal = x * dx - b;
      for (y = 0; y < h/num_proc - 1; y++){
        float dimag =  y * dy - b;
        complex double d = dreal + I * dimag;
        image[x][y] = cal_pixel(d, b, N);
        printf("%d\n",image[x][y]);
      }
    }
    */
    rc = MPI_Send(&image, w*(h/num_proc), MPI_INT, 0, 0, MPI_COMM_WORLD);  
  }
 /*
  if (rank == 0){//master
    for(i=1;i<num_proc;i++){
      rc = MPI_Recv(&image, w*(h/num_proc), MPI_INT, i, 0, MPI_COMM_WORLD, &status);
    }
  }
*/
  rc = MPI_Finalize();
  printf("%p",image);
  return 0;
}
