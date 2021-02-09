#include <stdio.h>
#include <complex.h>
#include "mpi.h"
#include <stdlib.h>
#include <string.h>
#define w 2028
#define h 2048
#define P 2
#define iteration 10

unsigned char cal_pixel(double complex d, double b,int N)
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
  int rank, dest, source, len, tag, rc, i;
  int num_proc = 2;
  int count = 0;
  int row = 0;
  MPI_Status status;
  rc = MPI_Init(&argc, &argv);
  rc = MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
  rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  double b=2;
  int dx = 2 * b / (w - 1);
  int dy = 2 * b / (h - 1);
  double cx = -b;
  double cy = -b;
  int x, y;
  int color[w][h];
  int image[w][h] = {0};
  int img_middle = w/2;
  int N = 256;
  if (rank == 0){//master
    /*
    if(argc<2){
      rc = MPI_Abort(MPI_COMM_WORLD, 1);
      exit(1);
    }
    */
    for(i=1;i<num_proc;i++){
    // rc = MPI_Recv(answer, 2, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD,&status);
    // answer[0] -- from whom
    // answer[1] -- '-1' or the X coordinate
      rc = MPI_Send(&image[i*(h/2)], 1, MPI_INT, i+1, 0, MPI_COMM_WORLD);
    }
    for(i=0;i<w*h;i++){
      rc = MPI_Recv(&image[i],1,MPI_INT,MPI_ANY_SOURCE, 0, MPI_COMM_WORLD,&status);
      printf("%p",color);
    }
      //rc = MPI_Send(&image[w/2], 2, MPI_DOUBLE, i, 4, MPI_COMM_WORLD);
    /*
    for(i=1;i<num_proc;i++){
      rc = MPI_Recv(answer, 2, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
      rc = MPI_Recv(&image[answer[1]][0],h,MPI_INT, answer[0],2, MPI_COMM_WORLD, &status);
      rc = MPI_Send(&term, 2, MPI_INT, answer[0], 3, MPI_COMM_WORLD);
    */
  }

  else{ //slave 
    rc = MPI_Recv(&image[0], 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
    for (x = 0; x < w - 1; x++){
      float dreal = x * dx - b;
      for (y = 0; y < h/num_proc - 1; y++){
        float dimag =  y * dy - b;
        complex double d = dreal + I * dimag;
        color[x][y] = cal_pixel(d, b, N);
      }
    }
  /*
    for(x=0; x<w; x++){
      for(y=row; y<row; y++){
        c.real=min_real + (float(x*scale_real));
        c.imag=min_imag + (float(y*scale_imag));
        color = cal_pixel(c);
        rc = MPI_Send(&c, &color,to master, MPI_COMM_WORLD);
      }
    }*/
  }
  
  rc = MPI_Finalize();
  printf("%p",image);
  return 0;
}
