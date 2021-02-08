#include <stdio.h>
#include <complex.h>
#include "mpi.h"
#include <stdlib.h>
#include <string.h>
#define w 2028
#define h 2048
#define N 256
#define b 2
#define P 2
#define iteration 10

int converges(double cx,double cy){
  int n=0;
  double zx=0;
  double new_zx=0;
  double zy=0;
  // we iterate until max. iteration count iter_n, or until z^2 (complex!) runs over 4 - this means, our series will run to infinity, so it's not part of the set
  while((n<iteration) && (zx*zx + zy*zy)<4){
  // z * z => new_zx = (zx*zx - zy*zy)  new_zy = (zx*zy + zx*zy)
  // we work with complex numbers
  // z*z + c = zx^2 - zy^2 +cx   +  i(zx*zy*2 + cy) 
    new_zx = zx*zx - zy*zy + cx;
    zy = 2*zx*zy + cy;
    zx = new_zx;
    n++;
  }
  return n;
}

int main(int argc, char **argv)
{
  int rank, dest, source, len, num_proc, tag, rc, i;
  MPI_Status status;
  rc = MPI_Init(&argc, &argv);
  rc = MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
  printf("%d",num_proc);
  rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int Tag1=6,Tag2=10;  /*Just Random int s. There is a max size consideration.*/
  int dx = 2 * b / (w - 1);
  int dy = 2 * b / (h - 1);
  double cx = -b;
  double cy = -b;
  int x, y;
  int color[w][h];
  int image[w][h] = {0};
  int img_middle = w/2;
  int img_line[h] = {0};
  int answer[2];
  double question;
  tag = 100;
  
  if (rank == 0){//master
    /*
    if(argc<2){
      rc = MPI_Abort(MPI_COMM_WORLD, 1);
      exit(1);
    }
    */
    for(i=0;i<num_proc;i++){
      rc = MPI_Recv(answer, 2, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD,&status);
    // answer[0] -- from whom
    // answer[1] -- '-1' or the X coordinate
      
      if(answer[1]>=0){ //not the first answer
      rc = MPI_Recv(&img_array[answer[1]][0], h, MPI_INT, answer[0],2, MPI_COMM_WORLD, &status);
      }
      
      rc = MPI_Send(&image[i*(w/2)], 2, MPI_INT, i, 3, MPI_COMM_WORLD);
      rc = MPI_Send(&image[w/2], 2, MPI_DOUBLE, i, 4, MPI_COMM_WORLD);
    }
    int term=-1;
    for(i=1;i<num_proc;i++){
      rc = MPI_Recv(answer, 2, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
      rc = MPI_Recv(&image[answer[1]][0],h,MPI_INT, answer[0],2, MPI_COMM_WORLD, &status);
      rc = MPI_Send(&term, 2, MPI_INT, answer[0], 3, MPI_COMM_WORLD);
    }
  }
  else{ //slave 
    int i;
    answer[0]=rank;
    answer[1]=-1;
    rc = MPI_Send(answer, 2, MPI_INT, 0, 1, MPI_COMM_WORLD);
    while(1){
      rc = MPI_Recv(&i, 2, MPI_INT, 0, 3, MPI_COMM_WORLD, &status);
      if(i<0) break; //got termination command!
        answer[1]=i;
      rc = MPI_Recv(&question, 2, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, &status);
      cy=-b;  // at every new step in X direction, we start at the first Y value 
      for(int j=0;j<h;j++){
        img_line[j]=converges(question,cy);
 	    cy=cy+dy;
      }
      rc = MPI_Send(answer, 2, MPI_INT, 0, 1, MPI_COMM_WORLD);
      rc = MPI_Send(img_line, h, MPI_INT, 0, 2, MPI_COMM_WORLD);
    }
  }
  rc = MPI_Finalize();
  printf("%p",image);
  return 0;
}
