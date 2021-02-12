#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <complex.h> //must be compiled with -lm flag
#include <unistd.h>
#include <stdlib.h>

unsigned char pixel_value(double complex d, double b, unsigned char N){
  /*
  forumla from lecture notes: z = 0 initially,
  z = z*z + d, where d is the complex coordinates of the point (re,im)
  N is the maximum number of allowed max iterations
  b = maximum radius of z.
  out: the number of iterations required for z > N
  */

  complex double z = 0;
  unsigned char count = 0;

  while (cabs(z)<b && count < N) {
    z = z*z + d;
    count++;
  }
  return count;
}

 int main(int argc, char **argv) {
   int rank, size, tag, rc, i, j, x, y;
   double b;
   MPI_Status status;
   char message[20];
   tag = 100;

   rc = MPI_Init(&argc, &argv);
   rc = MPI_Comm_size(MPI_COMM_WORLD, &size);
   rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

   //windowing parameters:
   long int h;
   long int w;
   unsigned char N = 255;

   if (argc < 2){
     h = 512;
     w = 512;
   }else{
     h = atoi(argv[1]); //pixel height
     w = atoi(argv[1]); //pixel width
   }
   int partition_width = h/size;

   unsigned char partition[h*partition_width];

   b = 2;

   //exit if the number of processors do not divides the number of columns evenly
   int mod = h%size;
   if (mod != 0) {
     printf("Uneven distribution of columns per process\n");
     exit(1);
   }

   printf("%d\n", mod);
   // double zoom = 0.00390625/4;
   double zoom = 2.0;
   //scaling and centering.
   double dx = zoom*b/(w-1); //step size
   double dy = zoom*b/(h-1);

   double x_offset = 0; // 2.057; // 1.25;//2.057; // offsets
   double y_offset = 0; // 2.656; // 1.85;//2.656;

   if (rank == 0) {
     printf("Number of nodes available: %d\n", size);
   }
   printf("Node with rank %d checking in.\n", rank);

   double complex d;
   double im;
   double re;
   int n = 0;
   for (x = rank*partition_width; x <partition_width*(rank+1); x++) {
     re = x*dx-b+ x_offset;
     for (y = 0; y < h; y++) {
       im = y*dy-b + y_offset;
       d = re + I*im;
       partition[n] = pixel_value(d, b, N);
       ++n;
     }
   }

   printf("Node %d finalized, gathering data. \n",rank);

   unsigned char map[h*w];
   rc = MPI_Gather(partition, partition_width*h, MPI_UNSIGNED_CHAR, &map, partition_width*h, MPI_UNSIGNED_CHAR, 0, MPI_COMM_WORLD);

if (rank == 0) {
   FILE *fp;
   int count = 0;
   fp = fopen(argv[2],"w");
   for (int x = 0; x < w; x++) {
     for (int y = 0; y < h; y++) {
       fprintf(fp, "%hhu ", map[x*w+y]);
     }
     fprintf(fp, "\n");
     //printf("LINE BREAK: %d\n", count);
     count++;
   }
   fclose(fp);
 }
   rc = MPI_Finalize();

   return 0;
}