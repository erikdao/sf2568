#include "stdio.h"
#include "complex.h"

unsigned char cal_pixel(double complex d, double b, unsigned char N)
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
  int N = 256;
  int b = 2;
  int w = 2048;
  int h = 2048;
  int dx = 2 * b / (w - 1);
  int dy = 2 * b / (h - 1);
  int x, y;
  int color[w][h];
  // p = 2; q = 3;
  int Q = 1;
  int P = 2;
  int p = 2;
  int wp = w / P;
  int hp = h;
  int xoff = p * w / P;
  int yoff = 0;

  for (x = 0; x < wp - 1; x++)
  {
    float dreal = (x + xoff) * dx - b;
    for (y = 0; y < h - 1; y++)
    {
      float dimag = (y + yoff) * dy - b;
      complex double d = dreal + I * dimag;
      color[x][y] = cal_pixel(d, b, N);
    }
  }
  return 0;
}
