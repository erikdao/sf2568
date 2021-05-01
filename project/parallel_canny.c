#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <mpi.h>
 
#define MAX_BRIGHTNESS 255

typedef struct {
    uint8_t magic[2];
} bmpfile_magic_t;
 
typedef struct {
    uint32_t filesz;
    uint16_t creator1;
    uint16_t creator2;
    uint32_t bmp_offset;
} bmpfile_header_t;
 
typedef struct {
    uint32_t header_sz;
    int32_t  width;
    int32_t  height;
    uint16_t nplanes;
    uint16_t bitspp;
    uint32_t compress_type;
    uint32_t bmp_bytesz;
    int32_t  hres;
    int32_t  vres;
    uint32_t ncolors;
    uint32_t nimpcolors;
} bitmap_info_header_t;
 
typedef struct {
    uint8_t r;
    uint8_t g;
    uint8_t b;
    uint8_t nothing;
} rgb_t;

int *read_bmp(const char *fname, bitmap_info_header_t *bmpInfoHeader, int num_procs);

bool save_bmp(const char *fname, const bitmap_info_header_t *bmpInfoHeader, const int *data);

int *canny_edge_detection(const int *in, const int width, int const height,
                          const int tmin, const int tmax, const float sigma);

void convolution(const int *in, int *out, const float *kernel,
                 const int nx, const int ny, const int kn,
                 const bool normalize);

void gaussian_filter(const int *in, int *out,
                     const int nx, const int ny, const float sigma);

int main(int argc, char *argv[]) {
    int size, rank, tag = 100;
    int num_pixels; // Number of pixels of the original image
    int my_count; // Number of pixels for each sub image

    int *in_image = NULL;  // Input image
    // int *recv_buf = NULL;  // Buffer image
    int im_width, im_height;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;

    if (rank == 0) {
        static bitmap_info_header_t ih;
        in_image = read_bmp(argv[1], &ih, size);
        if (in_image == NULL) {
            fprintf(stderr, "Main process: error while reading BMP\n");
            return -1;
        }
        im_width = ih.width;
        im_height = ih.height;
        num_pixels = im_width * im_height;
        my_count = num_pixels / size; 
        if (num_pixels % size != 0) {
            my_count += (num_pixels % size);
        }
        fprintf(stdout, "Master %d: my_count=%d, width=%d, height=%d\n", rank, my_count, im_width, im_height);
        for (int i = 0; i < 10; i++) {
            fprintf(stdout, "%d ", in_image[i]);
        }
        fprintf(stdout, "\n");
    }

    // Broadcast my_count to all processes
    MPI_Bcast(&my_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&im_width, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&im_height, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int *recv_buf = (int*) malloc(my_count * sizeof(int));

    MPI_Scatter(in_image, my_count, MPI_INT, recv_buf, my_count, MPI_INT, 0, MPI_COMM_WORLD);

    int sub_image_height = my_count / im_width; // Height of the sub-image
    int *edge_image = canny_edge_detection(recv_buf, im_width, sub_image_height, 45, 50, 1.0f);

    // Sequential write subimage
    int signal = 0;
    char *fname = "out.txt";
    // sprintf(fname, "outputs/out.", argv[1]);
    FILE *f; 

    if (rank == 0) {
        f = fopen(fname, "w");
        int count = 0;
        for (int h = 0; h < sub_image_height; h++) {
            for (int w = 0; w < im_width; w++) {
                fprintf(f, "%d ", edge_image[count++]);
            }
            fprintf(f, "\n");
        }
        fclose(f);
        fprintf(stdout, "Processor %d finished writing its subimage edges\n", rank);

        signal = 1;
        MPI_Send(&signal, 1, MPI_INT, rank+1, 100, MPI_COMM_WORLD);
    } else {
        MPI_Recv(&signal, 1, MPI_INT, rank-1, 100, MPI_COMM_WORLD, &status);

        if (signal == 1) {
            f = fopen(fname, "a");
            int count = 0;
            for (int h = 0; h < sub_image_height; h++) {
                for (int w = 0; w < im_width; w++) {
                    fprintf(f, "%d ", edge_image[count++]);
                }
                fprintf(f, "\n");
            }
            fclose(f);
            if (rank != size - 1) {
                MPI_Send(&signal, 1, MPI_INT, rank+1, 100, MPI_COMM_WORLD);
                fprintf(stdout, "Processor %d finished writing its subedges, send signal to %d\n", rank, rank+1);
            }
        }
    }

    // MPI_Gather(in_image, my_count, MPI_INT, recv_buf, my_count, MPI_INT, 0, MPI_COMM_WORLD);
    // free(edge_image);
    // free(in_image);
    // free(recv_buf);

    MPI_Finalize();
    return 0;
}

/**
 * Read BMP image
 * @param fname: 
 */
int *read_bmp(const char *fname, bitmap_info_header_t *bmpInfoHeader, int num_procs) {
    FILE *f = fopen(fname, "rb");
    if (f == NULL) {
        printf("Error while reading file");
        return NULL;
    }

    bmpfile_magic_t mag;
    if (fread(&mag, sizeof(bmpfile_magic_t), 1, f) != 1) {
        fclose(f);
        return NULL;
    }

    bmpfile_header_t bitmapFileHeader; // our bitmap file header
    // read the bitmap file header
    if (fread(&bitmapFileHeader, sizeof(bmpfile_header_t),
              1, f) != 1) {
        fclose(f);
        return NULL;
    }

    // read the bitmap info header
    if (fread(bmpInfoHeader, sizeof(bitmap_info_header_t),
              1, f) != 1) {
        fclose(f);
        return NULL;
    }

    // move file point to the beginning of bitmap data
    if (fseek(f, bitmapFileHeader.bmp_offset, SEEK_SET)) {
        fclose(f);
        return NULL;
    }
 
    int32_t original_height = bmpInfoHeader-> height;
    int added_lines = 0;
    if (original_height % num_procs != 0) {
        added_lines = num_procs - (original_height % num_procs);
    }

    int32_t new_height = original_height + added_lines;
    int32_t new_width = bmpInfoHeader-> width;

    // allocate enough memory for the bitmap image data
    // int *bitmapImage = malloc(bmpInfoHeader->bmp_bytesz * sizeof(int));
    int *bitmapImage = malloc(new_height * new_width * sizeof(int));

    int count = 0;
    // Pad the first line with 0
    for (size_t w = 0; w < new_width; w++) {
        bitmapImage[count++] = 0;
    }
    // read in the bitmap image data
    size_t pad;
    unsigned char c;
    pad = 4*ceil(bmpInfoHeader->bitspp*bmpInfoHeader->width/32.) - bmpInfoHeader->width;
    for(size_t i = 0; i < new_height; i++ ) {  // i<bmpInfoHeader->height; i++){
	    for(size_t j=0; j < new_width; j++) { // <bmpInfoHeader->width; j++){
		    if (fread(&c, sizeof(unsigned char), 1, f) != 1) {
			    fclose(f);
			    return NULL;
		    }
		    bitmapImage[count++] = (int) c;
	    }
	    fseek(f, pad, SEEK_CUR);
    }

    // Pad the rest of the height with 0
    for (size_t i = 0; i < added_lines - 1; i++) {
        for (size_t j = 0; j < new_width; j++) {
            bitmapImage[count++] = 0;
        }
    }

    fclose(f);
    return bitmapImage;
}

void convolution(const int *in, int *out, const float *kernel,
                 const int nx, const int ny, const int kn,
                 const bool normalize)
{
    assert(kn % 2 == 1);
    assert(nx > kn && ny > kn);
    const int khalf = kn / 2;
    float min = FLT_MAX, max = -FLT_MAX;
 
    if (normalize)
        for (int m = khalf; m < nx - khalf; m++)
            for (int n = khalf; n < ny - khalf; n++) {
                float pixel = 0.0;
                size_t c = 0;
                for (int j = -khalf; j <= khalf; j++)
                    for (int i = -khalf; i <= khalf; i++) {
                        pixel += in[(n - j) * nx + m - i] * kernel[c];
                        c++;
                    }
                if (pixel < min)
                    min = pixel;
                if (pixel > max)
                    max = pixel;
                }
 
    for (int m = khalf; m < nx - khalf; m++)
        for (int n = khalf; n < ny - khalf; n++) {
            float pixel = 0.0;
            size_t c = 0;
            for (int j = -khalf; j <= khalf; j++)
                for (int i = -khalf; i <= khalf; i++) {
                    pixel += in[(n - j) * nx + m - i] * kernel[c];
                    c++;
                }
 
            if (normalize)
                pixel = MAX_BRIGHTNESS * (pixel - min) / (max - min);
            out[n * nx + m] = (int)pixel;
        }
}

void gaussian_filter(const int *in, int *out,
                     const int nx, const int ny, const float sigma)
{
    const int n = 2 * (int)(2 * sigma) + 3;
    const float mean = (float)floor(n / 2.0);
    float kernel[n * n]; // variable length array
 
    size_t c = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            kernel[c] = exp(-0.5 * (pow((i - mean) / sigma, 2.0) +
                                    pow((j - mean) / sigma, 2.0)))
                        / (2 * M_PI * sigma * sigma);
            c++;
        }
 
    convolution(in, out, kernel, nx, ny, n, true);
}

/**
 * Canny algorithm for edge detection
 * @param in: input image
 * @param width: width of the input image
 * @param height: height of the input image
 * @param tmin: Minimum threshold
 */
int *canny_edge_detection(const int *in, const int width, int const height,
                          const int tmin, const int tmax, const float sigma)
{
    const int nx = width;
    const int ny = height;
 
    int *G = calloc(nx * ny * sizeof(int), 1);  // Gradient
    int *after_Gx = calloc(nx * ny * sizeof(int), 1);  // Gradient in x direction
    int *after_Gy = calloc(nx * ny * sizeof(int), 1);  // Gradient in y direction
    int *nms = calloc(nx * ny * sizeof(int), 1);
    int *out = malloc(nx * ny * sizeof(int));  // Output edge image
 
    if (G == NULL || after_Gx == NULL || after_Gy == NULL ||
        nms == NULL || out == NULL) {
        fprintf(stderr, "canny_edge_detection:"
                " Failed memory allocation(s).\n");
        exit(1);
    }
 
    gaussian_filter(in, out, nx, ny, sigma);
 
    const float Gx[] = {-1, 0, 1,
                        -2, 0, 2,
                        -1, 0, 1};
 
    convolution(out, after_Gx, Gx, nx, ny, 3, false);
 
    const float Gy[] = { 1, 2, 1,
                         0, 0, 0,
                        -1,-2,-1};
 
    convolution(out, after_Gy, Gy, nx, ny, 3, false);
 
    for (int i = 1; i < nx - 1; i++)
        for (int j = 1; j < ny - 1; j++) {
            const int c = i + nx * j;
            // G[c] = abs(after_Gx[c]) + abs(after_Gy[c]);
            G[c] = (int)hypot(after_Gx[c], after_Gy[c]);
        }
 
    // Non-maximum suppression, straightforward implementation.
    for (int i = 1; i < nx - 1; i++)
        for (int j = 1; j < ny - 1; j++) {
            const int c = i + nx * j;
            const int nn = c - nx;
            const int ss = c + nx;
            const int ww = c + 1;
            const int ee = c - 1;
            const int nw = nn + 1;
            const int ne = nn - 1;
            const int sw = ss + 1;
            const int se = ss - 1;
 
            const float dir = (float)(fmod(atan2(after_Gy[c],
                                                 after_Gx[c]) + M_PI,
                                           M_PI) / M_PI) * 8;
 
            if (((dir <= 1 || dir > 7) && G[c] > G[ee] &&
                 G[c] > G[ww]) || // 0 deg
                ((dir > 1 && dir <= 3) && G[c] > G[nw] &&
                 G[c] > G[se]) || // 45 deg
                ((dir > 3 && dir <= 5) && G[c] > G[nn] &&
                 G[c] > G[ss]) || // 90 deg
                ((dir > 5 && dir <= 7) && G[c] > G[ne] &&
                 G[c] > G[sw]))   // 135 deg
                nms[c] = G[c];
            else
                nms[c] = 0;
        }
 
    // Reuse array
    // used as a stack. nx*ny/2 elements should be enough.
    int *edges = (int*) after_Gy;
    memset(out, 0, sizeof(int) * nx * ny);
    memset(edges, 0, sizeof(int) * nx * ny);
 
    // Tracing edges with hysteresis . Non-recursive implementation.
    size_t c = 1;
    for (int j = 1; j < ny - 1; j++)
        for (int i = 1; i < nx - 1; i++) {
            if (nms[c] >= tmax && out[c] == 0) { // trace edges
                out[c] = MAX_BRIGHTNESS;
                int nedges = 1;
                edges[0] = c;
 
                do {
                    nedges--;
                    const int t = edges[nedges];
 
                    int nbs[8]; // neighbours
                    nbs[0] = t - nx;     // nn
                    nbs[1] = t + nx;     // ss
                    nbs[2] = t + 1;      // ww
                    nbs[3] = t - 1;      // ee
                    nbs[4] = nbs[0] + 1; // nw
                    nbs[5] = nbs[0] - 1; // ne
                    nbs[6] = nbs[1] + 1; // sw
                    nbs[7] = nbs[1] - 1; // se
 
                    for (int k = 0; k < 8; k++)
                        if (nms[nbs[k]] >= tmin && out[nbs[k]] == 0) {
                            out[nbs[k]] = MAX_BRIGHTNESS;
                            edges[nedges] = nbs[k];
                            nedges++;
                        }
                } while (nedges > 0);
            }
            c++;
        }
 
    free(after_Gx);
    free(after_Gy);
    free(G);
    free(nms);
 
    return out;
}
