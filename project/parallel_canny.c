/**
 * Paralleled Edge Detection for High-Resolution Images
 * Course Project for SF2568 - Parallel Computation for Large-Scale Problems
 * 
 * Cuong Dao (cuongdd@kth.se), Donggyun Park (donggyun@kth.se) 
 * 
 **/
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdbool.h>
#include <mpi.h>
#include <time.h>
 
#define MAX_BRIGHTNESS 255  // "White" pixel of the image

/*
 * Loading part taken from
 * http://www.vbforums.com/showthread.php?t=261522
 * BMP info:
 * http://en.wikipedia.org/wiki/BMP_file_format
 *
 * Note: the magic number has been removed from the bmpfile_header_t
 * structure since it causes alignment problems
 *     bmpfile_magic_t should be written/read first
 * followed by the
 *     bmpfile_header_t
 * [this avoids compiler-specific alignment pragmas etc.]
 */

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

/**
 * Function prototypes 
 */
int *read_bmp(const char *fname, bitmap_info_header_t *bmpInfoHeader);

bool save_bmp(const char *fname, const bitmap_info_header_t *bmpInfoHeader, const int *data);

int *canny_edge_detection(const int *in, const int kernel, const int width,
                          const int height, const int tmin, const int tmax,
                          const float sigma);

void convolution(const int *in, int *out, const float *kernel,
                 const int nx, const int ny, const int kn,
                 const bool normalize);

void gaussian_filter(const int *in, int *out,
                     const int nx, const int ny, const float sigma);

int *pad_image(const int *orig_image, int pad_one, int pad_two, int width, int height);

void write_output(const int rank, const int size, const char *fname, const int *image, int with, int height, int pad_one);

/**
 * Main entry to the program.
 * 
 * The command to run this program is
 * mpirun -np NUM_PROCS parallelcanny K input_image output_file
 * @param NUM_PROCS: number of processor
 * @param K: Convolution kernal size (3 or 5)
 * @param input_image: Path to the input BMP image
 * @param output_file: Path to save the output result 
 **/
int main(int argc, char *argv[]) {
    if (argc != 4) {
        fprintf(stderr, "Invalid parameters, run the program with: \n"
                        "parallelcanny K input_image output_file\n");
        exit(1);
    }
    int size, rank, tag = 100;

    int *image = NULL;
    int *sub_image = NULL;
    int sub_image_height;

    int orig_width, orig_height, im_width, im_height;
    int my_count; // How many element I will receive from master
    int my_sub_image_height; // The height of my sub image
    int startIndex, endIndex; // My local start and end indices

    int kernel = atoi(argv[1]); // Kernel size for the convolution

    /** We'll first pad the top and the bottom of the height with
     * `pad_one` lines of zero.
     * Then we'll pad the bottom of the image with `pad_two` lines 
     * of zeros so the image's height becomes divisble by `size
     */
    int pad_one = 0, pad_two = 0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;
    pad_one = (int) (kernel - 1)/2;

    /**
     * Read the input image, do some pre-processing such as blurring
     * and padding on the master processor before sending out sub-images to slaves
     **/
    if (rank == 0) {
        int *original_image = NULL; // Input image

        static bitmap_info_header_t ih;
        original_image = read_bmp(argv[2], &ih);
        if (original_image == NULL) {
            fprintf(stderr, "Main process: error while reading BMP\n");
            return -1;
        }
        orig_width = ih.width;
        orig_height = ih.height;
        fprintf(stdout, "Original height=%d, width=%d\n", orig_height, orig_width);

        fprintf(stdout, "Applying Gaussian filter to denoise image\n");
        // Apply Gaussian filter to denoise the image
        int *blurred_image = malloc(ih.width * ih.height * sizeof(int));
        gaussian_filter(original_image, blurred_image, ih.width, ih.height, 1.0f);
        im_width = orig_width;
        im_height = orig_height;
        im_width = orig_width;
        // Padding images
        im_height = orig_height + 2 * pad_one;
        if (im_height % size != 0) {
            pad_two = size - (im_height % size);
            im_height += pad_two;
        }

        fprintf(stdout, "Padding the images with zeros\n");
        image = pad_image(blurred_image, pad_one, pad_two, orig_width, orig_height);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    clock_t start = clock();
    MPI_Bcast(&im_height, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&im_width, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Calculate the local indices on each processor
    if (rank == 0) {
        startIndex = 0;
        endIndex = im_height / size + pad_one;
    } else if (rank > 0 && rank < size - 1) {
        startIndex = rank * im_height / size - pad_one;
        endIndex = (rank+1) * im_height / size + pad_one;
    } else { // Last partition
        startIndex = rank * im_height / size - pad_one;
        endIndex = (rank+1) * im_height / size;
    }
    my_sub_image_height = endIndex - startIndex;
    my_count = my_sub_image_height * im_width;

    // Allocate memory for the sub image to be sent
    sub_image = malloc(2 * my_count * sizeof(int));
    if (sub_image == NULL) {
        fprintf(stderr, "Error while allocating sub_image memory!");
    }

    if (rank == 0) {
        // Send subimage to each slave processor
        for (int p_i = 1; p_i < size; p_i++) {
            int pStartIndex, pEndIndex;
            // Copy the portion of the sub image from the `image`                
            if (p_i > 0 && p_i < size - 1) {
                pStartIndex = p_i * im_height / size - pad_one;
                pEndIndex = (p_i+1) * im_height / size + pad_one;
            } else { // Last partition
                pStartIndex = p_i * im_height / size - pad_one;
                pEndIndex = (p_i+1) * im_height / size;
            }
            int count = 0;;
            for (int i = pStartIndex; i < pEndIndex - 1; i++) {
                for (int j = 0; j < im_width; j++) {
                    sub_image[count] = image[i * im_width + j];
                    count++;
                }
            }

            MPI_Send(sub_image, my_count, MPI_INT, p_i, tag, MPI_COMM_WORLD);
        }

        free(sub_image);
        sub_image = malloc(my_count * sizeof(int));
        // For master processor, don't need to send its partition, just copy
        int count = 0;
        for (int i = startIndex; i < endIndex; i++) {
            for (int j = 0; j < im_width; j++) {
                sub_image[count] = image[i * im_width + j];
                count++;
            }
        }
    } else { // slave processes
        MPI_Recv(sub_image, my_count, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);
    }

    clock_t start_canny = clock();
    int *edge_image = canny_edge_detection(sub_image, kernel, im_width, my_sub_image_height, 45, 50, 1.0f);
    clock_t end_canny = clock();
    double edge_time = ((double) (end_canny - start_canny)) / CLOCKS_PER_SEC;
    fprintf(stdout, "Processor %d, Canny algorithm took %f seconds\n", rank, edge_time);

    write_output(rank, size, argv[3], edge_image, im_width, my_sub_image_height, pad_one);

    free(image);
    free(sub_image);
    MPI_Finalize();

    if (rank == size - 1) {
        clock_t final = clock();
        double total_time = ((double) (final - start)) / CLOCKS_PER_SEC;
        fprintf(stdout, "Processor %d, Total time %f seconds\n", rank, total_time);
    }
    return 0;
}

/**
 * Read BMP image
 * @param fname: 
 */
int *read_bmp(const char *fname, bitmap_info_header_t *bmpInfoHeader) {
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

    // allocate enough memory for the bitmap image data
    int *bitmapImage = malloc(bmpInfoHeader->bmp_bytesz * sizeof(int));

    // read in the bitmap image data
    unsigned int pad, count = 0;
    unsigned char c;
    pad = 4*ceil(bmpInfoHeader->bitspp*bmpInfoHeader->width/32.) - bmpInfoHeader->width;
    for(size_t i = 0; i < bmpInfoHeader->height; i++ ) {
	    for(size_t j=0; j < bmpInfoHeader->width; j++) {
		    if (fread(&c, sizeof(unsigned char), 1, f) != 1) {
			    fclose(f);
			    return NULL;
		    }
		    bitmapImage[count++] = (int) c;
	    }
	    fseek(f, pad, SEEK_CUR);
    }

    fclose(f);
    return bitmapImage;
}

/**
 * Pad the top and bottom rows of the image with zero so that the height
 * become divisible by the number of processor
 */
int *pad_image(const int *orig_image, int pad_one, int pad_two, int width, int height) {
    // Allocate memory for new image
    int *image = malloc(width * (height + 2 * pad_one + pad_two) * sizeof(int));

    // Pad the top lines
    int count = 0;
    for (int i = 0; i < pad_one; i++) {
        for (int j = 0; j < width; j++) {
            image[count++] = 0;
        }
    }

    // Copy original image to the new image
    int iterator = 0;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            image[count + iterator] = orig_image[iterator];
            iterator++;
        }
    }

    // Pad the bottom lines
    count += iterator;
    for (int i = 0; i < (pad_one + pad_two); i++) {
        for (int j = 0; j < width; j++) {
            image[count++] = 0;
        }
    }
    return image;
}

/**
 * Perform convolution operation
 *
 * @param in: input image
 * @param out: output image
 * @param kernel: the convolution kernel, a flatten 1D array
 * @param nx: width of the input image
 * @param ny: height of the input image
 * @param kn: kernel size
 **/ 
void convolution(const int *in, int *out, const float *kernel,
                 const int nx, const int ny, const int kn,
                 const bool normalize)
{
    const int khalf = kn / 2;
    float min = FLT_MAX, max = -FLT_MAX;
 
    if (normalize)
        for (int m = khalf; m < nx - khalf; m++)
            for (int n = khalf; n < ny - khalf; n++) {
                float pixel = 0.0;
                unsigned int c = 0;
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
 
    // Actual implementation of the convolution operation, i.e.
    // element-wise multiplication of the kernel and the image patch
    for (int m = khalf; m < nx - khalf; m++)
        for (int n = khalf; n < ny - khalf; n++) {
            float pixel = 0.0;
            unsigned int c = 0;
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

/**
 * Gaussian filter
 * @param in: input image
 * @param out: output image
 * @param nx: width of the input image
 * @param ny: height of the input image
 * @param sigma: the standard deviation of the gaussian distribution
 **/
void gaussian_filter(const int *in, int *out,
                     const int nx, const int ny, const float sigma)
{
    const int n = 2 * (int)(2 * sigma) + 3;
    const float mean = (float)floor(n / 2.0);
    float kernel[n * n]; // variable length array
 
    unsigned int c = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            kernel[c++] = exp(-0.5 * (pow((i - mean) / sigma, 2.0) + pow((j - mean) / sigma, 2.0))) / (2 * M_PI * sigma * sigma);
        }
 
    convolution(in, out, kernel, nx, ny, n, true);
}

/*
 * Reference:
 * http://en.wikipedia.org/wiki/Canny_edge_detector
 * http://www.tomgibara.com/computer-vision/CannyEdgeDetector.java
 * http://fourier.eng.hmc.edu/e161/lectures/canny/node1.html
 * http://www.songho.ca/dsp/cannyedge/cannyedge.html
 */

/**
 * Canny algorithm for edge detection
 * @param in: input image
 * @param width: width of the input image
 * @param height: height of the input image
 * @param tmin: Minimum threshold
 */
int *canny_edge_detection(const int *in, const int kernel, const int width, int const height, const int tmin, const int tmax, const float sigma)
{
    const int nx = width;
    const int ny = height;
 
    int *G = calloc(nx * ny * sizeof(int), 1);  // Gradient
    int *after_Gx = calloc(nx * ny * sizeof(int), 1);  // Gradient in x direction
    int *after_Gy = calloc(nx * ny * sizeof(int), 1);  // Gradient in y direction
    int *nms = calloc(nx * ny * sizeof(int), 1); // Output of non-maximum suppression
    int *out = malloc(nx * ny * sizeof(int));  // Output edge image
 
    // Sobel operators
    float Gx3[] = {1, 0, -1,
                  2, 0, -2,
                  1, 0, -1};
    float Gx5[] = { 2, 1, 0, -1, -2,
        2, 1, 0, -1, -2,
        4, 2, 0, -2, -4,
        2, 1, 0, -1, -2,
        2, 1, 0, -1, -2};
    float *Gx = Gx3;

    if (kernel == 5) {
        Gx = Gx5;
    }
    convolution(in, after_Gx, Gx, nx, ny, kernel, false);
 
    float Gy3[] = { 1, 2, 1,
                   0, 0, 0,
                   -1,-2,-1};
    float Gy5[] = {2,  2,  4,  2,  2,
              1,  1,  2,  1,  1,
              0,  0,  0,  0,  0,
             -1, -1, -2, -1, -1,
             -2, -2, -4, -2, -2};
    
    float *Gy = Gy3;
    if (kernel == 5) {
        Gy = Gy5;
    }
    convolution(in, after_Gy, Gy, nx, ny, kernel, false);
 
    for (int i = 1; i < nx - 1; i++)
        for (int j = 1; j < ny - 1; j++) {
            const int c = i + nx * j;
            G[c] = (int)hypot(after_Gx[c], after_Gy[c]);
        }
 
    // Non-maximum suppression
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
 
            const float dir = (float)(fmod(atan2(after_Gy[c], after_Gx[c]) + M_PI, M_PI) / M_PI) * 8;
 
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
 
    // Tracing edges with hysteresis
    unsigned int c = 1;
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

void write_output(const int rank, const int size, const char *fname, const int *image, int width, int height, int pad_one) {
    MPI_Status status;
    int signal = 0;
    FILE *f; 

    if (rank == 0) {
        f = fopen(fname, "w");
        int count = 0;
        for (int h = pad_one; h < height - 2 * pad_one; h++) {
            for (int w = 0; w < width; w++) {
                fprintf(f, "%d ", image[count++]);
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
            // Skip padded-top rows
            for (int i = 0; i < pad_one; i++) {
                for (int j = 0; j < width; j++) {
                    count++;
                }
            }
            if (rank > 0 && rank < size - 1) {
                for (int h = pad_one; h < height - 2 * pad_one; h++) {
                    for (int w = 0; w < width; w++) {
                        fprintf(f, "%d ", image[count++]);
                    }
                    fprintf(f, "\n");
                }
            } else {
                for (int h = pad_one; h < height; h++) {
                    for (int w = 0; w < width; w++) {
                        fprintf(f, "%d ", image[count++]);
                    }
                    fprintf(f, "\n");
                }
            }
            fclose(f);
            if (rank != size - 1) {
                MPI_Send(&signal, 1, MPI_INT, rank+1, 100, MPI_COMM_WORLD);
                fprintf(stdout, "Processor %d finished writing its subedges, send signal to %d\n", rank, rank+1);
            }
        }
    }
}