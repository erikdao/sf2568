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
#define KERNEL_SIZE 5

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

int *read_bmp(const char *fname, bitmap_info_header_t *bmpInfoHeader);

bool save_bmp(const char *fname, const bitmap_info_header_t *bmpInfoHeader, const int *data);

int *canny_edge_detection(const int *in, const int width, int const height,
                          const int tmin, const int tmax, const float sigma);

void convolution(const int *in, int *out, const float *kernel,
                 const int nx, const int ny, const int kn,
                 const bool normalize);

void gaussian_filter(const int *in, int *out,
                     const int nx, const int ny, const float sigma);

int *pad_image(const int *orig_image, int pad_one, int pad_two, int width, int height);

int main(int argc, char *argv[]) {
    int size, rank, tag = 100;

    int *image = NULL;
    int *sub_image = NULL;
    int sub_image_height;

    int orig_width, orig_height, im_width, im_height;
    int my_count; // How many element I will receive from master
    int my_sub_image_height; // The height of my sub image
    int startIndex, endIndex; // My local start and end indices

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
    pad_one = (int) (KERNEL_SIZE - 1)/2;

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        int *original_image = NULL; // Input image

        static bitmap_info_header_t ih;
        original_image = read_bmp(argv[1], &ih);
        if (original_image == NULL) {
            fprintf(stderr, "Main process: error while reading BMP\n");
            return -1;
        }
        orig_width = ih.width;
        orig_height = ih.height;
        fprintf(stdout, "Original height=%d, width=%d\n", orig_height, orig_width);
        im_width = orig_width;
        // Padding images
        im_height = orig_height + 2 * pad_one;
        if (im_height % size != 0) {
            pad_two = size - (im_height % size);
            im_height += pad_two;
        }

        image = pad_image(original_image, pad_one, pad_two, orig_width, orig_height);

        // Sanity check
        int non_zero = 0;
        int count = 0;
        for (int i = 0; i < im_height; i++) {
            for (int j = 0; j < im_width; j++) {
                if (image[count++] != 0) {
                    non_zero += 1;
                }
            }
        }

        fprintf(stdout, "Processor: %d, after padding non-zero elements: %d\n", rank, non_zero);
        // Broadcast sub images to slave processor
        MPI_Bcast(&im_height, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&im_width, 1, MPI_INT, 0, MPI_COMM_WORLD);
    } else {
        MPI_Bcast(&im_height, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&im_width, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    MPI_Barrier(MPI_COMM_WORLD);

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

    fprintf(stdout, "Processor: %d, im_height=%d, im_width=%d, pad_one=%d, my_count=%d, sub_image=%d\n",
                    rank, im_height, im_width, pad_one, my_count, my_sub_image_height);

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
            // printf("p_i = %d; end-start=%d, sub_image=%d, sub_image * im_width=%d, my_count = %d\n",
            //             p_i, endIndex - startIndex, my_sub_image_height * im_width, my_count);
            int count = 0;
            printf("p_i = %d; end-start=%d, end=%d, start=%d\n", p_i, pEndIndex - pStartIndex, pEndIndex, pStartIndex);
            for (int i = pStartIndex; i < pEndIndex - 1; i++) {
                for (int j = 0; j < im_width; j++) {
                    sub_image[count] = image[i * im_width + j];
                    count++;
                }
            }

            // Sanity check
            int non_zero = 0;
            count = 0;
            for (int i = 0; i < my_sub_image_height; i++) {
                for (int j = 0; j < im_width; j++) {
                    if (sub_image[count] != 0) {
                        non_zero += 1;
                    }
                    count++;
                }
            }

            fprintf(stdout, "Processor: %d, sub image before sending from master. non-zero elements: %d\n", rank, non_zero);

            MPI_Send(sub_image, my_count, MPI_INT, p_i, tag, MPI_COMM_WORLD);
        }

        free(sub_image);
        sub_image = malloc(my_count * sizeof(int));
        // For master processor, don't need to send its partition, just copy
        startIndex = 0;
        endIndex = im_height / size + pad_one;
        printf("startIndex %d, endIndex: %d\n", startIndex, endIndex);
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

    // Sanity check
    int non_zero = 0;
    int count = 0;
    for (int i = 0; i < my_sub_image_height; i++) {
        for (int j = 0; j < im_width; j++) {
            if (sub_image[count++] != 0) {
                non_zero += 1;
            }
        }
    }

    fprintf(stdout, "Processor: %d, non-zero elements: %d\n", rank, non_zero);

    int *edge_image = canny_edge_detection(sub_image, im_width, my_sub_image_height, 45, 50, 1.0f);

    // Sequential write subimage
    int signal = 0;
    char *fname = "out.txt";
    // sprintf(fname, "outputs/out.", argv[1]);
    FILE *f; 

    if (rank == 0) {
        f = fopen(fname, "w");
        int count = 0;
        for (int h = pad_one; h < my_sub_image_height - 2 * pad_one; h++) {
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
            // Skip padded-top rows
            for (int i = 0; i < pad_one; i++) {
                for (int j = 0; j < im_width; j++) {
                    count++;
                }
            }
            if (rank > 0 && rank < size - 1) {
                for (int h = pad_one; h < my_sub_image_height - 2 * pad_one; h++) {
                    for (int w = 0; w < im_width; w++) {
                        fprintf(f, "%d ", edge_image[count++]);
                    }
                    fprintf(f, "\n");
                }
            } else {
                for (int h = pad_one; h < my_sub_image_height; h++) {
                    for (int w = 0; w < im_width; w++) {
                        fprintf(f, "%d ", edge_image[count++]);
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

    MPI_Finalize();
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
    // int *bitmapImage = malloc(new_height * new_width * sizeof(int));

    // read in the bitmap image data
    size_t pad, count = 0;
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
 
    // const float Gx[] = {-1, 0, 1,
    //                     -2, 0, 2,
    //                     -1, 0, 1};

    // Sobel filtering
    // const float Gx[] = {1, 0, -1,
    //                     2, 0, -2,
    //                     1, 0, -1};
    const float Gx[] = {2, 2, 4, 2, 2,
                        1, 1, 2, 1, 1,
                        0, 0, 0, 0, 0,
                        -1, -1, -2, -1, -1,
                        -2, -2, -4, -2, -2};

    convolution(out, after_Gx, Gx, nx, ny, 3, false);
 
    // const float Gy[] = { 1, 2, 1,
    //                      0, 0, 0,
    //                     -1,-2,-1};
    const float Gy[] = { 2, 1, 0, -1, -2,
                         2, 1, 0, -1, -2,
                         4, 2, 0, -2, -4,
                         2, 1, 0, -1, -2,
                         2, 1, 0, -1, -2};
 
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
