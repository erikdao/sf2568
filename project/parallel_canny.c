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

// Define a type for image pixel values
typedef short int pixel_t;

pixel_t *read_bmp(const char *fname, bitmap_info_header_t *bmpInfoHeader);

bool save_bmp(const char *fname, const bitmap_info_header_t *bmpInfoHeader, const pixel_t *data);

int main(int argc, char *argv[]) {
    int size, rank,send_count, tag = 100;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;
    int width = 3360;
    int height = 5040;
    send_count = (width * height) / size;
    pixel_t *recv_data = malloc(send_count * sizeof(pixel_t));

    if (rank == 0) {
        static bitmap_info_header_t ih;
        pixel_t *in_image = read_bmp(argv[1], &ih);
        if (in_image == NULL) {
            fprintf(stderr, "Main process: error while reading BMP\n");
            return -1;
        }
        fprintf(stdout, "BMP image read!\n");
        fprintf(stdout, "Image size: w=%d, h=%d\n", ih.width, ih.height);
        
        //send_count = (ih.width * ih.height) / size;
        recv_data = malloc(send_count * sizeof(pixel_t));
        MPI_Scatter(in_image, send_count ,MPI_INT, recv_data, send_count, MPI_INT, 0, MPI_COMM_WORLD);


    } else {
        fprintf(stdout, "Processor: %d - No read\n", rank);

    }
    
    MPI_Finalize();
    return 0;
}

pixel_t *read_bmp(const char *fname, bitmap_info_header_t *bmpInfoHeader) {
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
    pixel_t *bitmapImage = malloc(bmpInfoHeader->bmp_bytesz *
                                  sizeof(pixel_t));

    // read in the bitmap image data
    size_t pad, count=0;
    unsigned char c;
    pad = 4*ceil(bmpInfoHeader->bitspp*bmpInfoHeader->width/32.) - bmpInfoHeader->width;
    for(size_t i=0; i<bmpInfoHeader->height; i++){
	    for(size_t j=0; j<bmpInfoHeader->width; j++){
		    if (fread(&c, sizeof(unsigned char), 1, f) != 1) {
			    fclose(f);
			    return NULL;
		    }
		    bitmapImage[count++] = (pixel_t) c;
	    }
	    fseek(f, pad, SEEK_CUR);
    }

    fclose(f);
    return bitmapImage;
}
