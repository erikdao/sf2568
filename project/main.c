/**
 * Parallel Edge Detection for High-Resolution Image
 * Course Project - SF2568
 * Author: Cuong Dao, Donggyun Park
 */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <jpeglib.h>

#define BYTES_PER_PIXEL = 3;  // Colored images
#define COLOR_SPACE = JCS_RGG;

// Array to hold image data
unsigned char *raw_image = NULL;

int read_jpege_file(char *fname);

int main(char *argc, char **argv) {
    // Main program
}

/**
 * Read an JPEG image file into a big 1D array
 */
int read_jpeg_file(char *fname) {
    struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;

    JSAMPROW row_pointer[1];

    FILE *infile = fopen(fname, "rb");

    unsigned long location = 0;
    int i = 0;

    if (!infile) {
        printf("Error while attempting to open JPEG file %s\n", fname);
        return -1;
    }

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&cinfo);
    jpeg_stdio_src(&cinfo, infile);
    jpeg_read_header(&cinfo, TRUE);

    jpeg_start_decompress(&cinfo);

	jpeg_start_decompress(&cinfo);

    // Allocate the image data pointer, the image is flatten into a 1D array
    // TODO: Consider using a separate 2D array for each of the R, G, B channels
    raw_image = (unsigned char*) malloc(cinfo.output_width * cinfo.output_height * cinfo.num_components);
    row_pointer[0] = (unsigned char *) malloc(cinfo.output_width * cinfo.num_components);

    while(cinfo.output_scanline < cinfo.image_height) {
        jpeg_read_scanlines(&cinfo, row_pointer, 1);
        for (i = 0; i < cinfo.image_width * cinfo.num_components; i++) 
            raw_image[location++] = row_pointer[0][i];
    }

    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    free(row_pointer[0]);
    fclose(infile);

    return 1;
}