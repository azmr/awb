/* AWB - A single-header library for colour-balancing photos */

/* Feature List / TODO:
 *  Illuminant Estimation Methods
 *  [ ] Gray World Theory
 *  	[ ] Vanilla
 *  	[ ] + Retinex
 *  	[ ] + Standard Deviation Weighting
 *  	[ ] + Standard Deviation & Luminance Weighting
 *  	[ ] Robust (Huo et al., 2005)
 *  [ ] Simplest Color Balance (Limare et al., 2011)
 *  [ ] Retinex/Perfect Reflector
 *  [ ] Sensor Correlation (?)
 *  [ ] ...
 *
 *  Chromatic Adaptation Transforms
 *  [ ] None
 *  [ ] von Kries
 *  [ ] Bradford
 *  [ ] Sharp
 *  [ ] W-CAT
 *  [ ] Custom
 *  [ ] ...
 *
 *  Illuminant
 *  [ ] D65?
 *
 *  Generality
 *  [ ] Bits per component (8/16/12?)
 *  [ ] Ordering of components (ARGB, BGR, RGBX, Greyscale?)
 *
 *  Parallelism
 *  [ ] OMP
 *  [ ] Implementation-generic threading
 *  [ ] SIMD
 *
 *  Benchmarks
 *  [ ] Against each other
 *  [ ] For common sizes of image (?) / fractions of image size
 *
 *  Terminology
 *  - Component/Channel?
 */

typedef struct awb_image
{
	unsigned char *Data;
	int W;
	int H;
	int Stride; /* 0 -> use Width * Bytes per pixel (BPC/8)* */
	int BitsPerComponent; /* 8/16 */
	char ComponentLayout[4]; /* e.g. "ARGB" */
} awb_image;

#if defined(AWB_IMPLEMENTATION) || defined(AWB_SELFTEST)
/* Input and output images - may be the same */
int awbSimplestColorBalance()
{

}
#endif/*AWB_IMPLEMENTATION*/

#ifdef AWB_SELFTEST
#include <stdio.h>
#define STBI_NO_SIMD
#define STB_IMAGE_IMPLEMENTATION
#define STBI_ONLY_JPEG
#include "../stb/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../stb/stb_image_write.h"

int main() {
	char *ImgName    = "test_images/GYM_9837.JPG";
	char *ImgOutName = "test_images/test.png";
	awb_image Img    = {0};
	int cComponents  = 0;
	int Result       = 0;

	Img.Data = stbi_load(ImgName, &Img.W, &Img.H, &cComponents, 0);
	if(! Img.Data) {
		fprintf(stderr, "Failed to open file %s", ImgName);
		Result = -1; goto end;
	}
	printf("Opened image %s for processing\n", ImgName);
	Img.Stride = Img.W * cComponents;
	printf("Image data:\n"
		   "   Width:  %d\n"
		   "   Height: %d\n"
		   "   Stride: %d\n"
		   "   Channels: %d\n"
		   , Img.W, Img.H, Img.Stride, cComponents);

	stbi_write_png(ImgOutName, Img.W, Img.H, cComponents, Img.Data, Img.Stride);

	end:
	return Result;
}
#endif/*AWB_SELFTEST*/
