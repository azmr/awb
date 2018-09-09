/* AWB - A single-header C library for colour-balancing photos */

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
 *  	[ ] Multiscale (IPOL - Petro et al., 2014)
 *  	[ ] Poisson/PDE (IPOL - Limare et al., 2011)
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
#ifdef AWB_SELFTEST
#include <stdio.h>

#ifndef awbAssert
#define awbAssert(expr) ((expr) || (fprintf(stderr, "\nAssertion Failed in %s (%s:%d):\n\t%s", __FUNCTION__, __FILE__, __LINE__, #expr), __debugbreak(), 1))
#endif/*awbAssert*/
void DebugHisto(char *Name, int N, int *Histo, int Step, int Denom)
{
	int i, j;
	printf("\n\n%s:", Name);
	for(i = Step; i < N; i+=Step) { /* i -> bucket */
		int Sum = 0;
		printf("\n%3d ", i);
		for(j = 0; j < Step; Sum += Histo[i - j++]);
		for(j = 0; j < Sum / Denom; ++j, putc('#', stdout));
		printf("  %d", Sum);
	}
	putc('\n', stdout);
}
#endif

typedef struct awb_image
{
	unsigned char *Data;
	int W;
	int H;
	int Stride; /* 0 -> use Width * Bytes per pixel (BPC/8)* */
	int BitsPerComponent; /* 8/16 */
	int NumComponents;
	char ComponentLayout[4]; /* e.g. "ARGB" */
} awb_image;

#if defined(AWB_IMPLEMENTATION) || defined(AWB_SELFTEST)
/* Input and output images may be the same */

/* Apply channel-wise affine transformation to stretch highest value pixels to white
 * and lowest to black */
/* Clips a small percentage at top and bottom to account for noise */
/* W,  and H should be the same between Src & Dst */
/* TODO (api): check that some pixels come under the sat values */
/* TODO (api): allow user-provided histo? */
void awbSimplestColorBalance(awb_image *Src, awb_image *Dst, float LowSaturate, float HighSaturate)
{
	int x, y, i;
	int W = Src->W, H = Src->H;
	int SrcStride = Src->Stride, DstStride = Dst->Stride;
	/* TODO: determine these from component layout */
	int SrcR = 0, SrcG = 1, SrcB = 2;
	int DstR = 0, DstG = 1, DstB = 2;

	int N          = W * H;
	float fN       = (float)N;
	int SatLoCount = (int)(LowSaturate  * fN + 0.5f);
	int SatHiCount = (int)(HighSaturate * fN + 0.5f);
	int SatLoPos   = SatLoCount;
	int SatHiPos   = W*H - SatHiCount - 1;
	
	/* TODO (opt): do I really need to know if more than 65K/255 pixels are in
	 * a single histo bucket? if not, I can save a lot of space (and maybe do some
	 * more efficent SIMD?), particularly for 16-bit channels
	 *
	 * Alternatively, can I get rid of the middle section, and only include the top
	 * and bottom X buckets of the histogram?
	 * > no, the image might be really dark so all its pixels are low */
	/* NOTE: 16-bit histos (with 65536 buckets) may be searched more efficiently
	 * with e.g. bisection */
	unsigned int HistoR[256] = {0};
	unsigned int HistoG[256] = {0};
	unsigned int HistoB[256] = {0};
	unsigned char MinR = 0, MaxR = 254;
	unsigned char MinG = 0, MaxG = 254;
	unsigned char MinB = 0, MaxB = 254;

	unsigned char *SrcRow = Src->Data, *DstRow = Dst->Data;
	for(y = 0; y < H; ++y)
	{ /* Make histogram of values for each channel */
		for(x = 0; x < W; ++x)
		{
			unsigned char *Px = SrcRow + (x*Src->NumComponents);
			++HistoR[Px[SrcR]];
			++HistoG[Px[SrcG]];
			++HistoB[Px[SrcB]];
		}
		SrcRow += SrcStride;
		DstRow += DstStride;
	}

	/* TODO (feat): make histogram available for user? */
	{ /* DEBUG: Histos */
		int HistScale = N/500;
		int HistStep = 16;
		DebugHisto("Red",   sizeof(HistoR)/sizeof(*HistoR), HistoR, HistStep, HistScale);
		DebugHisto("Green", sizeof(HistoG)/sizeof(*HistoR), HistoG, HistStep, HistScale);
		DebugHisto("Blue",  sizeof(HistoB)/sizeof(*HistoR), HistoB, HistStep, HistScale);
	}
	
	for(i = 1; i < 255; ++i)
	{ /* Make histograms cumulative */
		HistoR[i] += HistoR[i-1];
		HistoG[i] += HistoG[i-1];
		HistoB[i] += HistoB[i-1];
	}

	/* Find min and max sat points */
	for(MinR = 0; HistoR[MinR+1] <= SatLoPos; ++MinR);
	for(MinG = 0; HistoG[MinG+1] <= SatLoPos; ++MinG);
	for(MinB = 0; HistoB[MinB+1] <= SatLoPos; ++MinB);

	/* TODO: does this (254) mean that it will always scale? i.e. not idempotent? */
	for(MaxR = 254; HistoR[MaxR-1] > SatHiPos; --MaxR);
	for(MaxG = 254; HistoG[MaxG-1] > SatHiPos; --MaxG);
	for(MaxB = 254; HistoB[MaxB-1] > SatHiPos; --MaxB);

	if(MaxR < 254) { ++MaxR; }
	if(MaxG < 254) { ++MaxG; }
	if(MaxB < 254) { ++MaxB; }
	

	{ /* Saturate and rescale */
		unsigned char OutMax = 255, OutMin = 0;
		float ScaleR = (OutMax - OutMin) / (MaxR - MinR);
		float ScaleG = (OutMax - OutMin) / (MaxG - MinG);
		float ScaleB = (OutMax - OutMin) / (MaxB - MinB);

		SrcRow = Src->Data, DstRow = Dst->Data;
		for(y = 0; y < H; ++y)
		{ /* Apply scale, copying to each channel */
			for(x = 0; x < W; ++x)
			{
				unsigned char *SrcPx = SrcRow + (x*Src->NumComponents);
				unsigned char *DstPx = DstRow + (x*Dst->NumComponents);
					 if(SrcPx[SrcR] < MinR) { DstPx[DstR] = MinR; }
				else if(SrcPx[SrcR] > MaxR) { DstPx[DstR] = MaxR; }
				else                        { DstPx[DstR] = SrcPx[SrcR]; }

				     if(SrcPx[SrcG] < MinG) { DstPx[DstG] = MinG; }
				else if(SrcPx[SrcG] > MaxG) { DstPx[DstG] = MaxG; }
				else                        { DstPx[DstG] = SrcPx[SrcG]; }

				     if(SrcPx[SrcB] < MinB) { DstPx[DstB] = MinB; }
				else if(SrcPx[SrcB] > MaxB) { DstPx[DstB] = MaxB; }
				else                        { DstPx[DstB] = SrcPx[SrcB]; }

				DstPx[DstR] = (DstPx[DstR] - MinR) * ScaleR;
				DstPx[DstG] = (DstPx[DstG] - MinG) * ScaleG;
				DstPx[DstB] = (DstPx[DstB] - MinB) * ScaleB;
			}
			SrcRow += SrcStride;
			DstRow += DstStride;
		}
	}

}
#endif/*AWB_IMPLEMENTATION*/

#ifdef AWB_SELFTEST
#define STBI_NO_SIMD
#define STB_IMAGE_IMPLEMENTATION
#define STBI_ONLY_JPEG
#include "../stb/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../stb/stb_image_write.h"

int main() {
	char *ImgName    = "test_images/GYM_9837.JPG";
	awb_image Img    = {0};
	int Result       = 0;

	Img.Data = stbi_load(ImgName, &Img.W, &Img.H, &Img.NumComponents, 0);
	if(! Img.Data) {
		fprintf(stderr, "Failed to open file %s", ImgName);
		Result = -1; goto end;
	}
	printf("Opened image %s for processing\n", ImgName);
	Img.Stride = Img.W * Img.NumComponents;
	printf("Image data:\n"
		   "   Width:  %d\n"
		   "   Height: %d\n"
		   "   Stride: %d\n"
		   "   Channels: %d\n"
		   , Img.W, Img.H, Img.Stride, Img.NumComponents);

	/* puts("Writing first image"); */
	/* stbi_write_png("test_images/testOrig.png", Img.W, Img.H, Img.NumComponents, Img.Data, Img.Stride); */

	puts("Performing color balance");
	awbSimplestColorBalance(&Img, &Img, 0.01f, 0.02f);

	puts("Writing second image");
	stbi_write_png("test_images/testProc.png", Img.W, Img.H, Img.NumComponents, Img.Data, Img.Stride);

	/* NOTE: just for the new histogram */
	awbSimplestColorBalance(&Img, &Img, 0.01f, 0.02f);

	end:
	return Result;
}
#endif/*AWB_SELFTEST*/
