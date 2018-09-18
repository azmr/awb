/* AWB - A single-header C library for colour-balancing photos */

/* Feature List / TODO:
 *  Illuminant Estimation Methods
 *  [ ] Gray World Theory
 *  	[X] Vanilla
 *  	[ ] + Retinex
 *  	[ ] + Standard Deviation Weighting
 *  	[ ] + Standard Deviation & Luminance Weighting
 *  	[ ] Robust (Huo et al., 2005)
 *  [X] Simplest Color Balance (Limare et al., 2011) - colour and contrast
 *  [/] Retinex/Perfect Reflector
 *  	[ ] Multiscale (IPOL - Petro et al., 2014)
 *  	[ ] Poisson/PDE (IPOL - Limare et al., 2011)
 *  [ ] Sensor Correlation (?)
 *  [ ] Chiou's White Balance (white point detection, white balance judge, white balance adjustment)
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
#ifndef AWB_API
#define AWB_API static
#endif/*AWB_API*/

#ifdef AWB_SELFTEST
#include <stdio.h>

#ifndef awbAssert
#ifdef _MSC_VER
#define awbAssert(expr) ((expr) || (fprintf(stderr, "\nAssertion Failed in %s (%s:%d):\n\t%s", __FUNCTION__, __FILE__, __LINE__, #expr), __debugbreak, 1))
#else
#define awbAssert(expr) ((expr) || (fprintf(stderr, "\nAssertion Failed in %s (%s:%d):\n\t%s", __FUNCTION__, __FILE__, __LINE__, #expr), *(void*)0, 1))
#endif
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

typedef enum awb_result {
	AWB_Ok = 0,
	AWB_Error = -1,
	AWB_UnmatchedDimensions = -2,
} awb_result;

#if defined(AWB_IMPLEMENTATION) || defined(AWB_SELFTEST)
/* Input and output images may be the same */

AWB_API awb_result
awbGrayWorld(awb_image *Src, awb_image *Dst)
{
	awb_result Result = 0;
	int W = Src->W, H = Src->H, N = W*H;

	/* TODO: error checking */
	if(Dst->W != W || Dst->H != H)
	{ Result = AWB_UnmatchedDimensions; }

	if(Result == AWB_Ok)
	{
		int x, y, i;
		int SrcStride = Src->Stride, DstStride = Dst->Stride;
		int cSrcChannels = Src->NumComponents, cDstChannels = Dst->NumComponents;
		/* TODO: determine these from component layout */
		int SrcR = 0, SrcG = 1, SrcB = 2;
		int DstR = 0, DstG = 1, DstB = 2;

		double RSum = 0, GSum = 0, BSum = 0, RAvg, GAvg, BAvg, dN = (double)N;
		{ /* Find avg value of each */
			unsigned char *SrcRow = Src->Data;
			for(y = 0; y < H; ++y) {
				for(x = 0; x < W; ++x) {
					unsigned char *Px = SrcRow + (x*cSrcChannels);
					RSum += (double)Px[SrcR];
					GSum += (double)Px[SrcG];
					BSum += (double)Px[SrcB];
				}
				SrcRow += SrcStride;
			}

			RAvg = RSum/dN, GAvg = GSum/dN, BAvg = BSum/dN;
		}

		{ /* Do corrections */
			double RCorr = GAvg/RAvg, BCorr = GAvg/BAvg;
			unsigned char *SrcRow = Src->Data, *DstRow = Dst->Data;

#ifdef AWB_SELFTEST
			printf( "RSum = %lf,   GSum = %lf,   BSum = %lf\n"
					"RAvg = %lf,   GAvg = %lf,   BAvg = %lf\n"
					"RCorr = %lf,  BCorr = %lf\n",
					RSum, GSum, BSum,
					RAvg, GAvg, BAvg,
					RCorr, BCorr);
#endif/*AWB_SELFTEST*/

			for(y = 0; y < H; ++y) {
				for(x = 0; x < W; ++x) {
					unsigned char *SrcPx = SrcRow + (x*cSrcChannels);
					unsigned char *DstPx = DstRow + (x*cDstChannels);
					double R = RCorr * (double)SrcPx[SrcR] + 0.5;
					double B = BCorr * (double)SrcPx[SrcB] + 0.5;
					DstPx[DstR] = R < 255.0 ? (unsigned char)R : 255;
					DstPx[DstG] = SrcPx[SrcG];
					DstPx[DstB] = B < 255.0 ? (unsigned char)B : 255;
				}
				SrcRow += SrcStride;
				DstRow += DstStride;
			}
		}
	}

	return Result;
}

/* NOTE: does nothing if there are white pixels */
/* TODO: use histo & add saturation */
AWB_API awb_result
awbPerfectReflector(awb_image *Src, awb_image *Dst)
{
	awb_result Result = 0;
	int W = Src->W, H = Src->H, N = W*H;

	/* TODO: error checking */
	if(Dst->W != W || Dst->H != H)
	{ Result = AWB_UnmatchedDimensions; }

	if(Result == AWB_Ok)
	{
		int x, y, i;
		int SrcStride = Src->Stride, DstStride = Dst->Stride;
		int cSrcChannels = Src->NumComponents, cDstChannels = Dst->NumComponents;
		/* TODO: determine these from component layout */
		int SrcR = 0, SrcG = 1, SrcB = 2;
		int DstR = 0, DstG = 1, DstB = 2;

		unsigned char RMax = 0, GMax = 0, BMax = 0;
		{ /* Find max value of each */
			unsigned char *SrcRow = Src->Data;
			for(y = 0; y < H; ++y) {
				for(x = 0; x < W; ++x) {
					unsigned char *Px = SrcRow + (x*cSrcChannels);
					if(RMax < Px[SrcR]) { RMax = Px[SrcR]; }
					if(GMax < Px[SrcG]) { GMax = Px[SrcG]; }
					if(BMax < Px[SrcB]) { BMax = Px[SrcB]; }
				}
				SrcRow += SrcStride;
			}
		}

		{ /* Do corrections */
			double RCorr = GMax/(double)RMax, BCorr = GMax/(double)BMax;
			unsigned char *SrcRow = Src->Data, *DstRow = Dst->Data;

#ifdef AWB_SELFTEST
			printf( "RMax = %d,   GMax = %d,   BMax = %d\n"
					"RCorr = %lf,  BCorr = %lf\n",
					RMax, GMax, BMax,
					RCorr, BCorr);
#endif/*AWB_SELFTEST*/

			for(y = 0; y < H; ++y) {
				for(x = 0; x < W; ++x) {
					unsigned char *SrcPx = SrcRow + (x*cSrcChannels);
					unsigned char *DstPx = DstRow + (x*cDstChannels);
					double R = RCorr * (double)SrcPx[SrcR] + 0.5;
					double B = BCorr * (double)SrcPx[SrcB] + 0.5;
					DstPx[DstR] = R < 255.0 ? (unsigned char)R : 255;
					DstPx[DstG] = SrcPx[SrcG];
					DstPx[DstB] = B < 255.0 ? (unsigned char)B : 255;
				}
				SrcRow += SrcStride;
				DstRow += DstStride;
			}
		}
	}

	return Result;
}

AWB_API awb_result
awbGrayWorldRetinex(awb_image *Src, awb_image *Dst)
{
	awb_result Result = 0;
	int W = Src->W, H = Src->H, N = W*H;

	/* TODO: error checking */
	if(Dst->W != W || Dst->H != H)
	{ Result = AWB_UnmatchedDimensions; }

	if(Result == AWB_Ok)
	{
		int x, y, i;
		int SrcStride = Src->Stride, DstStride = Dst->Stride;
		int cSrcChannels = Src->NumComponents, cDstChannels = Dst->NumComponents;
		/* TODO: determine these from component layout */
		int SrcR = 0, SrcG = 1, SrcB = 2;
		int DstR = 0, DstG = 1, DstB = 2;
		double RSum = 0.0, GSum = 0.0, BSum = 0.0;
		double RSqSum = 0.0, BSqSum = 0.0;
		double ToNorm = 1.0/255.0;

		unsigned char uRMax = 0, uGMax = 0, uBMax = 0;
		{ /* Find max value of each */
			unsigned char *SrcRow = Src->Data;
			for(y = 0; y < H; ++y) {
				for(x = 0; x < W; ++x) {
					unsigned char *Px = SrcRow + (x*cSrcChannels);
					double R = ToNorm * (double)Px[SrcR],
					       G = ToNorm * (double)Px[SrcG],
					       B = ToNorm * (double)Px[SrcB];

					RSum += R;
					GSum += G;
					BSum += B;

					RSqSum += R*R;
					BSqSum += B*B;

					if(uRMax < Px[SrcR]) { uRMax = Px[SrcR]; }
					if(uGMax < Px[SrcG]) { uGMax = Px[SrcG]; }
					if(uBMax < Px[SrcB]) { uBMax = Px[SrcB]; }
				}
				SrcRow += SrcStride;
			}
		}

		/*
		 * | u |   | SumR^2   SumR |-1  | SumG |
		 * |   | = |               |  . |      |
		 * | v |   | MaxR^2   MaxR |    | MaxG |
		 *
		 * R' = u*R*R + v*R
		 */
		{ /* Do corrections */
			double RMax = ToNorm * (double)uRMax,                  /* d */
			       GMax = ToNorm * (double)uGMax,
			       BMax = ToNorm * (double)uBMax;
			double RSqMax = RMax*RMax, BSqMax = BMax*BMax;         /* c */
			double RDenom = 1.0 / (RSqSum * RMax - RSqMax * RSum), /* 1/(ad-bc) */
			       BDenom = 1.0 / (BSqSum * BMax - BSqMax * BSum);
			double Ru = (   RMax*GSum -   RSum*GMax) * RDenom,
			       Rv = (-RSqMax*GSum + RSqSum*GMax) * RDenom,
			       Bu = (   BMax*GSum -   BSum*GMax) * BDenom,
			       Bv = (-BSqMax*GSum + BSqSum*GMax) * BDenom;

			unsigned char *SrcRow = Src->Data, *DstRow = Dst->Data;
#ifdef AWB_SELFTEST
			printf( "RMax = %lf,   GMax = %lf,   BMax = %lf\n"
					"RSum = %lf,   GSum = %lf,   BSum = %lf\n"
					"RSqMax = %lf, BSqMax = %lf\n"
					"RSqSum = %lf, BSqSum = %lf\n\n"
					"RDenom: %.16lf, BDenom: %.16lf\n"
					"Ru = %lf,  Rv = %lf\n"
					"Bu = %lf,  Bv = %lf\n",
					RMax, GMax, BMax,
					RSum, GSum, BSum,
					RSqMax, BSqMax,
					RSqSum, BSqSum,
					RDenom, BDenom,
					Ru, Rv, Bu, Bv);
#endif/*AWB_SELFTEST*/

			for(y = 0; y < H; ++y) {
				for(x = 0; x < W; ++x) {
					unsigned char *SrcPx = SrcRow + (x*cSrcChannels);
					unsigned char *DstPx = DstRow + (x*cDstChannels);
					double R = ((double)SrcPx[SrcR] / 255.0);
					double B = ((double)SrcPx[SrcB] / 255.0);
					R = Ru * R*R + Rv * R;
					B = Bu * B*B + Bv * B;
					DstPx[DstR] = R < 1.0 ? (unsigned char)(255.0*R + 0.5) : 255;
					DstPx[DstG] = SrcPx[SrcG];
					DstPx[DstB] = B < 1.0 ? (unsigned char)(255.0*B + 0.5) : 255;
				}
				SrcRow += SrcStride;
				DstRow += DstStride;
			}
		}
	}

	return Result;
}



/* Apply channel-wise affine transformation to stretch highest value pixels to white
 * and lowest to black */
/* Clips a small percentage at top and bottom to account for noise */
/* W,  and H should be the same between Src & Dst */
/* TODO (api): check that some pixels come under the sat values */
/* TODO (api): allow user-provided histo? */
/* TODO: error checking */
AWB_API awb_result
awbSimplestColorBalance(awb_image *Src, awb_image *Dst, float LowSaturate, float HighSaturate)
{
	awb_result Result = 0;
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

	unsigned char *SrcRow = Src->Data, *DstRow = 0;
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

	return Result;
}
#endif/*AWB_IMPLEMENTATION*/

#ifdef AWB_SELFTEST
#define STBI_NO_SIMD
#define STB_IMAGE_IMPLEMENTATION
#define STBI_ONLY_JPEG
#include "../stb/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../stb/stb_image_write.h"
#include <stdlib.h>

int main() {
	int Result       = 0;
	char *ImgName    = "test_images/GYM_9837.JPG";
	unsigned char *OutData = 0;
	awb_image Img = {0};

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
	
	OutData = malloc(Img.Stride * Img.H);
	memset(OutData, 127, Img.Stride * Img.H); /* check that the data is being used) */ \

	/* puts("Writing first image"); */
	/* stbi_write_png("test_images/testOrig.png", Img.W, Img.H, Img.NumComponents, Img.Data, Img.Stride); */

#define TestMethod(method, ...) do { \
		awb_image Out = Img; \
		puts("\n\nPerforming " #method); \
		Out.Data = OutData; \
		awbAssert(! awb## method(&Img, __VA_ARGS__)); \
		puts("Writing image"); \
		stbi_write_png("test_images/"#method".png", Out.W, Out.H, Out.NumComponents, Out.Data, Out.Stride); \
	} while(0)

	/* TestMethod(SimplestColorBalance, &Out, 0.01f, 0.02f); */
	/* TestMethod(GrayWorld, &Out); */
	/* TestMethod(PerfectReflector, &Out); */
	TestMethod(GrayWorldRetinex, &Out);


	/* /1* NOTE: just for the new histogram *1/ */
	/* awbSimplestColorBalance(&Img, &Img, 0.01f, 0.02f); */

	end:
	return Result;
}
#endif/*AWB_SELFTEST*/
