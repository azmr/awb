/* AWB - A single-header C library for colour-balancing photos */

/* Feature List / TODO:
 *  Illuminant Estimation Methods
 *  [ ] Gray World Theory
 *  	[X] Vanilla
 *  	[X] + Retinex
 *  	[ ] + Standard Deviation Weighting
 *  	[ ] + Standard Deviation & Luminance Weighting
 *  	[X] Robust (Huo et al., 2005)
 *  [X] Simplest Color Balance (Limare et al., 2011) - colour and contrast
 *  [/] Retinex/Perfect Reflector
 *  	[ ] Multiscale (IPOL - Petro et al., 2014)
 *  	[ ] Poisson/PDE (IPOL - Limare et al., 2011)
 *  [ ] Sensor Correlation (?)
 *  [ ] Chiou's White Balance (white point detection, white balance judge, white balance adjustment)
 *  [ ] Edge detection based
 *  [ ] Combined methods/consensus
 *  [ ] ...
 *
 *  Chromatic Adaptation Transforms
 *  [ ] (Simple gain)
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
 *  [ ] Bits per channel (8/16/12?)
 *  [ ] Ordering of channels (ARGB, BGR, RGBX, Greyscale?)
 *
 *  Parallelism
 *  [ ] OMP
 *  [ ] Implementation-generic threading
 *  [ ] SIMD
 *
 *  Benchmarks
 *  [ ] Against each other
 *  [ ] For common sizes of image (?) / fractions of image size
 */
#ifndef AWB_API
#define AWB_API static
#endif/*AWB_API*/

#ifdef AWB_SELFTEST
#include <stdio.h>

#ifndef awbAssert
#if 0//#ifdef _MSC_VER
#define awbAssert(expr) ((expr) || (fprintf(stderr, "\nAssertion Failed in %s (%s:%d):\n\t%s\n", __FUNCTION__, __FILE__, __LINE__, #expr), __debugbreak, 1))
#else
#define awbAssert(expr) ((expr) || (fprintf(stderr, "\nAssertion Failed in %s (%s:%d):\n\t%s\n", __FUNCTION__, __FILE__, __LINE__, #expr), *(int*)0 == 1))
#endif
#endif/*awbAssert*/


AWB_API void DebugHisto(char *Name, int N, int *Histo, int Step, int Denom)
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
	int BitsPerChannel; /* 8/16 */
	int NumChannels;
	char ChannelLayout[4]; /* e.g. "ARGB" */
	unsigned char R, G, B, A; /* Offset into pixel - never assumed to be set */
} awb_image;

typedef enum awb_result {
	AWB_Ok = 0,
	AWB_Error = -1,
	AWB_UnmatchedDimensions = -2,
	AWB_InvalidConstant = -3,
} awb_result;

typedef struct awb_vector { double E[3]; } awb_vector;
typedef struct awb_matrix { double E[3][3]; } awb_matrix;

awb_matrix awbIdentity = { 1.0, 0.0, 0.0,
	                       0.0, 1.0, 0.0,
	                       0.0, 0.0, 1.0 };

#if defined(AWB_IMPLEMENTATION) || defined(AWB_SELFTEST)
/* Input and output images may be the same */

#define awbAbs(x) (((x) >= 0) ? (x) : -(x))
#define awbNorm(x) ((double)(x)/255.0)
#define awbDenorm(x) ((unsigned char)((x)*255.0 + 0.5))
#define awbClip(x, lo, hi) (
#define awbSign(x) \
((x)               \
	? ((x) > 0     \
		?  1       \
		: -1)      \
	: 0)
#define awbClamp(x, min, max) ( \
	((x) <= (min)) ? (min) : \
	((x) >= (max)) ? (max) : \
	(x))
		
#define awbClamp01(x) awbClamp(x, 0.0, 1.0)
#define awbMin(a,b) (((a) <= (b)) ? (a) : (b))
#define awbMax(a,b) (((a) >= (b)) ? (a) : (b))

AWB_API awb_vector
awbM3V3Mult(awb_matrix Mat, awb_vector Vec)
{ /* TODO (opt): SIMD: */
	awb_vector Result = {
		Mat.E[0][0]*Vec.E[0] + Mat.E[0][1]*Vec.E[1] + Mat.E[0][2]*Vec.E[2],
		Mat.E[1][0]*Vec.E[0] + Mat.E[1][1]*Vec.E[1] + Mat.E[1][2]*Vec.E[2],
		Mat.E[2][0]*Vec.E[0] + Mat.E[2][1]*Vec.E[1] + Mat.E[2][2]*Vec.E[2]
	};
	return Result;
}

AWB_API awb_vector
awbYUVFromRGB(awb_vector RGB)
{
	awb_matrix M = {  .299,  .587,  .114,
		             -.299, -.587,  .886,
		              .701, -.587, -.114 }; /* from Huo et al. */
	awb_vector YUV = awbM3V3Mult(M, RGB);
	return YUV;
}

AWB_API awb_vector
awbYUVFromRGB255(unsigned char R, unsigned char G, unsigned char B) {
	awb_vector RGB = { awbNorm(R), awbNorm(G), awbNorm(B) };
	return awbYUVFromRGB(RGB);
}

AWB_API awb_vector
awbRGBFromYUV(awb_vector YUV)
{
	awb_matrix M = { 1.0, 0.0,          1.0,
		             1.0, -.1942078365, -.5093696763,
		             1.0, 1.0,          0.0 }; /* inverse of above */
	awb_vector RGB = awbM3V3Mult(M, YUV);
	return RGB;
}

AWB_API void
awbRGB255FromYUV(awb_vector YUV, unsigned char *R, unsigned char *G, unsigned char *B) {
	awb_vector RGB = awbRGBFromYUV(YUV);
	*R = awbDenorm(RGB.E[0]), *G = awbDenorm(RGB.E[1]), *B = awbDenorm(RGB.E[2]);
}


AWB_API awb_result
awbRobustGrayWorld(awb_image *Src, awb_image *Dst, double GrayThreshold)
{
	awb_result Result = 0;
	int W = Src->W, H = Src->H, N = W*H;

	if(Dst->W != W || Dst->H != H)
	{ Result = AWB_UnmatchedDimensions; }
	else if(GrayThreshold > 1.0 || GrayThreshold < 0.0)
		//TODO: ErrorTolerance limits
	{ Result = AWB_InvalidConstant; }

	if(Result == AWB_Ok)
	{
		int SrcAdj = 0, DstAdj = 0,
		    SrcRem = 0, DstRem = 0;
		
		/* int iFeedback, U_V; */
		enum { gR, gG, gB, cRGB };

		Src->R = 0, Src->G = 1, Src->B = 2;
		Dst->R = 0, Dst->G = 1, Dst->B = 2;

		/* for(iFeedback = 0; ; ++iFeedback, Src = Dst) */
		{ /* iteratively approach 'optimal' white balance ***************************/
			int cSrcChannels = Src->NumChannels, cDstChannels = Dst->NumChannels;
			int SrcStride = Src->Stride, DstStride = Dst->Stride;
			int x, y, i;
			double GrayAvg[cRGB];
			{
				/* Find avg value of U and V for gray color points ***********************/
				int cGray = 0;
				double GraySum[cRGB] = {0}, dN;
				unsigned char *SrcRow = Src->Data;

				for(y = 0; y < H; ++y, SrcRow += SrcStride) {
					for(x = 0; x < W; ++x) {
						unsigned char *Px = SrcRow + (x*cSrcChannels),
									  r = Px[Src->R], g = Px[Src->G], b = Px[Src->B];

						awb_vector YUV = awbYUVFromRGB255(r, g, b);
						double Y_ = YUV.E[0], u = YUV.E[1], v = YUV.E[2];

						/* TODO (opt): would this be faster as 3 mults + add? */
						if(( (awbAbs(u) + awbAbs(v)) / Y_ ) < GrayThreshold) { /* Found a gray point */
							++cGray;
							GraySum[gR] += (double)r;
							GraySum[gG] += (double)g;
							GraySum[gB] += (double)b;
						}
					}
				}

				// TODO: what if no gray points found?
				dN = (double) cGray;
				GrayAvg[gR] = GraySum[gR] / dN,
				GrayAvg[gG] = GraySum[gG] / dN,
				GrayAvg[gB] = GraySum[gB] / dN;
			}

			{ /* Do corrections ********************************************************/
				double CorrR = GrayAvg[gG]/GrayAvg[gR], CorrB = GrayAvg[gG]/GrayAvg[gB];
				unsigned char *SrcRow = Src->Data, *DstRow = Dst->Data;


				for(y = 0; y < H; ++y, SrcRow += SrcStride, DstRow += DstStride) {
					for(x = 0; x < W; ++x) {
						unsigned char *SrcPx = SrcRow + (x*cSrcChannels);
						unsigned char *DstPx = DstRow + (x*cDstChannels);
						double fR = CorrR * (double)SrcPx[Src->R] + 0.5;
						double fB = CorrB * (double)SrcPx[Src->B] + 0.5;
						DstPx[Dst->R] = fR < 255.0 ? (unsigned char)fR : 255;
						DstPx[Dst->G] = SrcPx[Src->G];
						DstPx[Dst->B] = fB < 255.0 ? (unsigned char)fB : 255;
					}
				}
			}

		}
	}

	return Result;
}

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
		int x, y;
		int SrcStride = Src->Stride, DstStride = Dst->Stride;
		int cSrcChannels = Src->NumChannels, cDstChannels = Dst->NumChannels;
		/* TODO: determine these from channel layout */
		int SrcR = 0, SrcG = 1, SrcB = 2;
		int DstR = 0, DstG = 1, DstB = 2;

		double RSum = 0, GSum = 0, BSum = 0, RAvg, GAvg, BAvg, dN = (double)N;
		{ /* Find avg value of each */
			unsigned char *SrcRow = Src->Data;
			for(y = 0; y < H; ++y, SrcRow += SrcStride) {
				for(x = 0; x < W; ++x) {
					unsigned char *Px = SrcRow + (x*cSrcChannels);
					RSum += (double)Px[SrcR];
					GSum += (double)Px[SrcG];
					BSum += (double)Px[SrcB];
				}
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
	int W = Src->W, H = Src->H;

	/* TODO: error checking */
	if(Dst->W != W || Dst->H != H)
	{ Result = AWB_UnmatchedDimensions; }

	if(Result == AWB_Ok)
	{
		int x, y;
		int SrcStride = Src->Stride, DstStride = Dst->Stride;
		int cSrcChannels = Src->NumChannels, cDstChannels = Dst->NumChannels;
		/* TODO: determine these from channel layout */
		int SrcR = 0, SrcG = 1, SrcB = 2;
		int DstR = 0, DstG = 1, DstB = 2;

		unsigned char RMax = 0, GMax = 0, BMax = 0;
		{ /* Find max value of each */
			unsigned char *SrcRow = Src->Data;
			for(y = 0; y < H; ++y, SrcRow += SrcStride) {
				for(x = 0; x < W; ++x) {
					unsigned char *Px = SrcRow + (x*cSrcChannels);
					if(RMax < Px[SrcR]) { RMax = Px[SrcR]; }
					if(GMax < Px[SrcG]) { GMax = Px[SrcG]; }
					if(BMax < Px[SrcB]) { BMax = Px[SrcB]; }
				}
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

			for(y = 0; y < H; ++y, SrcRow += SrcStride, DstRow += DstStride) {
				for(x = 0; x < W; ++x) {
					unsigned char *SrcPx = SrcRow + (x*cSrcChannels);
					unsigned char *DstPx = DstRow + (x*cDstChannels);
					double R = RCorr * (double)SrcPx[SrcR] + 0.5;
					double B = BCorr * (double)SrcPx[SrcB] + 0.5;
					DstPx[DstR] = R < 255.0 ? (unsigned char)R : 255;
					DstPx[DstG] = SrcPx[SrcG];
					DstPx[DstB] = B < 255.0 ? (unsigned char)B : 255;
				}
			}
		}
	}

	return Result;
}

AWB_API awb_result
awbGrayWorldRetinex(awb_image *Src, awb_image *Dst)
{
	awb_result Result = 0;
	int W = Src->W, H = Src->H;

	/* TODO: error checking */
	if(Dst->W != W || Dst->H != H)
	{ Result = AWB_UnmatchedDimensions; }

	if(Result == AWB_Ok)
	{
		int x, y;
		int SrcStride = Src->Stride, DstStride = Dst->Stride;
		int cSrcChannels = Src->NumChannels, cDstChannels = Dst->NumChannels;
		/* TODO: determine these from channel layout */
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

			for(y = 0; y < H; ++y, SrcRow += SrcStride, DstRow += DstStride) {
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
			}
		}
	}

	return Result;
}


#ifdef AWB_SELFTEST
/* Exposure Adjustment ****************************************************** */
#include "../graph/graph.h"

typedef enum awbZone
{
	awbZONE_0 = 0,
	awbZONE_I,
	awbZONE_II,
	awbZONE_III,
	awbZONE_IV,
	awbZONE_V,
	awbZONE_VI,
	awbZONE_VII,
	awbZONE_VIII,
	awbZONE_IX,
	awbZONE_X,
} awbZone;
double awbZoneVal[] = { 0.0, .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0 };

typedef struct awb_segmentation {
	grf_ds Graph;
	double *InternalDiffs; // This should be kept parallel to vertices
	grf_i *Components; // the roots of each disjoint set (no parallel)
	// TODO: intrusive graphs?
} awb_segmentation;

static void
awbSegmentationFree(awb_segmentation Seg)
{
	arr_free(Seg.InternalDiffs);
	arr_free(Seg.Components);
}


static double
awbMinInternalDiffWithThreshold(awb_segmentation Seg, grf_i C1, grf_i C2, double ThresholdK, double *MinIntDiff)
{
	grf_ds G = Seg.Graph;
	grf_i RepC1 = grfDSFindRootCompress(G,C1),
	      RepC2 = grfDSFindRootCompress(G,C2);

	double thresh1 = ThresholdK / G.Vtxs[RepC1].Info.Size,
	       thresh2 = ThresholdK / G.Vtxs[RepC2].Info.Size,
	       diff1   = Seg.InternalDiffs[RepC1],
	       diff2   = Seg.InternalDiffs[RepC2],
	       d1      = diff1 + thresh1,
	       d2      = diff2 + thresh2;
		   
	double Result = awbMin(d1, d2);
	*MinIntDiff = awbMin(diff1, diff2);
	return Result;
}

/* int awbComponentComparisonPredicate(grf_ds G, grf_i C1, grf_i C2) { */
/* 	int Result = grfDSComponentDiff(G, C1, C2) > awbMinInternalDiff(G, C1, C2); */
/* 	return Result; */
/* } */

// Crude measure of brightness difference
static inline double
awbRGBDiff(awb_image *Src, unsigned char *PxA, unsigned char *PxB)
{
	// TODO: properly
	const int R = 0, G = 1, B = 2;
	double RGB_A = (double)PxA[R] + (double)PxA[G] + (double)PxA[B];
	double RGB_B = (double)PxB[R] + (double)PxB[G] + (double)PxB[B];
	double Result = RGB_A - RGB_B;
	Result = awbAbs(Result);
	return Result;
}

static inline double
awbRGBDist(awb_image *Src, unsigned char *PxA, unsigned char *PxB)
{
	// TODO: properly
	const int R = 0, G = 1, B = 2;
	double dR = awbNorm(PxA[R]) - awbNorm(PxB[R]),
	       dG = awbNorm(PxA[G]) - awbNorm(PxB[G]),
	       dB = awbNorm(PxA[B]) - awbNorm(PxB[B]);
	double Result = dR*dR + dG*dG + dB*dB;
	return sqrt(Result);
}

static inline double
awbRDiff(awb_image *Src, unsigned char *PxA, unsigned char *PxB) {
	double dR = awbNorm(*PxA) - awbNorm(*PxB);
	return awbAbs(dR);
}

static inline grf_ds_edge
awbAddWeightedEdge(awb_image *Src, grf_ds *G, grf_i Vi, grf_i Vj, char *Label) {
	awbAssert(Vi != Vj);
	grf_ds_edge Edge = G->Edges[grfDSAddEdge(G, Vi, Vj, 
		awbRDiff(Src, (unsigned char *)G->Vtxs[Vi].Data, (unsigned char *)G->Vtxs[Vj].Data),
		Label)];
	return Edge;
}

// TODO: may want to be able to specify edge creation
static awb_segmentation
awbSegmentImage(awb_image *Src, double GroupingThreshold)
{ /* Felzenszwalb & HuttenLocher, 2003: Efficient Graph-Based Image Segmentation */
	int x,y,cEdges,iEdge, W = Src->W, H = Src->H, Stride = Src->Stride, cChannels = Src->NumChannels;
	/* awb_edge *Edges = 0; // ? */
	awb_segmentation Seg = {0};

	/* for(int x = 0, y = 0, i = 0; x < ) */
	unsigned char *SrcRow = Src->Data;
	for(y = 0; y < H; ++y, SrcRow += Stride)
	{
		for(x = 0; x < W; ++x)
		{
			unsigned char *SrcPx = SrcRow + (x*cChannels);

			ahd_int iPx = grfDSAddSet(&Seg.Graph, SrcPx);
			arr_push(Seg.InternalDiffs, 0.0);

			// TODO: better edge weighting

			int AfterFirstCol = x >= 1, //iPx % W >= 1,
			    AfterFirstRow = y >= 1, //iPx >= W,
			    BeforeLastCol = x < W-1; //iPx % W != W-1;

			if(AfterFirstCol) {     awbAddWeightedEdge(Src, &Seg.Graph, iPx - 1,       iPx, "Horizontal"); }
			if(AfterFirstRow) {     awbAddWeightedEdge(Src, &Seg.Graph, iPx - W,       iPx, "Vertical");
				if(AfterFirstCol) { awbAddWeightedEdge(Src, &Seg.Graph, (iPx - W) - 1, iPx, "Diagonal Up-Left \\"); }
				if(BeforeLastCol) { awbAddWeightedEdge(Src, &Seg.Graph, (iPx - W) + 1, iPx, "Diagonal Up-Right /"); }
			}
		}
	}
	awbAssert(arr_len(Seg.Graph.Vtxs) == W * H);
	awbAssert(arr_len(Seg.Graph.Edges) == (W-1)*H + (H-1)*W + 2*(H-1)*(W-1)); 

	/* 0. Sort E into Pi = (o1, . . . , om), by non-decreasing (>=) edge weight. */
	grfDSEdgeSortWeightAsc(Seg.Graph.Edges);
	for(grf_i i = 1; i < arr_len(Seg.Graph.Edges); ++i) {
		if(Seg.Graph.Edges[i].Weight != Seg.Graph.Edges[i-1].Weight)
		{ int breakhere = 0; ++breakhere;}
	}

	/* 1. Start with a segmentation S0, where each vertex vi is in its own component. */
	/* awb_ds_cmp *Segmentation = 0; //? */
	/* 2. Repeat step 3 for q = 1, . . . ,m. (m = number of edges in graph) */
	cEdges = arr_len(Seg.Graph.Edges); // NOTE: otherwise keeps adding edges then comparing them!
	for(iEdge = 0; iEdge < cEdges; ++iEdge)
	/* arr_foreachv(Seg.Graph.Edges, iEdge, grf_ds_edge, Edge) */
	{
		/* 3. Construct S^q given S^(q−1) as follows. */
		/* a) Let vi and vj denote the vertices connected by the q-th edge
		   in the ordering, i.e., oq = (vi, vj ). */
		grf_ds_edge Edge = Seg.Graph.Edges[iEdge];
		grf_i Vi = Edge.Vi, Vj = Edge.Vj;

		/* b) If vi and vj are in disjoint components of Sq−1 and w(oq) is
		   small compared to the internal difference of both those components,
		   then merge the two components; otherwise do nothing.
		   More formally, let Ci^q−1 be the component of S^(q−1) containing vi and
		   Cj^(q−1) the component containing vj.
		   If Ci^(q−1) = Cj^(q−1) and
		   w(oq ) <= MInt(Ci^(q−1), Cj^(q−1))
		   then S^q is obtained from S^(q−1) by merging Ci^(q−1) and Cj^(q−1).
		   Otherwise S^q = S^(q−1).*/
		grf_i Rep_i = grfDSFindRootCompress(Seg.Graph, Vi);
		grf_i Rep_j = grfDSFindRootCompress(Seg.Graph, Vj);
		int InDisjointComponents = Rep_i != Rep_j;
		double MinIntDiff = 0;
		/* awb_component PrevComponent_i = awbComponentContaining(Segmentation, Vi); */
		/* awb_component PrevComponent_j = awbComponentContaining(Segmentation, Vj); */
		if(InDisjointComponents &&
				Edge.Weight <= awbMinInternalDiffWithThreshold(Seg, Rep_i, Rep_j, GroupingThreshold, &MinIntDiff) )
		{
			grf_i iRep = grfDSUnion(&Seg.Graph, Rep_i, Rep_j, Edge.Weight);
			// Edge.Weight will always be >= minIntDiff as the edges are sorted in increasing order
			awbAssert(Edge.Weight >= MinIntDiff);
			if(!strcmp(Edge.Label,  "Horizontal"))
			{ int b = 0; ++b; }
			if(!strcmp(Edge.Label,  "Vertical"))
			{ int b = 0; ++b; }
			if(!strcmp(Edge.Label,  "Diagonal Up-Left \\"))
			{ int b = 0; ++b; }
			if(!strcmp(Edge.Label,  "Diagonal Up-Right /"))
			{ int b = 0; ++b; }
			Seg.InternalDiffs[iRep] = Edge.Weight;
		}
	}

	// TODO(Opt): include in above loop
	arr_foreach(Seg.Graph.Vtxs, i) {
		if(grfDSFindRootCompress(Seg.Graph, i) == i)
		{ arr_push(Seg.Components, i); }
	}

	/* 4. Return S = Sm. */
	return Seg;
}

/* From: Yuan & Sun (date unknown): Automatic Exposure Correction of Consumer Photographs */
/* Pipeline:
           Auto-level stretch 
           /                \
Region Segmentation      High Level Features
            \              /
			Region Analysis
			       |
			Estimate Optimal Zone of Region
			       |
			Detail-preserving S-curve adjustment
*/
AWB_API awb_result
awbExposureCurveAdjustment(awb_image *Src, awb_image *Dst, double GroupingThreshold)
{
	awb_result Result = 0;
	int W = Src->W, H = Src->H, cSrcChannels = Src->NumChannels, SrcStride = Src->Stride, x, y;

	/* TODO: error checking */
	if(Dst->W != W || Dst->H != H)
	{ Result = AWB_UnmatchedDimensions; }

	if(Result == AWB_Ok)
	{
		/* Auto-level stretch */
		{ /* Region Segmentation */
			awb_segmentation Seg = awbSegmentImage(Src, GroupingThreshold);
			/* High Level Features */
		}
		/* Region Analysis */
		/* Estimate Optimal Zone of Region */
		/* Detail-preserving S-curve adjustment */
	}

	return Result;
}
#endif /*AWB_SELFTEST*/


#ifdef AWB_SELFTEST
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
	/* TODO: determine these from channel layout */
	int SrcR = 0, SrcG = 1, SrcB = 2;
	int DstR = 0, DstG = 1, DstB = 2;

	int N          = W * H;
	float fN       = (float)N;
	// TODO: ensure they're positive
	int SatLoCount = (int)(LowSaturate  * fN + 0.5f);
	int SatHiCount = (int)(HighSaturate * fN + 0.5f);
	unsigned int SatLoPos   = SatLoCount;
	unsigned int SatHiPos   = W*H - SatHiCount - 1;
	
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
			unsigned char *Px = SrcRow + (x*Src->NumChannels);
			++HistoR[Px[SrcR]];
			++HistoG[Px[SrcG]];
			++HistoB[Px[SrcB]];
		}
		SrcRow += SrcStride;
	}

	/* TODO (feat): make histogram available for user? */
#ifdef AWB_SELFTEST
	{ /* DEBUG: Histos */
		int HistScale = N/500;
		int HistStep = 16;
		DebugHisto("Red",   sizeof(HistoR)/sizeof(*HistoR), HistoR, HistStep, HistScale);
		DebugHisto("Green", sizeof(HistoG)/sizeof(*HistoR), HistoG, HistStep, HistScale);
		DebugHisto("Blue",  sizeof(HistoB)/sizeof(*HistoR), HistoB, HistStep, HistScale);
	}
#endif/*AWB_SELFTEST*/
	
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
				unsigned char *SrcPx = SrcRow + (x*Src->NumChannels);
				unsigned char *DstPx = DstRow + (x*Dst->NumChannels);
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
#endif/*SELFTEST*/
#endif/*AWB_IMPLEMENTATION*/

#ifdef AWB_SELFTEST
#define STBI_NO_SIMD
#define STB_IMAGE_IMPLEMENTATION
#define STBI_ONLY_JPEG
#include "../stb/stb_image.h"
#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include "../stb/stb_image_resize.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "../stb/stb_image_write.h"
#include <stdlib.h>
#define SWEET_NUM_TESTS 512
#define SWEET_NOCOLOUR
#include "../sweet/sweet.h"
typedef struct rgb {
	char r, g, b;
} rgb;
rgb rgbBlack = { 0, 0, 0 };
#define blit_pixel rgb
#include "..\blit-fonts\blit32.h"

#define PrintArrayD(var, N, Wrap) PrintArrayD_(#var": ", (double *)var, N, Wrap)
#define PrintElements(var, N, Wrap) PrintArrayD_(#var": ", (double *)var.E, N, Wrap)
void PrintArrayD_(char *Prefix, double *D, int N, int Wrap)
{
	int i;
	printf("%s{ ", Prefix);
	for(i = 0; i < N; ++i)
	{
		printf("%lf, ", D[i]);
		if(i%Wrap == Wrap-1 && i != N-1)
		{ printf("\n     ");}
	}
	printf("}\n");
}

void TestSegmentation(char *ImgName)
{// segmentation tests
	awb_image Img = {0};
	char buf[256];
	char buf2[256];
	snprintf(buf, sizeof(buf), "test_images/%s.jpg", ImgName);
	
	Img.Data = stbi_load(buf, &Img.W, &Img.H, &Img.NumChannels, 0);
	if(! Img.Data) {
		fprintf(stderr, "Failed to open file %s", buf);
		return;
	}
	unsigned char *GroupVisRand, *GroupVisOrig, *GroupOverlay;
	GroupVisRand = malloc(Img.W * Img.H * Img.NumChannels);
	GroupVisOrig = malloc(Img.W * Img.H * Img.NumChannels);
	GroupOverlay = malloc(Img.W * Img.H * Img.NumChannels);

	/* memcpy(GroupOverlay, Img.Data, Img.Stride * Img.H); */

	awb_segmentation Seg = awbSegmentImage(&Img, 6);
	grf_i PxSets = grfDSCountSets(Seg.Graph);
	awbAssert(PxSets < arr_len(Seg.Graph.Vtxs));

	arr_foreach(Seg.Graph.Vtxs, iVtx) {
		unsigned char *RandPx = GroupVisRand + 3 * iVtx; // assumes no gap between rows (bytes in W == Stride)
		unsigned char *OrigPx = GroupVisOrig + 3 * iVtx;
		unsigned char *OverPx = GroupOverlay + 3 * iVtx;
		double Alpha = 0.5, inv_Alpha = 1.0-Alpha;
		int R=0,G=1,B=2;
		grf_i iRoot = grfDSFindRootCompress(Seg.Graph, iVtx);

		unsigned char *RootPx = Img.Data + 3 * iRoot;
		OrigPx[R] = awbDenorm(awbClamp01(awbNorm(RootPx[R])*1.15));
		OrigPx[G] = awbDenorm(awbClamp01(awbNorm(RootPx[G])*1.15));
		OrigPx[B] = awbDenorm(awbClamp01(awbNorm(RootPx[B])*1.15));

		srand(iRoot);
		RandPx[R] = (unsigned char) (rand() % 200 + 55 );
		RandPx[G] = (unsigned char) (rand() % 200 + 55 );
		RandPx[B] = (unsigned char) (rand() % 200 + 55 );

		// alpha blend group colours over normal image
		OverPx[R] = RandPx[R]*Alpha + OverPx[R]*inv_Alpha;
		OverPx[G] = RandPx[G]*Alpha + OverPx[G]*inv_Alpha;
		OverPx[B] = RandPx[B]*Alpha + OverPx[B]*inv_Alpha;
	}

	snprintf(buf, sizeof(buf), "test_images/%s_GroupVisRand.jpg", ImgName);
	snprintf(buf2, sizeof(buf2), "%d Sets", arr_len(Seg.Components));
	blit32_TextExplicit((rgb*)GroupVisRand, rgbBlack, 4, Img.W, Img.H, blit_Clip, 30, 30, buf2);
	stbi_write_jpg(buf, Img.W, Img.H, Img.NumChannels, GroupVisRand, 80);

	snprintf(buf, sizeof(buf), "test_images/%s_GroupVisOrig.jpg", ImgName);
	stbi_write_jpg(buf, Img.W, Img.H, Img.NumChannels, GroupVisOrig, 80);

	snprintf(buf, sizeof(buf), "test_images/%s_GroupOverlay.jpg", ImgName);
	stbi_write_jpg(buf, Img.W, Img.H, Img.NumChannels, GroupOverlay, 80);


	awbSegmentationFree(Seg);
	free(GroupVisRand);
	free(GroupVisOrig);
	free(GroupOverlay);
	free(Img.Data);
}

int main() {
	int Result       = 0;
	{ /* UNIT TESTS */
		TestGroup("Matrices")
		{
			awb_vector Before = { 3.21, -.8, 48.8 };
			awb_vector After = awbM3V3Mult(awbIdentity, Before);
			TestEq(Before, After);
			PrintElements(Before, 3, 8);
			PrintElements(After, 3, 8);

			awb_matrix A = { 2, -5, 1,
	                         0,  3, 4,
	                        -7,  1, 8 };
	        awb_vector B = { 3, -1, 2 };
	        awb_vector ABExpected = { 13, 5, -6 };
			awb_vector ABActual = awbM3V3Mult(A, B);
	        TestEq(ABExpected, ABActual);
			PrintElements(A, 9, 3);
			PrintElements(B, 3, 3);
			PrintElements(ABExpected, 3, 3);
			PrintElements(ABActual, 3, 3);
		} EndTestGroup;
	}

	{ /* CONTENT TESTS */
#if 0
		char *ImgName	= "test_images/GYM_9837.JPG";
		/* char *ImgName	= "test_images/IMG_5041.JPG"; */
		unsigned char *OutData = 0;
		awb_image Img = {0};

		Img.Data = stbi_load(ImgName, &Img.W, &Img.H, &Img.NumChannels, 0);
		if(! Img.Data) {
			fprintf(stderr, "Failed to open file %s", ImgName);
			Result = -1; goto end;
		}
		printf("Opened image %s for processing\n", ImgName);
		Img.Stride = Img.W * Img.NumChannels;
		printf("Image data:\n"
				"   Width:  %d\n"
				"   Height: %d\n"
				"   Stride: %d\n"
				"   Channels: %d\n"
				, Img.W, Img.H, Img.Stride, Img.NumChannels);

		OutData = malloc(Img.Stride * Img.H);
		memset(OutData, 127, Img.Stride * Img.H); /* check that the data is being used) */ \
		/* memset(Img.Data, 127, Img.Stride * Img.H); /1* check that gray does nothing *1/ \ */

			/* puts("Writing first image"); */
			/* stbi_write_png("test_images/testOrig.png", Img.W, Img.H, Img.NumChannels, Img.Data, Img.Stride); */
#endif

#define TestMethod(method, ...) do { \
		awb_image Out = Img; \
		puts("\n\nPerforming " #method); \
		Out.Data = OutData; \
		awbAssert(! awb## method(&Img, __VA_ARGS__)); \
		puts("Writing image"); \
		stbi_write_jpg("test_images/"#method".jpg", Out.W, Out.H, Out.NumChannels, Out.Data, 150); \
		} while(0) \

		/* TestMethod(SimplestColorBalance, &Out, 0.01f, 0.02f); */
		/* TestMethod(GrayWorld, &Out); */
		/* TestMethod(PerfectReflector, &Out); */
		/* TestMethod(GrayWorldRetinex, &Out); */
		/* TestMethod(RobustGrayWorld, &Out, .25); */
		TestSegmentation("Abstract1");
		TestSegmentation("Abstract2");
		TestSegmentation("Abstract3");
		TestSegmentation("Abstract4");

		/* /1* NOTE: just for the new histogram *1/ */
		/* awbSimplestColorBalance(&Img, &Img, 0.01f, 0.02f); */
	}


	end:
	PrintTestResults(sweetCONTINUE);
	return Result;
}
#endif/*AWB_SELFTEST*/
