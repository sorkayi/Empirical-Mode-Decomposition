#include "EmpiricalModeDecomposition.h"

void emdSetup(emdData* emd, int order, int iterations, int locality) {
	emd->iterations = iterations;
	emd->order = order;
	emd->locality = locality;
	emd->size = 0;
	emd->imfs = NULL;
	emd->residue = NULL;
	emd->minPoints = NULL;
	emd->maxPoints = NULL;
	emd->min = NULL;
	emd->max = NULL;
}

void emdResize(emdData* emd, int size) {
	int i;

	emdClear(emd);

	emd->size = size;
	emd->imfs = cnew(float*, emd->order);
	for(i = 0; i < emd->order; i++)
		emd->imfs[i] = cnew(float, size);
	emd->residue = cnew(float, size);
	emd->minPoints = cnew(int, size / 2);
	emd->maxPoints = cnew(int, size / 2);
	emd->min = cnew(float, size);
	emd->max = cnew(float, size);
}

void emdCreate(emdData* emd, int size, int order, int iterations, int locality) {
	emdSetup(emd, order, iterations, locality);
	emdResize(emd, size);
}

void emdClear(emdData* emd) {
	int i;
	if(emd->imfs != NULL) {
		for(i = 0; i < emd->order; i++)
			cdelete(emd->imfs[i]);
		cdelete(emd->imfs);
		cdelete(emd->minPoints);
		cdelete(emd->maxPoints);
		cdelete(emd->min);
		cdelete(emd->max);
		cdelete(emd->residue);
	}
}

void emdDecompose(emdData* emd, const float* signal) {
	int i, j;
	memcpy(emd->imfs[0], signal, emd->size * sizeof(float));
	memcpy(emd->residue, signal, emd->size * sizeof(float));
	for(i = 0; i < emd->order - 1; i++) {
		float* curImf = emd->imfs[i];
		for(j = 0; j < emd->iterations; j++) {
			emdMakeExtrema(emd, curImf);
			if(emd->minSize < 4 || emd->maxSize < 4)
				break; // can't fit splines
			emdInterpolate(emd, curImf, emd->min, emd->minPoints, emd->minSize);
			emdInterpolate(emd, curImf, emd->max, emd->maxPoints, emd->maxSize);
			emdUpdateImf(emd, curImf);
		}
		emdMakeResidue(emd, curImf);
		memcpy(emd->imfs[i + 1], emd->residue, emd->size * sizeof(float));
	}
}

// Currently, extrema within (locality) of the boundaries are not allowed.
// A better algorithm might be to collect all the extrema, and then assume
// that extrema near the boundaries are valid, working toward the center.

void emdMakeExtrema(emdData* emd, const float* curImf) {
	int i, lastMin = 0, lastMax = 0;
	emd->minSize = 0;
	emd->maxSize = 0;
	for(i = 1; i < emd->size - 1; i++) {
		if(curImf[i - 1] < curImf[i]) {
			if(curImf[i] > curImf[i + 1] && (i - lastMax) > emd->locality) {
				emd->maxPoints[emd->maxSize++] = i;
				lastMax = i;
			}
		} else {
			if(curImf[i] < curImf[i + 1] && (i - lastMin) > emd->locality) {
				emd->minPoints[emd->minSize++] = i;
				lastMin = i;
			}
		}
	}
}

void emdInterpolate(emdData* emd, const float* in, float* out, int* points, int pointsSize) {
	int size = emd->size;
	int i, j, i0, i1, i2, i3, start, end;
	float a0, a1, a2, a3;
	float y0, y1, y2, y3, muScale, mu;
	for(i = -1; i < pointsSize; i++) {
		i0 = points[mirrorIndex(i - 1, pointsSize)];
		i1 = points[mirrorIndex(i, pointsSize)];
		i2 = points[mirrorIndex(i + 1, pointsSize)];
		i3 = points[mirrorIndex(i + 2, pointsSize)];

		y0 = in[i0];
		y1 = in[i1];
		y2 = in[i2];
		y3 = in[i3];

		a0 = y3 - y2 - y0 + y1;
		a1 = y0 - y1 - a0;
		a2 = y2 - y0;
		a3 = y1;

		// left boundary
		if(i == -1) {
			start = 0;
			i1 = -i1;
		} else
			start = i1;

		// right boundary
		if(i == pointsSize - 1) {
			end = size;
			i2 = size + size - i2;
		} else
			end = i2;

		muScale = 1.f / (i2 - i1);
		for(j = start; j < end; j++) {
			mu = (j - i1) * muScale;
			out[j] = ((a0 * mu + a1) * mu + a2) * mu + a3;
		}
	}
}

void emdUpdateImf(emdData* emd, float* imf) {
	int i;
	for(i = 0; i < emd->size; i++)
		imf[i] -= (emd->min[i] + emd->max[i]) * .5f;
}

void emdMakeResidue(emdData* emd, const float* cur) {
	int i;
	for(i = 0; i < emd->size; i++)
		emd->residue[i] -= cur[i];
}

inline int mirrorIndex(int i, int size) {
	if(i < size) {
		if(i < 0)
			return -i - 1;
		return i;
	}
	return (size - 1) + (size - i);
}
