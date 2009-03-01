#include <math.h>
#include <stdio.h>

void ofgm(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	int n = *size; double expBX; double AdivB = *a/ *b;
	while ( --n >= 0){
		expBX = exp(*b * x[n]);
		*ans = *ans + log(*c + *a * expBX) * censor[n];
		*ans = *ans -(*c * x[n]) - AdivB * (expBX - 1);
	}
}
void ofg(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	*c = 0;
	ofgm(a, b, c, s, x, size, censor, ans);
}
/* The objective function for the Logistic and Logistic-Makeham models (if c==0) */
void oflm(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{ 
	int n = *size; double expBX, add1toSAexpBXminus1divB;
	double prdASdivB = *a * *s/ *b;
	while ( --n >= 0){
		expBX = exp(*b * x[n]); 
		add1toSAexpBXminus1divB = 1 + prdASdivB * (expBX - 1);
		*ans = *ans + log((*c + *a * expBX/add1toSAexpBXminus1divB)) * censor[n];
		*ans = *ans + log(exp(- *c * x[n]) * pow(add1toSAexpBXminus1divB,(- 1/ *s)));
	}
}
void ofl(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	*c = 0;
	oflm(a, b, c, s, x, size, censor, ans);
}
void gg(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	int n = *size;
	double suma; double sumb;
	double expBX; double expBXminus1;
	double pow2B = pow(*b,2);
	while ( --n >= 0){
		expBX = exp(*b * x[n]); expBXminus1 = expBX -1;
		suma = suma + censor[n]/ *a;
		suma = suma - expBXminus1/ *b;
		sumb = sumb + x[n] * censor[n];
		sumb = sumb - *a * (expBX * x[n]/ *b - expBXminus1/pow2B);
	}
	ans[0] = suma; ans[1] = sumb;
}
// some values disagree; probably left out a term!
void ggm(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	int n = *size;
	double suma; double sumb; double sumc; 
	double expBX, expBXminus1divB, AexpBX;
	while ( --n >= 0){
		expBX = exp(*b * x[n]); expBXminus1divB = (expBX - 1)/ *b; 
		AexpBX = *a * expBX;
		suma = suma + censor[n]/(*a + *c/expBX); //ok
		suma = suma - expBXminus1divB; //ok
		sumb = sumb + censor[n] * x[n]/(1+ *c * expBX/ *a); //ok
		sumb = sumb + *a * expBXminus1divB/ *b - x[n] * AexpBX/ *b ;//ok
		sumc = sumc + censor[n]/(*c + AexpBX); //ok
		sumc = sumc - x[n]; //ok
	}
	ans[0] = suma; ans[1] = sumb; ans[2] = sumc;
}
void gl(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	int n = *size;
	double suma, sumb, sumc, sums;
	double expBX, expBXminus1, prdASexpBXminus1addB;
	double pow2B = pow(*b,2); double pow2S = pow(*s,2); 
	double prdAB = *a * *b; double prdAS = *a * *s; double prdABS = prdAB * *s;
	while ( --n >= 0){
		expBX = exp(*b * x[n]);
		expBXminus1 = expBX - 1;
		prdASexpBXminus1addB = prdAS * expBXminus1 + *b;
		suma = suma + censor[n] * *b/(*a * prdASexpBXminus1addB); //ok
		suma = suma - expBXminus1/(prdASexpBXminus1addB); //ok
		sumb = sumb + censor[n] * (prdAS * expBX + (pow2B - prdABS) *x[n] - prdAS)/(*b * prdASexpBXminus1addB); //ok
		sumb = sumb - ((prdAB * x[n] - *a) * expBX + *a)/(*b * prdASexpBXminus1addB); //ok
		sums = sums - censor[n] * *a * expBXminus1/prdASexpBXminus1addB; //ok
		sums = sums + (prdASexpBXminus1addB * log(prdASexpBXminus1addB/ *b) - prdAS * expBXminus1)/(pow2S * prdASexpBXminus1addB); //ok
	}
	ans[0] = suma; ans[1] = sumb; ans[2] = sums;
}

void glm(double *a, double *b, double *c, double *s, int *x, int *size, int *censor, double *ans)
{
	int n = *size;
	double suma; double sumb; double sumc; double sums;
	double exp1B = exp(1) * *b; 
	double pow2S = pow(*s,2); double pow2B = pow(*b,2); 
	double prdAS = *a * *s; double prdABS = prdAS * *b;
	double prdACS = prdAS * *c; double prdAB = *a * *b;
	double prdBC = *b * *c;
	double expBX, prdASexpBX, prdACSexpBX, prdABexpBX, BminusprdAS;
	while ( --n >= 0){
		expBX = exp(*b * x[n]); prdASexpBX = prdAS * expBX;
		prdACSexpBX = prdACS * expBX; prdABexpBX = prdAB * expBX;
		BminusprdAS = *b - prdAS;

		suma = suma + censor[n] * (pow2B * expBX)/((prdASexpBX + BminusprdAS)*(prdACSexpBX + prdABexpBX - prdACS + prdBC)); //ok
		suma = suma - expBX/(prdASexpBX + exp1B); //ok
		sumb = sumb + censor[n] * (*a * expBX * (prdASexpBX - prdABS * x[n] + pow2B * x[n] - prdAS))/((prdASexpBX + BminusprdAS)*(prdACSexpBX +prdABexpBX - prdACS + prdBC)); //ok
		sumb = sumb - ((prdAB * x[n] - *a) * expBX)/(prdABS * expBX + exp1B * *b); //ok
		sumc = sumc + censor[n] * (prdASexpBX + BminusprdAS)/(prdACSexpBX + prdABexpBX - prdACS + prdBC); //ok
		sumc = sumc - x[n];
		sums = sums - censor[n] * (*a * prdAB * (expBX-1) * expBX)/((prdASexpBX + BminusprdAS)*(prdACSexpBX + prdABexpBX - prdACS + prdBC)); //ok
		sums = sums + ((prdASexpBX + exp1B)*log(((prdASexpBX + exp1B))/exp1B)- prdASexpBX)/(pow2S * prdASexpBX + exp1B * pow2S); //ok
	}
	ans[0] = suma; ans[1] = sumb; ans[2] = sumc; ans[3] = sums;
}

