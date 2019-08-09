// no copyright
// placed in the public domain
// author intealls

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>
#include <unistd.h>

#ifndef RS_BLOCKSIZE
#define RS_BLOCKSIZE 4
#endif

typedef double T;
typedef double U;

typedef struct sinc_table {
	int order;
	int precision;
	double (*window)(double, double);
	T* kernel;
} sinc_table;

typedef struct sinc_state {
	double rate;
	double leftover;
	T* p;
	bool ready;
	sinc_table* tbl;
} sinc_state;

void     sr_process(sinc_state*, T*, int, U*, int*, int, double);
U        sr_intp(sinc_state*, double, T*, double);
inline U sr_mulsum(T* __restrict, T* __restrict);

static double
sinc(double x)
{
	if (x == 0)
		return 1;
	else
		return sin(M_PI * x) / (M_PI * x);
}

static double
blackman_nuttall(double x, double N)
{
	return 0.3635819 - 0.4891775 * cos((2 * M_PI * x) / N)
	                 + 0.1365995 * cos((4 * M_PI * x) / N)
	                 - 0.0106411 * cos((6 * M_PI * x) / N);
}

static double
kernel_func(sinc_state* st, double x, double cutoff)
{
	if (x > -st->tbl->order && x < st->tbl->order)
		return sinc(cutoff * x) *
		            st->tbl->window(x + 1 - st->tbl->order,
		                            st->tbl->order * 2 - 1);
	else
		return 0;
}

sinc_state*
sr_create(int order, int precision, double cutoff, double rate)
{
	int idx = 0;

	sinc_table* st = (sinc_table*) malloc(sizeof(sinc_table));
	assert(st);

	sinc_state* ss = (sinc_state*) malloc(sizeof(sinc_state));
	assert(ss);

	st->order = order;
	st->precision = precision;
	st->window = blackman_nuttall;

	st->kernel = (T*) malloc(st->order * 2 * (st->precision + 1) * sizeof(T));
	assert(st->kernel);

	ss->rate = rate;
	ss->leftover = 0.0;
	ss->ready = false;
	ss->tbl = st;

	ss->p = (T*) malloc(st->order * 2 * sizeof(T));
	assert(ss->p);

	for (int offset = 0; offset <= precision; offset++) {
		T scale = 0;

		int startidx = idx;

		for (int i = -order + 1; i <= order; i++) {
			st->kernel[idx] = kernel_func(ss,
			                              (T) offset / precision - i,
			                              cutoff);
			scale += st->kernel[idx];
			idx++;
		}

		idx = startidx;
		for (int i = -order + 1; i <= order; i++)
			st->kernel[idx++] /= scale;
	}

	return ss;
}

void
sr_destroy(sinc_state* ss)
{
	assert(ss);

	free(ss->tbl->kernel);
	free(ss->tbl);
	free(ss->p);
	free(ss);
}

inline U
sr_mulsum(T*__restrict samples, T*__restrict coeffs)
{
	T sum = 0;
	T tmp_sum[RS_BLOCKSIZE];

	for (int i = 0; i < RS_BLOCKSIZE; i++)
		tmp_sum[i] = samples[i] * coeffs[i];

	for (int i = 0; i < RS_BLOCKSIZE; i++)
		sum += tmp_sum[i];

	return sum;
}

U
sr_intp(sinc_state* ss, double x, T* input, double gain)
{
	int i;
	double sum = 0;

	double fraction = fabs(x - floor(x));
	int c_idx = (int) (fraction * ss->tbl->precision + 0.5f) * ss->tbl->order * 2;
	T* coeffs = &ss->tbl->kernel[c_idx];

	for (i = (int) floor(x) - ss->tbl->order + 1;
	         i <= (int) floor(x) + ss->tbl->order;
	             i += RS_BLOCKSIZE) {
		T* samp_ptr = 0;
		bool overlap = true;

		// sample index, [0..order*2] are in state->p, [order*2..n] are in input
		//                        state->p     input
		// total input vector = [0..order*2 order*2..n]
		int sample_index = i + ss->tbl->order * 2 - 1;

		// do we have blocksize samples in p?
		if (sample_index + RS_BLOCKSIZE < ss->tbl->order * 2) {
			samp_ptr = &ss->p[sample_index];
			overlap = false;
		} else if (sample_index >= ss->tbl->order * 2) {
			samp_ptr = &input[sample_index - ss->tbl->order * 2];
			overlap = false;
		}

		// if the samples saved from the previous iteration overlap
		// with the ones in this one, do them separately
		if (overlap) {
			for (int j = 0; j < RS_BLOCKSIZE; j++) {
				float sample;

				if (sample_index + j >= ss->tbl->order * 2)
					sample = input[sample_index + j - (ss->tbl->order * 2)];
				else
					sample = ss->p[sample_index + j];

				sum += sample * coeffs[j];
			}
		} else {
			sum += sr_mulsum(samp_ptr, coeffs);
		}

		coeffs += RS_BLOCKSIZE;
	}

	sum *= gain;

	return sum;
}

void
sr_process(sinc_state* ss,
           T* input,
           int n_input,
           U* output,
           int* n_output,
           int stride,
           double gain)
{
	int wr_pos = 0;
	double x = ss->leftover - (ss->ready ? ss->tbl->order : 0);

	if (n_input < ss->tbl->order * 2)
		return;

	while (x < n_input - ss->tbl->order) {
		output[wr_pos] = sr_intp(ss, x, input, gain);

		wr_pos += stride;

		x += ss->rate;
	}

	ss->leftover = x - (n_input - ss->tbl->order);

	memcpy(ss->p,
	       &input[n_input - ss->tbl->order * 2],
	       ss->tbl->order * 2 * sizeof(T));

	ss->ready = true;

	*n_output = wr_pos / stride;
}

static void
chirpd(T* buf, double i_phase, int nsamp, int fs, double f0, double f1)
{
	double k = (f1 - f0) / (nsamp / fs);

	for (int i = 0; i < nsamp; i++) {
		double t = (double) i / fs;
		buf[i] = sin(i_phase + 2 * M_PI * (k / 2 * t * t + f0 * t));
	}
}

static void
freqd(T* buf, double i_phase, int nsamp, int fs, double f0)
{
	for (int i = 0; i < nsamp; i++) {
		double t = (double) i / fs;
		buf[i] = sin(i_phase + 2 * M_PI * (f0 * t));
	}
}

useconds_t
us_passed(struct timeval* start, struct timeval* end)
{
	useconds_t diff;

	diff = (end->tv_sec - start->tv_sec) * 1e6;
	diff += (end->tv_usec - start->tv_usec);

	return diff;
}

int
main(int argc, char* argv[])
{
	FILE* out;

	int order = 64;
	double rt = 1.052348;
	double cutoff = 0.45;
	int fs = 48e3;
	int sec = 2;
	int precision = 10000;
	double gen_sine_f = 1000.f;

	int c;
	while ((c = getopt(argc, argv, "r:c:g:f:t:p:o:")) != -1) {
		switch (c) {
			case 'r':
				if (optarg)
					rt = atof(optarg);
				break;
			case 'c':
				if (optarg)
					cutoff = atof(optarg);
				break;
			case 'g':
				if (optarg)
					gen_sine_f = atof(optarg);
				break;
			case 'f':
				if (optarg)
					fs = atoi(optarg);
				break;
			case 't':
				if (optarg)
					sec = atoi(optarg);
				break;
			case 'p':
				if (optarg)
					precision = atoi(optarg);
				break;
			case 'o':
				if (optarg)
					order = atoi(optarg);
				break;
		}
	}

	int n_input = fs * sec;
	int n_output;

	const int ntrials = 20;

	printf("\nConfiguration\n");
	printf("--\n");
	printf("Sinc order (-o):       %d\n", order);
	printf("Generated sine (-g):   %f Hz\n", gen_sine_f);
	printf("Sample rate (-f):      %d Hz\n", fs);
	printf("Length (-t):           %d seconds\n", sec);
	printf("Resample rate (-r):    %f (%.2f Hz)\n", rt, fs / rt);
	printf("Frequency cutoff (-c): %f\n", cutoff);
	printf("Table precision (-p):  %.6f\n", 1.0 / precision);
	printf("--\n");

	sinc_state* ss = sr_create(order, precision, cutoff, rt);

	T* left = (T*) malloc(n_input * sizeof(T));
	T* right = (T*) malloc(n_input * sizeof(T));
	T* output = (T*) malloc(n_input * 2 * ceil(1 / rt) * sizeof(T));

	chirpd(left, 0, n_input, fs, 0, fs / 2);
	freqd(right, 0, n_input, fs, gen_sine_f);

	struct timeval start, end;

	gettimeofday(&start, NULL);

	for (int i = 0; i < ntrials; i++) {
		sr_process(ss, left, n_input, output, &n_output, 2, 1.0);
		sr_process(ss, right, n_input, output + 1, &n_output, 2, 1.0);
		ss->ready = false;
		ss->leftover = 0.0;
	}

	gettimeofday(&end, NULL);

	double ms = (double) (us_passed(&start, &end)) / 1e3 / ntrials / sec;
	printf("ms per second: %f\n", ms);
	printf("table size: %f MBytes\n",
	        (double) (order * 2 * (precision + 1) * sizeof(T)) / 1024 / 1024);

	out = fopen("data.csv", "w");

	printf("samples output (per channel): %d\n", n_output);

	fprintf(out, "%d,%f,", fs, 1.0);

	for (int i = 0; i < n_output * 2; i++) {
		double d = output[i];
		fprintf(out, "%.16f", d != d || isinf(d) || isnan(d) ? 0 : d);
		if (i != n_output * 2 - 1)
			fprintf(out, ",");
	}

	fclose(out);

	sr_destroy(ss);

	free(left);
	free(right);
	free(output);

	return 0;
}
