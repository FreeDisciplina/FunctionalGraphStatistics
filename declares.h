#ifndef DECLARES_H__
#define DECLARES_H__
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <time.h>
#include <math.h>
#include <string>
#include <cstring>
#include <malloc.h>
#include <memory.h>
#include <immintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <nmmintrin.h>
#include <wmmintrin.h> 
#include <ipp.h>
#include <ippcp.h>
#ifdef __GNUC__
#include <sched.h>
#include <sys/resource.h>
#include <x86intrin.h>
#else
#include <intrin.h>
#endif
#include "mpi.h"

#include <assert.h>

using namespace std;

typedef unsigned long long u64;
typedef unsigned int u32;
typedef unsigned short u16;
typedef unsigned char u8;
typedef long double ld64;

#ifndef nMin
#define nMin (12ULL)
#endif

#ifndef nMax
#define nMax (28ULL)
#endif

#define VERBOSE 0

#define N (1ULL << n)

#ifndef FG_SAMPLE_N
#define FG_SAMPLE_N (24*16*2)
#endif

#ifndef ThreadsN
#define ThreadsN (24*16)
#endif

#define TaskEachThread (FG_SAMPLE_N / ThreadsN)

#ifndef k
#define k ((n + 3ULL) / 4ULL)
#endif
#define mask ( (N==(1ULL << 32)) ? 0xffffffffUL : ((1ULL << n) - 1ULL))

#define data_t u32
#define insert  _mm_insert_epi32
#define extract _mm_extract_epi32

#define INF 0xffffffffffffffffULL

#if defined(_MSC_VER)
#define M_PI 3.1415926535897932384626433832795028841971l 
#endif

#endif