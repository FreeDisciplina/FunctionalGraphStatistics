#include "data_structures.h"

using namespace std;

int idx, size;

MPI_File fhcomponentN;
MPI_File fhcyclicNodeN;
MPI_File fhtailNodeN;
MPI_File fhterminalN;
MPI_File fhimageN;
MPI_File fhk_thNodeN;

ld64 componentN_theorem        = 0.0L;
ld64 cyclicNodeN_theorem       = 0.0L;
ld64 onlyCyclicNodeN_theorem   = 0.0L;
ld64 tailNodeN_theorem         = 0.0L;
ld64 terminalN_theorem         = 0.0L;
ld64 imageN_theorem            = 0.0L;
ld64 k_thNodeN_theorem         = 0.0L;

ld64 componentNsum_thread      = 0.0L;
ld64 cyclicNodeNsum_thread     = 0.0L;
ld64 onlyCyclicNodeNsum_thread = 0.0L;
ld64 tailNodeNsum_thread       = 0.0L;
ld64 terminalNsum_thread       = 0.0L;
ld64 imageNsum_thread          = 0.0L;
ld64 k_thNodeNsum_thread       = 0.0L;

ld64 componentNsum_total       = 0.0L;
ld64 cyclicNodeNsum_total      = 0.0L;
ld64 onlyCyclicNodeNsum_total  = 0.0L;
ld64 tailNodeNsum_total        = 0.0L;
ld64 terminalNsum_total        = 0.0L;
ld64 imageNsum_total           = 0.0L;
ld64 k_thNodeNsum_total        = 0.0L;

ld64 componentNsqu_thread      = 0.0L;
ld64 cyclicNodeNsqu_thread     = 0.0L;
ld64 onlyCyclicNodeNsqu_thread = 0.0L;
ld64 tailNodeNsqu_thread       = 0.0L;
ld64 terminalNsqu_thread       = 0.0L;
ld64 imageNsqu_thread          = 0.0L;
ld64 k_thNodeNsqu_thread       = 0.0L;

ld64 componentNsqu_total       = 0.0L;
ld64 cyclicNodeNsqu_total      = 0.0L;
ld64 onlyCyclicNodeNsqu_total  = 0.0L;
ld64 tailNodeNsqu_total        = 0.0L;
ld64 terminalNsqu_total        = 0.0L;
ld64 imageNsqu_total           = 0.0L;
ld64 k_thNodeNsqu_total        = 0.0L;

ld64 componentNmax_thread      = 0.0L;
ld64 cyclicNodeNmax_thread     = 0.0L;
ld64 onlyCyclicNodeNmax_thread = 0.0L;
ld64 tailNodeNmax_thread       = 0.0L;
ld64 terminalNmax_thread       = 0.0L;
ld64 imageNmax_thread          = 0.0L;
ld64 k_thNodeNmax_thread       = 0.0L;

ld64 componentNmax_total       = 0.0L;
ld64 cyclicNodeNmax_total      = 0.0L;
ld64 onlyCyclicNodeNmax_total  = 0.0L;
ld64 tailNodeNmax_total        = 0.0L;
ld64 terminalNmax_total        = 0.0L;
ld64 imageNmax_total           = 0.0L;
ld64 k_thNodeNmax_total        = 0.0L;

ld64 componentNmin_thread      = (ld64)INF;
ld64 cyclicNodeNmin_thread     = (ld64)INF;
ld64 onlyCyclicNodeNmin_thread = (ld64)INF;
ld64 tailNodeNmin_thread       = (ld64)INF;
ld64 terminalNmin_thread       = (ld64)INF;
ld64 imageNmin_thread          = (ld64)INF;
ld64 k_thNodeNmin_thread       = (ld64)INF;

ld64 componentNmin_total       = (ld64)INF;
ld64 cyclicNodeNmin_total      = (ld64)INF;
ld64 onlyCyclicNodeNmin_total  = (ld64)INF;
ld64 tailNodeNmin_total        = (ld64)INF;
ld64 terminalNmin_total        = (ld64)INF;
ld64 imageNmin_total           = (ld64)INF;
ld64 k_thNodeNmin_total        = (ld64)INF;

#define WIDTH 25

char componentN_log  [TaskEachThread][WIDTH + 2];
char cyclicNodeN_log [TaskEachThread][WIDTH + 2];
char tailNodeN_log   [TaskEachThread][WIDTH + 2];
char terminalN_log   [TaskEachThread][WIDTH + 2];
char imageN_log      [TaskEachThread][WIDTH + 2];
char k_thNodeN_log   [TaskEachThread][WIDTH + 2];

#if defined(__GNUC__)
template <typename T>
std::string to_string(T value)
{
	std::ostringstream os;
	os << value;
	return os.str();
}
typedef unsigned long mUL64;

inline unsigned long long read_tsc(void)
{
#if defined(__i386__)
	unsigned long long cycles;
	__asm__ volatile (".byte 0x0f, 0x31" : "=A"(cycles));
	return cycles;
#else
#if defined(__x86_64__)
	unsigned int hi, lo;
	__asm__ volatile ("rdtsc" : "a="(lo), "=d"(hi));
	return (((unsigned long long)lo) | ((unsigned long long)(hi) << 32));
#else
#error "Unsupported architecture for counting cycles"
#endif
#endif
}

#else
typedef unsigned long long mUL64;

#pragma intrinsic( __rdtsc )
__inline unsigned long long read_tsc(void)
{
	return __rdtsc();
}
#endif

#define RAND(a,b) (((a = 36969 * (a & 65535) + (a >> 16)) << 16) + \
	(b = 18000 * (b & 65535) + (b >> 16))  )

void block_rndfill(unsigned char *buf, const int len)
{
	static unsigned long a[2], mt = 1, count = 4;
	static unsigned char r[4];
	int                  i;

	if (mt) { mt = 0; *(unsigned long long*)a = read_tsc(); }

	for (i = 0; i < len; ++i)
	{
		if (count == 4)
		{
			*(unsigned long*)r = RAND(a[0], a[1]);
			count = 0;
		}

		buf[i] = r[count++];
	}
}

template <typename T>
std::string ld_to_string(T value, size_t digits = 20, size_t width = WIDTH)
{
	std::ostringstream os;
	os << std::fixed << std::setprecision(digits) << std::setfill(' ') << setw(width) << value;
	return os.str();
}

u64 n;
string fn;

IppsAESSpec * pCtx1;
IppsSMS4Spec * pCtx2;

void init_func()
{
	IppStatus status;
	int keylen;
	int ctxSize;
	Ipp8u pKey[16];

	/* Init Intel IPP library */
	ippInit();

	keylen = 16;
	block_rndfill(pKey, 16);
	//while (0 == _rdseed64_step((mUL64 *)&pKey[0])) {}
	//while (0 == _rdseed64_step((mUL64 *)&pKey[8])) {}
	initfunc(AES, pCtx1);
	initfunc(SMS4, pCtx2);
}

void end_func()
{
	endfunc(pCtx1);
	endfunc(pCtx2);
}

info_node * drawFunctionalG(info_node * FG, func_t func, int si)
{
#if (VERBOSE==1)
	ofstream fout;
	fout.open(fn.c_str(), ios::app);
	fout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	fout << " n = " << n << endl;
	fout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	clock_t t0, t1;
	t0 = clock();
#endif

	u64 componentN = 0ULL;
	u64 cyclicNodeN = 0ULL;
	u64 onlyCyclicNodeN = 0ULL;
	u64 tailNodeN = 0ULL;
	u64 terminalN = 0ULL;
	u64 imageN = 0ULL;
	u64 k_thNodeN = 0ULL;

	FG = (info_node *)malloc(sizeof(info_node) * N);
	assert(FG != NULL);

	u64 sum = 0ULL;
	data_t x = 0x0;
	u64 * indicator;
	indicator = (u64 *)malloc((N >> 6) * sizeof(u64));
	assert(indicator != NULL);

	for (u64 i = 0; i < (N >> 6); i++)
	{
		indicator[i] = 0xffffffffffffffffULL;
	}

#if (VERBOSE==1)
	cout << "Generating random mapping: " << endl;
#endif

	for (u64 i = 0; i < N; i++)
	{
		FG[i].image = func(i);
		FG[i].preimageN = 0ULL;
		FG[i].depth = 0ULL;
		setBit(indicator, FG[i].image);
	}

#if  (VERBOSE==1)
	cout << "Done." << endl;
#endif

	terminalN = getIndicatedNumber(indicator);
	imageN = N - terminalN;

	data_t * terminals;
	terminals = (data_t *)malloc(terminalN * sizeof(data_t));
	assert(terminals != NULL);
	for (u64 i = 0; i < (terminalN); i++)
	{
		terminals[i] = 0x0UL;
	}
	u64 getNumber = getIndicated(terminals, indicator);
	assert(getNumber == terminalN);

	if (indicator != NULL) { free(indicator); }

	u64 * indicator_cir;
	indicator_cir = (u64 *)malloc((N >> 6) * sizeof(u64));
	assert(indicator_cir != NULL);

	for (u64 i = 0; i < (N >> 6); i++)
	{
		indicator_cir[i] = 0xffffffffffffffffULL;
	}

	u64 * indicator_chain;
	indicator_chain = (u64 *)malloc((N >> 6ULL) * sizeof(u64));
	assert(indicator_chain != NULL);

	u64 count = 0ULL;
	while (count < terminalN)
	{
		for (u64 i = 0; i < (N >> 6); i++)
		{
			indicator_chain[i] = 0xffffffffffffffffULL;
		}
		u64 notInChain;
		data_t y = terminals[count];

		FG[y].depth = 0ULL;
		setBit(indicator_chain, y);
		setBit(indicator_cir, y);
		tailNodeN++;
		do
		{
			x = y;
			y = FG[x].image;

			FG[y].preimageN++;

			notInChain = getBit(indicator_chain, y);
			setBit(indicator_chain, y);

			if (FG[y].preimageN > 1)
			{
				u64 dold = FG[y].depth;
				if (notInChain == 1)
				{
					u64 dnew = FG[x].depth + 1;
					if (dnew > dold)
					{
						FG[y].depth = dnew;
						data_t z = FG[y].image;
						if (FG[z].depth != INF)
						{
							dnew = FG[y].depth + 1;
							while (dnew > FG[z].depth)
							{
								FG[z].depth = dnew;
								z = FG[z].image;
								if (FG[z].depth == INF) break;
								dnew = dnew + 1;
							}
						}
					}
				}
				else
				{
					componentN++;
					u64 dnew = FG[x].depth - dold + 1;
					cyclicNodeN = cyclicNodeN + dnew;
					tailNodeN = tailNodeN - dnew;
					FG[y].depth = INF; 
					data_t z = FG[y].image;
					while (z != y)
					{
						FG[z].depth = INF;
						z = FG[z].image;
					}
				}
				break;
			}
			else
			{
				setBit(indicator_cir, y);
				tailNodeN++;
				FG[y].depth = FG[x].depth + 1;
			}
		} while (FG[y].preimageN < 2);
		count++;
	}
	free(indicator_chain);
	if (terminals != NULL) { free(terminals); }

	onlyCyclicNodeN = N - (tailNodeN + cyclicNodeN);
	assert(onlyCyclicNodeN >= 0ULL);

	if (onlyCyclicNodeN > 0ULL)
	{
		BSTree * onlyCyclicNodes = NULL;
		u64 getNumber;
		onlyCyclicNodes = getIndicated_cir(&getNumber, onlyCyclicNodes, indicator_cir);
		assert(getNumber == onlyCyclicNodeN);

		while (onlyCyclicNodes != NULL)
		{
			data_t y = onlyCyclicNodes->val;
			data_t start = y;
			do
			{
				cyclicNodeN++;
				getNumber--;
				onlyCyclicNodes = extractBST(onlyCyclicNodes, y);
				FG[y].depth = INF;
				x = y;
				y = FG[x].image;
				FG[y].preimageN++;
			} while (y != start);
			componentN++;
		}
		assert(getNumber == 0ULL);
	}
	if (indicator_cir != NULL) { free(indicator_cir); }

	for (u64 i = 0; i < N; i++)
	{
		k_thNodeN += (FG[i].depth >= (1ULL << (k))) ? 1ULL : 0ULL;
	}
	ld64 kdt = (ld64)k_thNodeN / (ld64)N;

#if (VERBOSE==1)
	cout << setprecision(20) << setiosflags(ios::left) << endl;
	cout << setw(65) << "# Components: " << componentN << endl;
	cout << setw(65) << "# Theoretical Components : " << componentN_theorem << endl;
	cout << endl;

	cout << setw(65) << "# OnlyCyclic nodes: " << onlyCyclicNodeN << endl;
	cout << endl;

	cout << setw(65) << "# Cyclic nodes: " << cyclicNodeN << endl;
	cout << setw(65) << "# Cyclic nodes / # Total nodes^(1/2) : " << (ld64)cyclicNodeN / sqrt((ld64)N) << endl;
	cout << setw(65) << "# Theoretical Cyclic nodes / # Total nodes^(1/2) : " << sqrtl(M_PI / 2.0L) << endl;
	cout << endl;

	cout << setw(65) << "# Tail nodes: " << tailNodeN << endl;
	cout << setw(65) << "# Theoretical Tail nodes: " << tailNodeN_theorem << endl;
	cout << endl;

	cout << setw(65) << "# Terminal nodes: " << terminalN << endl;
	cout << setw(65) << "# Terminal nodes / # Total nodes : " << (ld64)terminalN / (ld64)N << endl;
	cout << setw(65) << "# Theoretical Terminal nodes / # Total nodes : " << expl(-1.0L) << endl;
	cout << endl;

	cout << setw(65) << "# Image nodes: " << imageN << endl;
	cout << setw(65) << "# Image nodes / # Total nodes : " << (ld64)imageN / (ld64)N << endl;
	cout << setw(65) << "# Theoretical Image nodes / # Total nodes : " << (1.0L - expl(-1.0L)) << endl;
	cout << endl;

	fout << setprecision(20) << setiosflags(ios::left) << endl;
	fout << setw(65) << "# Components: " << componentN << endl;
	fout << setw(65) << "# Theoretical Components : " << (n >> 1) << endl;
	fout << endl;

	fout << setw(65) << "# OnlyCyclic nodes: " << onlyCyclicNodeN << endl;
	fout << endl;

	fout << setw(65) << "# Cyclic nodes: " << cyclicNodeN << endl;
	fout << setw(65) << "# Cyclic nodes / # Total nodes^(1/2) : " << (ld64)cyclicNodeN / sqrt((ld64)N) << endl;
	fout << setw(65) << "# Theoretical Cyclic nodes / # Total nodes^(1/2) : " << sqrtl(M_PI / 2.0L) << endl;
	fout << endl;

	fout << setw(65) << "# Tail nodes: " << tailNodeN << endl;
	fout << setw(65) << "# Theoretical Tail nodes: " << tailNodeN_theorem << endl;
	fout << endl;

	fout << setw(65) << "# Terminal nodes: " << terminalN << endl;
	fout << setw(65) << "# Terminal nodes / # Total nodes : " << (ld64)terminalN / (ld64)N << endl;
	fout << setw(65) << "# Theoretical Terminal nodes / # Total nodes : " << expl(-1.0L) << endl;
	fout << endl;

	fout << setw(65) << "# Image nodes: " << imageN << endl;
	fout << setw(65) << "# Image nodes / # Total nodes : " << (ld64)imageN / (ld64)N << endl;
	fout << setw(65) << "# Theoretical Image nodes / # Total nodes : " << (1.0L - expl(-1.0L)) << endl;
	fout << endl;

	cout << setw(65) << "# 2^" + to_string((k)) + "-th iterate image nodes: " << k_thNodeN << endl;
	cout << setw(65) << "# 2^" + to_string((k)) + "-th iterate nodes / # Total nodes : " << kdt << " = 2^" << logl(kdt) / logl(2.0L) << endl;
	cout << setw(65) << "# Theoretical 2^" + to_string((k)) + "-th iterate nodes / # Total nodes : " << (k_thNodeN_theorem / (ld64)N) << " = 2^" << logl((k_thNodeN_theorem / (ld64)N)) / logl(2.0L) << endl;
	cout << endl;

	fout << setw(65) << "# 2^" + to_string((k)) + "-th iterate image nodes: " << k_thNodeN << endl;
	fout << setw(65) << "# 2^" + to_string((k)) + "-th iterate nodes / # Total nodes : " << kdt << " = 2^" << logl(kdt) / logl(2.0L) << endl;
	fout << setw(65) << "# Theoretical 2^" + to_string((k)) + "-th iterate nodes / # Total nodes : " << (k_thNodeN_theorem / (ld64)N) << " = 2^" << logl((k_thNodeN_theorem / (ld64)N)) / logl(2.0L) << endl;
	fout << endl;

	t1 = clock();
	cout << "Takes time: " << (double)(t1 - t0) / ((double)CLOCKS_PER_SEC * 60.0) << " mins." << endl;
	fout << "Takes time: " << (double)(t1 - t0) / ((double)CLOCKS_PER_SEC * 60.0) << " mins." << endl;

	cout << endl;
	fout << endl;

	fout.close();
#endif

	componentNsum_thread      += componentN ;
	cyclicNodeNsum_thread     += cyclicNodeN;
	tailNodeNsum_thread       += tailNodeN  ;
	terminalNsum_thread       += terminalN  ;
	imageNsum_thread          += imageN     ;
	k_thNodeNsum_thread       += k_thNodeN  ;

	componentNsqu_thread      += powl(componentN , 2.0L);
	cyclicNodeNsqu_thread     += powl(cyclicNodeN, 2.0L);
	tailNodeNsqu_thread       += powl(tailNodeN  , 2.0L);
	terminalNsqu_thread       += powl(terminalN  , 2.0L);
	imageNsqu_thread          += powl(imageN     , 2.0L);
	k_thNodeNsqu_thread       += powl(k_thNodeN  , 2.0L);

	componentNmax_thread      = componentNmax_thread  >= componentN  ? componentNmax_thread  : componentN ;
	cyclicNodeNmax_thread     = cyclicNodeNmax_thread >= cyclicNodeN ? cyclicNodeNmax_thread : cyclicNodeN;
	tailNodeNmax_thread       = tailNodeNmax_thread   >= tailNodeN   ? tailNodeNmax_thread   : tailNodeN  ;
	terminalNmax_thread       = terminalNmax_thread   >= terminalN   ? terminalNmax_thread   : terminalN  ;
	imageNmax_thread          = imageNmax_thread      >= imageN      ? imageNmax_thread      : imageN     ;
	k_thNodeNmax_thread       = k_thNodeNmax_thread   >= k_thNodeN   ? k_thNodeNmax_thread   : k_thNodeN  ;

	componentNmin_thread      = componentNmin_thread  <= componentN  ? componentNmin_thread  : componentN ;
	cyclicNodeNmin_thread     = cyclicNodeNmin_thread <= cyclicNodeN ? cyclicNodeNmin_thread : cyclicNodeN;
	tailNodeNmin_thread       = tailNodeNmin_thread   <= tailNodeN   ? tailNodeNmin_thread   : tailNodeN  ;
	terminalNmin_thread       = terminalNmin_thread   <= terminalN   ? terminalNmin_thread   : terminalN  ;
	imageNmin_thread          = imageNmin_thread      <= imageN      ? imageNmin_thread      : imageN     ;
	k_thNodeNmin_thread       = k_thNodeNmin_thread   <= k_thNodeN   ? k_thNodeNmin_thread   : k_thNodeN  ;

	memcpy(componentN_log  [si], (ld_to_string(log2l((ld64) componentN )) + ", ").c_str(), WIDTH + 2);
	memcpy(cyclicNodeN_log [si], (ld_to_string(log2l((ld64) cyclicNodeN)) + ", ").c_str(), WIDTH + 2);
	memcpy(tailNodeN_log   [si], (ld_to_string(log2l((ld64) tailNodeN  )) + ", ").c_str(), WIDTH + 2);
	memcpy(terminalN_log   [si], (ld_to_string(log2l((ld64) terminalN  )) + ", ").c_str(), WIDTH + 2);
	memcpy(imageN_log      [si], (ld_to_string(log2l((ld64) imageN     )) + ", ").c_str(), WIDTH + 2);
	memcpy(k_thNodeN_log   [si], (ld_to_string(log2l((ld64) k_thNodeN  )) + ", ").c_str(), WIDTH + 2);

	return FG;
}

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &idx);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	string nstr = "n" + to_string(nMin) + "_n" + to_string(nMax) + "_";
	string cpstr = nstr + "componentN.txt";
	string cnstr = nstr + "cyclicNodeN.txt";
	string tnstr = nstr + "tailNodeN.txt";
	string tmstr = nstr + "terminalN.txt";
	string instr = nstr + "imageN.txt";
	string knstr = nstr + "k_thNodeN.txt";
	MPI_File_open(MPI_COMM_WORLD, cpstr.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fhcomponentN );
	MPI_File_open(MPI_COMM_WORLD, cnstr.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fhcyclicNodeN);
	MPI_File_open(MPI_COMM_WORLD, tnstr.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fhtailNodeN  );
	MPI_File_open(MPI_COMM_WORLD, tmstr.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fhterminalN  );
	MPI_File_open(MPI_COMM_WORLD, instr.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fhimageN     );
	MPI_File_open(MPI_COMM_WORLD, knstr.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fhk_thNodeN  );

	clock_t t0;
	ofstream fout;

	string pre;
	string suf = "]\n";

	ld64 componentN_theorems[nMax - nMin + 1];
	ld64 cyclicNodeN_theorems[nMax - nMin + 1];
	ld64 tailNodeN_theorems[nMax - nMin + 1];
	ld64 terminalN_theorems[nMax - nMin + 1];
	ld64 imageN_theorems[nMax - nMin + 1];
	ld64 k_thNodeN_theorems[nMax - nMin + 1];

	for (n = nMin; n <= nMax; n++)
	{
		componentN_theorem        = 0.0L;
		cyclicNodeN_theorem       = 0.0L;
		onlyCyclicNodeN_theorem   = 0.0L;
		tailNodeN_theorem         = 0.0L;
		terminalN_theorem         = 0.0L;
		imageN_theorem            = 0.0L;
		k_thNodeN_theorem         = 0.0L;
		
		componentNsum_thread      = 0.0L;
		cyclicNodeNsum_thread     = 0.0L;
		onlyCyclicNodeNsum_thread = 0.0L;
		tailNodeNsum_thread       = 0.0L;
		terminalNsum_thread       = 0.0L;
		imageNsum_thread          = 0.0L;
		k_thNodeNsum_thread       = 0.0L;
		
		componentNsqu_thread      = 0.0L;
		cyclicNodeNsqu_thread     = 0.0L;
		onlyCyclicNodeNsqu_thread = 0.0L;
		tailNodeNsqu_thread       = 0.0L;
		terminalNsqu_thread       = 0.0L;
		imageNsqu_thread          = 0.0L;
		k_thNodeNsqu_thread       = 0.0L;
		
		componentNmax_thread      = 0.0L;
		cyclicNodeNmax_thread     = 0.0L;
		onlyCyclicNodeNmax_thread = 0.0L;
		tailNodeNmax_thread       = 0.0L;
		terminalNmax_thread       = 0.0L;
		imageNmax_thread          = 0.0L;
		k_thNodeNmax_thread       = 0.0L;
			
		componentNmin_thread      = (ld64)INF;
		cyclicNodeNmin_thread     = (ld64)INF;
		onlyCyclicNodeNmin_thread = (ld64)INF;
		tailNodeNmin_thread       = (ld64)INF;
		terminalNmin_thread       = (ld64)INF;
		imageNmin_thread          = (ld64)INF;
		k_thNodeNmin_thread       = (ld64)INF;
			
		if (idx == 0)
		{
			fn = "FG_AES128_n" + to_string(n) + "_samples" + to_string(FG_SAMPLE_N) + "_statistic.txt";
			fout.open(fn.c_str(), ios::app);
			fout << "=====================================================================================" << endl;
			fout << " Functional Graphs Statistics (chopped AES-128): n = " << n << endl;
			fout << "=====================================================================================" << endl;
			cout << "=====================================================================================" << endl;
			cout << " Functional Graphs Statistics (chopped AES-128): n = " << n << endl;
			cout << "=====================================================================================" << endl;
			fout.close();

			pre = "n" + to_string(n) + " = [";
			MPI_File_write_shared(fhcomponentN , pre.c_str(), pre.size(), MPI_CHAR, MPI_STATUS_IGNORE);
			MPI_File_write_shared(fhcyclicNodeN, pre.c_str(), pre.size(), MPI_CHAR, MPI_STATUS_IGNORE);
			MPI_File_write_shared(fhtailNodeN  , pre.c_str(), pre.size(), MPI_CHAR, MPI_STATUS_IGNORE);
			MPI_File_write_shared(fhterminalN  , pre.c_str(), pre.size(), MPI_CHAR, MPI_STATUS_IGNORE);
			MPI_File_write_shared(fhimageN     , pre.c_str(), pre.size(), MPI_CHAR, MPI_STATUS_IGNORE);
			MPI_File_write_shared(fhk_thNodeN  , pre.c_str(), pre.size(), MPI_CHAR, MPI_STATUS_IGNORE);

			t0 = clock();
		
			componentNsum_total       = 0.0L;
			cyclicNodeNsum_total      = 0.0L;
			onlyCyclicNodeNsum_total  = 0.0L;
			tailNodeNsum_total        = 0.0L;
			terminalNsum_total        = 0.0L;
			imageNsum_total           = 0.0L;
			k_thNodeNsum_total        = 0.0L;
			
			componentNsqu_total       = 0.0L;
			cyclicNodeNsqu_total      = 0.0L;
			onlyCyclicNodeNsqu_total  = 0.0L;
			tailNodeNsqu_total        = 0.0L;
			terminalNsqu_total        = 0.0L;
			imageNsqu_total           = 0.0L;
			k_thNodeNsqu_total        = 0.0L;
			
			componentNmax_total       = 0.0L;
			cyclicNodeNmax_total      = 0.0L;
			onlyCyclicNodeNmax_total  = 0.0L;
			tailNodeNmax_total        = 0.0L;
			terminalNmax_total        = 0.0L;
			imageNmax_total           = 0.0L;
			k_thNodeNmax_total        = 0.0L;
			
			componentNmin_total       = (ld64)INF;
			cyclicNodeNmin_total      = (ld64)INF;
			onlyCyclicNodeNmin_total  = (ld64)INF;
			tailNodeNmin_total        = (ld64)INF;
			terminalNmin_total        = (ld64)INF;
			imageNmin_total           = (ld64)INF;
			k_thNodeNmin_total        = (ld64)INF;
		}

		ld64 t = 0.0L;
		for (u64 i = 1; i <= (1ULL << (k)); i++)
		{
			t = expl(-1.0L + t);
		}
		t = 1.0L - t;

		componentN_theorem = (ld64)n / 2.0L;
		cyclicNodeN_theorem = sqrtl((M_PI * (ld64)N) / 2.0L);
		tailNodeN_theorem = (ld64)N - cyclicNodeN_theorem;
		terminalN_theorem = expl(-1.0L) * (ld64)N;
		imageN_theorem = (1.0L - expl(-1.0L)) * (ld64)N;
		k_thNodeN_theorem = t * (ld64)N;

		if (idx == 0)
		{
			componentN_theorems[n - nMin]  = log2l(componentN_theorem);
			cyclicNodeN_theorems[n - nMin] = log2l(cyclicNodeN_theorem);
			tailNodeN_theorems[n - nMin]   = log2l(tailNodeN_theorem);
			terminalN_theorems[n - nMin]   = log2l(terminalN_theorem);
			imageN_theorems[n - nMin]      = log2l(imageN_theorem);
			k_thNodeN_theorems[n - nMin]   = log2l(k_thNodeN_theorem);
		}

		MPI_Barrier(MPI_COMM_WORLD);

		for (u64 threadSamplei = 0; threadSamplei < TaskEachThread; threadSamplei++)
		{
			info_node * FG1 = NULL;

			init_func();
			FG1 = drawFunctionalG(FG1, func1, threadSamplei);
			end_func();

			if (FG1 != NULL) { free(FG1); }
		}

		MPI_Barrier(MPI_COMM_WORLD);

		MPI_File_write_ordered(fhcomponentN , componentN_log , TaskEachThread * (WIDTH + 2), MPI_CHAR, MPI_STATUS_IGNORE);
		MPI_File_write_ordered(fhcyclicNodeN, cyclicNodeN_log, TaskEachThread * (WIDTH + 2), MPI_CHAR, MPI_STATUS_IGNORE);
		MPI_File_write_ordered(fhtailNodeN  , tailNodeN_log  , TaskEachThread * (WIDTH + 2), MPI_CHAR, MPI_STATUS_IGNORE);
		MPI_File_write_ordered(fhterminalN  , terminalN_log  , TaskEachThread * (WIDTH + 2), MPI_CHAR, MPI_STATUS_IGNORE);
		MPI_File_write_ordered(fhimageN     , imageN_log     , TaskEachThread * (WIDTH + 2), MPI_CHAR, MPI_STATUS_IGNORE);
		MPI_File_write_ordered(fhk_thNodeN  , k_thNodeN_log  , TaskEachThread * (WIDTH + 2), MPI_CHAR, MPI_STATUS_IGNORE);

		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Reduce(&componentNsum_thread     , &componentNsum_total     , 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&cyclicNodeNsum_thread    , &cyclicNodeNsum_total    , 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&tailNodeNsum_thread      , &tailNodeNsum_total      , 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&terminalNsum_thread      , &terminalNsum_total      , 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&imageNsum_thread         , &imageNsum_total         , 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&k_thNodeNsum_thread      , &k_thNodeNsum_total      , 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		MPI_Reduce(&componentNsqu_thread     , &componentNsqu_total     , 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&cyclicNodeNsqu_thread    , &cyclicNodeNsqu_total    , 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&tailNodeNsqu_thread      , &tailNodeNsqu_total      , 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&terminalNsqu_thread      , &terminalNsqu_total      , 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&imageNsqu_thread         , &imageNsqu_total         , 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		MPI_Reduce(&k_thNodeNsqu_thread      , &k_thNodeNsqu_total      , 1, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		MPI_Reduce(&componentNmax_thread     , &componentNmax_total     , 1, MPI_LONG_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(&cyclicNodeNmax_thread    , &cyclicNodeNmax_total    , 1, MPI_LONG_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(&tailNodeNmax_thread      , &tailNodeNmax_total      , 1, MPI_LONG_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(&terminalNmax_thread      , &terminalNmax_total      , 1, MPI_LONG_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(&imageNmax_thread         , &imageNmax_total         , 1, MPI_LONG_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
		MPI_Reduce(&k_thNodeNmax_thread      , &k_thNodeNmax_total      , 1, MPI_LONG_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
	
		MPI_Reduce(&componentNmin_thread     , &componentNmin_total     , 1, MPI_LONG_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
		MPI_Reduce(&cyclicNodeNmin_thread    , &cyclicNodeNmin_total    , 1, MPI_LONG_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
		MPI_Reduce(&tailNodeNmin_thread      , &tailNodeNmin_total      , 1, MPI_LONG_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
		MPI_Reduce(&terminalNmin_thread      , &terminalNmin_total      , 1, MPI_LONG_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
		MPI_Reduce(&imageNmin_thread         , &imageNmin_total         , 1, MPI_LONG_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
		MPI_Reduce(&k_thNodeNmin_thread      , &k_thNodeNmin_total      , 1, MPI_LONG_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

		if (idx == 0)
		{
			ld64 avg;
			ld64 dev;
			ld64 mav;
			ld64 miv;

			fout.open(fn.c_str(), ios::app);
			cout << setprecision(20) << setiosflags(ios::left) << endl;
			fout << setprecision(20) << setiosflags(ios::left) << endl;

			avg = componentNsum_total / (ld64)FG_SAMPLE_N;
			dev = sqrt(componentNsqu_total / (ld64)FG_SAMPLE_N - avg * avg);

			cout << setw(55) << "Theoretical  Average # Components: " << setw(30) << componentN_theorem << " = 2^" << log2l(componentN_theorem) << endl;
			cout << setw(55) << "Experimental Average # Components: " << setw(30) << avg << " = 2^" << log2l(avg) << endl;
			cout << setw(55) << "Maximum: " << setw(30) << componentNmax_total << " = 2^" << log2l(componentNmax_total) << endl;
			cout << setw(55) << "Minimum: " << setw(30) << componentNmin_total << " = 2^" << log2l(componentNmin_total) << endl;
			cout << setw(55) << "Standard Deviation: " << setw(30) << dev << " = 2^" << log2l(dev) << endl;
			cout << setw(55) << "Standard Deviation / Average: " << setw(30) << dev / avg << " = 2^" << log2l(dev/ avg) << endl;
			cout << endl;
			fout << setw(55) << "Theoretical  Average # Components: " << setw(30) << componentN_theorem << " = 2^" << log2l(componentN_theorem) << endl;
			fout << setw(55) << "Experimental Average # Components: " << setw(30) << avg << " = 2^" << log2l(avg) << endl;
			fout << setw(55) << "Maximum: " << setw(30) << componentNmax_total << " = 2^" << log2l(componentNmax_total) << endl;
			fout << setw(55) << "Minimum: " << setw(30) << componentNmin_total << " = 2^" << log2l(componentNmin_total) << endl;
			fout << setw(55) << "Standard Deviation: " << setw(30) << dev << " = 2^" << log2l(dev) << endl;
			fout << setw(55) << "Standard Deviation / Average: " << setw(30) << dev / avg << " = 2^" << log2l(dev / avg) << endl;
			fout << endl;

			avg = cyclicNodeNsum_total / (ld64)FG_SAMPLE_N;
			dev = sqrt(cyclicNodeNsqu_total / (ld64)FG_SAMPLE_N - avg * avg);

			cout << setw(55) << "Theoretical  Average # Cyclic nodes: " << setw(30) << cyclicNodeN_theorem << " = 2^" << log2l(cyclicNodeN_theorem) << endl;
			cout << setw(55) << "Experimental Average # Cyclic nodes: " << setw(30) << avg << " = 2^" << log2l(avg) << endl;
			cout << setw(55) << "Maximum: " << setw(30) << cyclicNodeNmax_total << " = 2^" << log2l(cyclicNodeNmax_total) << endl;
			cout << setw(55) << "Minimum: " << setw(30) << cyclicNodeNmin_total << " = 2^" << log2l(cyclicNodeNmin_total) << endl;
			cout << setw(55) << "Standard Deviation: " << setw(30) << dev << " = 2^" << log2l(dev) << endl;
			cout << setw(55) << "Standard Deviation / Average: " << setw(30) << dev / avg << " = 2^" << log2l(dev / avg) << endl;
			cout << endl;
			fout << setw(55) << "Theoretical  Average # Cyclic nodes: " << setw(30) << cyclicNodeN_theorem << " = 2^" << log2l(cyclicNodeN_theorem) << endl;
			fout << setw(55) << "Experimental Average # Cyclic nodes: " << setw(30) << avg << " = 2^" << log2l(avg) << endl;
			fout << setw(55) << "Maximum: " << setw(30) << cyclicNodeNmax_total << " = 2^" << log2l(cyclicNodeNmax_total) << endl;
			fout << setw(55) << "Minimum: " << setw(30) << cyclicNodeNmin_total << " = 2^" << log2l(cyclicNodeNmin_total) << endl;
			fout << setw(55) << "Standard Deviation: " << setw(30) << dev << " = 2^" << log2l(dev) << endl;
			fout << setw(55) << "Standard Deviation / Average: " << setw(30) << dev / avg << " = 2^" << log2l(dev / avg) << endl;
			fout << endl;

			avg = tailNodeNsum_total / (ld64)FG_SAMPLE_N;
			dev = sqrt(tailNodeNsqu_total / (ld64)FG_SAMPLE_N - avg * avg);

			cout << setw(55) << "Theoretical  Average # Tail nodes: " << setw(30) << tailNodeN_theorem << " = 2^" << log2l(tailNodeN_theorem) << endl;
			cout << setw(55) << "Experimental Average # Tail nodes: " << setw(30) << avg << " = 2^" << log2l(avg) << endl;
			cout << setw(55) << "Maximum: " << setw(30) << tailNodeNmax_total << " = 2^" << log2l(tailNodeNmax_total) << endl;
			cout << setw(55) << "Minimum: " << setw(30) << tailNodeNmin_total << " = 2^" << log2l(tailNodeNmin_total) << endl;
			cout << setw(55) << "Standard Deviation: " << setw(30) << dev << " = 2^" << log2l(dev) << endl;
			cout << setw(55) << "Standard Deviation / Average: " << setw(30) << dev / avg << " = 2^" << log2l(dev / avg) << endl;
			cout << endl;
			fout << setw(55) << "Theoretical  Average # Tail nodes: " << setw(30) << tailNodeN_theorem << " = 2^" << log2l(tailNodeN_theorem) << endl;
			fout << setw(55) << "Experimental Average # Tail nodes: " << setw(30) << avg << " = 2^" << log2l(avg) << endl;
			fout << setw(55) << "Maximum: " << setw(30) << tailNodeNmax_total << " = 2^" << log2l(tailNodeNmax_total) << endl;
			fout << setw(55) << "Minimum: " << setw(30) << tailNodeNmin_total << " = 2^" << log2l(tailNodeNmin_total) << endl;
			fout << setw(55) << "Standard Deviation: " << setw(30) << dev << " = 2^" << log2l(dev) << endl;
			fout << setw(55) << "Standard Deviation / Average: " << setw(30) << dev / avg << " = 2^" << log2l(dev / avg) << endl;
			fout << endl;

			avg = terminalNsum_total / (ld64)FG_SAMPLE_N;
			dev = sqrt(terminalNsqu_total / (ld64)FG_SAMPLE_N - avg * avg);

			cout << setw(55) << "Theoretical  Average # Terminal nodes: " << setw(30) << terminalN_theorem << " = 2^" << log2l(terminalN_theorem) << endl;
			cout << setw(55) << "Experimental Average # Terminal nodes: " << setw(30) << avg << " = 2^" << log2l(avg) << endl;
			cout << setw(55) << "Maximum: " << setw(30) << terminalNmax_total << " = 2^" << log2l(terminalNmax_total) << endl;
			cout << setw(55) << "Minimum: " << setw(30) << terminalNmin_total << " = 2^" << log2l(terminalNmin_total) << endl;
			cout << setw(55) << "Standard Deviation: " << setw(30) << dev << " = 2^" << log2l(dev) << endl;
			cout << setw(55) << "Standard Deviation / Average: " << setw(30) << dev / avg << " = 2^" << log2l(dev / avg) << endl;
			cout << endl;
			fout << setw(55) << "Theoretical  Average # Terminal nodes: " << setw(30) << terminalN_theorem << " = 2^" << log2l(terminalN_theorem) << endl;
			fout << setw(55) << "Experimental Average # Terminal nodes: " << setw(30) << avg << " = 2^" << log2l(avg) << endl;
			fout << setw(55) << "Maximum: " << setw(30) << terminalNmax_total << " = 2^" << log2l(terminalNmax_total) << endl;
			fout << setw(55) << "Minimum: " << setw(30) << terminalNmin_total << " = 2^" << log2l(terminalNmin_total) << endl;
			fout << setw(55) << "Standard Deviation: " << setw(30) << dev << " = 2^" << log2l(dev) << endl;
			fout << setw(55) << "Standard Deviation / Average: " << setw(30) << dev / avg << " = 2^" << log2l(dev / avg) << endl;
			fout << endl;

			avg = imageNsum_total / (ld64)FG_SAMPLE_N;
			dev = sqrt(imageNsqu_total / (ld64)FG_SAMPLE_N - avg * avg);

			cout << setw(55) << "Theoretical  Average # Image nodes: " << setw(30) << imageN_theorem << " = 2^" << log2l(imageN_theorem) << endl;
			cout << setw(55) << "Experimental Average # Image nodes: " << setw(30) << avg << " = 2^" << log2l(avg)<< endl;
			cout << setw(55) << "Maximum: " << setw(30) << imageNmax_total << " = 2^" << log2l(imageNmax_total) << endl;
			cout << setw(55) << "Minimum: " << setw(30) << imageNmin_total << " = 2^" << log2l(imageNmin_total) << endl;
			cout << setw(55) << "Standard Deviation: " << setw(30) << dev << " = 2^" << log2l(dev) << endl;
			cout << setw(55) << "Standard Deviation / Average: " << setw(30) << dev / avg << " = 2^" << log2l(dev / avg) << endl;
			cout << endl;
			fout << setw(55) << "Theoretical  Average # Image nodes: " << setw(30) << imageN_theorem << " = 2^" << log2l(imageN_theorem) << endl;
			fout << setw(55) << "Experimental Average # Image nodes: " << setw(30) << avg << " = 2^" << log2l(avg) << endl;
			fout << setw(55) << "Maximum: " << setw(30) << imageNmax_total << " = 2^" << log2l(imageNmax_total) << endl;
			fout << setw(55) << "Minimum: " << setw(30) << imageNmin_total << " = 2^" << log2l(imageNmin_total) << endl;
			fout << setw(55) << "Standard Deviation: " << setw(30) << dev << " = 2^" << log2l(dev) << endl;
			fout << setw(55) << "Standard Deviation / Average: " << setw(30) << dev / avg << " = 2^" << log2l(dev / avg) << endl;
			fout << endl;

			avg = k_thNodeNsum_total / (ld64)FG_SAMPLE_N;
			dev = sqrt(k_thNodeNsqu_total / (ld64)FG_SAMPLE_N - avg * avg);

			cout << setw(55) << "Theoretical  Average # 2^" + to_string((k)) + "-th iterates image nodes: " << setw(30) << k_thNodeN_theorem << " = 2^" << log2l(k_thNodeN_theorem) << endl;
			cout << setw(55) << "Experimental Average # 2^" + to_string((k)) + "-th iterates image nodes: " << setw(30) << avg << " = 2^" << log2l(avg) << endl;
			cout << setw(55) << "Maximum: " << setw(30) << k_thNodeNmax_total << " = 2^" << log2l(k_thNodeNmax_total) << endl;
			cout << setw(55) << "Minimum: " << setw(30) << k_thNodeNmin_total << " = 2^" << log2l(k_thNodeNmin_total) << endl;
			cout << setw(55) << "Standard Deviation: " << setw(30) << dev << " = 2^" << log2l(dev) << endl;
			cout << setw(55) << "Standard Deviation / Average: " << setw(30) << dev / avg << " = 2^" << log2l(dev / avg) << endl;
			cout << endl;
			fout << setw(55) << "Theoretical  Average # 2^" + to_string((k)) + "-th iterates image nodes: " << setw(30) << k_thNodeN_theorem << " = 2^" << log2l(k_thNodeN_theorem) << endl;
			fout << setw(55) << "Experimental Average # 2^" + to_string((k)) + "-th iterates image nodes: " << setw(30) << avg << " = 2^" << log2l(avg) << endl;
			fout << setw(55) << "Maximum: " << setw(30) << k_thNodeNmax_total << " = 2^" << log2l(k_thNodeNmax_total) << endl;
			fout << setw(55) << "Minimum: " << setw(30) << k_thNodeNmin_total << " = 2^" << log2l(k_thNodeNmin_total) << endl;
			fout << setw(55) << "Standard Deviation: " << setw(30) << dev << " = 2^" << log2l(dev) << endl;
			fout << setw(55) << "Standard Deviation / Average: " << setw(30) << dev / avg << " = 2^" << log2l(dev / avg) << endl;
			fout << endl;

			t0 = clock() - t0;
			cout << "Takes time: " << (double)(t0) / ((double)CLOCKS_PER_SEC * 60.0) << " mins." << endl;
			fout << "Takes time: " << (double)(t0) / ((double)CLOCKS_PER_SEC * 60.0) << " mins." << endl;
			cout << endl;
			fout << endl;
			fout.close();

			MPI_File_write_shared(fhcomponentN , suf.c_str(), suf.size(), MPI_CHAR, MPI_STATUS_IGNORE);
			MPI_File_write_shared(fhcyclicNodeN, suf.c_str(), suf.size(), MPI_CHAR, MPI_STATUS_IGNORE);
			MPI_File_write_shared(fhtailNodeN  , suf.c_str(), suf.size(), MPI_CHAR, MPI_STATUS_IGNORE);
			MPI_File_write_shared(fhterminalN  , suf.c_str(), suf.size(), MPI_CHAR, MPI_STATUS_IGNORE);
			MPI_File_write_shared(fhimageN     , suf.c_str(), suf.size(), MPI_CHAR, MPI_STATUS_IGNORE);
			MPI_File_write_shared(fhk_thNodeN  , suf.c_str(), suf.size(), MPI_CHAR, MPI_STATUS_IGNORE);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}

	if (idx == 0)
	{
		string nt = "xt = [";
		string co = "yt = [";
		string cy = "yt = [";
		string ta = "yt = [";
		string te = "yt = [";
		string im = "yt = [";
		string kt = "yt = [";

		for (n = nMin; n <= nMax; n++)
		{
			nt += to_string(n - nMin + 1) + ", ";
			co += ld_to_string(componentN_theorems[n - nMin] ) + ", ";
			cy += ld_to_string(cyclicNodeN_theorems[n - nMin]) + ", ";
			ta += ld_to_string(tailNodeN_theorems[n - nMin]  ) + ", ";
			te += ld_to_string(terminalN_theorems[n - nMin]  ) + ", ";
			im += ld_to_string(imageN_theorems[n - nMin]     ) + ", ";
			kt += ld_to_string(k_thNodeN_theorems[n - nMin]  ) + ", ";
		}
		nt += "]\n";
		co += "]\n" + nt;
		cy += "]\n" + nt;
		ta += "]\n" + nt;
		te += "]\n" + nt;
		im += "]\n" + nt;
		kt += "]\n" + nt;

		MPI_File_write_shared(fhcomponentN , co.c_str(), co.size(), MPI_CHAR, MPI_STATUS_IGNORE);
		MPI_File_write_shared(fhcyclicNodeN, cy.c_str(), cy.size(), MPI_CHAR, MPI_STATUS_IGNORE);
		MPI_File_write_shared(fhtailNodeN  , ta.c_str(), ta.size(), MPI_CHAR, MPI_STATUS_IGNORE);
		MPI_File_write_shared(fhterminalN  , te.c_str(), te.size(), MPI_CHAR, MPI_STATUS_IGNORE);
		MPI_File_write_shared(fhimageN     , im.c_str(), im.size(), MPI_CHAR, MPI_STATUS_IGNORE);
		MPI_File_write_shared(fhk_thNodeN  , kt.c_str(), kt.size(), MPI_CHAR, MPI_STATUS_IGNORE);
	}

	MPI_File_close(&fhcomponentN);
	MPI_File_close(&fhcyclicNodeN);
	MPI_File_close(&fhtailNodeN);
	MPI_File_close(&fhterminalN);
	MPI_File_close(&fhimageN);
	MPI_File_close(&fhk_thNodeN);

	MPI_Finalize();

#if defined(_MSC_VER)
	system("Pause");
#endif
	return 0;
}

