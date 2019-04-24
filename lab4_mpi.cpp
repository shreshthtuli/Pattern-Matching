#include "lab4_mpi.h"

#include <malloc.h>
#include <math.h>
#include "mpi.h"

#define MASTER 0

void calcPi(int* pi_set, int* m_set, int* p_set, int num_patterns)
{
	for(int i = 0; i < num_patterns; i++)
		pi_set[i] = min(p_set[i], ceil(float(m_set[i])/2));
}

void pi(int period, int m){
	return min(period, ceil(float(m)/2));
}

void calcWitness(int** witness, int* pi_set, char** pattern_set, int num_patterns)
{
	int i, j, k;
	for(i = 0; i < num_patterns; i++){
		int* phi = new int[pi_set[i]];
		phi[0] = 0;
		for(j = 1; j < pi_set[i]; j++){
			for(k = 0; k < pi_set[i] - j; k++){
				if(pattern_set[i][k] != pattern_set[i][j+k])
					break;
			}
			phi[j] = k;
		}
		witness[i] = *phi;
	}
}

void witn(int* witness, int pi, char* pattern)
{
	int* phi = new int[pi];
	phi[0] = 0;
	for(j = 1; j < pi; j++){
		for(k = 0; k < pi - j; k++){
			if(pattern[k] != pattern[j+k])
				break;
		}
		phi[j] = k;
	}
	witness = *phi;
}

int duel(char* Z, int n; char* pattern, int* phi, int i, int j)
{
	int k = phi[j-i];
	if((j+k < 0 || j+k > n) || Z[j+k] != Y[k])
		return i;
	return j;
}

bool match(char* T, int i, char* p, int m)
{
	for(int j = 0; j < m; j++){
		if(T[i+j] != p[j])
			return false;
	}
	return true;
}

void npTextAnalysis(char* T, int n, char* p, int m, int* phi, int procs, int* match_count, int** matched_pos)
{
	int blocks = 12 / ceil(double(m)/2);
	int bi, first, last, j, i;
	int* potential_pos = new int[blocks];
	int* matched_pos_temp = new int[blocks];
	int count = 0;

	// Parallelise across procs
	for(bi = 0; bi < b; bi++)
	{
		first = n*bi/b; last = n*(bi+1)/b;
		i = first;
		for(j = i+1; j < last; j++){
			i = duel(T, p, phi, i, j);
		}
		potential_pos[bi] = i;
	}

	
	for(i = 0; i  blocks; i++){
		if(matched(T, i, p, m)){
			matched_pos[count] = i;
			count++;
		}
	}

	*match_counts = count;
	*matched_pos = matched_pos_temp;
}

void pTextAnalysis(char* T, int n, char* p, int m, int period, int procs)
{
	char* pPrime = new char[2*p-1];
	for(int i = 0; i < 2*p-1; i++)
		pPrime[i] = p[i];
	
	int ppi = pi(period, m);
	int* witness = new int[ppi];
	witn(witness, ppi, p);
	int count, *pos;

	npTextAnalysis(T, n, p, m, witness, procs, &count, &pos);

	char* u = new char[period];
	for(int i = 0; i < p; i++)
		u[i] = p[i];

	int k = floor(float(m)/period);
	char* v = new char[m - k*period];
	int i;
	for(int j = 0; j < count; j++){
		i = pos[j];
		if()
	}
}

void periodic_pattern_matching (
		int n, 
		char *text, 
		int num_patterns, 
		int *m_set, 
		int *p_set, 
		char **pattern_set, 
		int **match_counts, 
		int **matches)
{
	int procs, rank;
	
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	std::cout << "Numprocs = " << procs << " rank = " << rank << std::endl;
	std::cout << "first = " << first << " last = " << last-1 << std::endl;

	int* counts = new int[num_patterns];
	int* pi_set = new int[num_patterns];
	int** witness = new int*[num_patterns];

	calcPi(pi_set, m_set, p_set, num_patterns);
	calcWitness(witness, pi_set, pattern_set, num_patterns);
	
	int i, j, k; bool matched;
	
	MPI_Finalize();
	*match_counts = counts;
}
