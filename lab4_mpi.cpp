#include "lab4_mpi.h"

#include <malloc.h>
#include <math.h>
#include "mpi.h"

#define MASTER 0

using namespace std;

void calcPi(int* pi_set, int* m_set, int* p_set, int num_patterns)
{
	for(int i = 0; i < num_patterns; i++)
		pi_set[i] = min(p_set[i], (int)ceil(float(m_set[i])/2));
}

int pi(int period, int m){
	return min(period, (int)ceil(float(m)/2));
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
		witness[i] = phi;
	}
}

void witn(int** witness, int pi, char* pattern)
{
	int* phi = new int[pi];
	phi[0] = 0; int j, k;
	for(j = 1; j < pi; j++){
		for(k = 0; k < pi - j; k++){
			if(pattern[k] != pattern[j+k])
				break;
		}
		phi[j] = k;
	}
	*witness = phi;
}

int duel(char* Z, int n, char* pattern, int* phi, int i, int j)
{
	int k = phi[j-i];
	if((j+k < 0 || j+k > n) || Z[j+k] != pattern[k])
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
	int b = n / ceil(float(m)/2);
	int bi, first, last, j, i;
	int* potential_pos = new int[b];
	int* matched_pos_temp = new int[b];
	int count = 0;

	cout << "Text=" << T << ", p=" << p << endl;

	// Parallelise across procs
	for(bi = 0; bi < b; bi++)
	{
		first = n*bi/b; last = n*(bi+1)/b;
		i = first;
		for(j = i+1; j < last; j++){
			i = duel(T, n, p, phi, i, j);
			cout << "duel(" << i << "," << j << ") = " << i << endl;
		}
		potential_pos[bi] = i;
	}

	for(int i = 0; i < n; i++)
		cout << potential_pos[i] << " ";
	cout << endl;

	for(i = 0; i < b; i++){
		if(match(T, potential_pos[i], p, m)){
			matched_pos_temp[count] = potential_pos[i];
			count++;
		}
	}

	*match_count = count;
	*matched_pos = matched_pos_temp;
}

void pTextAnalysis(char* T, int n, char* p, int m, int period, int procs)
{
	char* pPrime = new char[2*period-1];
	for(int i = 0; i < 2*period-1; i++)
		pPrime[i] = p[i];
	
	cout << p << " " << pPrime << " ";
	for(int i = 0; i < 2*period-1; i++)
		cout << pPrime[i];
	cout << endl;
	int ppi = pi(period, m);
	int* witness = new int[ppi];
	witn(&witness, ppi, pPrime);

	cout << witness[0] << " " << witness[1] << " - " << ppi << endl;
	int count, *pos;

	npTextAnalysis(T, n, pPrime, 2*period-1, witness, procs, &count, &pos);

	for(int i = 0; i < count; i++)
		cout << pos[i] << " ";
	cout << "- " << count << endl;

	char* u = new char[period];
	for(int i = 0; i < period; i++)
		u[i] = p[i];

	int k = floor(float(m)/period);
	char* v = new char[m - k*period];
	for(int i = 0; i < (m-k*period); i++)
		v[i] = p[i];

	char* u2v = new char[2*period + (m - k*period)];
	for(int i = 0; i < 2*period; i++)
		u2v[i] = u[i%period];
	for(int i = 0; i < m-k*period; i++)
		u2v[i+2*period] = v[i];

	int* M = new int[n];
	for(int i = 0; i < n; i++)
		M[i] = 0;
		
	int i;
	for(int j = 0; j < count; j++){
		i = pos[j];
		if(match(T, i, u2v, 2*period+(m-k*period)))
			M[i] = 1;
	}

	int* S = new int[n/period];
	int** C = new int*[period];
	for(int i = 0; i < n; i+=period)
		S[i/period] = M[i];

	for(int i=0; i < period; i++){
		C[i] = new int[n/period];
		for(int j = 0; j < n/period; j++){
			C[i][j] = 1;
			for(int p = j; p < j+k; p++){
				if(S[p] == 0){
					C[i][j] = 0;
					break;
				}
			}
		}
	}

	int* Match = new int[n-m];
	for(int i = 0; i < period; i++){
		for(int l = 0; l < n/period; l++){
			if(i+l*period < n-m+1)
				Match[i+l*period] = C[i][l];
		}
	}

	for(int i = 0; i < n-m+1; i++)
		cout << "M[" << i << "] = " << M[i] << endl;
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
	// std::cout << "first = " << first << " last = " << last-1 << std::endl;

	int* counts = new int[num_patterns];
	int* pi_set = new int[num_patterns];
	int** witness = new int*[num_patterns];

	// calcPi(pi_set, m_set, p_set, num_patterns);
	// calcWitness(witness, pi_set, pattern_set, num_patterns);
	
	int i, j, k; bool matched;

	for(int i = 0; i < num_patterns; i++){
		cout << "n = " << n << " m = " << m_set[i] << " period = " << p_set[i] << endl;
		pTextAnalysis(text, n, pattern_set[i], m_set[i], p_set[i], procs);
	}
	
	MPI_Finalize();
	*match_counts = counts;
}
