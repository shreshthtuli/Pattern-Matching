#include "lab4_mpi.h"

#include <malloc.h>
#include "mpi.h"

#define MASTER 0

using namespace std;

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

	MPI_Comm_size(MPI_COMM_WORLD, &procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// std::cout << "Numprocs = " << procs << " rank = " << rank << std::endl;

	int first = num_patterns*rank/procs;
	int last = num_patterns*(rank+1)/procs;

	// std::cout << "first = " << first << " last = " << last-1 << std::endl;

	int* counts = new int[num_patterns];
	int* pos_temp = new int[num_patterns*n];
	int** pos = new int*[num_patterns];

	int i, j, k; bool matched;

	int idx = 0;
	// Loop over all patterns
	for(i = first; i < last; i++){
		counts[i] = 0;
		// Loop over text
		for(j = 0; j < n; j++){
			// Check pattern i and position j of text
			matched = true;
			for(k = 0; k < m_set[i]; k++){
				if(pattern_set[i][k] != text[j+k]){
					matched = false; break;
				}
			}
			if(!matched)
				continue;
			counts[i] += 1;
			pos_temp[idx] = j;
			idx++;
		}
	}

	idx = 0;
	for(int i = first; i < last; i++){
		pos[i] = new int[counts[i]];
		for(int j = 0; j < counts[i]; j++){
			pos[i][j] = pos_temp[idx];
			idx++;
		}
	}

	// std::cout <<"Rank = " << rank << std::endl;

	if(rank == MASTER){
		MPI_Status s;
		for(int i = 1; i < procs; i++){
			first = (num_patterns*i/procs); last = num_patterns*(i+1)/procs;
			MPI_Recv(counts+first, (last-first), MPI_INT, i, 0, MPI_COMM_WORLD, &s);
			// cout << "Recv from proc " << i << endl;
			for(int j = first; j < last; j++){
				pos[j] = new int[counts[j]];
				MPI_Recv(pos[j], counts[j], MPI_INT, i, 0, MPI_COMM_WORLD, &s);
			}
		}
		// for(int i = 0; i < num_patterns; i++)
		// 	std::cout << "Count of pattern " << i << " is = " << counts[i] << std::endl;
		int* matched_pos_res = new int[num_patterns*n]; int p = 0;
		for(int i = 0; i < num_patterns; i++){
			for(int j = 0; j < counts[i]; j++){
				matched_pos_res[p] = pos[i][j];
				p++;
			}
		}
		*match_counts = counts;
		*matches = matched_pos_res;
		MPI_Barrier(MPI_COMM_WORLD);
	}
	else{
		MPI_Send(counts+first, (last-first), MPI_INT, MASTER, 0, MPI_COMM_WORLD);
		for(int j = first; j < last; j++)
			MPI_Send(pos[j], counts[j], MPI_INT, MASTER, 0, MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
}
