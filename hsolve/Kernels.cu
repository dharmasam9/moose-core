#include "HSolveActive.h"

__global__
void get_lookup_rows_and_fractions_cuda(
		double* lookups,
		double* table,
		double min, double max, double dx,
		int* rows, double* fracs,
		unsigned int nColumns, unsigned int size){

	int tid = threadIdx.x + blockIdx.x * blockDim.x;

	if(tid < size){
		double x = lookups[tid];

		if ( x < min )
			x = min;
		else if ( x > max )
			x = max;

		double div = ( x - min ) / dx;
		unsigned int integer = ( unsigned int )( div );

		rows[tid] = integer*nColumns;
		fracs[tid] = div-integer;
	}
}

void HSolveActive::get_lookup_rows_and_fractions_cuda_wrapper(int gpu_load_count){
	int num_comps = V_.size();

	int THREADS_PER_BLOCK = 512;
	int BLOCKS = gpu_load_count/THREADS_PER_BLOCK;
	BLOCKS = (gpu_load_count + THREADS_PER_BLOCK-1)/THREADS_PER_BLOCK;


	// Getting lookup metadata for Vm
	get_lookup_rows_and_fractions_cuda<<<BLOCKS,THREADS_PER_BLOCK>>>(u_V,
    		d_V_table,
    		vTable_.get_min(), vTable_.get_max(), vTable_.get_dx(),
    		u_V_rows, u_V_fracs,
    		vTable_.get_num_of_columns(), gpu_load_count);
}

/*
 * Based on the near lookup value and fraction value, the function
 * interpolates the value and uses it to update appropriate state variables.
 * "indices" array is a subset of compartment id's which are
 * voltage dependent gate indices or Calcium dependent gate indices
 */
__global__
void advance_channels_opt_cuda(
		int* rows,
		double* fracs,
		double* table,
		int* indices,
		int* gate_to_comp,
		double* gate_values,
		int* gate_columns,
		int* state2chanId,
		int* chan_instants,
		unsigned int nColumns,
		double dt,
		int size
		){
	int tid = threadIdx.x + blockIdx.x * blockDim.x;
	if(tid < size){
		double a,b,C1,C2;
		int index, lookup_index, row_start_index, column;

		index = indices[tid];
		lookup_index = gate_to_comp[tid];
		row_start_index = rows[lookup_index];
		column = gate_columns[index];

		a = table[row_start_index + column];
		b = table[row_start_index + column + nColumns];

		C1 = a + (b-a)*fracs[lookup_index];

		a = table[row_start_index + column + 1];
		b = table[row_start_index + column + 1 + nColumns];

		C2 = a + (b-a)*fracs[lookup_index];

		if(!chan_instants[state2chanId[tid]]){
			a = 1.0 + dt/2.0 * C2; // reusing a
			gate_values[index] = ( gate_values[index] * ( 2.0 - a ) + dt * C1 ) / a;
		}
		else{
			gate_values[index] = C1/C2;
		}
	}
}

/*
 * Advance Channels performing both lookup and state update
 */
__global__
void advance_channels_for_externalCalcium(
		double* d_externalCalcium,
		int* d_exCalgate_indices,
		int* state2chanId,
		double* gate_values,
		int* state2Column,
		int* chan_instants,
		double* table,
		double min, double max, double dx,
		unsigned int nColumns, double dt, unsigned int size){

	int tid = threadIdx.x + blockIdx.x * blockDim.x;

	if(tid < size){
		int index = d_exCalgate_indices[tid]; // Index in the state_ array
		int chan_id = state2chanId[index];
		// Update state only if there is a contribution
		if(d_externalCalcium[chan_id] != 0){
			double a,b,C1,C2;
			double x = d_externalCalcium[chan_id];

			if ( x < min )
				x = min;
			else if ( x > max )
				x = max;

			double div = ( x - min ) / dx;
			unsigned int integer = ( unsigned int )( div );

			int row_start_index = integer*nColumns;
			double frac = div-integer;

			// Perform the update
			int column = state2Column[index];

			a = table[row_start_index + column];
			b = table[row_start_index + column + nColumns];

			C1 = a + (b-a)*frac;

			a = table[row_start_index + column + 1];
			b = table[row_start_index + column + 1 + nColumns];

			C2 = a + (b-a)*frac;

			if(!chan_instants[chan_id]){
				a = 1.0 + dt/2.0 * C2; // reusing a
				gate_values[index] = ( gate_values[index] * ( 2.0 - a ) + dt * C1 ) / a;
			}
			else{
				gate_values[index] = C1/C2;
			}

		}
	}
}

