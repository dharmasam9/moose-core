/**********************************************************************
** This program is part of 'MOOSE', the
** Messaging Object Oriented Simulation Environment.
**   copyright (C) 2003-2007 Upinder S. Bhalla, Niraj Dudani and NCBS
** It is made available under the terms of the
** GNU Lesser General Public License version 2.1
** See the file COPYING.LIB for the full notice.
**********************************************************************/

#include "header.h"
#include <queue>
#include "HSolveStruct.h"
#include "HinesMatrix.h"
#include "HSolvePassive.h"
#include "RateLookup.h"
#include "HSolveActive.h"
#include "HSolve.h"
#include "../biophysics/CompartmentBase.h"
#include "../biophysics/Compartment.h"
#include "../biophysics/CaConcBase.h"
#include "../biophysics/ChanBase.h"
#include "ZombieCaConc.h"
using namespace moose;
//~ #include "ZombieCompartment.h"
//~ #include "ZombieCaConc.h"

#include <omp.h>

extern ostream& operator <<( ostream& s, const HinesMatrix& m );

const int HSolveActive::INSTANT_X = 1;
const int HSolveActive::INSTANT_Y = 2;
const int HSolveActive::INSTANT_Z = 4;

HSolveActive::HSolveActive()
{
    caAdvance_ = 1;

    // Default lookup table size
    //~ vDiv_ = 3000;    // for voltage
    //~ caDiv_ = 3000;   // for calcium
}

//////////////////////////////////////////////////////////////////////
// Solving differential equations
//////////////////////////////////////////////////////////////////////
void HSolveActive::step( ProcPtr info )
{
    if ( nCompt_ <= 0 )
        return;

    if ( !current_.size() )
    {
        current_.resize( channel_.size() );
    }

    advanceChannels( info->dt );
    calculateChannelCurrents();
    updateMatrix();
    HSolvePassive::forwardEliminate();
    HSolvePassive::backwardSubstitute();
    advanceCalcium();
    advanceSynChans( info );

    sendValues( info );
    sendSpikes( info );

    externalCurrent_.assign( externalCurrent_.size(), 0.0 );
}

void HSolveActive::calculateChannelCurrents()
{
    vector< ChannelStruct >::iterator ichan;
    vector< CurrentStruct >::iterator icurrent = current_.begin();

    if ( state_.size() != 0 )
    {
        double* istate = &state_[ 0 ];

        for ( ichan = channel_.begin(); ichan != channel_.end(); ++ichan )
        {
            ichan->process( istate, *icurrent );
            ++icurrent;
        }
    }
}

void HSolveActive::updateMatrix()
{
    /*
     * Copy contents of HJCopy_ into HJ_. Cannot do a vector assign() because
     * iterators to HJ_ get invalidated in MS VC++
     */
    if ( HJ_.size() != 0 )
        memcpy( &HJ_[ 0 ], &HJCopy_[ 0 ], sizeof( double ) * HJ_.size() );

    double GkSum, GkEkSum; vector< CurrentStruct >::iterator icurrent = current_.begin();
    vector< currentVecIter >::iterator iboundary = currentBoundary_.begin();
    vector< double >::iterator ihs = HS_.begin();
    vector< double >::iterator iv = V_.begin();

    vector< CompartmentStruct >::iterator ic;
    for ( ic = compartment_.begin(); ic != compartment_.end(); ++ic )
    {
        GkSum   = 0.0;
        GkEkSum = 0.0;
        for ( ; icurrent < *iboundary; ++icurrent )
        {
            GkSum   += icurrent->Gk;
            GkEkSum += icurrent->Gk * icurrent->Ek;
        }

        *ihs = *( 2 + ihs ) + GkSum;
        *( 3 + ihs ) = *iv * ic->CmByDt + ic->EmByRm + GkEkSum;

        ++iboundary, ihs += 4, ++iv;
    }

    map< unsigned int, InjectStruct >::iterator inject;
    for ( inject = inject_.begin(); inject != inject_.end(); ++inject )
    {
        unsigned int ic = inject->first;
        InjectStruct& value = inject->second;

        HS_[ 4 * ic + 3 ] += value.injectVarying + value.injectBasal;

        value.injectVarying = 0.0;
    }

    // Synapses are being handled as external channels.
    //~ double Gk, Ek;
    //~ vector< SynChanStruct >::iterator isyn;
    //~ for ( isyn = synchan_.begin(); isyn != synchan_.end(); ++isyn ) {
    //~ get< double >( isyn->elm_, synGkFinfo, Gk );
    //~ get< double >( isyn->elm_, synEkFinfo, Ek );
    //~
    //~ unsigned int ic = isyn->compt_;
    //~ HS_[ 4 * ic ] += Gk;
    //~ HS_[ 4 * ic + 3 ] += Gk * Ek;
    //~ }

    ihs = HS_.begin();
    vector< double >::iterator iec;
    for ( iec = externalCurrent_.begin(); iec != externalCurrent_.end(); iec += 2 )
    {
        *ihs += *iec;
        *( 3 + ihs ) += *( iec + 1 );

        ihs += 4;
    }

    stage_ = 0;    // Update done.
}

void HSolveActive::advanceCalcium()
{
    vector< double* >::iterator icatarget = caTarget_.begin();
    vector< double >::iterator ivmid = VMid_.begin();
    vector< CurrentStruct >::iterator icurrent = current_.begin();
    vector< currentVecIter >::iterator iboundary = currentBoundary_.begin();

    /*
     * caAdvance_: This flag determines how current flowing into a calcium pool
     * is computed. A value of 0 means that the membrane potential at the
     * beginning of the time-step is used for the calculation. This is how
     * GENESIS does its computations. A value of 1 means the membrane potential
     * at the middle of the time-step is used. This is the correct way of
     * integration, and is the default way.
     */
    if ( caAdvance_ == 1 )
    {
        for ( ; iboundary != currentBoundary_.end(); ++iboundary )
        {
            for ( ; icurrent < *iboundary; ++icurrent )
            {
                if ( *icatarget )
                    **icatarget += icurrent->Gk * ( icurrent->Ek - *ivmid );

                ++icatarget;
            }

            ++ivmid;
        }
    }
    else if ( caAdvance_ == 0 )
    {
        vector< double >::iterator iv = V_.begin();
        double v0;

        for ( ; iboundary != currentBoundary_.end(); ++iboundary )
        {
            for ( ; icurrent < *iboundary; ++icurrent )
            {
                if ( *icatarget )
                {
                    v0 = ( 2 * *ivmid - *iv );

                    **icatarget += icurrent->Gk * ( icurrent->Ek - v0 );
                }

                ++icatarget;
            }

            ++ivmid, ++iv;
        }
    }

    vector< CaConcStruct >::iterator icaconc;
    vector< double >::iterator icaactivation = caActivation_.begin();
    vector< double >::iterator ica = ca_.begin();
    for ( icaconc = caConc_.begin(); icaconc != caConc_.end(); ++icaconc )
    {
        *ica = icaconc->process( *icaactivation );
        ++ica, ++icaactivation;
    }

    caActivation_.assign( caActivation_.size(), 0.0 );
}

void HSolveActive::advanceChannels( double dt )
{
    // Required variables
    vector<double> table = vTable_.get_table();
    double min = vTable_.get_min();
    double max = vTable_.get_max();
    double dx = vTable_.get_dx();
    unsigned int nColumns = vTable_.get_num_of_columns();

    // Looking up values for Voltages
    for(unsigned int tid = 0; tid < V_.size(); ++tid)
    {
        double x = V_[tid];

        if ( x < min )
            x = min;
        else if ( x > max )
            x = max;

        double div = ( x - min ) / dx;
        unsigned int integer = ( unsigned int )( div );

        h_V_rows[tid] = integer*nColumns;
        h_V_fracs[tid] = div-integer;
    }

    // Updating state variable of Voltage dependent gates
    for(unsigned int tid = 0; tid < h_vgate_indices.size(); ++tid) {
        double a,b,C1,C2;
        int index, lookup_index, row_start_index, column;

        index = h_vgate_indices[tid];
        lookup_index = h_vgate_compIds[tid];
        row_start_index = h_V_rows[lookup_index];
        column = h_state2column[index];

        a = table[row_start_index + column];
        b = table[row_start_index + column + nColumns];

        C1 = a + (b-a)*h_V_fracs[lookup_index];

        a = table[row_start_index + column + 1];
        b = table[row_start_index + column + 1 + nColumns];

        C2 = a + (b-a)*h_V_fracs[lookup_index];

        if(!h_chan_instant[h_state2chanId[index]]){
            a = 1.0 + dt/2.0 * C2; // reusing a
            state_[index] = ( state_[index] * ( 2.0 - a ) + dt * C1 ) / a;
        }
        else{
            state_[index] = C1/C2;
        }
    }

	if(h_cagate_indices.size() > 0){
        cout << "Number of Calcium dep gates " << h_cagate_indices.size() << endl;
        for (unsigned int tid = 0; tid < h_cagate_indices.size(); ++tid) {
                int index, lookup_index, column;

                index = h_cagate_indices[tid];
                lookup_index = h_cagate_capoolIds[tid];
                column = h_state2column[index];

                LookupRow caRow;
                LookupColumn caCol;
                caTable_.row(ca_[lookup_index], caRow);
                caCol.column = column;
                double C1,C2;

                caTable_.lookup( caCol, caRow, C1, C2 );

                if(!h_chan_instant[h_state2chanId[index]]){
                    double temp  = 1.0 + dt/2.0 * C2; // reusing a
                    state_[index] = ( state_[index] * ( 2.0 - temp ) + dt * C1 ) / temp;
                }else{
                    state_[index] = C1/C2;
                }
            }
    }
	


#if 0
    vector< double >::iterator iv;
    vector< double >::iterator istate = state_.begin();
    vector< int >::iterator ichannelcount = channelCount_.begin();
    vector< ChannelStruct >::iterator ichan = channel_.begin();
    vector< ChannelStruct >::iterator chanBoundary;
    vector< unsigned int >::iterator icacount = caCount_.begin();
    vector< double >::iterator ica = ca_.begin();
    vector< double >::iterator caBoundary;
    vector< LookupColumn >::iterator icolumn = column_.begin();
    vector< LookupRow >::iterator icarowcompt;
    vector< LookupRow* >::iterator icarow = caRow_.begin();
    vector< double >::iterator iextca = externalCalcium_.begin();

    LookupRow vRow;
    LookupRow dRow;
    double C1, C2;

    for ( iv = V_.begin(); iv != V_.end(); ++iv )
    {
        vTable_.row( *iv, vRow );
        icarowcompt = caRowCompt_.begin();
        caBoundary = ica + *icacount;
        for ( ; ica < caBoundary; ++ica )
        {
            caTable_.row( *ica, *icarowcompt );

            ++icarowcompt;
        }

        /*
         * Optimize by moving "if ( instant )" outside the loop, because it is
         * rarely used. May also be able to avoid "if ( power )".
         *
         * Or not: excellent branch predictors these days.
         *
         * Will be nice to test these optimizations.
         */
        chanBoundary = ichan + *ichannelcount;
        for ( ; ichan < chanBoundary; ++ichan )
        {

	  caTable_.row( *iextca, dRow );

            if ( ichan->Xpower_ > 0.0 )
            {
                vTable_.lookup( *icolumn, vRow, C1, C2 );
                //~ *istate = *istate * C1 + C2;
                //~ *istate = ( C1 + ( 2 - C2 ) * *istate ) / C2;
                if ( ichan->instant_ & INSTANT_X )
                    *istate = C1 / C2;
                else
                {
                    double temp = 1.0 + dt / 2.0 * C2;
                    *istate = ( *istate * ( 2.0 - temp ) + dt * C1 ) / temp;
                }

                ++icolumn, ++istate;
            }

            if ( ichan->Ypower_ > 0.0 )
            {
                vTable_.lookup( *icolumn, vRow, C1, C2 );
                //~ *istate = *istate * C1 + C2;
                //~ *istate = ( C1 + ( 2 - C2 ) * *istate ) / C2;
                if ( ichan->instant_ & INSTANT_Y )
                    *istate = C1 / C2;
                else
                {
                    double temp = 1.0 + dt / 2.0 * C2;
                    *istate = ( *istate * ( 2.0 - temp ) + dt * C1 ) / temp;

}
                ++icolumn, ++istate;
            }

            if ( ichan->Zpower_ > 0.0 )
            {
                LookupRow* caRow = *icarow;

                if ( caRow )
                {
                    caTable_.lookup( *icolumn, *caRow, C1, C2 );

                }
                 else if (*iextca >0)

		   {
		     caTable_.lookup( *icolumn, dRow, C1, C2 );
		   }
		else
                {
		  vTable_.lookup( *icolumn, vRow, C1, C2 );

                }

                //~ *istate = *istate * C1 + C2;
                //~ *istate = ( C1 + ( 2 - C2 ) * *istate ) / C2;
                if ( ichan->instant_ & INSTANT_Z )
                    *istate = C1 / C2;
                else
                {
                    double temp = 1.0 + dt / 2.0 * C2;
                    *istate = ( *istate * ( 2.0 - temp ) + dt * C1 ) / temp;
                }

                ++icolumn, ++istate, ++icarow;

            }
	    ++iextca;
        }

        ++ichannelcount, ++icacount;
    }
#endif
}

/**
 * SynChans are currently not under solver's control
 */
void HSolveActive::advanceSynChans( ProcPtr info )
{
    return;
}

void HSolveActive::sendSpikes( ProcPtr info )
{
    vector< SpikeGenStruct >::iterator ispike;
    for ( ispike = spikegen_.begin(); ispike != spikegen_.end(); ++ispike )
        ispike->send( info );
}

/**
 * This function dispatches state values via any source messages on biophysical
 * objects which have been taken over.
 *
 */
void HSolveActive::sendValues( ProcPtr info )
{
    vector< unsigned int >::iterator i;

    for ( i = outVm_.begin(); i != outVm_.end(); ++i )
        Compartment::VmOut()->send(
            //~ ZombieCompartment::VmOut()->send(
            compartmentId_[ *i ].eref(),
            V_[ *i ]
        );


    for ( i = outIk_.begin(); i != outIk_.end(); ++i ){

        unsigned int comptIndex = chan2compt_[ *i ];

        assert( comptIndex < V_.size() );

        ChanBase::IkOut()->send(channelId_[*i].eref(),
				(current_[ *i ].Ek - V_[ comptIndex ]) * current_[ *i ].Gk);

    }

    for ( i = outCa_.begin(); i != outCa_.end(); ++i )
        //~ CaConc::concOut()->send(
        CaConcBase::concOut()->send(
            caConcId_[ *i ].eref(),
            ca_[ *i ]
        );
}

void HSolveActive::preProcess(){
	int num_compts = V_.size();
	int num_channels = channel_.size();
	int num_cmprsd_gates = state_.size();
	int num_Ca_pools = ca_.size();

	// Channel variables
	h_chan_Gbar = new double[num_channels];
	h_chan_instant = new int[num_channels];
	h_chan_modulation= new double[num_channels];
	h_chan_to_comp = new int[num_channels];

	// Gathering data for each channel
	for(unsigned int i=0;i<channel_.size();i++){
		h_chan_Gbar[i] = channel_[i].Gbar_;

		h_chan_instant[i] = channel_[i].instant_;
		h_chan_modulation[i] = channel_[i].modulation_;

		// Channel to Compartment Info
		h_chan_to_comp[i] = chan2compt_[i];
	}

	// Constructing ca and channel row ptrs with nCompt as rows.
	int ca_rowPtr[V_.size()+1];
	int chan_rowPtr[V_.size()+1];
	int sum1 = 0, sum2 = 0;
	for(unsigned int i=0;i<=V_.size();i++){
		ca_rowPtr[i] = sum1;
		chan_rowPtr[i] = sum2;

		if(i < V_.size()){
			// Last one should be just set.
			sum1 += caCount_[i];
			sum2 += channelCount_[i];
		}
	}

	// Optimized version
	h_state_powers = new double[num_cmprsd_gates];
	h_state2chanId = new int[num_cmprsd_gates];
	h_state2column = new int[num_cmprsd_gates];
	int* h_state_rowPtr = new int[num_channels+1]();

	// Gathering gate information and separating gates (with power < 0) , (with vm dependent) , (with ca dependent)
	int cmprsd_gate_index = 0; // If the logic is true cmprsd_gate_index value at the end of for loop = # of gates with powers > 0
	double h_gate_powers[3]; // Reusable array for holding powers of a channel
	for(unsigned int i=0;i<V_.size();i++){
		for(int j=chan_rowPtr[i]; j<chan_rowPtr[i+1]; j++){

			// Setting powers
			h_gate_powers[0] = channel_[j].Xpower_;
			h_gate_powers[1] = channel_[j].Ypower_;
			h_gate_powers[2] = channel_[j].Zpower_;

			for(int k=0;k<3;k++){
				if(h_gate_powers[k] > 0){

					// Collecting power of valid gate
					switch(k){
						case 0:
							h_state_powers[cmprsd_gate_index] = channel_[j].Xpower_;
							break;
						case 1:
							h_state_powers[cmprsd_gate_index] = channel_[j].Ypower_;
							break;
						case 2:
							h_state_powers[cmprsd_gate_index] = channel_[j].Zpower_;
							break;
					}

					// Collecting channel and column of valid gate
					h_state2chanId[cmprsd_gate_index] = j;
					h_state2column[cmprsd_gate_index] = column_[cmprsd_gate_index].column;

					// Partitioning of vm and ca dependent gates.
					if(k == 2 && caDependIndex_[j] != -1){
						h_cagate_indices.push_back(cmprsd_gate_index);
						h_cagate_capoolIds.push_back(ca_rowPtr[i] + caDependIndex_[j]);
					}else{
						h_vgate_compIds.push_back((int)chan2compt_[j]);
						h_vgate_indices.push_back(cmprsd_gate_index);

						// k=2 gate might depend on externalCalcium
						if(k == 2)
							h_exCalgate_indices.push_back(cmprsd_gate_index);
					}
					h_state_rowPtr[j] += 1;
					cmprsd_gate_index++; // cmprsd_gate_index is incremented only if power > 0 is found.
				}
			}
		}
	}
	assert(cmprsd_gate_index == num_cmprsd_gates);

	// Converting rowCounts to rowptr
	int csum = 0, ctemp;
	int zero_count = 0;
	for (int i = 0; i < num_channels+1; ++i) {
		ctemp = h_state_rowPtr[i];
		if(i < num_channels && h_state_rowPtr[i] == 0) zero_count++;
		h_state_rowPtr[i] = csum;
		csum += ctemp;
	}

    // Advance Channels
    h_V_rows = new int[num_compts];
    h_V_fracs = new double[num_compts];
}

