#ifndef _INT_FIRE_H
#define _INT_FIRE_H


class IntFire: public Data
{
	friend void testStandaloneIntFire();
	friend void testSynapse();
	public: 
		IntFire();
		IntFire( double thresh, double tau );
		void process( const ProcInfo* p, const Eref& e );
		void reinit( Eref& e );

		/**
 		 * Inserts an event into the pendingEvents queue for spikes.
 		 */
		void addSpike( DataId synIndex, const double& time );
		
		////////////////////////////////////////////////////////////////
		// Field assignment stuff.
		////////////////////////////////////////////////////////////////
		
		void setVm( double v );
		double getVm() const;
		void setTau( double v );
		double getTau() const;
		void setThresh( double v );
		double getThresh() const;
		void setNumSynapses( unsigned int v );
		unsigned int getNumSynapses() const;

		static const Cinfo* initCinfo();

		Synapse* synapse( unsigned int i );
	private:
		double Vm_; // State variable: Membrane potential. Resting pot is 0.
		double thresh_;	// Firing threshold
		double tau_; // Time course of membrane settling.
		vector< Synapse > synapses_;
		priority_queue< Synapse > pendingEvents_;
};

#endif // _INT_FIRE_H
