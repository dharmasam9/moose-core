/**********************************************************************
** This program is part of 'MOOSE', the
** Messaging Object Oriented Simulation Environment.
**           Copyright (C) 2003-2015 Upinder S. Bhalla. and NCBS
** It is made available under the terms of the
** GNU Lesser General Public License version 2.1
** See the file COPYING.LIB for the full notice.
**********************************************************************/

#include "header.h"
#include "../shell/Shell.h"
#include "../utility/Vec.h"
#include "ReadSwc.h"
#include "CompartmentBase.h"
#include "Compartment.h"
#include "SymCompartment.h"
#include <fstream>


/*
#include <cstdlib>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <cassert>
*/

unsigned short SwcSegment::BadSegment = 8;

// Minimum allowed radius of segment, in microns
static const double MinRadius = 0.05; 

SwcSegment::SwcSegment( const string& line )
{
	vector< string > args;
	stringstream ss( line );
	string temp;
	while (ss >> temp ) {
		args.push_back( temp );
	}
	if ( args.size() == 7 ) {
		myIndex_ = atoi( args[0].c_str() );
		type_ = atoi( args[1].c_str() );
		double x = atof( args[2].c_str() );
		double y = atof( args[3].c_str() );
		double z = atof( args[4].c_str() );
		v_ = Vec( x, y, z );
		radius_ = atof( args[5].c_str() );
		int pa = atoi( args[6].c_str() );
		if ( pa > 0 )
			parent_ = pa;
		else
			parent_ = ~0U;
	} else {
		type_ = BadSegment;
	}
}

double SwcSegment::L( const SwcSegment& other ) const 
{
	double r = ( radius_ < other.radius_ ) ? radius_:other.radius_;
	if ( r > 0.0 )
		return v_.distance( other.v_ ) / sqrt( r );
	return 1.0;
}

void SwcSegment::figureOutType()
{
	if ( type_ == 1 ) // already defined as soma
		return;
	if ( kids_.size() > 1 )
		type_ = 5; // fork point
	else if ( kids_.size() == 0 )
		type_ = 6; // end point 
	else if ( kids_.size() == 1 && ( type_ == 5 || type_ == 6 ) )
		type_ = 3; // Fix it with the most generic option.
}

//////////////////////////////////////////////////////////////////////

SwcBranch::SwcBranch( int i,  const SwcSegment& start, double len, double L,
				const vector< int >& cable )
				: SwcSegment( start ),
				r0( start.radius() ),
				r1( start.radius() ),
				geomLength( len ),
				electroLength( L )
{
	myIndex_ = i;
	parent_ = 0;
	kids_.resize( 0 );
	segs_.resize( cable.size() );
	// Put the contents of cable into segs, in reverse order.
	vector< int >::const_reverse_iterator j = cable.rbegin();
	vector< int >::iterator k = segs_.begin();
	for ( k = segs_.begin(); k != segs_.end(); ++k )
		*k = *j++;
}

void SwcBranch::printDiagnostics() const
{
	cout << myIndex() << ":  " << segs_[0] << " -> " << segs_.back() <<
			" = " << segs_.size() << 
			" :	pa = " << parent() << " ,	length=( " << 
			geomLength << ", " << electroLength << " )\n";
}

//////////////////////////////////////////////////////////////////////

ReadSwc::ReadSwc( const string& fname )
{
	ifstream fin( fname.c_str() );
	if ( !fin ) {
		cerr << "ReadSwc:: could not open file " << fname << endl;
		return;
	}

	string temp;
	int badSegs = 0;
	while( getline( fin, temp ) ) {
		if ( temp.length() == 0 ) 
			continue;
		string::size_type pos = temp.find_first_not_of( "\t " );
		if ( pos == string::npos )
			continue;
		if ( temp[pos] == '#' )
			continue;

		SwcSegment t( temp );
		if ( t.OK() )
			segs_.push_back( SwcSegment( temp ) );
		else
			badSegs++;
	}
	bool valid = validate();
	if ( valid ) {
		assignKids();
		cleanZeroLength();
		parseBranches();
	}
	cout << fname << "	: NumSegs = " << segs_.size() << 
			", bad = " << badSegs <<
			", Validated = " << valid << 
			", numBranches = " << branches_.size() << 
			endl;
	diagnostics();
}

bool ReadSwc::validate() const
{
	int numStart = 0;
	int numOrphans = 0;
	int badIndex = 0;
	int badRadius = 0;
	for ( unsigned int i = 0; i < segs_.size(); ++i ) {
		const SwcSegment& s = segs_[i];
		if ( s.myIndex() != i + 1 )
			badIndex++;
		if ( s.parent() == ~0U ) {
			numStart++;
		} else {
			if ( s.parent() > i ) {
				numOrphans++;
			}
		}
		if ( s.radius() < MinRadius ) {
			badRadius++;
		}
	}
	bool valid = ( numStart == 1 && numOrphans == 0 && badRadius == 0 );
	if ( !valid ) {
		cout << "NumSegs = " << segs_.size() << 
				", numStart = " << numStart <<
				", orphans = " << numOrphans << 
				", badIndex = " << badIndex <<
				", badRadius = " << badRadius <<
				", numBranches = " << branches_.size() <<
				endl;
	}
	return valid;
}

void ReadSwc::assignKids()
{
	for ( unsigned int i = 0; i < segs_.size(); ++i ) {
		const SwcSegment& s = segs_[i];
		assert ( s.parent() != s.myIndex() );
		if ( s.parent() != ~0U ) {
			segs_[s.parent() - 1].addChild( i + 1 );
		}
	}
	for ( unsigned int i = 0; i < segs_.size(); ++i ) {
		segs_[i].figureOutType();
	}
}

void ReadSwc::cleanZeroLength()
{
	static double EPSILON = 1e-2; // Assume units in microns.
	for ( unsigned int i = 1; i < segs_.size(); ++i ) {
		SwcSegment& s = segs_[i];
		SwcSegment& pa = segs_[ s.parent() - 1 ];
		if ( s.length( pa ) < EPSILON ) {
			// Remove the zero length child from pa.kids_
			vector< int > temp;
			for ( unsigned int j = 0; j < pa.kids().size(); ++j ) {
				if ( static_cast< unsigned int >( pa.kids()[j] ) != s.myIndex() )
					temp.push_back( pa.kids()[j] );
			}
			// Go through all kids of s and reparent them.
			for ( unsigned int j = 0; j < s.kids().size(); ++j ) {
				SwcSegment& kid = segs_[ s.kids()[j] - 1 ];
				kid.setParent( pa.myIndex() );
				temp.push_back( kid.myIndex() );
			}
			pa.replaceKids( temp );
			s.setBad();
			cout << "Cleaned zero length " << s.myIndex() << endl;
		}
	}
}

void ReadSwc::traverseBranch( const SwcSegment& s, 
		double& len, double& L, vector< int >& cable ) const
{
	const SwcSegment* prev = &s;
	cable.resize( 1, s.myIndex() ); // Always include the starting seg.
	// Note that the cable is filled up with entries in reverse order.

	if ( s.parent() == ~0U ) {
		len = s.radius();
		L = sqrt( len );
		return ;
	}

	do {
		// Beware the indexing!
		const SwcSegment& pa = segs_[prev->parent() - 1];
		len += pa.length( *prev );
		L += pa.L( *prev );
		cable.push_back( pa.myIndex() );
		prev = &pa;
	} while ( (prev->parent() != ~0U) && (prev->kids().size() == 1) );
	cable.pop_back(); // Get rid of the last entry, it is on the parent.
}

void ReadSwc::parseBranches()
{
	// Fill vector of all branches.
	for ( unsigned int i = 0; i < segs_.size(); ++i ) {
		const SwcSegment& s = segs_[i];
		if ( s.OK() && s.kids().size() != 1 ) { // Either use a fork or an end.
			vector< int > cable;
			// int branchIndex = branches_.
			// branches_.push_back( i + 1 );
			double len = 0;
			double L = 0;
			traverseBranch( s, len, L, cable );
			// branchGeomLength_.push_back( len );
			// branchElectroLength_.push_back( L );
			SwcBranch br( branches_.size(), s, len, L, cable );
			branches_.push_back( br );
		}
	}
	// Assign the parent of each branch. This is known because the
	// parent of the first segment in the branch is the last segment
	// in the parent branch. I construct a reverse lookup table to find
	// the branch # from its last segment number.
	vector< int > reverseSeg ( segs_.size() + 1, 0 );
	for ( unsigned int i = 0; i < branches_.size(); ++i )
		reverseSeg[ branches_[i].segs_.back() ] = i;
	for ( unsigned int i = 0; i < branches_.size(); ++i ) {
		int parentSeg = segs_[ branches_[i].segs_[0] - 1 ].parent();
		assert( parentSeg != 0 ); // Note that segment indices start from 1
		branches_[i].setParent( reverseSeg[ parentSeg ] );
	}
}

void ReadSwc::diagnostics() const
{
	vector< int > diag( 8 );
	static string diagName[] = {"undef", "soma", "axon", "dend", 
			"apical dend", "fork", "end", "custom" };
	for ( unsigned int i = 0; i < segs_.size(); ++i ) {
		const SwcSegment& s = segs_[i];
		if ( s.type() < 8 )
			diag[s.type()]++;
	}
	for ( int i = 0; i < 8; ++i )
		cout << diagName[i] << " :	" << diag[i] << endl;

	for ( unsigned int i = 0; i < branches_.size(); ++i )
		branches_[i].printDiagnostics();
}

static Id makeCompt( Id parent, 
		const SwcSegment& seg, const SwcSegment& pa,
		double RM, double RA, double CM,
		unsigned int i, unsigned int j	)
{
	Shell* shell = reinterpret_cast< Shell* >( Id().eref().data() );
	double len = seg.radius() * 2.0;
	string name = "soma";
	Id compt;
	double x0, y0, z0;
	if ( seg.parent() != ~0U ) {
		len = seg.length( pa );
		stringstream ss;
		ss << "br" << i << "_" << j;
		name = ss.str();
		x0 = pa.vec().a0();
		y0 = pa.vec().a1();
		z0 = pa.vec().a2();
	} else {
		x0 = seg.vec().a0() - len;
		y0 = seg.vec().a1();
		z0 = seg.vec().a2();
	}
	assert( len > 0.0 );
	compt = shell->doCreate( "Compartment", parent, name, 1 );
	Eref er = compt.eref();
	moose::CompartmentBase *cptr = reinterpret_cast< moose::CompartmentBase* >(
					compt.eref().data() );
	double xa = seg.radius() * seg.radius() * PI * 1e-12;
	len *= 1e-6;
	cptr->setRm( er, RM / ( len * seg.radius() * PI ) );
	cptr->setRa( er, RA * len / xa );
	cptr->setCm( er, CM * ( len * seg.radius() * PI ) );
	cptr->setDiameter( seg.radius() * 2.0e-6 );
	cptr->setLength( len * 1e-6 );
	cptr->setX0( x0 * 1e-6 );
	cptr->setY0( y0 * 1e-6 );
	cptr->setZ0( z0 * 1e-6 );
	cptr->setX( seg.vec().a0() * 1e-6 );
	cptr->setY( seg.vec().a1() * 1e-6 );
	cptr->setZ( seg.vec().a2() * 1e-6 );
	return compt;
}

bool ReadSwc::build( Id parent, 
				double lambda, double RM, double RA, double CM )
{
	Shell* shell = reinterpret_cast< Shell* >( Id().eref().data() );
	vector< Id > compts( segs_.size() );
	for ( unsigned int i = 0; i < branches_.size(); ++i ) {
		SwcBranch& br = branches_[i];
		for ( unsigned int j = 0; j < br.segs_.size(); ++j ) {
			Id compt;
			SwcSegment& seg = segs_[ br.segs_[j] -1 ];
			unsigned int paIndex = seg.parent();
			if ( paIndex == ~0U ) { // soma
				compt = makeCompt( parent, seg, seg, RM, RA, CM, i, j );
			} else {
				SwcSegment& pa = segs_[ paIndex - 1 ];
				compt = makeCompt( parent, seg, pa, RM, RA, CM, i, j );
				assert( compt != Id() );
				assert( compts[ paIndex -1 ] != Id() );
				shell->doAddMsg( "Single", 
					compts[paIndex-1], "axial", compt, "raxial" );
			}
			assert( compt != Id() );
			compts[ seg.myIndex() -1 ] = compt;
		}
	}
	return true;
}