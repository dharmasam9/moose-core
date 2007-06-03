/**********************************************************************
** This program is part of 'MOOSE', the
** Messaging Object Oriented Simulation Environment,
** also known as GENESIS 3 base code.
**           copyright (C) 2003-2005 Upinder S. Bhalla. and NCBS
** It is made available under the terms of the
** GNU General Public License version 2
** See the file COPYING.LIB for the full notice.
**********************************************************************/

#include "moose.h"
#include "MsgDest.h"
#include "ParFinfo.h"
#include <mpi.h>
#include "PostMaster.h"
#include <sstream>


ParFinfo::ParFinfo( const string& name ) 
	: Finfo( name, Ftype0::global() )
{
	;
}


/**
 * This adds a message from a postmaster to a destination object.
 * It is only called when a message is sent between nodes.
 * It fills up the MsgSrc slot and incomingFunc_ array starting
 * at msgIndex.
 */
bool ParFinfo::add(
	Element* e, Element* destElm, const Finfo* destFinfo,
	unsigned int msgIndex
	// This extra argument is the slot in which the msg goes.
) const
{
	FuncList srcFl;
	FuncList destFl;
	unsigned int destIndex;
	unsigned int numDest;

	// destFinfo->ftype()->appendIncomingFunc( returnFl );
	
	const Ftype* dftype = destFinfo->ftype()->makeMatchingType();
	// Here we put the srcFl entries into the srcFl list
	// in case it is a SharedFinfo.
	// These are recvFuncs generated from the destFinfo Ftype
	// to handle return function calls.

	if ( destFinfo->respondToAdd( destElm, e, dftype,
							srcFl, destFl,
							destIndex, numDest ) )
	{
		unsigned int originatingConn;
		unsigned int targetConn;

		// First we decide where to put the originating Conn.
		if ( destFl.size() == 0 ) {  // The msg'src' has only dests
			originatingConn = e->insertConnOnDest( msgIndex_, 1);
		} else { // The usual case: put it on MsgSrc.
			originatingConn = 
					e->insertConnOnSrc( msgIndex, destFl, 0, 0 );
		}

		// Now the target Conn
		if ( srcFl.size() == 0 ) { // Target has only dests.
			targetConn = destElm->insertConnOnDest( destIndex, 1 );
		} else { // Here we need to put it on target MsgSrc.
			targetConn = 
					destElm->insertConnOnSrc( destIndex, srcFl, 0, 0 );
		}

		// Finally, we're ready to do the connection.
		e->connect( originatingConn, destElm, targetConn );

		// Now we put in the incomingFuncs to
		// invoke the recvFuncs of the destination elements.

		assert( e->className() == "PostMaster" );
		PostMaster* post = static_cast< PostMaster* >( e->data() );
		// vector< IncomingFunc > inFl;
		// destFinfo->ftype()->appendIncomingFunc( inFl );
		// post->placeIncomingFuncs( inFl, msgIndex );

		return 1;
	}
	return 0;
}


/**
 * This responds to a message request on a postmaster directed to
 * a remote node. It sends off a request for the remote node to
 * complete the message, creates the local node message and returns 
 * immediately. Later the remote node will complete the request and
 * if all is well, nothing happens. If the remote node fails to make
 * the message it will trigger a cleanup function which deletes this
 * message and issues an error call.
 *
 * The Element* e is a dummy Element with extra info
 * to indicate the id of the true target, and provision to store the
 * target field name. It is made on the fly by
 * the Element::element() function when the lookup returns a
 * node index rather than the element.
 * Here we basically package the entire argument list into a string
 * for sending abroad.
 * - It is a comma separated list.
 * - The target id is obtained from the dummy Element e.
 * - The originating id is obtained from the element src.
 * - The kind of message ( scheduled vs sporadic ) is obtained by
 *   traversing up through the message
 * - The type info is converted to strings using RTTI features and
 *   appended to the msg info.
 * - The srcFL and returnFl are immediately dealt with and are not
 *   transmitted.
 * - destIndex and numDest are immediately filled in for the local
 *   postmaster messaging info, and are not transmitted.
 */
bool ParFinfo::respondToAdd(
					Element* e, Element* src, const Ftype *srcType,
					FuncList& srcFl, FuncList& returnFl,
					unsigned int& destIndex, unsigned int& numDest
) const
{
	assert( e != 0 );
	assert( e->className() == "PostMaster" );
	// PostMaster* post = static_cast< PostMaster* >( e->data());
	
	// const Ftype* srcType = srcFinfo->ftype();

	unsigned int targetId;
	get< unsigned int >( e, "targetId", targetId );
	string targetField;
	get< string >( e, "targetField", targetField );
//	string respondString;
	ostringstream oss;
// 	oss.str( respondString );
	oss << src->id() << " " << targetId << " " << 
			targetField << " " << ftype2str( srcType );

	// Post irecv for return status of this message request.
	// Send out this message request.
	// Create local node message.
	
	// First, fill up the returnFl with functions provided by the
	// srcType.
	returnFl.resize( 0 );
	/// \todo: work out how to decide if a message should be sync/async.
	returnFl.push_back( srcType->asyncFunc() );
	numDest = returnFl.size();

	// Second, refer this to the postmaster which does the following:
	// - post irecv for return status of this message
	// - send out message request
	// - create local message
	destIndex = PostMaster::respondToAdd( e, oss.str(), numDest );
	
	// destIndex = post->outgoingSlotNum_;
	// post->outgoingSlotNum_ += numDest;

	// Still need to figure out how to deal with stuff in the srcFl.
	return 1;
}



////////////////////////////////////////////////////////////////////

void ParFinfo::dropAll( Element* e ) const
{
	;
}

bool ParFinfo::drop( Element* e, unsigned int i ) const
{
		return 0;
}

unsigned int ParFinfo::numIncoming( const Element* e ) const
{
		return 0;
}

unsigned int ParFinfo::numOutgoing( const Element* e ) const
{
		return 0;
}

unsigned int ParFinfo::incomingConns(
					const Element* e, vector< Conn >& list ) const
{
		return 0;
}

unsigned int ParFinfo::outgoingConns(
					const Element* e, vector< Conn >& list ) const
{
		return 0;
}

/// This Finfo does not support recvFuncs.
RecvFunc ParFinfo::recvFunc() const
{
		return &dummyFunc;
}

void ParFinfo::countMessages( 
					unsigned int& srcNum, unsigned int& destNum )
{
		return;
}
