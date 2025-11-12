/***************************************************************************
 *   This file is part of the MagiC++ library.                             *
 *                                                                         *
 *   Copyright (C) 1998-2005 Marko Grönroos <magi@iki.fi>                  *
 *                                                                         *
 ***************************************************************************
 *                                                                         *
 *  This library is free software; you can redistribute it and/or          *
 *  modify it under the terms of the GNU Library General Public            *
 *  License as published by the Free Software Foundation; either           *
 *  version 2 of the License, or (at your option) any later version.       *
 *                                                                         *
 *  This library is distributed in the hope that it will be useful,        *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU      *
 *  Library General Public License for more details.                       *
 *                                                                         *
 *  You should have received a copy of the GNU Library General Public      *
 *  License along with this library; see the file COPYING.LIB.  If         *
 *  not, write to the Free Software Foundation, Inc., 59 Temple Place      *
 *  - Suite 330, Boston, MA 02111-1307, USA.                               *
 *                                                                         *
 ***************************************************************************/

#include <magic/mclass.h>
#include <nhp/individual.h>
#include <inanna/annetwork.h>
#include <inanna/initializer.h>
#include "annalee/anngenes.h"
#include "annalee/miller.h"

impl_dynamic (MillerEncoding, {ANNEncoding});


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//    ----                  ___   o                     -----                //
//    |   )            ___  |  \         ___   ___   |  |       _    ___     //
//    |---  |   | |/\ /   ) |   | | |/\ /   ) |   \ -+- |---  |/ \  |   \    //
//    |     |   | |   |---  |   | | |   |---  |      |  |     |   | |        //
//    |      \__! |    \__  |__/  | |    \__   \__/   \ |____ |   |  \__/    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/*******************************************************************************
 * Standard constructor.
 *
 * @param params["pruneInputs"] Should the input neurons be encoded in the
 * genome also. [Default=1]
 *
 * @param params["pcVariance"] Variance distribution of the
 * connection probability of the initial neural networks. [Default=0]
 *
 * @param params["pcAverage"] Average connection probability for
 * initial neural networks. [Default=0.5]
 *
 * @param name Gene name, usually "brainplan".
 ******************************************************************************/
MillerEncoding::MillerEncoding (
	const GeneticID& name,   //< Name of the gene
	const StringMap& params) //< Parameters
		: ANNEncoding (name, params)
{
	mPruneInputs	= params["MillerEncoding.pruneInputs"].toInt ();
	mPcVariance		= params["MillerEncoding.pcVariance"].toDouble ();
	mPcAverage		= params["MillerEncoding.pcAverage"].toDouble ();
}

MillerEncoding::MillerEncoding (const MillerEncoding& other) : ANNEncoding (other) {
	mPruneInputs = other.mPruneInputs;
	mPcVariance  = other.mPcVariance;
	mPcAverage   = other.mPcAverage;
}

/*******************************************************************************
 * Implementation for @ref Genstruct.
 ******************************************************************************/
void MillerEncoding::copy (const Genstruct& o) {
	ANNEncoding::copy (o);
	const MillerEncoding& other = static_cast<const MillerEncoding&>(o);
	mPruneInputs = other.mPruneInputs;
	mPcVariance = other.mPcVariance;
	mPcAverage = other.mPcAverage;
}

/*******************************************************************************
 * Implementation for @ref Genstruct.
 ******************************************************************************/
void MillerEncoding::addPrivateGenes (Gentainer& g, const StringMap& params) {
	Gentainer::addPrivateGenes (g, params);

	// Create genes for the units
	int totalUnits = mInputs+mMaxHidden+mOutputs;
	for (int i=0; i<totalUnits-mOutputs; i++) {

		// Existance of a neuron. Not encoded for input units if input
		// pruning is not enabled, nor output units which always exist
		if (i>=mInputs || mPruneInputs)
			add (new BinaryGene (format ("E%d", i)));
		
		// Connect input units and hidden units to all successive neurons
		for (int j=i+1; j<totalUnits; j++)
				add (&(new BinaryGene (format ("W%d-%d", i, j)))->hide());
	}
}

/*******************************************************************************
 * Implementation for @ref Genstruct.
 ******************************************************************************/
void MillerEncoding::init () {
	// Find random p_c (connection probability) for the genome
	double pc=0.5;
	do {
		pc = mPcAverage+gaussrnd(mPcVariance);
	} while (pc<0.0 || pc>1.0); // Not before it's in a range of a probability value

	// Change the pc of all binary genes we own
	for (int i=0; i<size(); i++)
		if (BinaryGene* gene = dynamic_cast<BinaryGene*> (&(*this)[i]))
			gene->setInitP(pc);

	// Let the superclass implement the initialization
	Gentainer::init ();
}

/*******************************************************************************
 * Implementation for @ref Genstruct.
 ******************************************************************************/
bool MillerEncoding::execute (const GeneticMsg& msg) const {
	// Take pictures only if this is a picture-taking recreation
	bool takePics = dynamic_cast<const TakeBrainPicsMsg*>(&msg) != NULL;

	ANNetwork* net = new ANNetwork (format ("%d-%d-%d", mInputs, mMaxHidden, mOutputs));
	int totalUnits = mInputs+mMaxHidden+mOutputs;
	
	// Create and zero the connection matrix. Using this matrix is
	// useful mostly just for informative purposes (in logging)
	PackTable<int> cmatrix (totalUnits,totalUnits);
	for (int i=0; i<cmatrix.rows; i++)
		for (int j=0; j<cmatrix.rows; j++)
			cmatrix.get(i,j) = 0;

	// Go trough each input and hidden unit and check if it exists
	for (int i=0; i<totalUnits-mOutputs; i++) {
		if (i>=mInputs || mPruneInputs)
			// Enable or disable the unit from the network
			cmatrix.get(i,i) =
				static_cast<const BinaryGene&> ((*this)[(CONSTR)format ("E%d", i)]).getvalue();
	}


	// Fill the connection matrix
	for (int i=0; i<totalUnits-mOutputs; i++)
		for (int j=i+1; j<totalUnits; j++)
			if (j>=mInputs)
				if (static_cast<const BinaryGene&> (
					(*this)[(CONSTR)format ("W%d-%d", i, j)]).getvalue())
					cmatrix.get(i,j) = 1;

	// Calculate some statistics (probability of connection)
	int conns=0, totconns=0;
	for (int i=0; i<totalUnits-mOutputs; i++)
		for (int j=(i>=mInputs)?i+1:mInputs; j<totalUnits; j++) {
			totconns++;
			conns += cmatrix.get(i,j);
		}
	double pConn = double(conns) / double(totconns);
	msg.mrHost.set ("pConn", new String (format ("%f", pConn)));

	// Enable outputs
	for (int i=totalUnits-mOutputs; i<totalUnits; i++)
		cmatrix.get(i,i) = 1;
	
	// Enable and connect according to the connection matrix (if both
	// source and target units exist and also the connection)
	for (int i=0; i<cmatrix.rows; i++) {
		(*net)[i].enable(cmatrix.get(i,i));
		if (cmatrix.get(i,i))
			for (int j=i+1; j<cmatrix.rows; j++)
				if (cmatrix.get(j,j) && cmatrix.get(i,j))
					net->connect (i,j);
	}

	//for (int i=totalUnits-mOutputs; i<totalUnits; i++)
	//	(*net)[i].setTFunc (FreeNeuron::LINEAR_TF);

	// "Print" network connectivity matrix
	String desc; // Printed connection matrix
	if (takePics) {
		desc.reserve(totalUnits*(totalUnits+1)+10);
		for (int i=0; i<cmatrix.rows; i++) {
			for (int j=0; j<cmatrix.cols; j++)
				desc += cmatrix.get(i,j)? '1':'0';
			desc += '\n';
		}
	}
	
	if (net) {
		net->setInitializer (new GaussianInitializer (0.5));

		// Take some nice photos
		net->cleanup ();
		if (takePics)
			msg.mrHost.set ("brainpic1", new String (net->drawEPS()));
		net->cleanup (true, mPrunePassthroughs);
		if (takePics) {
			msg.mrHost.set ("brainpic2", new String (net->drawEPS()));
			net->drawFeedForward();
			msg.mrHost.set ("brainpic3", new String (net->drawEPS()));
			msg.mrHost.set ("braindesc1", new String (desc));
			delete net; // The net was created only for taking babypics
		} else {
			// Place the brain description into host
			msg.mrHost.set ("brainplan", net);
		}
	}

	return true;
}

