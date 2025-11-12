/***************************************************************************
 *   This file is part of the Annalee library.                             *
 *                                                                         *
 *   Copyright (C) 1998-2002 Marko Grönroos <magi@iki.fi>                  *
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

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  ___   o                     -----                      | o              //
//  |  \         ___   ___   |  |       _    ___           |     _          //
//  |   | | |/\ /   ) |   \ -+- |---  |/ \  |   \  __   ---| | |/ \   ___   //
//  |   | | |   |---  |      |  |     |   | |     /  \ (   | | |   | (   \  //
//  |__/  | |    \__   \__/   \ |____ |   |  \__/ \__/  ---| | |   |  ---/  //
//                                                                    __/   //
//////////////////////////////////////////////////////////////////////////////

DirectEncoding::DirectEncoding (const GeneticID& name, const StringMap& params) : ANNEncoding (name, params) {
	mPruneInputs   = isnull(params["prune_inputs"])? true  : params["prune_inputs"].toInt ();
	mPruneWeights  = isnull(params["prune_weights"])? false : params["prune_weights"].toInt ();
	mEncodeWeights = isnull(params["encode_weights"])? false : params["encode_weights"].toInt ();
}

DirectEncoding::DirectEncoding (const DirectEncoding& other) : ANNEncoding (other) {
	mPruneInputs = other.mPruneInputs;
	mPruneWeights = other.mPruneWeights;
	mEncodeWeights = other.mEncodeWeights;
}

void DirectEncoding::copy (const Genstruct& o) {
	ANNEncoding::copy (o);
	const DirectEncoding& other = static_cast<const DirectEncoding&>(o);
	mPruneInputs = other.mPruneInputs;
	mPruneWeights = other.mPruneWeights;
	mEncodeWeights = other.mEncodeWeights;
}

void DirectEncoding::addPrivateGenes (Gentainer& g, const StringMap& params) {
	Gentainer::addPrivateGenes (g, params);

	if (mPruneInputs)
		for (int i=0; i<mInputs; i++)
			add (new BinaryGene (format ("R%d", i), 1.0));

	// Inputs can be ...what?
	// TODO
	for (int i=0; i < mInputs + mMaxHidden; i++)
		for (int j=0; j < mMaxHidden + mOutputs; j++) {
			// Encode weight existence.
			if (mPruneWeights)
				add (new BinaryGene (format ("WX%d-%d", i, j), 1.0));

			// Encode weight value.
			if (mEncodeWeights)
				add (new BinaryGene (format ("WV%d-%d", i, j), 1.0));
		}
		
	// Prune hidden
	for (int i=0; i<mMaxHidden; i++)
		add (new BinaryGene (format ("H%d", i), 1.0));
}

bool DirectEncoding::execute (const GeneticMsg& msg) const {
	FreeNetwork* net = new FreeNetwork (format ("%d-%d-%d", mInputs, mMaxHidden, mOutputs));
	
	// Go trough each hidden unit and check if it exists
	bool hidexists [mMaxHidden];
	for (int h=0; h<mMaxHidden; h++) {
		hidexists[h] = static_cast<const BinaryGene&> (
			(*this)[(CONSTR)format ("H%d", h)]).getvalue();
		// TRACE2 ("%d=%d", h, int(hidexists[h]));
		
		// Enable or disable it from the network
		(*net)[h+mInputs].enable(hidexists[h]);
	}


	bool input_exists;
	bool w_exists;
	
	// Connect each input to every hidden neuron
	for (int i=0; i < mInputs; i++) {
		
		// See if this input unit "exists" (if input pruning is enabled)
		input_exists = true;
		if (mPruneInputs)
			input_exists = static_cast<const BinaryGene&> (
				(*this)[(CONSTR)format ("R%d", i)]).getvalue();

		(*net)[i].enable (input_exists);
		
		// If it exists...
		if (input_exists) {
			
			// Go trough each (existing) hidden unit
			for (int h=0; h<mMaxHidden; h++) {
				w_exists = hidexists[h];
				
				// If the hidden unit exists
				if (mPruneWeights) {
					w_exists = static_cast<const BinaryGene&> (
						(*this)[(CONSTR)format ("W%d-%d", i, h)]).getvalue();
				}
				
				// Create the connection if it exists
				if (w_exists)
					net->connect (i, h+mInputs);
			}
		}
	}

	//
	// Connect hiddens to outputs
	//

	for (int o=0; o<mOutputs; o++) {
	
		// Go trough each (existing) hidden unit
		for (int h=0; h<mMaxHidden; h++) {
			w_exists = hidexists[h];
			
			// If the hidden unit exists
			if (mPruneWeights) {
				w_exists = static_cast<const BinaryGene&> (
					(*this)[(CONSTR)format ("W%d-%d", o, h)]).getvalue();
			}
			
			// Create the connection if it exists
			if (w_exists)
				net->connect (h+mInputs, o+mInputs+mMaxHidden);
		}
	}

	// OStream out2;
	// *net >> out2;

	net->init (0.5);
	msg.host.set ("brainplan", net);

	return true;
}
