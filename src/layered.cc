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

#include <inanna/annetwork.h>
#include <inanna/initializer.h>
#include <nhp/individual.h>
#include <magic/mclass.h>

#include "annalee/layered.h"

impl_dynamic (LayeredEncoding, {ANNEncoding});


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  |                                     | -----                      |    //
//  |      ___         ___       ___      | |       _    ___           |    //
//  |      ___| \   | /   ) |/\ /   )  ---| |---  |/ \  |   \  __   ---|    //
//  |     (   |  \  | |---  |   |---  (   | |     |   | |     /  \ (   |    //
//  |____  \__|   \_/  \__  |    \__   ---| |____ |   |  \__/ \__/  ---| O  //
//               \_/                                                        //
//////////////////////////////////////////////////////////////////////////////

LayeredEncoding::LayeredEncoding (const GeneticID& name,
								  const StringMap& params) : ANNEncoding (name, params)
{
	mPruneInputs   = isnull(params["prune_inputs"])? true  : params["prune_inputs"].toInt ();
	mPruneWeights  = isnull(params["prune_weights"])? false : params["prune_weights"].toInt ();
	mEncodeWeights = isnull(params["encode_weights"])? false : params["encode_weights"].toInt ();

	// Parse layering description.
	String layering = params["layering"];
	Array<String> layers;
	layering.split(layers, '-');
	mLayering.make(layers.size());
	for (int i=0; i<layers.size(); i++)
		mLayering[i] = layers[i].toInt();
}

LayeredEncoding::LayeredEncoding (const LayeredEncoding& other) : ANNEncoding (other)
{
	mPruneInputs  = other.mPruneInputs;
	mPruneWeights = other.mPruneWeights;
	mLayering     = other.mLayering;
}

void LayeredEncoding::copy (const Genstruct& o)
{
	ANNEncoding::copy (o);
	const LayeredEncoding& other = static_cast<const LayeredEncoding&>(o);
	mPruneInputs  = other.mPruneInputs;
	mPruneWeights = other.mPruneWeights;
	mLayering     = other.mLayering;
}

void LayeredEncoding::addPrivateGenes (Gentainer& g, const StringMap& params)
{
	Gentainer::addPrivateGenes (g, params);

	if (mPruneInputs)
		for (int i=0; i<mInputs; i++)
			add (new BinaryGene (format ("R%d", i), 1.0));

	for (int layer=0; layer<mLayering.size(); layer++)
		for (int i=0; i<mInputs; i++)
			for (int j=0; j<mMaxHidden; j++) {
				if (mPruneWeights)
					add (new BinaryGene (format ("WX%d:%d-%d", layer, i, j), 1.0));
				if (mEncodeWeights)
					add (new FloatGene (format ("W%d:%d-%d", layer, i, j), 0.0, 1.0, 1.0));
			}

	// Prune hidden
	for (int i=0; i<mMaxHidden; i++)
		add (new BinaryGene (format ("H%d", i), 1.0));
}

bool LayeredEncoding::execute (const GeneticMsg& msg) const
{
	ANNetwork* net = new ANNetwork (format ("%d-%d-%d", mInputs, mMaxHidden, mOutputs));
	
	// Go trough each hidden unit and check if it exists
	bool hidexists [mMaxHidden];
	for (int h=0; h<mMaxHidden; h++) {
		hidexists[h] = static_cast<const BinaryGene&> (
			(*this)[(CONSTR)format ("H%d", h)]).getvalue();
		// TRACE2 ("%d=%d", h, int(hidexists[h]));
		
		// Enable or disable it from the network
		(*net)[h+mInputs].enable(hidexists[h]);
	}

	// Connect inputs to hiddens
	bool input_exists; // Does an input unit exist?
	bool w_exists;     // Does a weight exist?
	for (int i=0; i<mInputs; i++) {
		
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

	// Connect hiddens to outputs

	// To each output unit...
	for (int o=0; o<mOutputs; o++) {
	
		// ...connect every (existing) hidden unit
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

	msg.mrHost.set ("brainplan", net);

	return true;
}
