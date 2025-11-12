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

impl_dynamic (PureLayeredEncoding, {LayeredEncoding});


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  |                                     | -----                      |    //
//  |      ___         ___       ___      | |       _    ___           |    //
//  |      ___| \   | /   ) |/\ /   )  ---| |---  |/ \  |   \  __   ---|    //
//  |     (   |  \  | |---  |   |---  (   | |     |   | |     /  \ (   |    //
//  |____  \__|   \_/  \__  |    \__   ---| |____ |   |  \__/ \__/  ---| O  //
//               \_/                                                        //
//////////////////////////////////////////////////////////////////////////////

PureLayeredEncoding::PureLayeredEncoding (const GeneticID& name,
										  const StringMap& params)
		: LayeredEncoding (name, params)
{
}

PureLayeredEncoding::PureLayeredEncoding (const PureLayeredEncoding& other)
		: LayeredEncoding (other)
{
}

void PureLayeredEncoding::copy (const Genstruct& o)
{
	LayeredEncoding::copy (o);
}

void PureLayeredEncoding::addPrivateGenes (Gentainer& g, const StringMap& params)
{
	LayeredEncoding::addPrivateGenes(g, params);

	// Encode weights for each layer.
	for (int layer=0; layer<mLayers-1; i++)
		for (int i=0; i<mInputs; i++)
			for (int j=0; j<mMaxHidden; j++)
				add (new FloatGene (format ("W%d:%d-%d", layer, i, j)));
}

bool PureLayeredEncoding::execute (const GeneticMsg& msg) const
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
