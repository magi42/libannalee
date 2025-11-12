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

#include <nhp/individual.h>
#include <inanna/annetwork.h>
#include <inanna/topology.h>
#include <inanna/rprop.h>
#include <magic/mclass.h>

#include "annalee/anngenes.h"

impl_abstract (ANNGene, {Gentainer});
impl_abstract (ANNEncoding, {Gentainer});

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                   _   |   | |   |  ----                                   //
//                  / \  |\  | |\  | |      ___    _    ___                  //
//                 /   \ | \ | | \ | | --- /   ) |/ \  /   )                 //
//                 |---| |  \| |  \| |   \ |---  |   | |---                  //
//                 |   | |   | |   | |___/  \__  |   |  \__                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

void ANNGene::addPrivateGenes (Gentainer& parent, const StringMap& params) {
	Gentainer::addPrivateGenes (parent, params);
}

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//         _   |   | |   | -----                      | o                   //
//        / \  |\  | |\  | |       _    ___           |     _               //
//       /   \ | \ | | \ | |---  |/ \  |   \  __   ---| | |/ \   ___        //
//       |---| |  \| |  \| |     |   | |     /  \ (   | | |   | (   \       //
//       |   | |   | |   | |____ |   |  \__/ \__/  ---| | |   |  ---/       //
//                                                               __/        //
//////////////////////////////////////////////////////////////////////////////

ANNEncoding::ANNEncoding (const GeneticID& name,
						  const StringMap& params)
		: Gentainer (name)
{
	ASSERT (!isnull(params["inputs"]) && !isnull(params["outputs"]));

	mInputs            = params["inputs"].toInt ();
	mMaxHidden         = params["ANNEncoding.maxHidden"].toInt ();
	mPrunePassthroughs = params["ANNEncoding.prunePassthroughs"].toInt ();
	mOutputs           = params["outputs"].toInt ();
}

ANNEncoding::ANNEncoding (const ANNEncoding& other) : Gentainer (other) {
	mInputs            = other.mInputs;
	mMaxHidden         = other.mMaxHidden;
	mOutputs           = other.mOutputs;
	mPrunePassthroughs = other.mPrunePassthroughs;
}

void ANNEncoding::copy (const Genstruct& o) {
	Gentainer::copy (o);
	const ANNEncoding& other = static_cast<const ANNEncoding&>(o);

	mInputs            = other.mInputs;
	mMaxHidden         = other.mMaxHidden;
	mOutputs           = other.mOutputs;
	mPrunePassthroughs = other.mPrunePassthroughs;
}

void ANNEncoding::check () const {
	Gentainer::check ();

	ASSERT (mInputs    >  0);
	ASSERT (mInputs    <  1000);   // Sensible upper limit
	ASSERT (mMaxHidden >= 0);
	ASSERT (mMaxHidden <  100000); // Sensible upper limit
	ASSERT (mOutputs   >  0);
	ASSERT (mOutputs   <  1000);   // Sensible upper limit
}



