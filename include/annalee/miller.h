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

#ifndef __ANNALEE_MILLER_H__
#define __ANNALEE_MILLER_H__

#include "anngenes.h"

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//     |   | o | |           -----                      | o                 //
//     |\ /|   | |  ___      |       _    ___           |     _             //
//     | V | | | | /   ) |/\ |---  |/ \  |   \  __   ---| | |/ \   ___      //
//     | | | | | | |---  |   |     |   | |     /  \ (   | | |   | (   \     //
//     |   | | | |  \__  |   |____ |   |  \__/ \__/  ---| | |   |  ---/     //
//                                                                 __/      //
//////////////////////////////////////////////////////////////////////////////

/*******************************************************************************
 * This is a rather generic version of the direct encoding proposed by
 * Miller, Todd and Hedge (1989). It works on freely feed-forward
 * topologies.
 *
 * The basic idea in direct encoding is to have one binary gene for
 * each potential connection in the neural network that tells whether
 * or not that connection actually exists.
 *
 * In most implementations the connection genes are initialized (in
 * the initial generation) with a certain probability of connection
 * p(c). Often p(c) is simply 0.5, but can vary.
 *
 * Our implementation of the encoding method is somewhat different
 * from the original method by Miller et al.
 *
 * Most important modification in our implementation is the
 * initialization of the genomes in the first generation. Instead of
 * using a fixed probability of connection p(c), we initialize the
 * different genomes with different p(c). The p(c) of a certain genome
 * is determined from a normally distributed random value, controlled
 * by the pcVariance and pcAverage parameters (see below).
 ******************************************************************************/
class MillerEncoding : public ANNEncoding {
	decl_dynamic (MillerEncoding);
	bool	mPruneInputs;	// This is stored just for easy access
	double	mPcVariance;
	double	mPcAverage;

  public:
						MillerEncoding		() {FORBIDDEN}
						MillerEncoding		(const GeneticID& name,
											 const StringMap& params);
						MillerEncoding		(const MillerEncoding& other);
	
	// Implementations

	virtual Genstruct*	replicate			() const {return new MillerEncoding (*this);}
	virtual void		copy				(const Genstruct& other);
	virtual bool		execute				(const GeneticMsg& msg) const;
	virtual void		addPrivateGenes		(Gentainer& g, const StringMap& params);
	virtual void		init				();
};



#endif
