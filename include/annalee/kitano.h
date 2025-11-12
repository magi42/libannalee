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

#ifndef __KITANO_H__
#define __KITANO_H__

#include <magic/mtable.h>
#include <annalee/anngenes.h>

// Externals
class ANNetwork;


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  |  / o                      -----                      | o              //
//  | /     |   ___    _        |       _    ___           |     _          //
//  |/   | -+-  ___| |/ \   __  |---  |/ \  |   \  __   ---| | |/ \   ___   //
//  | \  |  |  (   | |   | /  \ |     |   | |     /  \ (   | | |   | (   \  //
//  |  \ |   \  \__| |   | \__/ |____ |   |  \__/ \__/  ---| | |   |  ---/  //
//                                                                    __/   //
//////////////////////////////////////////////////////////////////////////////

/** The gene for the Graph Generation Grammar encoding method by
 * Kitano (1990).
 **/
class KitanoEncoding : public ANNEncoding {
	decl_dynamic (KitanoEncoding);
  public:

						KitanoEncoding () {FORBIDDEN}
	
	/** Standard constructor.
	 *
	 * @param name Name of the gene; always "brainplan".
	 * @param params Dynamic parameters in a @ref String @ref Map.
	 * @param params["rewrites"] Number of rewriting iterations. Typically 3-10. The maximum number of neurons is 2^rewrites.
	 * @param params["rules"] Number of rewriting rules in the genome. Typically 32 or 64.
	 * @param params["nonTerminals"] Number of nonterminals, for example A-Z. Typical value is 26.
	 **/
						KitanoEncoding		(const GeneticID& name,
											 const StringMap& params);
						KitanoEncoding		(const KitanoEncoding& orig);
	
	// Implementations

	/** Implementation for @ref Genstruct. */
	virtual Genstruct*	replicate			() const {return new KitanoEncoding (*this);}
	/** Implementation for @ref Genstruct. */
	virtual void		copy				(const Genstruct& other);
	/** Implementation for @ref Genstruct. */
	virtual bool		execute				(const GeneticMsg& msg) const;
	/** Implementation for @ref Genstruct. */
	virtual void		addPrivateGenes		(Gentainer& g, const StringMap& params);
	/** Implementation for @ref Object. */
	virtual void		check				() const;
	
  private:
	int	mIters;
	int	mNonTerminals;
	int mRules;

	/** Exceptional values in the rewriting matrix */
	enum strangevalues {VOIDAREA=-1, FINALZERO=-2, FINALONE=-3, UNRESOLVED=-4};

	static Gentainer*	makeGenes			();

	/** Internal recursive decoding function for constructing a
	 * connection matrix.
	 *
	 * @return A connection matrix.
	 **/
	PackTable<int>*		decodeMatrix		(const PackTable<int>& string,
											 const PackTable<int>& rules,
											 int l) const;

	/** Eats a connection matrix and returns a corresponding freenetwork.
	 **/
	ANNetwork*		makeNet				(const PackTable<int>& connmat) const;
};

#endif
