/***************************************************************************
 *   This file is part of the Annalee library.                             *
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

/*******************************************************************************
   Strategiassa henkisen suuntasi ei tule olla erilainen
   normaalista. Sekä taistelussa että jokapäiväisessä elämässä sinun
   tulisi olla päättäväinen vaikkakin tyyni. Kohtaa tilanne ilman
   jännittyneisyyttä, ei kuitenkaan kärsimättömästi, henkesi ollessa
   asettunut muttei herpaantunut. Edes kun henkesi on tyyni älä anna
   ruumiisi rentoutua, ja kun ruumiisi on rentoutunut, älä anna
   henkesi herpaantua. Älä anna ruumiisi vaikuttaa henkeesi.  Älä ole
   alihenkinen äläkä ylihenkinen. Ylennetty henki on heikko ja
   alhainen henki on heikko. Älä anna vihollisesi nähdä henkeäsi.

   Pienten ihmisten täytyy tuntea täydellisesti isojen ihmisten henki
   ja isojen ihmisten täytyy tuntea pienikokoisten ihmisten
   henki. Mikä hyvänsä on kokosi, älä anna oman ruumiisi reaktioiden
   hämätä sinua.
*******************************************************************************/

#ifndef __ANNALEE_ANNGENES_H__
#define __ANNALEE_ANNGENES_H__

#include <nhp/genetics.h>
#include <nhp/genes.h>
#include <magic/mpackarray.h>

// Externals
class PatternSet;

// An active gene is the best gene



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                   _   |   | |   |  ----                                   //
//                  / \  |\  | |\  | |      ___    _    ___                  //
//                 /   \ | \ | | \ | | --- /   ) |/ \  /   )                 //
//                 |---| |  \| |  \| |   \ |---  |   | |---                  //
//                 |   | |   | |   | |___/  \__  |   |  \__                  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

/*******************************************************************************
 * Gene for constructing an artificial neural network for an individual.
 *
 * This is the training object, not the network itself, which is
 * constructed by ANNEncoding. The gene is usually given id "brain",
 * and its phenotypic feature in the Individual has the same name.
 ******************************************************************************/
class ANNGene : public Gentainer { // Abstract
	decl_dynamic (ANNGene);
  public:
						ANNGene				(const GeneticID& name=NULL) : Gentainer (name) {;}
						ANNGene				(const ANNGene& orig) : Gentainer (orig) {;}
	
	/** Implementation for @ref Genstruct. */
	virtual void		init				() {Gentainer::init ();}
	/** Implementation for @ref Genstruct. */
	virtual void		copy				(const Genstruct& other) {Gentainer::copy (other);}
	/** Implementation for @ref Genstruct. */
	virtual void		addPrivateGenes		(Gentainer& g, const StringMap& params);
};



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//         _   |   | |   | -----                      | o                   //
//        / \  |\  | |\  | |       _    ___           |     _               //
//       /   \ | \ | | \ | |---  |/ \  |   \  __   ---| | |/ \   ___        //
//       |---| |  \| |  \| |     |   | |     /  \ (   | | |   | (   \       //
//       |   | |   | |   | |____ |   |  \__/ \__/  ---| | |   |  ---/       //
//                                                               __/        //
//////////////////////////////////////////////////////////////////////////////

/*******************************************************************************
 * Baseclass for the neural network encoding method genes.
 *
 * Any ANN encodings should be built under this abstract class. The
 * encoding genes should be instantiated with ID "brainplan".
 *
 * During the ontogenesis, the "brainplan" gene constructs the topology of
 * the neural network of the Individual that is going to be evaluated with
 * the evaluation data in the learning environment.
 ******************************************************************************/
class ANNEncoding : public Gentainer { // Abstract
  public:
	
						ANNEncoding	(const GeneticID& name, const StringMap& params);
						ANNEncoding	(const ANNEncoding& other);

	// Implementations

	/** Implementation for @ref Genstruct. */
	virtual void		copy				(const Genstruct& other);
	/** Implementation for @ref Genstruct. */
	virtual void		addPrivateGenes		(Gentainer& g, const StringMap& params) {MUST_OVERLOAD}
	/** Implementation for @ref Object. */
	virtual void		check				() const;
	
  protected:
						ANNEncoding	() {FORBIDDEN}

	int		mInputs, mMaxHidden, mOutputs;
	bool	mPrunePassthroughs;
	decl_dynamic (ANNEncoding);
};



// Any man who wants to master the essence of my strategy must
// research diligently, training morning and evening. Thus can he
// polish his skill, become free from self, and realise extraordinary
// ability. He will come to possess miraculous power. -- Miyamoto Musashi

// Macros are the biggest enemy of the objectkind

/*******************************************************************************
* A genetic message for instructing the individual to save a
* snapshot of its "brain" to a log file.
*******************************************************************************/
class TakeBrainPicsMsg : public GeneticMsg {
  public:
	TakeBrainPicsMsg (const GeneticID& rcvr, Individual& ind) : GeneticMsg (rcvr, ind) {}
};

#endif


