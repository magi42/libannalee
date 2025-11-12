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

#include "annalee/cangelosi.h"
#include "annalee/cangelosinet.h"

impl_dynamic (CangelosiEncoding, {NolfiEncoding});



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//      ___                          |            o -----                   //
//     /   \  ___    _          ___  |       ____   |       _    ___        //
//     |      ___| |/ \   ___  /   ) |  __  (     | |---  |/ \  |   \       //
//     |     (   | |   | (   \ |---  | /  \  \__  | |     |   | |           //
//     \___/  \__| |   |  ---/  \__  | \__/ ____) | |____ |   |  \__/ O     //
//                        __/                                               //
//////////////////////////////////////////////////////////////////////////////

CangelosiEncoding::CangelosiEncoding (const GeneticID& name, const StringMap& params) : NolfiEncoding (name, params) {
}

CangelosiEncoding::CangelosiEncoding (const CangelosiEncoding& other) : NolfiEncoding (other) {
}

void CangelosiEncoding::copy (const Genstruct& o) {
	NolfiEncoding::copy (o);
	// const CangelosiEncoding& other = static_cast<const CangelosiEncoding&>(o);
}

void CangelosiEncoding::addPrivateGenes (Gentainer& g, const StringMap& params) {
	Gentainer::addPrivateGenes (g, params);

	// The rewriting-rules
	CangCellDescr::addGenesTo (*this, params);

	// This is exactly as in NolfiEncoding
	if (params["NolfiEncoding.tipRadius"] == "auto-network")
		add (&(new BitFloatGene	("tipr", 1, 10, 8, params))->hide());
}

bool CangelosiEncoding::execute (const GeneticMsg& msg) const {

	// Read rules from the genome
	Array<CangCellDescr> rules (16*2);
	for (int rule=0; rule<16; rule++)
		for (int daughter=0; daughter<2; daughter++)
			rules.put (new CangCellDescr (*this, rule, daughter), rule*2+daughter);

	// Fetch genome-global tip radius, if it is encoded
	double tipRadius = mTipRadius;
	if (getGene("tipr"))
		tipRadius	= ((const AnyFloatGene*) getGene("tipr"))->getvalue();
	
	// Rewrite for some cycles
	CangelosiNet cnet (mInputs, mMaxHidden, mOutputs, mXSize, mYSize, tipRadius, mAxonScale);
	cnet.rewrite (rules, int(log(mMaxHidden*1.0)/log(2.0)+0.99));

	// Build the network
	ANNetwork* net = cnet.growNet ();
		
	// Add the plan to the host
	if (net) {
		net->setInitializer (new GaussianInitializer ());
		
		// Take pictures only if this is a picture-taking recreation
		bool takePics = dynamic_cast<const TakeBrainPicsMsg*>(&msg) != NULL;
		if (takePics) {
			msg.mrHost.set ("brainpic1", new String (cnet.drawEPS()));
			msg.mrHost.set ("brainpic2", new String (net->drawEPS()));
		}
		net->cleanup (true, mPrunePassthroughs);
		if (takePics) {
			net->drawFeedForward();
			msg.mrHost.set ("brainpic3", new String (net->drawEPS()));
			delete net; // The net was created only for taking babypics
		} else {
			// Place the brain description into host
			msg.mrHost.set ("brainplan", net);
		}
	}

	return true;
}



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//     ___                     ___        | | ___                            //
//    /   \  ___    _         /   \  ___  | | |  \   ___   ____  ___         //
//    |      ___| |/ \   ___  |     /   ) | | |   | /   ) (     |   \ |/\    //
//    |     (   | |   | (   \ |     |---  | | |   | |---   \__  |     |      //
//    \___/  \__| |   |  ---/ \___/  \__  | | |__/   \__  ____)  \__/ |      //
//                       __/                                                 //
///////////////////////////////////////////////////////////////////////////////

void CangCellDescr::addGenesTo (Gentainer& g, const StringMap& params) {
	Array<String> slrange;
	params["CangelosiEncoding.segLenMulRange"].split (slrange, ',');
	double sMin = slrange[0].toDouble ();
	double sMax = slrange[1].toDouble ();

	ASSERT (sMin>=-2 && sMin<=1 && sMax>0 && sMax<=5);

	// Insert the genes for the rules
	for (int rule=0; rule<16; rule++)
		for (int daughter=0; daughter<2; daughter++) {
			String num = format ("%d%c", rule, daughter+'a');
			g.add (&(new BitIntGene (String("T")+num, 0, 15, 4, params))->hide());
			g.add (&(new BitFloatGene (String("b")+num, -1, 1, 10, params))->hide());
			g.add (&(new BitFloatGene (String("w")+num, -1, 1, 10, params))->hide());
			g.add (&(new BitIntGene (String("d")+num, 0, 7, 3, params))->hide());
			if (params["NolfiEncoding.faceGene"].toInt ())
				g.add (&(new BitFloatGene (String("f")+num, 0, 1, 4, params))->hide());
			g.add (&(new BitFloatGene (String("s")+num, sMin, sMax, 4, params))->hide());
			g.add (&(new BitFloatGene (String("a")+num, -1, 1, 6, params))->hide());
			if (params["NolfiEncoding.tipRadius"] == "auto-cell")
				g.add (&(new BitFloatGene (String("r")+num, 0, 2, 8, params))->hide());
		}
}

void CangCellDescr::decodeFrom (const Gentainer& g, int r, int d) {
	char c = d + 'a';
	mNewType		= ((const AnyIntGene&)		g[format ("T%d%c", r, c)]).getvalue();
	mBiasVar		= ((const AnyFloatGene&)	g[format ("b%d%c", r, c)]).getvalue();
	mWeightVar		= ((const AnyFloatGene&)	g[format ("w%d%c", r, c)]).getvalue();
	mDaughterLoc	= ((const AnyIntGene&)		g[format ("d%d%c", r, c)]).getvalue();
	mSegLengthVar	= ((const AnyFloatGene&)	g[format ("s%d%c", r, c)]).getvalue();
	mSegAngleVar	= ((const AnyFloatGene&)	g[format ("a%d%c", r, c)]).getvalue();
	if (!isnull(g[format ("f%d%c", r, c)]))
		mFaceVar	= ((const AnyFloatGene&)	g[format ("f%d%c", r, c)]).getvalue();
	else
		mFaceVar	= 0;
	if (!isnull(g[format ("r%d%c", r, c)])) {
		mTipRadiusMul = ((const AnyFloatGene&) g[format ("r%d%c", r, c)]).getvalue();
		TRACELINE;
	} else
		mTipRadiusMul = -666;
}



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//        ___                          |            o  ___        | |       //
//       /   \  ___    _          ___  |       ____   /   \  ___  | |       //
//       |      ___| |/ \   ___  /   ) |  __  (     | |     /   ) | |       //
//       |     (   | |   | (   \ |---  | /  \  \__  | |     |---  | |       //
//       \___/  \__| |   |  ---/  \__  | \__/ ____) | \___/  \__  | |       //
//                          __/                                             //
//////////////////////////////////////////////////////////////////////////////

CangelosiCell::CangelosiCell (const CangelosiCell& mother, const CangCellDescr& rule) : NolfiCell (mother) {
	// This cell is now a clone of its mother. Now we make
	// modifications according to the given rule
	add (rule);
	mFace = 0.0;
	if (mTipRadius<0.5)
		mTipRadius = 1.0;
}

void CangelosiCell::make () {
	NolfiCell::make ();

	// This method is called only for constructing the mother, so we
	// initialize the mother cell a bit differently from the NolfiCell.
	mFace = 0.0;
	mSegmentLength = 0.0;
	mSegmentAngle = 0.0;
	if (mTipRadius<0.5)
		mTipRadius = 1.0;
}

void CangelosiCell::copy (const CangelosiCell& o) {
	NolfiCell::copy (o);
	mFace = o.mFace;
}

void CangelosiCell::add (const CangCellDescr& rule) {
	static double daughterlocx[] = {0,1,1,1,0,-1,-1,-1};
	static double daughterlocy[] = {1,1,0,-1,-1,-1,0,1};

	mTypeID = rule.mNewType;
	mBias += rule.mBiasVar;
	mWeight += rule.mWeightVar;
	mFace += rule.mFaceVar;
	mSegmentLength += rule.mSegLengthVar;
	mSegmentAngle += 0.1*rule.mSegAngleVar;
	mCoord.x += daughterlocx[rule.mDaughterLoc];
	mCoord.y += daughterlocy[rule.mDaughterLoc];
	if (rule.mTipRadiusMul>0)
		mTipRadius *= rule.mTipRadiusMul;
}

OStream& CangelosiCell::operator>> (OStream& out) const {
	out.printf ("face=%+02.2f, ", mFace);
	NolfiCell::operator>> (out);
	return out;
}



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//        ___                          |            o |   |                 //
//       /   \  ___    _          ___  |       ____   |\  |  ___   |        //
//       |      ___| |/ \   ___  /   ) |  __  (     | | \ | /   ) -+-       //
//       |     (   | |   | (   \ |---  | /  \  \__  | |  \| |---   |        //
//       \___/  \__| |   |  ---/  \__  | \__/ ____) | |   |  \__    \       //
//                          __/                                             //
//////////////////////////////////////////////////////////////////////////////

CangelosiNet::CangelosiNet (int inputs, int maxhidden, int outputs, double xsize, double ysize, double tipRadius, double axonScale) : NolfiNet (inputs, maxhidden, outputs, xsize, ysize, tipRadius, axonScale) {
	// mSize = ((mInputs>mOutputs)? mInputs+1: mOutputs+1);
	makeMother ();
}

void CangelosiNet::makeMother () {
	cells.make (0);
	cells.add (new CangelosiCell ());
	//cells[0].setPos (mXSize/2, mXSize/2);
	cells[0].setPos (0.0, 0.0);
	if (mTipRadius>0)
		cells[0].setTipRadius(mTipRadius);
}

void CangelosiNet::rewrite (const Array<CangCellDescr>& rules, int cycles) {
	for (int cycle=0; cycle<cycles; cycle++) {
		Array<CangelosiCell> daughters (2 * cells.size());

		// Rewrite once
		for (int mother=0; mother < cells.size(); mother++) {
			const CangelosiCell& mom = static_cast<const CangelosiCell&>(cells[mother]);
			for (int daughter=0; daughter<2; daughter++)
				daughters.put (new CangelosiCell (mom, rules[mom.mTypeID*2+daughter]),
							   mother*2+daughter);
		}

		// Move the daughters to the cell array
		cells.empty ();
		cells.make (daughters.size());
		for (int i=0; i < daughters.size(); i++) {
			cells.put (&daughters[i], i);
			daughters.cut (i);
		}
	}

	// Normalize the coordinates to range [0,1]
	for (int i=0; i<cells.size(); i++) {
		const Coord2D& coord = cells[i].pos ();
		cells[i].setPos (Coord2D(coord.x/(cycles*2+1)+0.5, coord.y/(cycles*2+1)+0.5));
	}
}
