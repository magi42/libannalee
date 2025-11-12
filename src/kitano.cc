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

#include <magic/mclass.h>
#include <inanna/annetwork.h>
#include <inanna/initializer.h>
#include <nhp/individual.h>
#include "annalee/kitano.h"

impl_dynamic (KitanoEncoding, {Gentainer});

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  |  / o                      -----                      | o              //
//  | /     |   ___    _        |       _    ___           |     _          //
//  |/   | -+-  ___| |/ \   __  |---  |/ \  |   \  __   ---| | |/ \   ___   //
//  | \  |  |  (   | |   | /  \ |     |   | |     /  \ (   | | |   | (   \  //
//  |  \ |   \  \__| |   | \__/ |____ |   |  \__/ \__/  ---| | |   |  ---/  //
//                                                                    __/   //
//////////////////////////////////////////////////////////////////////////////

KitanoEncoding::KitanoEncoding (const GeneticID& name, const StringMap& params)
		: ANNEncoding (name, params)
{
	mIters        = getOrDefault (params, "KitanoEncoding.rewrites", String(5)).toInt ();
	mNonTerminals = getOrDefault (params, "KitanoEncoding.nonTerminals", String(26)).toInt ();
	mRules        = getOrDefault (params, "KitanoEncoding.rules", String(64)).toInt ();
}

KitanoEncoding::KitanoEncoding (const KitanoEncoding& orig) : ANNEncoding (orig) {
	mIters        = orig.mIters;
	mNonTerminals = orig.mNonTerminals;
	mRules        = orig.mRules;
}

void KitanoEncoding::copy (const Genstruct& o) {
	ANNEncoding::copy (o);
	const KitanoEncoding& orig = dynamic_cast<const KitanoEncoding&> (o);
	mIters        = orig.mIters;
	mNonTerminals = orig.mNonTerminals;
	mRules        = orig.mRules;
}

void KitanoEncoding::addPrivateGenes (Gentainer& p, const StringMap& params) {
	Gentainer::addPrivateGenes (p, params);

	for (int i=0; i<mRules; i++) {
		// Add nonterminal rule (N->NNNN, where N is a nont. symbol)
		// Left-hand-side:
		add (&(new IntGene (format ("R%d-0", i), 16, 16+mNonTerminals))->set(i).hide());
		for (int j=1; j<5; j++) // Right-hand-side:
			add (&(new IntGene (format ("R%d-%d", i, j), 0, 15+mNonTerminals))->hide());
	}
}

bool KitanoEncoding::execute (const GeneticMsg& msg) const
{
	// Cache the grammar into tables

	// One rule for each possible leftside.
	PackTable<int> rules (16+mNonTerminals,4);

	// Set terminal rules
	for (int i=0; i<16; i++)
		for (int j=0; j<4; j++)
			rules.get (i,j) = (i&(1<<j))? FINALONE : FINALZERO;
	
	// Set all nonterminal rules default to undefined value
	for (int i=0; i<mNonTerminals; i++)
		for (int j=0; j<4; j++)
			rules.get (16+i,j) = VOIDAREA;

	// Read the rules from the genome
	for (int i=16; i<rules.rows; i++) {
		int leftside=-1, symb=-1;
		for (int j=0; j<5; j++) {
			if (i==16 && j==0)
				symb = 16;	// The first rule in the chromosome is fixed
			else
				symb = static_cast<const IntGene&> (*getGene (
					format ("R%d-%d",i,j))).getvalue ();

			if (j==0) // Read the LHS of the production
				leftside = symb;
			else
				rules.get (leftside, j-1) = symb;
		}
	}

	//
	// Decode the grammar
	//
	
	// Make start rule
	PackTable<int> axiom (1,1);
	axiom.get (0,0) = 16; // 16 means the first nonterminal

	// Decode nonterminals recursively
	PackTable<int>* connmat = decodeMatrix (axiom, rules, mIters);

	// Calculate some statistics (probability of connection)
	int conns=0, totconns=0;
	for (int i=0; i<connmat->rows-mOutputs; i++)
		for (int j=(i>=mInputs)?i+1:mInputs; j<connmat->rows; j++) {
			totconns++;
			conns += (connmat->get(i,j)==1);
		}

	double pConn = double(conns) / double(totconns);
	msg.mrHost.set ("pConn", new String (format ("%f", pConn)));

	// Create a network from the connection matrix
	ANNetwork* net = makeNet (*connmat);

	if (net) {
		net->setInitializer (new GaussianInitializer (0.5));
		
		// Take baby pictures only if this is a picture-taking recreation
		bool takePics = dynamic_cast<const TakeBrainPicsMsg*>(&msg) != NULL;

		if (takePics) {
			net->cleanup (false);
			msg.mrHost.set ("brainpic1", new String (net->drawEPS()));
		}
		net->cleanup (true, mPrunePassthroughs);
		if (takePics) {
			msg.mrHost.set ("brainpic2", new String (net->drawEPS()));
			net->drawFeedForward();
			msg.mrHost.set ("brainpic3", new String (net->drawEPS()));

			const String& matStr = dynamic_cast<const String&> (net->getAttribute("matrix"));
			msg.mrHost.set ("braindesc1", new String (matStr));
			delete net;
		} else {
			// Place the brain description into host
			msg.mrHost.set ("brainplan", net);
		}
	}
	delete connmat;

	return true;
}

PackTable<int>* KitanoEncoding::decodeMatrix (const PackTable<int>& string, const PackTable<int>& rules, int l) const
{

	// Create a table to place the results of the next production
	PackTable<int>* nextString = new PackTable<int> (string.rows*2, string.cols*2);

	// Iterate through nonterminal matrix
	for (int i=0; i<string.rows; i++)
		for (int j=0; j<string.cols; j++) {
			// Decode a value in the matrix
			int k = string.get (i,j);
			
			// For nonterminals
			nextString->get (i*2  ,j*2  ) = (k<0)? k : rules.get (k, 0);
			nextString->get (i*2+1,j*2  ) = (k<0)? k : rules.get (k, 1);
			nextString->get (i*2  ,j*2+1) = (k<0)? k : rules.get (k, 2);
			nextString->get (i*2+1,j*2+1) = (k<0)? k : rules.get (k, 3);
		}
	
	// Recurse
	if (l>1) {
		PackTable<int>* retval = decodeMatrix (*nextString, rules, l-1);
		delete nextString;
		return retval;
	} else  {
		// Convert the temporary values to final
		for (int i=0; i<nextString->rows; i++)
			for (int j=0; j<nextString->cols; j++) {
				int& value = nextString->get(i,j);
				switch (value) {
				  case FINALONE:	value=1; break;
				  case FINALZERO:	value=0; break;
				  case VOIDAREA:	value=VOIDAREA; break;
				  default: 			value=UNRESOLVED; // Unresolved
				};
			}
		
		return nextString;
	}
}

ANNetwork* KitanoEncoding::makeNet (const PackTable<int>& connmatPar) const {
	ASSERT (connmatPar.rows == connmatPar.cols);
	ASSERT (connmatPar.rows != 0);
	ASSERT (connmatPar.rows > mInputs+mOutputs);

	PackTable<int> connmat = connmatPar;
	int hiddens = connmat.rows-(mInputs+mOutputs);
	
	// Outputs always enabled
	//for (int i=0; i<connmat.rows; i++)
	for (int i=mInputs+hiddens; i<connmat.rows; i++)
		connmat.get(i,i) = 1;

	String pic;
	pic.reserve(connmat.rows*(connmat.cols+1)+10);
 	for (int i=0; i<connmat.rows; i++) {
	 	for (int j=0; j<connmat.rows; j++) {
			switch (connmat.get(i,j)) {
			  case 0:			pic += '0';		break;
			  case 1:  			pic += '1';		break;
			  case VOIDAREA:	pic += ' ';		break;
			  case UNRESOLVED:	pic += 'x';		break;
			  default:
				  ASSERT (false);
			};
		}
	 	pic += '\n';
	} 
	
	ANNetwork* net = new ANNetwork (format("%dl-%d-%dl", mInputs, hiddens, mOutputs));
	net->setAttribute ("matrix", new String(pic));

	// Connect the network
	for (int i=0; i<connmat.rows; i++) {

		// Set the position of the unit
		if (i>=mInputs)
			(*net)[i].moveTo (double(i-mInputs)/hiddens*10+5, 20*frnd(), 20*frnd());

		// Disable unit if 0 at diagonal
		if (connmat.get(i,i)!=1) {
			(*net)[i].enable(false);
			(*net)[i].moveTo (-1,-1,-1);
		} else	// Connect the unit
			for (int j=mInputs; j<connmat.cols; j++)
				if (i<j && connmat.get(i,j)==1) // && i>=mInputs 
					net->connect (i,j);
	}

	net->check ();
	
	return net;
}

void KitanoEncoding::check () const {
	ANNEncoding::check ();
	ASSERT (mIters>0);
	ASSERT (mIters<10);
	ASSERT (mNonTerminals>0);
	ASSERT (mNonTerminals<100);
	ASSERT (mRules>=16);
	ASSERT (mRules<=1000);
}

///////////////////////////////////////////////////////////////////////////////
// Here is some code for statistical analysis of neural
// structures. This is not nice code and it exists for purely
// "debugging" purposes. More precisely, validation of propabilistic
// analysis of Kitano encoding.

void kitanoStat () {
}
