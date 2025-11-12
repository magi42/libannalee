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

#include <magic/mgdev-eps.h>
#include <magic/mlsystem.h>
#include <magic/mturtle.h>
#include <magic/mtextstream.h>

#include <inanna/annetwork.h>
#include "nhp/individual.h"

#include "annalee/nolfi.h"
#include "annalee/nolfinet.h"




//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//              |   |      |     o -----               |                    //
//              |\  |      |  __     |              |  |  ___               //
//              | \ |  __  | /   |   |   |   | |/\ -+- | /   )              //
//              |  \| /  \ | +-- |   |   |   | |    |  | |---               //
//              |   | \__/ | |   |   |    \__! |     \ |  \__               //
//                           |                                              //
//////////////////////////////////////////////////////////////////////////////

// This turtle device calculates the axon tip positions for neurons in
// the Nolfi encoding
class NolfiTurtle : public TurtleDevice {
	Array<Coord2D>	mTips;
	int				mTipPos;
  public:

			NolfiTurtle	(int tipEstimate=64) {
				mTips.make (tipEstimate);
				mTipPos=0;
			}

	virtual void	forwardLine	(const Coord2D& s, const Coord2D& e) {}

	void	tip			(const Coord2D& tp) {
		if (mTipPos>=mTips.size())
			ASSERT (false);
		// mTips.resize (int(mTips.size*1.5));
		mTips[mTipPos++] = tp;
	}

	void	getTips		(Array<Coord2D>& tips) const {
		tips.make (mTipPos);
		for (int i=0; i<tips.size(); i++)
			tips[i] = mTips[i];
	}

};



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//               ----- ----   ---- -----               |                     //
//               |     |   ) (       |              |  |  ___                //
//               |---  |---   ---    |   |   | |/\ -+- | /   )               //
//               |     |         )   |   |   | |    |  | |---                //
//               |____ |     ___/    |    \__! |     \ |  \__                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// This turtle is used for making nice pictures of the Nolfi networks
class EPSTurtle : public TurtleDevice {
	EPSDevice&	mDevice;
	double		mXScale, mYScale;
	double		mTipRadius;
  public:
				EPSTurtle	() : mDevice ((EPSDevice&)*((EPSDevice*)NULL)) {}
				EPSTurtle	(EPSDevice& dev, double tipR) : mDevice(dev), mTipRadius(tipR) {}
	
	void		forwardLine	(const Coord2D& start, const Coord2D& end) {
		mDevice.line (start,end);
	}

	void		tip			(const Coord2D& tipPoint) {
		mDevice.saveState ().setGray(0.5);
		mDevice.circle (tipPoint, mTipRadius);
		mDevice.restoreState ();
	}
};



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                    |   |      |     o  ___        | |                    //
//                    |\  |      |  __   /   \  ___  | |                    //
//                    | \ |  __  | /   | |     /   ) | |                    //
//                    |  \| /  \ | +-- | |     |---  | |                    //
//                    |   | \__/ | |   | \___/  \__  | |                    //
//                                 |                                        //
//////////////////////////////////////////////////////////////////////////////

// Creates the axon description string (see below) for NolfiCells
// using an L-System.
String makeAxonString () {
	String result;

	LGrammar grammar;
	grammar.addRule ("X", "F[-X][+X]");
	result = "X";
	grammar.applyTo (result, 4);

	// Change the Xs to Fs
	LGrammar grammar2;
	grammar2.addRule ("X", "F");
	grammar2.applyTo (result, 1);

	return result;
}

// The axon description string for NolfiCells.  The letters in the
// string are used for guiding the NolfiTurtle (or EPSTurtle) that
// determines the locations of axon tips.
String NolfiCell::smAxonString = makeAxonString ();

/*******************************************************************************
 * Initializes the cell.
 ******************************************************************************/
void NolfiCell::make ()
{
	mExpression = true;
	mCoord.x = 0.0;
	mCoord.y = 0.0;
	mSegmentAngle = 0.5;
	mSegmentLength = 0.5;
	mWeight = 0.5;
	mBias = 0.5;
	mTypeID = 0;
	mTipRadius = 0.5;
	
	mFinalID = EMPTYID;
	mFinalType = CT_NONE;
}

void NolfiCell::copy (const NolfiCell& o) {
	mExpression = o.mExpression;
	mCoord.x = o.mCoord.x;
	mCoord.y = o.mCoord.y;
	mSegmentAngle = o.mSegmentAngle;
	mSegmentLength = o.mSegmentLength;
	mWeight = o.mWeight;
	mBias = o.mBias;
	mTypeID = o.mTypeID;
	mTipRadius = o.mTipRadius;

	mFinalID = o.mFinalID;
	mFinalType = o.mFinalType;
}


/*******************************************************************************
 * Add the genes for a single cell (with given index) to the given
 * genome.
 *
 * @param g Genetic container where the genes are placed. This is
 * typically @ref NolfiEncoding.
 *
 * @param i Index number of the gene. It will be used during
 * ontogenesis to access the gene.
 *
 * @param types The number of different neuron types in the
 * genome. The values are always a power of two; typical values
 * are 8, 16, 32, 64.
 *
 * @param xsize X-dimension of the cellular space of the encoding
 * method. A power of two. Typical value is 16.
 *
 * @param ysize Y-dimension of the cellular space of the encoding
 * method. Typical value is 16.
 *
 * @param params Other params in a @ref String @ref Map.
 **/
void NolfiCell::addGenesTo (
	Gentainer&       g,
	int              i,
	int              types,
	int              xsize,
	int              ysize,
	const StringMap& params)
{
	if (params["NolfiEncoding.existenceGene"].toInt())
		g.add (&(new BinaryGene	(format ("e%d", i)))->hide());
	g.add (&(new BitFloatGene	(format ("x%d", i), 0, 1, 3, params))->hide());
	g.add (&(new BitFloatGene	(format ("y%d", i), 0, 1, 5, params))->hide());
	g.add (&(new BitFloatGene	(format ("a%d", i), -1, 1, 6, params))->hide());
	g.add (&(new BitFloatGene	(format ("s%d", i), 0, 1, 4, params))->hide());
	g.add (&(new BitFloatGene	(format ("w%d", i), -1, 1, 10, params))->hide());
	g.add (&(new BitFloatGene	(format ("b%d", i), -1, 1, 10, params))->hide());
	int typebits = int(log(types*1.0)/log(2.0)+.999);
	g.add (&(new BitIntGene		(format ("t%d", i), 0, (1<<typebits)-1, typebits, params))->hide());
	if (params["NolfiEncoding.tipRadius"]=="auto-cell")
		g.add (&(new BitFloatGene	(format ("r%d", i), 1, 10, 8, params))->hide());
}

/*******************************************************************************
 * Read the cell with given index description from given genetic
 * container.
 *
 * @param g Genetic container that contains the genes to individualize the cell.
 *
 * @param i The identifier number for the cell. Same as in @ref addGenesTo.
 ******************************************************************************/
void NolfiCell::decodeFrom (const Gentainer& g, int i) {
	if (!isnull(g[(CONSTR)format ("e%d", i)]))
		mExpression	= ((const BinaryGene&)		g[(CONSTR)format ("e%d", i)]).getvalue();
	else
		mExpression	= true;
	mCoord.x		= ((const AnyFloatGene&)	g[(CONSTR)format ("x%d", i)]).getvalue();
	mCoord.y		= ((const AnyFloatGene&)	g[(CONSTR)format ("y%d", i)]).getvalue();
	mBias			= ((const AnyFloatGene&)	g[(CONSTR)format ("b%d", i)]).getvalue();
	mWeight			= ((const AnyFloatGene&)	g[(CONSTR)format ("w%d", i)]).getvalue();
	mSegmentLength	= ((const AnyFloatGene&)	g[(CONSTR)format ("s%d", i)]).getvalue();
	mSegmentAngle	= ((const AnyFloatGene&)	g[(CONSTR)format ("a%d", i)]).getvalue();
	mTypeID			= ((const AnyIntGene&)		g[(CONSTR)format ("t%d", i)]).getvalue();
	if (!isnull(g[(CONSTR)format ("r%d", i)]))
		mTipRadius	= ((const AnyFloatGene&)	g[(CONSTR)format ("r%d", i)]).getvalue();

	mFinalID		= EMPTYID;
	mFinalType		= CT_NONE;
	//sout << *this; sout.print("\n");
}

/*******************************************************************************
 *
 ******************************************************************************/
OStream& NolfiCell::operator>> (OStream& out) const
{
	out.printf ("e=%01d, b=%+02.2f, w=%+02.2f, seglen=%02.2f, segang=%+02.2f, "
				"x=%+02.2f, y=%+02.2f, t=%02d, n=%03d, nt=%02d, r=%02.2f",
				int(mExpression), mBias, mWeight, mSegmentLength, mSegmentAngle,
				mCoord.x, mCoord.y, mTypeID, mFinalID, mFinalType, mTipRadius);
	return out;
}

/*******************************************************************************
 * Grows an "axon tree" L-System from the cell.
 ******************************************************************************/
void NolfiCell::developAxon (
	Array<Coord2D>& result, //< An array where the coordinates of the axon tips are stored.
	double          scale   //< parameter tells usually the size of the neural space.
	) const
{
	if (true) {
		NolfiTurtle turtleDevice;
		Turtle turtle (turtleDevice, scale*mSegmentLength, mSegmentAngle*180/M_PI);
		turtle.jumpTo (mCoord+Coord2D(0.5,0));
		turtle.drawLSystem (smAxonString);
		turtleDevice.getTips (result);
	} else {
		result.make (11);
		for (int d=0; d<=10; d++) {
			double angle = mSegmentAngle*(d-5)*M_PI/10.0;
			result[d].x = mCoord.x+scale*fabs(mSegmentLength)*cos(angle);
			result[d].y = mCoord.y+scale*fabs(mSegmentLength)*sin(angle);
		}
	}
}

void NolfiCell::drawEPS (EPSDevice& devcon, double scale) const {
	// Cell body
	devcon.circle (mCoord, 0.5, mExpression);
	
	// Axon tree
	if (mExpression) {
		EPSTurtle turtleDevice (devcon, mTipRadius);
		Turtle turtle (turtleDevice, scale*mSegmentLength, mSegmentAngle*180/M_PI);
		turtle.jumpTo (mCoord+Coord2D(0.5,0));
		turtle.drawLSystem (smAxonString);
	}
}

void NolfiCell::check () const {
	ASSERTWITH (mTipRadius>=0.5, "Internal error");
}



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                    |   |      |     o |   |                              //
//                    |\  |      |  __   |\  |  ___   |                     //
//                    | \ |  __  | /   | | \ | /   ) -+-                    //
//                    |  \| /  \ | +-- | |  \| |---   |                     //
//                    |   | \__/ | |   | |   |  \__    \                    //
//                                 |                                        //
//////////////////////////////////////////////////////////////////////////////

/*******************************************************************************
 * Standard constructor.
 *
 * @param tipRadius Global axon tip radius value. This value is not
 * used if the tip radius is encoded in the genome for each separate
 * cell ("auto-cell").
 ******************************************************************************/
NolfiNet::NolfiNet (
	int    inputs,     //< Number of inputs in the final network.
	int    maxhidden,  //< Maximum number of hidden neurons in the final network.
	int    outputs,    //< Number of output neurons in the final network.
	double xsize,      //< Length of the x-dimension of the cell space.
	double ysize,      //< Length of the y-dimension of the cell space.
	double tipRadius,
	double axonScale)  //< Axon length scaling factor. Typically 1.0.
{
	mInputs = inputs;
	mMaxHidden = maxhidden;
	mHiddens = 0;
	mOutputs = outputs;
	mXSize = xsize;
	mYSize = ysize;
	mSize = ((mInputs>mOutputs)? mInputs+1: mOutputs+1);
	mInputBorder = 0.3;
	mOutputBorder = 0.7;
	mTipRadius = tipRadius;
	mAxonScale = axonScale;
}

/*******************************************************************************
 * Decodes the rewriting rules from the given genome.
 ******************************************************************************/
void NolfiNet::decodeFrom (const Gentainer& g)
{
	// Fetch genome-global tip radius, if it is encoded
	if (!isnull(g["tipr"]))
		mTipRadius	= ((const AnyFloatGene&)	g["tipr"]).getvalue();

	cells.make (mMaxHidden);
	for (int i=0; i<cells.size(); i++) {
		cells[i].setTipRadius (mTipRadius); // This can be changed by the decodeFrom below
		cells[i].decodeFrom (g, i);
		cells[i].check ();
	}
	
}

OStream& NolfiNet::operator>> (OStream& out) const
{
	out.printf ("{%d-%d(%d)-%d (=%d), %fx%f}\n",
				mInputs, mMaxHidden, mHiddens, mOutputs,
				cells.size(), mXSize, mYSize);

	for (int i=0; i<cells.size(); i++)
		out << cells[i] << "\n";
	return out;
}

bool NolfiNet::indexCells () {

	// Resolve the type for each cell and count the types
	int inputs, hiddens, outputs;
	resolveTypes (inputs, hiddens, outputs);
	// TRACE3 ("%d-%d-%d", inputs, hiddens, outputs);
	
	// If there are no neurons of some neuron type, abort
	if (inputs==0 || hiddens==0 || outputs==0 || outputs<mOutputs)
		return false;

	indexInputs ();
	indexHiddens (hiddens);
	indexOutputs (hiddens);

	// Remove the cells with no index
	//for (int i=0; i<cells.size; i++)
	//	if (cells[i].mFinalID==EMPTYID)
	//		ASSERT (false);
	
	mHiddens = hiddens;
	return true;
}

void NolfiNet::resolveTypes (int& inputs, int& hiddens, int& outputs) {
	inputs = hiddens = outputs = 0;

	// Resolve the type of the cell according to it's X coordinate

	for (int i=0; i<cells.size(); i++)
		if (cells[i].mExpression) {
			if (cells[i].mCoord.x<mInputBorder) {
				// Input unit.
				if (cells[i].mCoord.x>=0 && cells[i].mCoord.y>=0 && cells[i].mCoord.y<=1) {
					cells[i].mFinalType = CT_INPUT;
					inputs++;
				}
			} else if (cells[i].mCoord.x>mOutputBorder) {
				// Output unit.
				if (cells[i].mCoord.x<=1 && cells[i].mCoord.y>=0 && cells[i].mCoord.y<=1) {
					cells[i].mFinalType = CT_OUTPUT;
					outputs++;
				}
			} else {
				// Hidden unit.
				if (cells[i].mCoord.y>=0 && cells[i].mCoord.y<=1) {
					cells[i].mFinalType = CT_HIDDEN;
					hiddens++;
				} else {
					sout << cells[i] << "\n";
					ASSERT (false);
				}
			}
		}
}

void NolfiNet::indexInputs () {
}

void NolfiNet::indexHiddens (int hiddens) {
	// Index ALL units according to their X-position
	for (int i=0; i<cells.size(); i++) {
		// Find the unindexed hidden cell with the smallest X coordinate
		double minX = 666.0;
		int minNeuron = -1;
		
		for (int j=0; j<cells.size(); j++)
			if (cells[j].mCoord.x<=minX && cells[j].mFinalID==EMPTYID) {
				minX = cells[j].mCoord.x;
				minNeuron = j;
			}
		
		// Set the index of the smallest found unindexed cell
		cells[minNeuron].mFinalID = i+mInputs;
	}
}

void NolfiNet::indexOutputs (int hiddens) {
}

void NolfiNet::removeDuplicates (int& rOutputs)
{
	// Remove input and output cells with duplicate indices
	for (int i=0; i<cells.size(); i++)
		for (int j=i+1; j<cells.size(); j++)
			if (cells[i].mFinalID!=EMPTYID && cells[i].mFinalID == cells[j].mFinalID) {
				cells[j].mFinalID = EMPTYID;
				// We want to calculate the number of output units exactly
				if (cells[i].mFinalType==CT_OUTPUT)
					rOutputs--;
			}
}


/*******************************************************************************
 * Executes the ontogeny of the cell space; starts rewriting from a
 * single mother cell, rewrites a few times, grows the axons, etc.
 *
 * @return The resulting phenotypic network topology, ready to be
 * inserted as the individual's "brainplan" attribute.
 ******************************************************************************/
ANNetwork* NolfiNet::growNet ()
{
	if (!indexCells ())
		return NULL;

	// Do some scaling
	for (int i=0; i<cells.size(); i++) {
		// Scale the coordinates
		//cells[i].mCoord.x = (mYSize/mXSize)*(0.25+int(mXSize*cells[i].mCoord.x));
		cells[i].mCoord.x *= mYSize;
		cells[i].mCoord.y *= mYSize;
		
		// Scale the axon segment length. Since there are five
		// segments, this leads to axons about ½ the size of the neural
		// space
		cells[i].mSegmentLength *= 0.1*mYSize;
	}

	// :DEBUG:
	/*
	sout << *this;
	String pic = this->drawEPS();
	FILE* fout = fopen ("cangelosipic.eps", "w");
	fprintf (fout, "%s", (CONSTR) pic);
	fclose (fout);
	*/
	
	// Create the network object
	ANNetwork* result = new ANNetwork (format("%d-%d-%d", mInputs, cells.size(), mOutputs));

	int connections = connect (*result);
	
	if (connections==0) {
		delete result;
		return NULL;
	}

	// Set final input units to be linear (propably unnecessary)
	for (int i=0; i<mInputs; i++)
		(*result)[i].setTFunc (Neuron::LINEAR_TF);

	// Set the coordinates for output units
	double spacing = mYSize/mOutputs;
	double startY = spacing/2+(mSize-mYSize)/2.0;
	for (int i=0; i<mOutputs; i++)
		(*result)[cells.size()-mOutputs+i].moveTo (Coord2D(9+mYSize, startY+i*spacing));
	
	return result;
}

int NolfiNet::connect (ANNetwork& network) const
{
	int connections=0;
	Coord2D scaling ((mYSize+1)/mYSize, mSize/mYSize);
	for (int i=0; i<cells.size(); i++) {
		//(sout << cells[i]).print("\n");

		Neuron& neuron = network[cells[i].mFinalID];

		// Set the neuron attributes
		neuron.setBias (cells[i].mBias);
		neuron.moveTo (cells[i].mCoord*scaling+Coord2D(3,0.5));

		if (!cells[i].mExpression) {
			neuron.enable (false);
			continue;
		}

		// If an "input unit", make linear connection from a final
		// input unit
		if (cells[i].mFinalType==CT_INPUT) {
			// neuron.setTFunc (FreeNeuron::LINEAR_TF);
			network.connect (cells[i].mTypeID % mInputs, cells[i].mFinalID);
			//TRACE1 ("%d", cells[i].mTypeID);
		}

		// If an "output unit", make linear connection to a final
		// output unit
		if (cells[i].mFinalType==CT_OUTPUT) {
			// neuron.setTFunc (FreeNeuron::LINEAR_TF);
			network.connect (cells[i].mFinalID,
							 network.size()-mOutputs+(cells[i].mTypeID % mOutputs));
		}
		
		// Generate branch tip points
		Array<Coord2D> tips;
		cells[i].developAxon (tips, mAxonScale); //mXSize

		// Now find other cells that lie near these points
		for (int j=0; j<cells.size()-mOutputs; j++)
			if (cells[j].mFinalID > cells[i].mFinalID && cells[j].mFinalType!=CT_INPUT
				&& !(cells[i].mFinalType==CT_OUTPUT && cells[j].mFinalType==CT_OUTPUT))
				for (int k=0; k<tips.size(); k++) {
					double d = sqrt(cells[j].mCoord.sqdist (tips[k]));
					if (d < cells[i].mTipRadius) {
						// TRACE4 ("%d->%d: d(%d)=%f", i, j, k, d);
						if (!neuron.connectedFrom (network[cells[j].mFinalID])) {
							// Found a new target, add it to the network
							network.connect (cells[i].mFinalID, cells[j].mFinalID);
							connections++;
							break;	// Break from the k loop
						}
					}
				}
		
		// Set the initial connection weights to the encoded weight
		for (int j=0; j<neuron.incomings(); j++)
			neuron.incoming(j).setWeight (cells[i].mWeight);
	}

	return connections;
}


/*******************************************************************************
 * Draws an EPS image of the cell space, with the axon trees and all.
 *
 * @return The EPS image in a string.
 ******************************************************************************/
String NolfiNet::drawEPS () const
{
	Coord2D picSize (175, 175);
	EPSDevice epsdevice (picSize);					// Graphics device
	epsdevice.framedStyle (mYSize, mYSize);			// Clipping with frame
	
	// Draw input/output region borders
	epsdevice.lineWidth (0);
	epsdevice.lineStyle ("dashed", mYSize/picSize.y);
	epsdevice.line (Coord2D(mYSize*mInputBorder, 0),
					Coord2D(mYSize*mInputBorder, mYSize));
	epsdevice.line (Coord2D(mYSize*mOutputBorder, 0),
					Coord2D(mYSize*mOutputBorder, mYSize));
	epsdevice.lineStyle ("solid");

	// Draw cells
	for (int i=0; i<cells.size(); i++)
		cells[i].drawEPS (epsdevice, mAxonScale);

	epsdevice.printFooter ();
	return epsdevice.getBuffer ();
}

/*
CString NolfiNet::drawLatex () const {
	CString result;
	result.reserve (10000);
	double picXSize=176.0, picYSize=152;
	result = format ("\\psset{xunit=1pt,yunit=1pt}\n"
					 "\\psclip{\\psframe(0,0)(%f,%f)}\n"
					 "\\psset{fillstyle=solid,fillcolor=white}\n"
					 "\\psline[linestyle=dashed](%f,0)(%f,%f)\n"
					 "\\psline[linestyle=dashed](%f,0)(%f,%f)\n",
					 picXSize, picYSize,
					 picXSize*mInputBorder, picXSize*mInputBorder, picYSize,
					 picXSize*0.75, picXSize*0.75, picYSize
		);
	
		double xscale = picXSize/mXSize; 
	double yscale = picYSize/mYSize;
	
	for (int i=0; i<cells.size; i++) {
		// Cell body
		result += format ("\\cnode(%f,%f){5pt}{N%d}\n",
						  cells[i].mCoord.x*xscale, picYSize-cells[i].mCoord.y*yscale, i);

		result += cells[i].drawLatex (picYSize, yscale);

		// Generate branch tip points
		PackArray<Coord2D> tips;
		cells[i].developAxon (tips, mXSize);

		for (int j=0; j<tips.size; j++)
			result += format ("\\psline[linestyle=dotted](%f,%f)(%f,%f)\n",
							  cells[i].mCoord.x*xscale, cells[i].mCoord.y*yscale,
							  tips[j].x*xscale, tips[j].y*yscale);
						
	}
	
	result += format ("\\endpsclip\n");

	return result;
}

*/
