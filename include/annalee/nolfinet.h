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

#ifndef __NOLFINET_H__
#define __NOLFINET_H__

#include <magic/mcoord.h>
#include "annalee/nolfi.h"

class NolfiNet;

// Externals
namespace MagiC {
	class GDevice;
	class EPSDevice;
}


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                    |   |      |     o  ___        | |                    //
//                    |\  |      |  __   /   \  ___  | |                    //
//                    | \ |  __  | /   | |     /   ) | |                    //
//                    |  \| /  \ | +-- | |     |---  | |                    //
//                    |   | \__/ | |   | \___/  \__  | |                    //
//                                 |                                        //
//////////////////////////////////////////////////////////////////////////////

/** A temporary cell object used in the @ref NolfiEncoding encoding
 * method. The NolfiCells are always contained in a @ref NolfiNet. The
 * cell has an ability to sprout a L-System axon tree from itself,
 * where the shape of the tree is determined by parameters encoded in
 * the genome.
 **/
class NolfiCell : public Object {
  public:

							NolfiCell		() {make();}
							NolfiCell		(const NolfiCell& other) {copy (other);}

	
	static void				addGenesTo		(Gentainer& g, int i, int types,
											 int xsize, int ysize, const StringMap& params);
	void					make			();
	void					decodeFrom		(const Gentainer& g, int i);
	void					developAxon		(Array<Coord2D>& result, double scale) const;

	// Access functions

	/**
	 * @return The 2-dimensional coordinates of the cell in the cell
	 * space.
	**/
	const Coord2D&			pos				() const {return mCoord;}

	/** Sets the coordinates of the cell in the cell space.
	 * @see setPos(double x, double y)
	 **/
	void					setPos			(const Coord2D& npos) {mCoord=npos;}

	/** Sets the coordinates of the cell in the cell space.
	 * @see setPos(const Coord2D& npos)
	 **/
	void					setPos			(double x, double y) {mCoord.moveTo (x, y);}

	/** Sets the connection radius of the axon tips. This radius
	 * determines how far away the tip can be from another neuron to
	 * create an actual connection.
	 *
	 * @param r Typical values for the radius are 0.5...10.
	 **/
	void					setTipRadius	(double r) {mTipRadius = r;}

	/** The NolfiCell is an temporary object. The internal neuron
	 * index tells it's identity in the final network. If the value is
	 * EMPTYID (-1), the neuron will not exist in the final network.
	 **/
	void					mapTo			(int neuron) {mFinalID = neuron;}

	// Object handling and IO
	
	void					copy			(const NolfiCell& o);
	OStream&				operator>>		(OStream& out) const;

	/** Draws the image of the cell and it's axon tree into the given
	 * graphics driver.
	 *
	 * @param dc Encapsulates PostScript graphics driver.
	 * @param scale Drawing scale. The size of the image is multiplied by this factor.
	 **/
	void					drawEPS			(EPSDevice& dc, double scale) const;

	/** Implementation for @ref Object. */
	void					check			() const;

  protected:

	// Basic model attributes

	bool	mExpression;	//< Existence of the neuron.
	Coord2D	mCoord;			//< Coordinates.
	double	mSegmentAngle;	//< Axon segment angle.
	double	mSegmentLength;	//< Axon segment length.
	double	mWeight;		//< Initial weight for all the connections.
	double	mBias;			//< Initial bias.
	int		mTypeID;		//< I/O neuron semantical identifier.
	
	// Additional attributes

	double	mTipRadius;		//< Connection radius from axon tips and branching points.
	int		mFinalID;		//< Neuron index in the final brain.
	int		mFinalType;		//< input, hidden, output or none.

	friend class NolfiNet;

  private:
	static String	smAxonString; //< Axon description string generated with a L-System
};

enum celltypes {CT_NONE=-1, CT_INPUT=0, CT_HIDDEN=1, CT_OUTPUT=2};
enum cellids {EMPTYID=-1};


//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                    |   |      |     o |   |                              //
//                    |\  |      |  __   |\  |  ___   |                     //
//                    | \ |  __  | /   | | \ | /   ) -+-                    //
//                    |  \| /  \ | +-- | |  \| |---   |                     //
//                    |   | \__/ | |   | |   |  \__    \                    //
//                                 |                                        //
//////////////////////////////////////////////////////////////////////////////

/** Temporary network object representing the 2-dimensional cell space
 * used by @ref NolfiEncoding.
 **/
class NolfiNet : public Object {
  public:

							NolfiNet			(int inputs, int hiddens, int outputs,
												 double xsize, double ysize, double tipRadius,
												 double axonScale);

	virtual void			decodeFrom			(const Gentainer& g);
	virtual ANNetwork*		growNet				();

	// Implementations
	
	virtual OStream&		operator>>			(OStream& out) const;
	String					drawEPS				() const;
	
  protected:
	Array<NolfiCell>	cells;			//< The cells in the cell space
	double				mSize;			//< The size of the 2d-space
	int					mInputs;		//< Number of input cells in the final network
	int					mHiddens;		//< Number of hidden cells in the final network
	int					mMaxHidden;		//< Maximum number of hidden cells in the final network
	int					mOutputs;		//< Number of output cells in the final network
	double				mXSize;			//< X-size of the cell space.
	double				mYSize;			//< Y-size of the cell space.
	double				mInputBorder;	//< Size of the input-portion of the cell space. Default is 0.3.
	double				mOutputBorder;	//< 1 - size of the output-portion of the cell space. Default is 0.75.
	double				mTipRadius;		//< Global axon tip radius. Default is 0.5.
	double				mAxonScale;		//< Axon scaling factor. Default is 1.0.
	
  protected:
	/** Gives the cells a neuron index number and a neuron type
	 * specifier according to their location on the 2D-space.
	 *
	 * @return TRUE return value tells that the indexing has resulted in
	 * a network with at least one input unit and all output units.
	 **/
	virtual bool			indexCells			();

	/** Resolves types (input, hidden or output) for each cell
	 * according to their location in the cell space. The behaviour of
	 * this method is controlled by the mInputBorder and mOutputBorder
	 * variables.
	**/
	virtual void			resolveTypes		(int& inputs, int& hiddens, int& outputs);

	/** Gives input neurons an index value in the final network */
	virtual void			indexInputs			();

	/** Gives hidden neurons an index value in the final network */
	virtual void			indexHiddens		(int hiddens);

	/** Gives output neurons an index value in the final network */
	virtual void			indexOutputs		(int hiddens);

	/** Removes input and output cells that have duplicate ID values */
	virtual void			removeDuplicates	(int& rOutputs);

	/** Connects the units in a FreeNetwork according to the
	 * connections between the temporary cells.
	 *
	 * @return The number of connections in the network.
	 **/
	virtual int				connect				(ANNetwork& net) const;
};

#endif
