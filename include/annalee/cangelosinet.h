#ifndef __CANGELOSINET_H__
#define __CANGELOSINET_H__

#include "cangelosi.h"
#include "nolfinet.h"

class CangCellDescr;
class CangelosiCell;
class CangelosiNet;



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//     ___                     ___        | | ___                            //
//    /   \  ___    _         /   \  ___  | | |  \   ___   ____  ___         //
//    |      ___| |/ \   ___  |     /   ) | | |   | /   ) (     |   \ |/\    //
//    |     (   | |   | (   \ |     |---  | | |   | |---   \__  |     |      //
//    \___/  \__| |   |  ---/ \___/  \__  | | |__/   \__  ____)  \__/ |      //
//                       __/                                                 //
///////////////////////////////////////////////////////////////////////////////

/** Decoded daughter cell descriptor for Cangelosi-type neuron used in
 * @ref CangelosiEncoding. This is an internal, temporary class; it is
 * used during the ontogenetic process.
 *
 * The rewriting rules have form C -> CC, where C is a cell. This
 * class implements the definition of C, but not it's instance.
 *
 * Notice that the rewriting replaces some parameters of the mother
 * cell totally, while some parameters it only modifies (such as the
 * x,y-location).
 **/
class CangCellDescr : public Object {
  public:

						CangCellDescr	() {}

	/** Standard constructor.
	 *
	 * See decodeFrom() for the description of the parameters.
	 **/
						CangCellDescr	(const Gentainer& g, int rule, int daughter) {
							decodeFrom (g, rule, daughter);
						}

	/** Decodes the cell description from the given genome.
	 *
	 * @param g Genome that contains the genetic code for the cell
	 * @param rule Tells the type-specifier for the mother of this cell.
	 * @param daughter Tells which daughter this cell is.
	 *                 With this and the above rule-number, the genetic code for this
	 *                 cell can be fetched from the genome.
	 **/
	void				decodeFrom		(const Gentainer& g, int rule, int daughter);
	
	/** Add the genes that this class requires to the given genome.
	 **/
	static void			addGenesTo		(Gentainer& g, const StringMap& params);

  protected:
	int		mNewType;
	double	mBiasVar;
	double	mWeightVar;
	int		mDaughterLoc;
	double	mFaceVar;
	double	mSegLengthVar;
	double	mSegAngleVar;
	double	mTipRadiusMul;

	friend class CangelosiCell;
};



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//        ___                          |            o  ___        | |       //
//       /   \  ___    _          ___  |       ____   /   \  ___  | |       //
//       |      ___| |/ \   ___  /   ) |  __  (     | |     /   ) | |       //
//       |     (   | |   | (   \ |---  | /  \  \__  | |     |---  | |       //
//       \___/  \__| |   |  ---/  \__  | \__/ ____) | \___/  \__  | |       //
//                          __/                                             //
//////////////////////////////////////////////////////////////////////////////

/** A temporary cell used in the ontogeny process of the @ref
 * CangelosiEncoding. These cells have ability to divide, sprout axon
 * trees, and then mature as neurons.
 **/
class CangelosiCell : public NolfiCell {
  public:

							CangelosiCell	() {make ();}

	/** Creates the cell as a daughter cell of the given mother which
		will be modified by the given rule.
	**/
							CangelosiCell	(const CangelosiCell& mother,
											 const CangCellDescr& rule);

	/** Makes this cell a mother cell.
	 **/
	void					make			();
	
	/** Modifies the cell according to the given rule.
	 **/
	void					add				(const CangCellDescr& rule);

	/** Build connections for the neuron by modifying the final network object.
	 **/
	void					connect			(const Array<CangelosiCell>& cells,
											 ANNetwork& ann) const;

	// Implementations
	
	/** Implementation */
	virtual void			copy			(const CangelosiCell& other);

	/** Implementation */
	virtual OStream&		operator>>		(OStream& out) const;

  protected:
	/** Direction angle of axon */
	double	mFace;
	
	friend class CangelosiNet;
};



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//        ___                          |            o |   |                 //
//       /   \  ___    _          ___  |       ____   |\  |  ___   |        //
//       |      ___| |/ \   ___  /   ) |  __  (     | | \ | /   ) -+-       //
//       |     (   | |   | (   \ |---  | /  \  \__  | |  \| |---   |        //
//       \___/  \__| |   |  ---/  \__  | \__/ ____) | |   |  \__    \       //
//                          __/                                             //
//////////////////////////////////////////////////////////////////////////////

/** Temporary network object representing the 2-dimensional cell space
 * used by @ref CangelosiEncoding.
 *
 * This is otherwise same as the parent class @ref NolfiNet, except
 * that the cells in the cell space are generated using an iterated
 * rewriting process.
 **/
class CangelosiNet : public NolfiNet {
  public:

	/** Constructor for the cell space.
	 *
	 * See the constructor of the parent class @ref NolfiNet for
	 * description of the parameters.
	 **/
							CangelosiNet	(int inputs, int maxhidden, int outputs,
											 double xsize, double ysize, double tipRadius,
											 double axonScale);

	/** Rewrites all the cells in the cell space for a specified
	 * number of iterations. Typically just the mother cell exists in
	 * the cell space when this is called.
	 *
	 * @param rules Rules, indexed by the left-hand-side value *
	 * 2. The right-hand-side value for rule #n is pair (#n*2, #n*2+1).
	 *
	 * @param cycles Number of rewriting cycles.
	 **/
	void					rewrite			(const Array<CangCellDescr>& rules,
											 int cycles=5);

  protected:

	/** Creates the mother cell in the middle of the cell space. The
     *  cell space should be empty when this is called.
	**/
	void					makeMother		();
};

#endif

