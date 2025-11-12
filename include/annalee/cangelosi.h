#ifndef __CANGELOSI_H__
#define __CANGELOSI_H__

#include "anngenes.h"
#include "nolfi.h"

// Internals
class CangelosiEncoding;



//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//      ___                          |            o -----                   //
//     /   \  ___    _          ___  |       ____   |       _    ___        //
//     |      ___| |/ \   ___  /   ) |  __  (     | |---  |/ \  |   \       //
//     |     (   | |   | (   \ |---  | /  \  \__  | |     |   | |           //
//     \___/  \__| |   |  ---/  \__  | \__/ ____) | |____ |   |  \__/ O     //
//                        __/                                               //
//////////////////////////////////////////////////////////////////////////////

/** The encoding method by Cangelosi, Nolfi and Parisi (1994).
 *
 * The encoding is base on the @ref NolfiEncoding, except that the
 * cells in the cell space are not encoded directly. Instead, they are
 * generated using cell-rewriting.
 *
 * The cell rewriting starts from a single "mother" cell in the cell
 * space. It is rewritten into two daughters according to rewriting
 * rules stored in the genome. The rewriting is continued for certain
 * number of iterations. After that, the axon trees are generated just
 * like in @ref NolfiEncoding.
 **/
class CangelosiEncoding : public NolfiEncoding {
	decl_dynamic (CangelosiEncoding);
  public:
						CangelosiEncoding () {FORBIDDEN}
	
	/** Standard constructor, called by LearningGAEnv
	 *
	 * @param name Name of the gene.
	 * @param params Dynamic parameter @ref String @ref Map.
	 *               See @ref NolfiEncoding for parameters common with that class.
	 * @param params["segLenMulRange"] Segment length multiplier range for axon trees. Typical value: "0,0.5".
	 **/
						CangelosiEncoding	(const GeneticID& name,
											 const StringMap& params);
						CangelosiEncoding	(const CangelosiEncoding& other);
	
	// Implementations

	/** Implementation for @ref Genstruct. */
	virtual Genstruct*	replicate			() const {return new CangelosiEncoding (*this);}
	/** Implementation for @ref Genstruct. */
	virtual void		copy				(const Genstruct& other);
	/** Implementation for @ref Genstruct. */
	virtual bool		execute				(const GeneticMsg& msg) const;
	/** Implementation for @ref Genstruct. */
	virtual void		addPrivateGenes		(Gentainer& g, const StringMap& params);
};






#endif
