#ifndef __LAYERED_H__
#define __LAYERED_H__

#include "anngenes.h"

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//  |                                     | -----                      |    //
//  |      ___         ___       ___      | |       _    ___           |    //
//  |      ___| \   | /   ) |/\ /   )  ---| |---  |/ \  |   \  __   ---|    //
//  |     (   |  \  | |---  |   |---  (   | |     |   | |     /  \ (   |    //
//  |____  \__|   \_/  \__  |    \__   ---| |____ |   |  \__/ \__/  ---| O  //
//               \_/                                                        //
//////////////////////////////////////////////////////////////////////////////

/** A simple direct encoding method with three-layered feedforward MLP
 * structure.
 *
 * This straightforward gene allows encoding weights, biases, and the existences of
 * inputs, weights and hidden units of an MLP using direct encoding scheme. It is truly
 * marvellous.
 **/
class LayeredEncoding : public ANNEncoding {
	decl_dynamic (LayeredEncoding);
	bool		mPruneInputs;	// This is stored just for easy access
	bool		mPruneWeights;	// This is stored just for easy access
	bool		mEncodeWeights;
	Array<int>	mLayering;

  public:
						LayeredEncoding		() {FORBIDDEN}
						LayeredEncoding		(const GeneticID& name,
											 const StringMap& params);
						LayeredEncoding		(const LayeredEncoding& other);
	
	// Implementations

	/** Implementation for @ref Genstruct. */
	virtual Genstruct*	replicate			() const {return new LayeredEncoding (*this);}
	/** Implementation for @ref Genstruct. */
	virtual void		copy				(const Genstruct& other);
	/** Implementation for @ref Genstruct. */
	virtual bool		execute				(const GeneticMsg& msg) const;
	/** Implementation for @ref Genstruct. */
	virtual void		addPrivateGenes		(Gentainer& g, const StringMap& params);
};

#endif
