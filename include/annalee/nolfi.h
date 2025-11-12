#ifndef __NOLFI_H__
#define __NOLFI_H__

#include "anngenes.h"

// Internals
class NolfiEncoding;


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//       |   |      |     o -----                      | o                   //
//       |\  |      |  __   |       _    ___           |     _               //
//       | \ |  __  | /   | |---  |/ \  |   \  __   ---| | |/ \   ___        //
//       |  \| /  \ | +-- | |     |   | |     /  \ (   | | |   | (   \       //
//       |   | \__/ | |   | |____ |   |  \__/ \__/  ---| | |   |  ---/       //
//                    |                                           __/        //
///////////////////////////////////////////////////////////////////////////////

/** Encoding method by Nolfi and Parisi (1992). It uses a cell space
 * where the potential neurons are located. The neurons grow axon
 * trees, and if the tips of the axon trees touch other neurons,
 * connections are made.
 *
 * Each of the neurons is encoded separately in the genome. Parameters
 * include position in the cell space as well as the parameters
 * dictating the shape of the axon tree.
 *
 * Our implementation encodes only the topology of the network; the
 * weights are not encoded, but are learned by a separate neural
 * training algorithm.
 **/
class NolfiEncoding : public ANNEncoding {
	decl_dynamic (NolfiEncoding);
  public:
						NolfiEncoding () {FORBIDDEN}
	
	/** Standard constructor.
	 * @param name The gene name. Usually "brainplan".
	 * @param params Dynamic parameter @ref String @ref Map.
	 * @param params["existenceGene"] Boolean parameter that controls whether or not each neuron has an additional expression gene that tells if the neuron is enabled or not. This gene was used in the original method, but it may cause unnecessary ...
	 * @param params["faceGene"] Should the "face" gene be used, that allows axon trees to sprout to any direction. [Default: 0]
	 * @param params["neurons"] Number of encoded cells in the genome; the maximum number of neurons in the phenotypic neural network. Typical values are 32 or 64.
	 * @param params["types"] Number of neuron types. Must be a power of 2; typically 8,16,32,64. This value MUST exceed the number of input and output neurons.
	 * @param params["tipRadius"] Connection radius of axon tips. Must be >=0.5, typical values are 0.5..10. Special values "auto-cell" or "auto-network" cause the tip radius to be encoded in the genome. The value "auto-network" uses a common gene for all the cells in the genome; "auto-cell" uses a separate gene for every cell.
	 * @param params["axonScale"] Axon size scale coefficient. [Default=1.0]
	 * @param params["xSize"] Size of x-dimension of the cell space. [Default=9]
	 * @param params["ySize"] Size of y-dimension of the cell space. [Default=22]
	 **/
						NolfiEncoding		(const GeneticID& name,
											 const StringMap& params);
						NolfiEncoding		(const NolfiEncoding& other);
	
	/** Implementation for @ref Genstruct. */
	virtual Genstruct*	replicate			() const {return new NolfiEncoding (*this);}
	/** Implementation for @ref Genstruct. */
	virtual void		copy				(const Genstruct& other);
	/** Implementation for @ref Genstruct. */
	virtual bool		execute				(const GeneticMsg& msg) const;
	/** Implementation for @ref Genstruct. */
	virtual void		addPrivateGenes		(Gentainer& g, const StringMap& params);

  protected:
	int		mTypes;			// Number of cell identities. Originally 16
	int		mXSize;			// X-size of the cell grid. Originally 7
	int		mYSize;			// Y-size of the cell grid. Originally 21
	double	mTipRadius;		// Tip radius, default: 1.0
	double	mAxonScale;		// Axon size scaling. Default: 1.0
	
};


#endif
