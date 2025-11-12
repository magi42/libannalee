#include <magic/mgdev-eps.h>
#include <magic/mclass.h>
#include <magic/mlsystem.h>
#include <magic/mturtle.h>

#include <inanna/annetwork.h>
#include <inanna/initializer.h>
#include <nhp/individual.h>

#include "annalee/nolfi.h"
#include "annalee/nolfinet.h"

impl_dynamic (NolfiEncoding, {ANNEncoding});



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//       |   |      |     o -----                      | o                   //
//       |\  |      |  __   |       _    ___           |     _               //
//       | \ |  __  | /   | |---  |/ \  |   \  __   ---| | |/ \   ___        //
//       |  \| /  \ | +-- | |     |   | |     /  \ (   | | |   | (   \       //
//       |   | \__/ | |   | |____ |   |  \__/ \__/  ---| | |   |  ---/       //
//                    |                                           __/        //
///////////////////////////////////////////////////////////////////////////////

NolfiEncoding::NolfiEncoding (const GeneticID& name, const StringMap& params) : ANNEncoding (name, params) {
	/*if (isnull(params["NolfiEncoding.types"])) {
		mTypes = (mInputs>mOutputs)? mInputs:mOutputs;
		if (mTypes<16)
			mTypes = 16;
	} else
	*/
	mTypes = params["NolfiEncoding.types"].toInt ();

	mXSize		= getOrDefault (params, "NolfiEncoding.xSize", String(9)).toInt ();
	mYSize		= getOrDefault (params, "NolfiEncoding.ySize", String(22)).toInt ();
	mAxonScale	= getOrDefault (params, "NolfiEncoding.axonScale", String(0.5)).toDouble ();
	mTipRadius	= getOrDefault (params, "NolfiEncoding.tipRadius", String(0.5)).toDouble ();
	mMaxHidden	= getOrDefault (params, "NolfiEncoding.neurons", String(0.5)).toInt ();

	ASSERTWITH (mTipRadius>=0.5 || params["NolfiEncoding.tipRadius"]=="auto-network"
				|| params["NolfiEncoding.tipRadius"]=="auto-cell",
				format ("NolfiEncoding.tipRadius must be a real-valued radius >=0.5, "
						"\042auto-genome\042, or \042auto-cell\042; "
						"not \042%s\042", (CONSTR) params["NolfiEncoding.tipRadius"]));
}

NolfiEncoding::NolfiEncoding (const NolfiEncoding& other) : ANNEncoding (other) {
	mTypes		= other.mTypes;
	mXSize		= other.mXSize;
	mYSize		= other.mYSize;
	mAxonScale	= other.mAxonScale;
	mTipRadius	= other.mTipRadius;
}

void NolfiEncoding::copy (const Genstruct& o) {
	ANNEncoding::copy (o);
	const NolfiEncoding& other = static_cast<const NolfiEncoding&>(o);
	mTypes		= other.mTypes;
	mXSize		= other.mXSize;
	mYSize		= other.mYSize;
	mAxonScale	= other.mAxonScale;
	mTipRadius	= other.mTipRadius;
}

void NolfiEncoding::addPrivateGenes (Gentainer& g, const StringMap& params) {
	Gentainer::addPrivateGenes (g, params);
	
	//StringMap params;
	//params.set("graycoding","0");
	for (int i=0; i<mMaxHidden; i++)
		NolfiCell::addGenesTo (*this, i, mTypes, mXSize, mYSize, params);

	if (params["NolfiEncoding.tipRadius"] == "auto-network")
		add (&(new BitFloatGene	("tipr", 1, 10, 8, params))->hide());
}

bool NolfiEncoding::execute (const GeneticMsg& msg) const {

	// Read the cells from the genome
	NolfiNet nnet (mInputs, mMaxHidden, mOutputs, mXSize, mYSize, mTipRadius, mAxonScale);
	nnet.decodeFrom (*this);

	// Build the network
	ANNetwork* net = nnet.growNet ();

	//TRACE1 ("%s", (CONSTR) (net->getLayering().toString()));
	
	// Add the plan to the host
	if (net) {
		net->setInitializer (new GaussianInitializer ());

		// Take pictures only if this is a picture-taking recreation
		bool takePics = dynamic_cast<const TakeBrainPicsMsg*>(&msg) != NULL;

		if (takePics) {
			msg.mrHost.set ("brainpic1", new String (nnet.drawEPS()));
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
