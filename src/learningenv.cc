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

#include <magic/mmap.h>
#include <magic/mclass.h>
#include <magic/mtextstream.h>

#include <nhp/individual.h>
#include <inanna/annetwork.h>
#include <inanna/annfilef.h>
#include <inanna/patternset.h>
#include <inanna/rprop.h>

#include "annalee/learningenv.h"
#include "annalee/anngenes.h"
#include "annalee/layered.h"
#include "annalee/miller.h"
#include "annalee/cangelosi.h"
#include "annalee/kitano.h"
//#include "chaosenc.h"

impl_dynamic (LearningEAEnv, {EAEnvironment});



///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  |                           o              ----   _   -----              //
//  |      ___   ___        _       _         |      / \  |       _          //
//  |     /   )  ___| |/\ |/ \  | |/ \   ___  | --- /   \ |---  |/ \  |   |  //
//  |     |---  (   | |   |   | | |   | (   \ |   \ |---| |     |   |  \ /   //
//  |____  \__   \__| |   |   | | |   |  ---/ |___/ |   | |____ |   |   V    //
//                                       __/                                 //
///////////////////////////////////////////////////////////////////////////////

/*******************************************************************************
 * Dummy constructor, shouldn't be used (results in runtime error).
 *
 * Dummy. For some reason, we have to have a dummy constructor, so
 * this is what we get...
 ******************************************************************************/
LearningEAEnv::LearningEAEnv () :
		mTrainData ((PatternSet&) *new PatternSet()),
		mTrainSet ((PatternSet&) *new PatternSet()),
		mEvaluationSet ((PatternSet&) *new PatternSet()),
		mReportSet ((PatternSet&) *new PatternSet()),
		mParams((StringMap&) *new StringMap())
{
	FORBIDDEN; // But IT'S NEVER CALLED! (whew...)
}

/*******************************************************************************
 * Primary constructor. Most of the parameters are passed in a string map.
 *
 * @param trainset Used for local training of the ANN weights with a backprop algorithm.
 * @param evaluationset Used for evaluating the individuals to measure their fitness.
 *
 * @param reportset - Used for "final tests" in generation reports.
 * Although we take the training set and evaluation set separately,
 * the sets may be recombined and permutated later if so dictated by
 * the permutate-parameter (see below).
 *
 *  @param params Dynamic parameter map.
 *	@param params["evals"] - Minimum number of evaluations per individual per generation [Default=1]
 *	@param params["encoding"] - The name of the encoding method to be used: layered, miller, kitano, nolfi, cangelosi [No default - required]
 *  @param params["noise"] - Amount of artificial noise to be added [Default=0]
 *	@param params["permutate"] - Should we permutate the training and evaluation sets during evolution? [Default=0 (no)]
 *	@param params["evalPart"] - Portion of EA evaluation set as a fraction [Default=0.333]
 *	@param params["maxTrainCycles"] - Number of maximum training cycles. [Default=3000]
 *	@param params["stripLen"] - Interval between validations for early termination. [Default=10]
 *	@param params["terminator"] - Termination method and parameters. See @ref Terminator for defaults. [Default="GL5"]
 *	@param params["termPart"] - Termination set portion of training set as a fraction [Default=0.25]
 *	@param params["logDir"] - Logging directory [Default="log"]
 *	@param params["optParams"] - Should the learning parameters be optimized by evolution? [Default=0 (no)]
 ******************************************************************************/
LearningEAEnv::LearningEAEnv (const PatternSet& trainSet,
							  const PatternSet& evalSet,
							  const PatternSet& testSet,
							  StringMap& params)
		: mTrainSet      (dynamic_cast<const PatternSet&>(trainSet)),
		  mEvaluationSet (dynamic_cast<const PatternSet&>(evalSet)),
		  mReportSet     (dynamic_cast<const PatternSet&>(testSet)),
		  mParams        (params)
{
	ASSERTWITH (trainSet.patterns>0, "Must have train patterns");
	ASSERTWITH (evalSet.patterns>0, "Must have evaluation patterns");
	
	// Set some default values if not given explicitly
	mNEvals			= getOrDefault (mParams, "LearningEAEnv.evals", String(1)).toInt ();
	mNoise			= getOrDefault (mParams, "LearningEAEnv.noise", String(0.0)).toDouble ();
	mPermutate		= getOrDefault (mParams, "LearningEAEnv.permutate", String(0)).toInt ();
	mMaxTrainCycles	= getOrDefault (mParams, "LearningEAEnv.maxTrainCycles", String(3000)).toInt ();
	mReportCycles	= mMaxTrainCycles;
	mValidInterval	= getOrDefault (mParams, "LearningEAEnv.stripLen", String(10)).toInt ();
	mTermMethod		= getOrDefault (mParams, "LearningEAEnv.terminator", String("GL5"));
	mTermPart		= getOrDefault (mParams, "LearningEAEnv.termPart", String(0.25)).toInt ();
	logDir (getOrDefault (params, "logdir", String("log")));

	mEvalPart		= evalSet.patterns/double(evalSet.patterns+trainSet.patterns);
	mProblemType	= (trainSet.outputs>1)? CLASSIFICATION : CLASSIFICATION2;
	
	if (mTermMethod=="none")
		mTermPart = 0;

	// Join the given training and evaluation sets
	mTrainData.join (trainSet, evalSet);

	// Then split them again, just to be sure that the splitting is
	// done exactly the same way if the training data is permutated or
	// something
	splitTrainData ();

	ASSERTWITH (trainSet.patterns == mTrainSet.patterns,
			"Bug in algorithm: Uneven splitting of training data");

	ASSERT (mTrainSet.patterns>0);
	ASSERT (mEvaluationSet.patterns>0);
	ASSERT (mReportSet.patterns>0);
}

/*******************************************************************************
 * Implementation for @ref EAEnvironment.
 ******************************************************************************/
void LearningEAEnv::addFeaturesTo (Genome& genome) const
{
	// TODO: bool optimize_parameters = mParams["LearningEAEnv.optimizeParams"].toInt ();
	
	mParams.set ("inputs", String(mTrainData.inputs));
	mParams.set ("outputs", String(mTrainData.outputs));
	
	// Brain, according to the selected encoding scheme
	String encoding = mParams["LearningEAEnv.encoding"];
	if (encoding=="layered")
		genome.add (new LayeredEncoding ("brainplan", mParams));
	else if (encoding == "miller")
		genome.add (new MillerEncoding ("brainplan", mParams));
	else if (encoding == "nolfi")
		genome.add (new NolfiEncoding ("brainplan", mParams));
	else if (encoding == "cangelosi")
		genome.add (new CangelosiEncoding ("brainplan", mParams));
	else if (encoding == "kitano")
		genome.add (new KitanoEncoding ("brainplan", mParams));
	//	else if (encoding == "chaos")
	//		genome.add (new ChaosEncoding ("brainplan", mParams));
	else
		ASSERTWITH (false, format ("Unknown EANN encoding '%s'", (CONSTR) encoding));

	// At initialization of an individual, invoke the brainplan
	genome.add (new InterGene ("init", "brainplan"));
}

void LearningEAEnv::permutate () {
	FORBIDDEN; // Temporarily, we don't want anyone to use this
	
	// First permutate the whole dataset
	if (mProblemType==CLASSIFICATION2)
		mTrainData.recombine2 ();
	else
		mTrainData.recombine ();

	splitTrainData ();
}

/*******************************************************************************
 * Splits training data into a training set and an evaluation set.
 *
 * The training set is used for training the neural network, and evaluation
 * set for evaluating the fitness after training.
 ******************************************************************************/
void LearningEAEnv::splitTrainData ()
{
	mTrainSet.copy (mTrainData, 0, int(mTrainData.patterns*(1-mEvalPart))-1);
	mEvaluationSet.copy (mTrainData, int(mTrainData.patterns*(1-mEvalPart)),
						 mTrainData.patterns-1);
}

/*******************************************************************************
 * Evaluates an individual in the learning environment.
 *
 * Implementation for @ref EAEnvironment.
 *
 * Trains the individual with the training set and then evaluates the
 * it with the evaluation set.
 *
 * If early stopping is enabled, the training set is further divided
 * into an actual training set and termination set. If permutation of
 * training patterns is enabled, the patterns are shuffled before
 * division into actual training set and termination set.
 *
 * Prints statistics if they are enabled.
 *
 * @return Measured fitness of the individual.
 ******************************************************************************/
double LearningEAEnv::evaluateg (const Individual& ind)
{
	// Shuffle patterns if such is enabled
	if (mPermutate)
		permutate ();
	
	// Create a separate training set and GA evaluation set
	PatternSet trainSet, terminSet;
	trainSet.copy (mTrainSet, 0, int(mTrainSet.patterns*(1-mTermPart))-1);
	terminSet.copy (mTrainSet, int(mTrainSet.patterns*(1-mTermPart)), mTrainSet.patterns-1);

	// Get the I/O interface of the individual and set the parameters
	// which it doesn't know yet
	const ANNetwork* brainplan = dynamic_cast<const ANNetwork*> (ind.getFeature ("brainplan"));
	ASSERT (brainplan);
	//io.logDir (mLogDir);

	ANNetwork* brain = new ANNetwork (*brainplan);

	Trainer* pTrainer = createTrainer ();

	// Train the individual for a while
	double trainMSE = pTrainer->train (*brain,
									  trainSet,
									  mMaxTrainCycles,
									  &terminSet,
									  mValidInterval);
	delete pTrainer;
	trainMSE = 0.0; // Dispose. WARNING

	// Measure the fitness of the network with several criteria
	
	// Test with evaluation set
	double fitn_MSE = brain->test (mEvaluationSet);
	delete brain;

	double fitn_conns	= 0;
	double fitn_hiddens	= 0;
	double fitn_inputs	= 0;

	// Collect the factors together
	double fitness = fitn_MSE*1.0 + fitn_conns*0.0 + fitn_hiddens*0.0 + fitn_inputs*0.0;

	// Print some stats
	const Object& stats = ind["stats"];
	if (!isnull(stats))
		sout.printf (", stats=%s", (CONSTR) dynamic_cast<const String&>(stats));
	else
		sout.printf (", stats=0 0");
	const Object& pConn = ind["pConn"];
	if (!isnull(pConn))
		sout.printf (", pConn=%s", (CONSTR) dynamic_cast<const String&>(pConn));
	
	return fitness;
}

/*******************************************************************************
 *
 ******************************************************************************/
Trainer* LearningEAEnv::createTrainer () const
{
	StringMap trainParams;
	trainParams.set ("RPropTrainer.delta0", "1.0");
	RPropTrainer* pTrainer = new RPropTrainer();
	pTrainer->init (trainParams);
	pTrainer->setTerminator (mTermMethod);
	return pTrainer;
}

/*******************************************************************************
 *
 ******************************************************************************/
void LearningEAEnv::cycle_report (OStream& log, OStream& out) {
	
	// Get the I/O interface of the individual
	// best->execute (GeneticMsg ("IO", *best));
	ASSERT (mpBest);
	/*
	LearningIO& io = static_cast<LearningIO&> ((*best)["IO"]);

	// Determine where to log the picture of the King
	String cycleLogDir = mLogDir; // No generational logging
	if (0) { // TODO: Make this configurable if we want more logging
		cycleLogDir += format("/gens/gen%04d", mCycles);
		system (String("mkdir -p ")+cycleLogDir); // Create the directory
	}
	io.logDir (cycleLogDir);
	*/
	
	// Create a separate training set and GA evaluation set
	PatternSet trainSet, terminSet;
	trainSet.copy (mTrainSet, 0, int(mTrainSet.patterns*(1-mTermPart))-1);
	terminSet.copy (mTrainSet, int(mTrainSet.patterns*(1-mTermPart)),
					mTrainSet.patterns-1);
	
	ANNetwork& brain = dynamic_cast <ANNetwork&> ((*mpBest)["brainplan"]);
	
	Trainer* pTrainer = createTrainer ();
	
	// Train the individual for a while
	pTrainer->train (brain, trainSet, mMaxTrainCycles, &terminSet, mValidInterval);
	
	// Save this to a file
	ANNFileFormatLib::save (mLogDir + "/einstein.net", brain);
	//brain.saveBrain (mLogDir + "/einstein.net", "Best brain found by Annalee");
	
	// Save the pattern sets _only_ for the first cycle
	/*
	  if (mCycles==1) {
	  trainSet.save (mLogDir+"/train.pat");
	  terminSet.save (mLogDir+"/termin.pat");
	  mEvaluationSet.save (mLogDir+"/eval.pat");
	  ((PatternSet&)mReportSet).save (mLogDir+"/report.pat");
	  }
	*/
	
	// Excuse me, I'm looking for the Einstein Brain, can you tell me
	// where I can find the Einstein Brain, please?
	
	// Test with each evaluation pattern while counting the correct
	// predictions
	switch (mProblemType) {
	  case CLASSIFICATION:
	  case CLASSIFICATION2: {
		  //const ANNetwork& net = io.getNet ();
		  if (!isnull (brain)) {
			  ClassifResults* clsresults = brain.testClassify (mReportSet);
			  double perc = double(clsresults->failures) / double(mReportSet.patterns);
			  out.printf ("Number of incorrect predictions: "
						  "%4d out of %4d (%0.2f%%), mse=%f\n",
						  clsresults->failures, mReportSet.patterns, perc*100,
						  clsresults->mse);
			  log.printf ("%f %f", clsresults->mse, perc);
			  delete clsresults;
		  } else {
			  out.printf ("Einstein is out of his mind. "
						  "Propably his brain doesn't exist at all\n");
			  log.printf ("666.0 666.0");
		  }
	  } break;
	  
	  case APPROXIMATION: {
		  double mse = brain.test (mReportSet);
		  out.printf ("MSE=%f", mse);
	  } break;
	  
	  default:
		  throw generic_exception (format ("Unknown problem type %d for LearningEAEnv",
										   mProblemType));
	}
	log.flush ();
	
	// Print any network pictures to corresponding log files
	
	// Invoke the brainplan to make the network pictures
	mpBest->execute (TakeBrainPicsMsg ("brainplan", *mpBest));
	
	// Check if any of these exists in the individual's properties
	char pnames[][20]={"brainpic1","brainpic2","brainpic3","braindesc1"};
	char fnames[][20]={"/einstein-pic1.eps","/einstein-pic2.eps",
					   "/einstein-pic3.eps","/einstein-desc1.txt"};
	for (int p=0; p<4; p++) {
		const String& pic = static_cast<const String&> ((*mpBest)[pnames[p]]);
		if (!isnull(pic)) {
			// A property exists -> save it
			FILE* fout = fopen (mLogDir + fnames[p], "w");
			ASSERT (fout);
			fprintf (fout, "%s", (CONSTR) pic);
			fclose (fout);
		}
	}
}

/*******************************************************************************
 *
 ******************************************************************************/
DataOStream& LearningEAEnv::operator>> (DataOStream& out) const {
	out.name("mProblemType") << mProblemType;
	out.name("mEvalPart") << mEvalPart;
	out.name("mTermPart") << mTermPart;
	out.name("mMaxTrainCycles") << mMaxTrainCycles;
	out.name("mReportCycles") << mReportCycles;
	out.name("mValidInterval") << mValidInterval;
	out.name("mTermMethod") << mTermMethod;
	out.name("mPermutate") << mPermutate;
	return out;
}

/*******************************************************************************
 *
 ******************************************************************************/
void LearningEAEnv::check () const {
	EAEnvironment::check ();
	mTrainData.check ();
	mTrainSet.check ();
	mEvaluationSet.check ();
	mReportSet.check ();
	mParams.check ();

	ASSERT (mProblemType>=0 && mProblemType<=2);
	ASSERT (mEvalPart>=0 && mEvalPart<1);
	ASSERT (mTermPart>=0 && mTermPart<1);
	ASSERT (mMaxTrainCycles>=0 && mMaxTrainCycles<100000);
	ASSERT (mReportCycles>=0 && mReportCycles<100000);
	ASSERT (mValidInterval>=0 && mValidInterval<100000);
}
