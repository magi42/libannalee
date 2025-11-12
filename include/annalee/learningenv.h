/***************************************************************************
 *   This file is part of the MagiC++ library.                             *
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

#ifndef __NEURALENV__
#define __NEURALENV__

#include <nhp/gaenvrnmt.h>
#include <inanna/patternset.h>

// Externals
class PatternSet;
template<class K, class V> class Map;
//typedef Map<String,String> StringMap;

/*******************************************************************************
 * An EA environment that evaluates individuals by their learning
 * capabilities. The environment's main attribute is the dataset used
 * for training.
 *
 * This class is currently designed mainly for using a local search in
 * addition to the evolutionary learning of the neural network
 * topology. Whether or not the environment uses local training
 * depends *on the parameters.
 ******************************************************************************/
class LearningEAEnv : public EAEnvironment {
	decl_dynamic (LearningEAEnv);
  public:

						LearningEAEnv	();
						LearningEAEnv	(const PatternSet& trainset,
										 const PatternSet& evaluationset,
										 const PatternSet& reportset,
										 StringMap& params);

	/** Sets problem type to be a classification task or a function
     *    approximation task.
	 *
	 * Default is CLASSIFICATION unless there is only one output in
	 * which case the default is CLASSIFICATION2.
	 *
	 * @see problemtypes
	**/
	void				setProblemType	(int pt) {mProblemType=pt;}

	/** Problem types.
	 *
	 * @param CLASSIFICATION many classes that are represented by one output each
	 * @param CLASSIFICATION2 two classes that are represented by a single output
	 * @param APPROXIMATION many outputs, error always calculated as MSE.
	 *
	 **/	 
	enum problemtypes {CLASSIFICATION=0, CLASSIFICATION2, APPROXIMATION};
	
	// Implementations

	virtual void		addFeaturesTo	(Genome& genome) const;

	/** Implementation for @ref EAEnvironment. */
	virtual void		cycle_report	(OStream& log, OStream& out);

	/** Implementation for @ref EAEnvironment. */
	virtual double		evaluateg		(const Individual& genome);

	/** Implementation for @ref Object. */
	virtual DataOStream& operator>>		(DataOStream& out) const;

	/** Implementation for @ref Object. */
	virtual void		check			() const;

  protected:

	/** Permutates and redivides the mTrainData into mTrainSet and
		mEvaluationSet.
	**/
	void permutate ();

	/** Resplit the dataset into training set and evaluation set */
	void				splitTrainData	();
	Trainer*			createTrainer	() const;
	
  private:
	PatternSet			mTrainData;		// Full training data
	PatternSet			mTrainSet;		// Training part extracted from mTrainData
	PatternSet			mEvaluationSet;	// Evaluation part extracted from mTrainData
	const PatternSource& mReportSet;	// Set for "final" testing. Should not even be here.
	StringMap&			mParams;		// Stored application-level parameters, hmmmmm....
	int					mProblemType;	// classif./2-classif./approx.
	double				mEvalPart;		// Part of training data used for evaluation
	double				mTermPart;		// Part of training data used for termination
	int					mMaxTrainCycles;// Max number of cycles to train
	int					mReportCycles;	// Max number of cycles to train for report set
	int					mValidInterval;	// Training termination check interval
	String				mTermMethod;	// Termination method name (default=UP2)
	bool				mPermutate;		// Permutate training data during evolution
};

#endif
