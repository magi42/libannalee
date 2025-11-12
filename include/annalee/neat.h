/***************************************************************************
 *   This file is part of the Annalee library.                             *
 *                                                                         *
 *   Copyright (C) 1998-2008 Marko Grönroos <magi@iki.fi>                  *
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

#ifndef __ANNALEE_NEAT_H__
#define __ANNALEE_NEAT_H__

#include "anngenes.h"

/*******************************************************************************
 * NEAT encoding for neural networks.
 *
 * The NEAT encoding method is described in, for example,
 * - Kenneth O. Stanley and Risto Miikkulainen, Evolving Neural Networks
 *   through Augmenting Topologies, 2002.
 * - Kenneth O. Stanley, Efficient Evolution of Neural Networks
 *   through Complexification, Report AI-TR-04-314 August 2004.
 ******************************************************************************/
class NEATEncoding : public ANNEncoding {
	decl_dynamic (NEATEncoding);
  public:

	/**
	 * Gene for a single NEAT node, that is, a neuron.
	 *
	 * Nodes can be input, output, or hidden nodes.
	 **/
	class NodeGene : public Gentainer {
		decl_dynamic (NodeGene);
	  public:
							NodeGene			(const GeneticID& name=NULL);
							NodeGene			(const NodeGene& orig);
	
		/** Implementation for @ref Genstruct. */
		virtual void		init				();
	
		/** Implementation for @ref Genstruct. */
		virtual void		copy				(const Genstruct& other) {Gentainer::copy (other);}
	
		/** Implementation for @ref Genstruct. */
		virtual void		addPrivateGenes		(Gentainer& g, const StringMap& params);
	};

	/**
	 * Gene for a NEAT node connection.
	 *
	 * Connections are between two neurons and have a weight. A
	 * connection is expressed only if it is enabled. The special
	 * feature of NEAT connections is that they have an "innovation
	 * number", which is used to find corresponding genes during
	 * crossover.
	 **/
	class ConnectionGene : public Gentainer {
		decl_dynamic (ConnectionGene);
	  private:
		const int mInnovation;
		int       mSource;
		int       mTarget;
	  public:
							ConnectionGene		(const GeneticID& name=NULL, int innovation=0, int source=0, int target=0);
							ConnectionGene		(const ConnectionGene& orig);
	
		virtual void		init				();
		virtual void		copy				(const Genstruct& other) {Gentainer::copy (other);}
		virtual void		addPrivateGenes		(Gentainer& g, const StringMap& params);
	};
	
  public:
						NEATEncoding 		() {FORBIDDEN}
						NEATEncoding		(const GeneticID& name,
											 const StringMap& params);
						NEATEncoding		(const NEATEncoding& other);
	
	// Implementations

	/** Implementation for @ref Genstruct. */
	virtual Genstruct*	replicate			() const;

	/** Implementation for @ref Genstruct. */
	virtual void		copy				(const Genstruct& other);

	/** Implementation for @ref Genstruct. */
	virtual bool		execute				(const GeneticMsg& msg) const;

	/** Implementation for @ref Genstruct. */
	virtual void		addPrivateGenes		(Gentainer& g, const StringMap& params);

	/** Implementation for @ref Genstruct. */
	virtual void		init				();

	/** Implementation for @ref Genstruct. */
	virtual bool		pointMutate			(const MutationRate& r);

  protected:
	int mInnovationCounter;
};

#endif
