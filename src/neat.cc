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

#include <magic/mclass.h>
#include <annalee/neat.h>

impl_inner_dynamic (NEATEncoding, NodeGene, {Gentainer});
impl_inner_dynamic (NEATEncoding, ConnectionGene, {Gentainer});
impl_dynamic (NEATEncoding, {ANNEncoding});

//////////////////////////////////////////////////////////////////////////////
//              |   |          |        ----                                //
//              |\  |          |  ___  |      ___    _    ___               //
//              | \ |  __   ---| /   ) | --- /   ) |/ \  /   )              //
//              |  \| /  \ (   | |---  |   \ |---  |   | |---               //
//              |   | \__/  ---|  \__  |___/  \__  |   |  \__               //
//////////////////////////////////////////////////////////////////////////////

NEATEncoding::NodeGene::NodeGene (const GeneticID& name)
		: Gentainer (name)
{
}

NEATEncoding::NodeGene::NodeGene (const NodeGene& orig)
		: Gentainer (orig)
{
}

void NEATEncoding::NodeGene::init ()
{
	Gentainer::init ();
}

/** Add genes needed by a NEAT node.
 *
 * A basic NEAT node does not contain any genes.
 */
void NEATEncoding::NodeGene::addPrivateGenes (
	Gentainer& g,
	const StringMap& params)
{
}

/////////////////////////////////////////////////////////////////////////////////
//  ___                                   o             ----                   //
// /   \        _     _    ___   ___   |           _   |      ___    _    ___  //
// |      __  |/ \  |/ \  /   ) |   \ -+- |  __  |/ \  | --- /   ) |/ \  /   ) //
// |     /  \ |   | |   | |---  |      |  | /  \ |   | |   \ |---  |   | |---  //
// \___/ \__/ |   | |   |  \__   \__/   \ | \__/ |   | |___/  \__  |   |  \__  //
//                                                                             //
/////////////////////////////////////////////////////////////////////////////////

/*******************************************************************************
 * Create a connection gene.
 *
 * A basic NEAT connection specifies two nodes: the source and the
 * target node. These can be changed during the lifetime of a
 * connection.
 *
 * A special feature of connections in NEAT is the "innovation
 * number", which allows finding corresponding genes during
 * crossover. The innovation number is enumerated globally in
 * NEATEncoding and assigned to connections on creation.
 *
 * Stanley & Miikkulainen 2002, p. 107.
 */ 
NEATEncoding::ConnectionGene::ConnectionGene (
	const GeneticID& name,       //< Name of the gene.
	int              innovation, //< Innovation number.
	int              source,     //< Source node.
	int              target)     //< Target node
		: Gentainer (name),
		  mInnovation (innovation),
		  mSource (source),
		  mTarget (target)
{
}

NEATEncoding::ConnectionGene::ConnectionGene (const ConnectionGene& orig)
		: Gentainer (orig),
		  mInnovation (orig.mInnovation),
		  mSource (orig.mSource),
		  mTarget (orig.mTarget)
{
}

/** Implementation for @ref Genstruct. */
void NEATEncoding::ConnectionGene::init ()
{
	Gentainer::init ();
}

/*******************************************************************************
 * Add genes needed by a NEAT connection.
 *
 * A connection has an evolvable weight and an enabled
 * bit that says whether or not the gene is expressed.
 *
 * Stanley & Miikkulainen 2002, p. 107.
 */
void NEATEncoding::ConnectionGene::addPrivateGenes (
	Gentainer& g,
	const StringMap& params)
{
	Gentainer::addPrivateGenes(g, params);

	// The expression gene - is the connection enabled?
	add (new BinaryGene ("E"));

	// Weight
	add (new FloatGene ("W", -1.0, +1.0));
}


//////////////////////////////////////////////////////////////////////////////
//    |   | -----   _   ----- -----                      | o                //
//    |\  | |      / \    |   |       _    ___           |     _            //
//    | \ | |---  /   \   |   |---  |/ \  |   \  __   ---| | |/ \   ___     //
//    |  \| |     |---|   |   |     |   | |     /  \ (   | | |   | (   \    //
//    |   | |____ |   |   |   |____ |   |  \__/ \__/  ---| | |   |  ---/    //
//                                                                  __/     //
//////////////////////////////////////////////////////////////////////////////

/*******************************************************************************
 * Constructor.
 ******************************************************************************/
NEATEncoding::NEATEncoding (const GeneticID& name, const StringMap& params)
		: ANNEncoding (name, params),
		  mInnovationCounter(0)
{
}

/*******************************************************************************
 *
 ******************************************************************************/
NEATEncoding::NEATEncoding (const NEATEncoding& other)
		: ANNEncoding (other),
		  mInnovationCounter(0)
{
}
	
/*******************************************************************************
 *
 ******************************************************************************/
Genstruct* NEATEncoding::replicate () const
{
	return new NEATEncoding (*this);
}

/*******************************************************************************
 *
 ******************************************************************************/
void NEATEncoding::copy (const Genstruct& other)
{
	ANNEncoding::copy (other);
}

/*******************************************************************************
 * Creates the NEAT genome.
 *
 * The genome consists of two lists: node list and connection list.
 * The node list will include all input and output nodes. The params
 * argument must include the following parameters:
 *
 * NEATEncoding.inputs   - Number of input units
 * NEATEncoding.outputs  - Number of output units
 ******************************************************************************/
void NEATEncoding::addPrivateGenes (Gentainer& g, const StringMap& params)
{
	Gentainer::addPrivateGenes (g, params);

	mInputs  = params["NEATEncoding.inputs"].toInt();
	mOutputs = params["NEATEncoding.outputs"].toInt();
	
	// Add input nodes
	for (int i=0; i<mInputs; i++, mInnovationCounter++) {
		NodeGene* gene = new NodeGene (strformat ("U%d", mInnovationCounter));
		gene->hide();
		add (gene);
	}
	
	// Add output nodes
	for (int i=0; i<mOutputs; i++, mInnovationCounter++) {
		NodeGene* gene = new NodeGene (strformat ("U%d", mInnovationCounter));
		gene->hide();
		add (gene);
	}
}

/*******************************************************************************
 * Implementation for Genstruct.
 ******************************************************************************/
bool NEATEncoding::pointMutate (const MutationRate& r)
{
	// Normal point mutations

	// Add connection mutation

	// Find two previously unconnected nodes.

	// Add node mutation

	// Find a random connection to replace with a node and two connections

	return false;
}

/*******************************************************************************
 *
 ******************************************************************************/
void NEATEncoding::init ()
{
}

/*******************************************************************************
 *
 ******************************************************************************/
bool NEATEncoding::execute (const GeneticMsg& msg) const
{
	return false;
}

