
void NolfiNet::indexInputs () {
	for (int i=0; i<cells.size; i++)
		if (cells[i].mFinalType==CT_OUTPUT)
			cells[i].mFinalID = int (cells[i].mCoord.y+0.5);
}

void NolfiNet::indexHiddens (int hiddens) {
	// Index the hidden units according to their X-position. This is
	// one easy way to solve the problem
	int ind=0;
	for (int i=0; i<hiddens; i++) {
		// Find the unindexed hidden cell with the smallest X coordinate
		double minX = 666.0;
		int minNeuron = -1;
		
		for (int j=0; j<cells.size; j++)
			if (cells[j].mFinalType==CT_HIDDEN && cells[j].mCoord.x<=minX &&
				cells[j].mFinalID==EMPTYID) {
				minX = cells[j].mCoord.x;
				minNeuron = j;
			}
		
		ASSERT (minNeuron!=-1);
		
		// Set the index of the smallest found unindexed cell
		cells[minNeuron].mFinalID = mInputs+ind;
		ind++;
	}
}

void NolfiNet::indexOutputs (int hiddens) {
	int totalUnits = mInputs+hiddens+mOutputs;
	for (int i=0; i<cells.size; i++)
		if (cells[i].mFinalType==CT_OUTPUT) {
			cells[i].mFinalID = mInputs+hiddens;//int ((mSize/mOutputs)/cells[i].mCoord.y);
			
			ASSERT (cells[i].mFinalID>=totalUnits-mOutputs);
			ASSERT (cells[i].mFinalID<totalUnits);
		}
}

void NolfiNet::removeDuplicates (int& rOutputs) {
	// Remove input and output cells with duplicate indices
	for (int i=0; i<cells.size; i++)
		for (int j=i+1; j<cells.size; j++)
			if (cells[i].mFinalID!=EMPTYID && cells[i].mFinalID == cells[j].mFinalID) {
				cells[j].mFinalID = EMPTYID;
				// We want to calculate the number of output units exactly
				if (cells[i].mFinalType==CT_OUTPUT)
					rOutputs--;
			}
}

int NolfiNet::connect (FreeNetwork& network) {
	int connections=0;
	for (int i=0; i<cells.size; i++) {
		// (sout << cells[i]).print("\n");

		// Only the cells that have been given a neuron index
		if (cells[i].mFinalID!=EMPTYID) {
			FreeNeuron& neuron = network[cells[i].mFinalID];

			// Set the neuron attributes
			neuron.setBias (cells[i].mBias);
			neuron.moveTo (cells[i].mCoord);
			
			// Generate the branch tip points
			PackArray<Coord2D> tips;
			cells[i].developAxon (tips, mSize);
					
			// Now find other cells that lie near these points
			for (int j=0; j<cells.size; j++)
				if (cells[j].mFinalID > cells[i].mFinalID && cells[j].mFinalType>CT_INPUT)
					for (int k=0; k<tips.size; k++) {
						double d = cells[j].mCoord.sqdist (tips[k]);
						if (d < 25.0) {
							// TRACE4 ("%d->%d: d(%d)=%f", i, j, k, d);
							if (!neuron.connectedTo (cells[j].mFinalID)) {
								// Found a new target, add it to the network
								network.connect (cells[i].mFinalID, cells[j].mFinalID);
								connections++;
								break;	// Break from the k loop
							}
						}
					}

			// Set the initial connection weights to the encoded weight
			for (int j=0; j<neuron.conns(); j++)
				neuron[j].setWeight (cells[i].mWeight);
		}
	}

	return connections;
}
