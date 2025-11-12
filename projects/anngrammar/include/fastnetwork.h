#ifndef __FASTNETWORK_H__
#define __FASTNETWORK_H__

class FastNetwork {
  public:
					FastNetwork		(int inputs, int h1, int h2, int outputs);
					FastNetwork		(const FastNetwork& other);
					~FastNetwork	();

	void			setInputs		(const double* pInputs);
	double*			getInputs		() {return mpActivation;}
	double*			getWeights		() {return mpWeights;}
	const double*	getOutputs		() const {return mpActivation + mUnits -mLayerSizes[3];}
	inline int		size			() const {return mUnits;}
	inline int		weights			() const {return mWeights;}
	inline int		layerSize		(int i) const {return mLayerSizes[i];}
	void			update			();

  private:
	int mLayerSizes[4];
	int mUnits;
	int mWeights;
	
	double* mpActivation;
	double* mpWeights;
};

#endif
