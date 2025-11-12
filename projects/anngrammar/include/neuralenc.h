#ifndef __NEURALENC_H__
#define __NEURALENC_H__

#include <magic/mmatrix.h>

#include "matrixenc.h"

// external
class FastNetwork;

class NeuralMatrixEnc : public MatrixEnc {
  public:
						NeuralMatrixEnc		(const GeneticID& name, const StringMap& paramMap);
						NeuralMatrixEnc		(const NeuralMatrixEnc& orig);

	// Implementations
	Genstruct*			replicate			() const {return new NeuralMatrixEnc (*this);}
	void				copy				(const Genstruct& other);
	bool				execute				(const GeneticMsg& msg) const;
	void				addPrivateGenes		(Gentainer& g, const StringMap& params);
	void				check				() const;

  private:
	Matrix*				decodeMatrix		(const Matrix& matrix,
											 FastNetwork& network, int l) const;
};

#endif
