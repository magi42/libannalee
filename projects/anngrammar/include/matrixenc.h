#ifndef __MATRIXENC_H__
#define __MATRIXENC_H__

#include <nhp/genetics.h>
#include <magic/mtable.h>

class MatrixEnc : public Gentainer {
  public:
						MatrixEnc			(const GeneticID& name, const StringMap& paramMap);
						MatrixEnc			(const MatrixEnc& orig);
	
	// Implementations
	void				addPrivateGenes		(Gentainer& g, const StringMap& params){MUST_OVERLOAD}
	void				copy				(const Genstruct& other);

  protected:
	int	mIterations;	// Rewriting iterations needed
};

#endif
