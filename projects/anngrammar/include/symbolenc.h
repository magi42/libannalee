#ifndef __SYMBOLENC_H__
#define __SYMBOLENC_H__

#include "matrixenc.h"

class SymbolMatrixEnc : public MatrixEnc {
  public:
						SymbolMatrixEnc		(const GeneticID& name, const StringMap& paramMap);
						SymbolMatrixEnc		(const SymbolMatrixEnc& orig);

	// Implementations
	Genstruct*			replicate			() const {return new SymbolMatrixEnc (*this);}
	void				copy				(const Genstruct& other);
	bool				execute				(const GeneticMsg& msg) const;
	void				addPrivateGenes		(Gentainer& g, const StringMap& params);
	void				check				() const;

  protected:
	PackTable<uchar>*	decodeMatrix		(const PackTable<uchar>& matrix,
											 const PackTable<uchar>& rules, int l) const;
};

#endif
