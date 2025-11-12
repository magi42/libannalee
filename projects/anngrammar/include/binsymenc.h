#include "symbolenc.h"

class BinSymMatrixEnc : public SymbolMatrixEnc {
  public:
						BinSymMatrixEnc		(const GeneticID& name, const StringMap& paramMap);
						BinSymMatrixEnc		(const BinSymMatrixEnc& orig);

	// Implementations
	Genstruct*			replicate			() const {return new BinSymMatrixEnc (*this);}
	void				copy				(const Genstruct& other);
	bool				execute				(const GeneticMsg& msg) const;
	void				addPrivateGenes		(Gentainer& g, const StringMap& params);
	void				check				() const;

  private:
};
