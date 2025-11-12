#ifndef __MATRIXENV__
#define __MATRIXENV__

#include <magic/mmatrix.h>
#include <magic/mmap.h>
#include <magic/mclass.h>
#include <nhp/gaenvrnmt.h>

class MatrixEnv : public EAEnvironment {
	decl_dynamic (MatrixEnv);
  public:
						MatrixEnv		() {FORBIDDEN;}
						MatrixEnv		(const StringMap& params);

	// Implementations
	
	virtual void		addFeaturesTo	(Genome& genome) const;
	virtual void		cycle_report	(OStream& log, OStream& out);
	virtual double		evaluateg		(const Individual& genome);
	TextOStream&		operator>>		(TextOStream& out) const;
	void				check			() const;

  private:
	Matrix				mMatrix;
	int					mPowerOf2;
	const StringMap*	mrpParamMap;
};

#endif
