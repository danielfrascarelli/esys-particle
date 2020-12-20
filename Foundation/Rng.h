/////////////////////////////////////////////////////////////
//                                                         //
// Copyright (c) 2003-2017 by The University of Queensland //
// Centre for Geoscience Computing                         //
// http://earth.uq.edu.au/centre-geoscience-computing      //
//                                                         //
// Primary Business: Brisbane, Queensland, Australia       //
// Licensed under the Open Software License version 3.0    //
// http://www.apache.org/licenses/LICENSE-2.0              //
//                                                         //
/////////////////////////////////////////////////////////////


#ifndef ESYS_LSMRNG_H
#define ESYS_LSMRNG_H

#include <boost/random.hpp>

namespace esys
{
  namespace lsm
  {
    /**
     * Uniform distribution random number generator.
     */
    template <typename TmplRng=boost::mt19937>
    class UniformRng
    {
    public:
      typedef TmplRng Rng;
      typedef boost::uniform_real<> UniformReal;
      typedef boost::variate_generator<Rng &, UniformReal> Generator;

      /**
       * Construct RNG, note that (minRn < maxRn), ie strictly less than.
       * @param minRn minimum value in uniform dist range.
       * @param maxRn maximum value in uniform dist range.
       */
      UniformRng(double minRn, double maxRn);

      /**
       * Generates a random number in this generators
       * range.
       */
      double operator()();

      /**
       * Seeds the RNG with default seed.
       */
      void seed();

      /**
       * Seeds the RNG with specifed seed.
       * @param s Unsigned int.
       */
      template <typename Tmpl>
      void seed(Tmpl &s);

      /**
       * Seeds the RNG with specifed sequence of values.
       * @param begin *begin is the first seed-value in sequence.
       * @param end Signifies termination of sequence.
       */
      template <typename TmplIt>
      void seed(TmplIt begin, TmplIt end);

    private:
      Rng         m_rng;
      UniformReal m_uniform;
      Generator   m_generator;
    };

    typedef UniformRng<> DefaultUniformRng;

    namespace rng
    {
      extern DefaultUniformRng s_zeroOneUniform;
    }
  }
}

#include "Foundation/Rng.hpp"

#endif
