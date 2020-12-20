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


namespace esys
{
  namespace lsm
  {
    template <typename TmplParticle>
    ParticleComparer<TmplParticle>::ParticleComparer(const Particle &p)
      : m_pParticle(&p)
    {
    }

    template <typename TmplParticle>
    bool ParticleComparer<TmplParticle>::operator()(
      const Particle &p1,
      const Particle &p2
    ) const
    {
      return
        (
          ((p1.getPos() - m_pParticle->getPos()).norm() - p1.getRad())
          <
          ((p2.getPos() - m_pParticle->getPos()).norm() - p2.getRad())
        );
    }

    template <typename TmplParticle>
    bool ParticleComparer<TmplParticle>::operator()(
      const Particle *p1,
      const Particle *p2
    ) const
    {
      return this->operator()(*p1,*p2);
    }
  }
}
