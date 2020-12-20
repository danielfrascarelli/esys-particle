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
    template <typename TmplParticleCollection>
    Packer<TmplParticleCollection>::Packer(NTablePtr nTablePtr)
      : m_nTablePtr(nTablePtr),
        m_particlePoolPtr(new ParticlePool),
        m_particleCollectionPtr(
          new ParticleCollection(m_particlePoolPtr)
        ),
        m_idSet()
    {
    }

    template <typename TmplParticleCollection>
    Packer<TmplParticleCollection>::Packer(
      ParticlePoolPtr particlePoolPtr,
      NTablePtr       nTablePtr
    )
      : m_nTablePtr(nTablePtr),
        m_particlePoolPtr(particlePoolPtr),
        m_particleCollectionPtr(
          new ParticleCollection(m_particlePoolPtr)
        ),
        m_idSet()
    {
    }

    template <typename TmplParticleCollection>
    Packer<TmplParticleCollection>::~Packer()
    {
    }

    template <typename TmplParticleCollection>
    void Packer<TmplParticleCollection>::setNTablePtr(NTablePtr nTablePtr)
    {
      m_nTablePtr = nTablePtr;
    }

    template <typename TmplParticleCollection>
    typename Packer<TmplParticleCollection>::NTable &Packer<TmplParticleCollection>::getNTable()
    {
      return *(m_nTablePtr);
    }

    template <typename TmplParticleCollection>
    const typename Packer<TmplParticleCollection>::NTable &
    Packer<TmplParticleCollection>::getNTable() const
    {
      return *m_nTablePtr;
    }

    template <typename TmplParticleCollection>
    typename Packer<TmplParticleCollection>::ParticlePool &
    Packer<TmplParticleCollection>::getParticlePool()
    {
      return *m_particlePoolPtr;
    }

    template <typename TmplParticleCollection>
    typename Packer<TmplParticleCollection>::ParticlePoolPtr
    Packer<TmplParticleCollection>::getParticlePoolPtr()
    {
      return m_particlePoolPtr;
    }

    template <typename TmplParticleCollection>
    const typename Packer<TmplParticleCollection>::ParticlePool &
    Packer<TmplParticleCollection>::getParticlePool() const
    {
      return *m_particlePoolPtr;
    }

    template <typename TmplParticleCollection>
    typename Packer<TmplParticleCollection>::ParticleCollection &
    Packer<TmplParticleCollection>::getParticleCollection()
    {
      return *m_particleCollectionPtr;
    }

    template <typename TmplParticleCollection>
    const typename Packer<TmplParticleCollection>::ParticleCollection &
    Packer<TmplParticleCollection>::getParticleCollection() const
    {
      return *m_particleCollectionPtr;
    }

    template <typename TmplParticleCollection>
    typename Packer<TmplParticleCollection>::Particle &
    Packer<TmplParticleCollection>::constructParticle(const Particle &particle)
    {
      return getParticleCollection().createParticle(particle);
    }

    template <typename TmplParticleCollection>
    typename Packer<TmplParticleCollection>::ParticleIterator
    Packer<TmplParticleCollection>::getParticleIterator()
    {
      return getParticleCollection().getParticleIterator();
    }

    template <typename TmplParticleCollection>
    typename Packer<TmplParticleCollection>::ParticleConstIterator
    Packer<TmplParticleCollection>::getParticleIterator() const
    {
      return getParticleCollection().getParticleIterator();
    }

    template <typename TmplParticleCollection>
    int
    Packer<TmplParticleCollection>::getNumParticles() const
    {
      return getParticleCollection().getNumParticles();
    }

    template <typename TmplParticleCollection>
    int Packer<TmplParticleCollection>::getNextParticleId()
    {
      return static_cast<int>(getNTable().getNumParticles());
    }


    template <typename TmplParticleCollection>
    typename Packer<TmplParticleCollection>::Particle &
    Packer<TmplParticleCollection>::createAndInsertParticle(
      const Particle &particle
    )
    {
      Particle *pParticle = &(constructParticle(particle));
      pParticle->setId(getNextParticleId());
      m_idSet.insert(pParticle->getId());
      getNTable().insert(pParticle);
      return *pParticle;
    }

    template <typename TmplParticleCollection>
    bool Packer<TmplParticleCollection>::contains(const Particle &particle) const
    {
      return (m_idSet.find(particle.getID()) != m_idSet.end());
    }
  }
}
