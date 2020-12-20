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


#ifndef ESYS_LSMITERATIVEREADER_H
#define ESYS_LSMITERATIVEREADER_H

#include <sstream>
#include <memory>

namespace esys
{
  namespace lsm
  {
    /**
     *
     */
    template <class TmplData>
    class IStreamIterator
    {
    public:
      typedef TmplData value_type;
      
      IStreamIterator(std::istream &iStream, int numElements);

      virtual ~IStreamIterator();

      /**
       * Returns true if there are any elements remaining in
       * the stream.
       */
      bool hasNext() const;

      /**
       * Returns the next element in the stream.
       */
      const TmplData &next();

      /**
       * Returns the number of elements remaining in the stream.
       */
      int getNumRemaining() const;

    protected:
      IStreamIterator(const IStreamIterator &it);
      IStreamIterator &operator=(const IStreamIterator &it);

      virtual void readDataFromStream();
      
      std::istream *m_pIStream;
      TmplData     m_data;
      int          m_numRemaining;
    };

    /**
     * Template class which provides an iterator for reading
     * multiple data-items from a stream.
     *
     * @param TmplData data class.
     *        The <code>operator<<(std::istream &iStream, TmplData &data)</code>
     *        operator is used to assign stream data inside the
     *        methods of the IterativeReader<TmplData>::Iterator class.
     */
    template <class TmplIterator>
    class IterativeReader
    {
    public:
      typedef TmplIterator Iterator;
      
      IterativeReader(std::istream &iStream);

      virtual ~IterativeReader();

      /**
       * Creates the iterator using the istream and
       * using the value returned by getNumElements.
       */
      virtual void initialise();

      /**
       * Returns the number of elements to be read from the stream.
       */
      int getNumElements() const;

      /**
       * Returns whether this reader is initialised, that is,
       * whether an iterator has been created.
       */
      bool isInitialised() const;

      Iterator &getIterator();

      protected:
        typedef std::auto_ptr<Iterator> IteratorAutoPtr;

        void setNumElements(int numElements);

        std::istream &getIStream();

        const std::istream &getIStream() const;

        /**
         * Returns a new Iterator object. Caller of this method is
         * to take ownership for the returned object.
         */
        virtual Iterator *createNewIterator();

      private:
        
        int             m_numElements;
        std::istream    *m_pIStream;
        IteratorAutoPtr m_iteratorPtr;
    };
  }
}

#include "Parallel/IterativeReader.hpp"

#endif
