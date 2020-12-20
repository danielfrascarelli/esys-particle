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


#ifndef ESYS_LSMSTRINGUTIL_H
#define ESYS_LSMSTRINGUTIL_H

#include <string>
#include <vector>
#include <sstream>

namespace esys
{
  namespace lsm
  {
    typedef std::vector<std::string> StringVector;
    /**
     * Convenience functions for string manipulation.
     */
    namespace StringUtil
    {
      typedef esys::lsm::StringVector StringVector;
      template <class TmplIterator>
      class StdOStreamOp
      {
      public:
        std::string operator()(const TmplIterator &it) const
        {
          std::stringstream sStream;
          sStream << *it;

          return sStream.str();
        }
      };

      template<class TmplIterator, class TmplStringOperator>
      inline
      std::string join(
        TmplIterator begin,
        TmplIterator end,
        const std::string &delim,
        TmplStringOperator toStringOp = StdOStreamOp<TmplIterator>())
      {
        std::string str;
        if (begin != end)
        {
          TmplIterator it = begin;
          str = toStringOp(begin);
          for (it++; it != end; it++)
          {
            str += delim + toStringOp(it);
          }
        }
        return str;
      }

      template<class TmplContainer, class TmplStringOperator>
      inline
      std::string join(
        const TmplContainer &container,
        const std::string &delim,
        TmplStringOperator toStringOp = StdOStreamOp<typename TmplContainer::const_iterator>())
      {
        typename TmplContainer::const_iterator it  = container.begin();
        typename TmplContainer::const_iterator end = container.end();
        std::string str;
        if (it != end)
        {
          str = toStringOp(it);
          for (it++; it != end; it++)
          {
            str += delim + toStringOp(it);
          }
        }
        return str;
      }

      inline
      std::string joinStringVector(
        const StringVector &container,
        const std::string &delim
      )
      {
        StringVector::const_iterator it  = container.begin();
        StringVector::const_iterator end = container.end();
        std::string str;
        if (it != end)
        {
          str = *it;
          for (it++; it != end; it++)
          {
            str += delim + *it;
          }
        }
        return str;
      }

      template <class TmplData>
      class StdIStreamOp
      {
      public:
        TmplData operator()(const std::string &str) const
        {
          std::stringstream sStream(str);
          TmplData data;
          sStream >> data;

          return data;
        }
      };
      
      template <typename TmplData>
      TmplData to(const std::string &str)
      {
        return StdIStreamOp<TmplData>()(str);
      }
      
      template <typename TmplData>
      std::string toString(const TmplData &data)
      {
        std::stringstream sStream;
        sStream << data;
        return sStream.str();
      }

      template<class TmplData, class TmplStdStreamOp>
      inline      
      std::vector<TmplData> split(
        const std::string &str,
        const std::string &delim,
        TmplStdStreamOp fromStringOp = StdIStreamOp<TmplData>()
      )
      {
        std::vector<TmplData> vec;
        std::string::size_type prevPos = 0;
        std::string::size_type currPos = str.find(delim);
        while (currPos != std::string::npos)
        {
          vec.push_back(
            fromStringOp(str.substr(prevPos, currPos - prevPos))
          );
          prevPos = currPos + delim.length();
          currPos = str.find(delim, prevPos);
        }
        if (prevPos < str.length())
        {
          vec.push_back(
            fromStringOp(str.substr(prevPos, str.length()))
          );
        }
        else if (prevPos == str.length())
        {
          vec.push_back(fromStringOp(""));
        }

        return vec;
      }

      inline
      StringVector splitStrings(const std::string &str, const std::string &delim)
      {
        return split<std::string, StdIStreamOp<std::string> >(str, delim, StdIStreamOp<std::string>());
      }
      
      inline
      std::string trim(const std::string &str)
      {
        size_t b = 0;
        for (; (b < str.size()) && std::isspace(str.at(b)); b++)
        {
          // continue
        }

        size_t e = 0;
        if (str.size() > 0)
        {
          e = str.size()-1;
          for (; (e > b) && std::isspace(str.at(e)); e--)
          {
            // continue
          }
        }        
        return ((e >= b) ? str.substr(b, e-b+1) : "");
      }
    }
  }
}

#endif
