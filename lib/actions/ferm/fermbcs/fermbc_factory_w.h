// -*- C++ -*-
// $Id: fermbc_factory_w.h,v 2.2 2006-03-16 03:00:11 edwards Exp $
/*! \file
 *  \brief Fermion Boundary Condition factories
 */

#ifndef __fermbc_factory_w_h__
#define __fermbc_factory_w_h__

#include "singleton.h"
#include "objfactory.h"

#include "fermbc.h"


namespace Chroma
{
  //! FermBC factory
  /*! \ingroup fermbcs */
  typedef SingletonHolder< 
    ObjectFactory<FermBC<LatticeFermion>, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  FermBC<LatticeFermion>* (*)(XMLReader&, const std::string&), 
		  StringFactoryError> >
  TheWilsonTypeFermBCFactory;


  //! FermBC array factory
  /*! \ingroup fermbcs */
  typedef SingletonHolder< 
    ObjectFactory<FermBC< multi1d<LatticeFermion> >, 
		  std::string,
		  TYPELIST_2(XMLReader&, const std::string&),
		  FermBC< multi1d<LatticeFermion> >* (*)(XMLReader&, const std::string&), 
		  StringFactoryError> >
  TheWilsonTypeFermBCArrayFactory;

}


#endif
