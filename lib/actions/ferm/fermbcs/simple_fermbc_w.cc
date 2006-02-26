// $Id: simple_fermbc_w.cc,v 2.3 2006-02-26 03:47:52 edwards Exp $
/*! \file
 *  \brief Simple fermionic BC
 */

#include "actions/ferm/fermbcs/simple_fermbc_w.h"
#include "actions/ferm/fermbcs/fermbc_factory_w.h"

namespace Chroma
{

  //! Name and registration
  namespace WilsonTypeSimpleFermBCEnv
  {
    //! Callback function
    FermBC<LatticeFermion>* createFermBC(XMLReader& xml_in, const std::string& path)
    {
      return new SimpleFermBC<LatticeFermion>(SimpleFermBCParams(xml_in, path));
    }

    //! Callback function
    FermBC< multi1d<LatticeFermion> >* createFermBCArray(XMLReader& xml_in, const std::string& path)
    {
      return new SimpleFermBC< multi1d<LatticeFermion> >(SimpleFermBCParams(xml_in, path));
    }

    //! Name to be used
    const std::string name = "SIMPLE_FERMBC";

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= Chroma::TheWilsonTypeFermBCFactory::Instance().registerObject(name, createFermBC);
      foo &= Chroma::TheWilsonTypeFermBCArrayFactory::Instance().registerObject(name, createFermBCArray);
      return foo;
    }

    const bool registered = registerAll();
  }

}
