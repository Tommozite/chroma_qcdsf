// $Id: unprec_hamberwu_fermact_w.cc,v 2.3 2006-02-26 03:47:51 edwards Exp $
/*! \file
 *  \brief Unpreconditioned Hamber-Wu fermion action
 */

#include "chromabase.h"
#include "actions/ferm/fermacts/unprec_hamberwu_fermact_w.h"
#include "actions/ferm/linop/unprec_hamberwu_linop_w.h"

#include "actions/ferm/fermacts/fermact_factory_w.h"
#include "actions/ferm/fermbcs/fermbcs_reader_w.h"

#include "io/param_io.h"

namespace Chroma
{

  //! Hooks to register the class with the fermact factory
  namespace UnprecHamberWuFermActEnv
  {
    //! Callback function
    WilsonTypeFermAct<LatticeFermion,multi1d<LatticeColorMatrix> >* createFermAct4D(XMLReader& xml_in,
										    const std::string& path)
    {
      return new UnprecHamberWuFermAct(WilsonTypeFermBCEnv::reader(xml_in, path), 
				       UnprecHamberWuFermActParams(xml_in, path));
    }

    //! Callback function
    /*! Differs in return type */
    FermionAction<LatticeFermion>* createFermAct(XMLReader& xml_in,
						 const std::string& path)
    {
      return createFermAct4D(xml_in, path);
    }

    //! Name to be used
    const std::string name = "UNPRECONDITIONED_HAMBER-WU";

    //! Register all the factories
    bool registerAll()
    {
      bool foo = true;
      foo &= Chroma::TheFermionActionFactory::Instance().registerObject(name, createFermAct);
      foo &= Chroma::TheWilsonTypeFermActFactory::Instance().registerObject(name, createFermAct4D);
      return foo;
    }

    //! Register the fermact
    const bool registered = registerAll();
  }


  //! Read parameters
  UnprecHamberWuFermActParams::UnprecHamberWuFermActParams(XMLReader& xml, const string& path)
  {
    XMLReader paramtop(xml, path);

    // Read the stuff for the action
    if (paramtop.count("Mass") != 0) 
    {
      read(paramtop, "Mass", Mass);
      if (paramtop.count("Kappa") != 0) 
      {
	QDPIO::cerr << "Error: found both a Kappa and a Mass tag" << endl;
	QDP_abort(1);
      }
    }
    else if (paramtop.count("Kappa") != 0)
    {
      Real Kappa;
      read(paramtop, "Kappa", Kappa);
      Mass = kappaToMass(Kappa);    // Convert Kappa to Mass
    }
    else
    {
      QDPIO::cerr << "Error: neither Mass or Kappa found" << endl;
      QDP_abort(1);
    }

    read(paramtop, "u0", u0);

    //  Read optional aniso
    if (paramtop.count("AnisoParam") != 0) 
      read(paramtop, "AnisoParam", anisoParam);
  }


  //! Read parameters
  void read(XMLReader& xml, const string& path, UnprecHamberWuFermActParams& param)
  {
    UnprecHamberWuFermActParams tmp(xml, path);
    param = tmp;
  }



  //! Produce a linear operator for this action
  /*!
   * The operator acts on the entire lattice
   *
   * \param state	    gauge field     	       (Read)
   */
  const UnprecLinearOperator< LatticeFermion, multi1d<LatticeColorMatrix> >*
  UnprecHamberWuFermAct::linOp(Handle<const ConnectState> state) const
  {
    if (param.anisoParam.anisoP)
    {
      QDPIO::cerr << "UnprecHamberWuFermAct::linOp - currently no aniso support" << endl;
      QDP_abort(1);
    }

    return new UnprecHamberWuLinOp(state->getLinks(), param.Mass, param.u0); 
  }

}
