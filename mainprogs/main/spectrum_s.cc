/*! \file
 *  \brief Main program to run all measurement codes.
 */

#include "chroma.h"

using namespace Chroma;
extern "C" { 
 void _mcleanup();
};

/*
 * Input 
 */

 

bool linkageHack(void)
{
  bool foo = true;

  // Inline Measurements
  foo &= InlineAggregateEnv::registerAll();
//  foo &= GaugeInitEnv::registerAll();

  return foo;
}

//! Main program to run all measurement codes
/*! \defgroup chromamain Main program to run all measurement codes.
 *  \ingroup main
 *
 * Main program to run all measurement codes.
 */

int main(int argc, char *argv[]) 
{
  // Chroma Init stuff
  Chroma::initialize(&argc, &argv);
  
  START_CODE();

  QDPIO::cout << "Linkage = " << linkageHack() << std::endl;

  StopWatch snoop;
  snoop.reset();
  snoop.start();

  XMLReader xml_in;

  // Input parameter structure
  //Inline_input_t  input;
  try
  {
    xml_in.open(Chroma::getXMLInputFileName());
//    read(xml_in, "/chroma", input);
  }
  catch(const std::string& e) 
  {
    QDPIO::cerr << "CHROMA: Caught Exception reading XML: " << e << std::endl;
    QDP_abort(1);
  }
  catch(std::exception& e) 
  {
    QDPIO::cerr << "CHROMA: Caught standard library exception: " << e.what() << std::endl;
    QDP_abort(1);
  }
  catch(...)
  {
    QDPIO::cerr << "CHROMA: caught generic exception reading XML" << std::endl;
    QDP_abort(1);
  }

  QDPIO::cout << "CHROMA: ran successfully" << std::endl;

  END_CODE();

  Chroma::finalize();
  exit(0);
}

