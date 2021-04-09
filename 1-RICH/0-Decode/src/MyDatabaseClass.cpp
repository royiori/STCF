#include <stdlib.h>

#include <TROOT.h>
#include <TRandom.h>
#include "MyDatabaseClass.h"

//______________________________________________________________________________
//
// MyDatabaseClass class implementation
//______________________________________________________________________________

ClassImp(MyDatabaseClass);

////////////////////////////////////////////////////////////////////////////////
/// Create an Event object.
/// When the constructor is invoked for the first time, the
/// class static variables fgParticles and fgTracks is 0 and
/// the TClonesArray fgParticles is created.

MyDatabaseClass::MyDatabaseClass()
{
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor

MyDatabaseClass::~MyDatabaseClass()
{
   reset();
}
