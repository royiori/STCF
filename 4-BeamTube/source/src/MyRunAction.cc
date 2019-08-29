//*********************************************
//  This is auto generated by G4gen 0.6
//                                  author:Qian


#include "G4Run.hh"
#include "G4RunManager.hh"
#include "Verbose.hh"
#include "g4root.hh"

#include "MyRunAction.hh"
#include "MyGunAction.hh"
#include "MyDetectorConstruction.hh"

MyRunAction::MyRunAction()
    : G4UserRunAction()
{
  if (verbose)
    G4cout << "====>MyRunAction::MyRunAction()" << G4endl;
}

MyRunAction::~MyRunAction()
{
}

void MyRunAction::BeginOfRunAction(const G4Run *)
{
  if (verbose)
    G4cout << "====>MyRunAction::BeginOfRunAction()" << G4endl;
}

void MyRunAction::EndOfRunAction(const G4Run *)
{
  if (verbose)
    G4cout << "====>MyRunAction::EndOfRunAction()" << G4endl;
}