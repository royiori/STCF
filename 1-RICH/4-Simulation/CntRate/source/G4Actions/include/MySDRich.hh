#ifndef MySDRich_h
#define MySDRich_h 1

#include "G4VSensitiveDetector.hh"

class G4Step;

/// Sensitive detector to be attached to the GDML geometry

class MySDRich : public G4VSensitiveDetector
{
  public:
      MySDRich(const G4String&);
     ~MySDRich();

      virtual void Initialize(G4HCofThisEvent*);
      virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
      virtual void EndOfEvent(G4HCofThisEvent*);

  private:

};

#endif

