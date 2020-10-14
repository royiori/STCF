#ifndef MyPhysListEM_h
#define MyPhysListEM_h 1

#include "G4VPhysicsConstructor.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MyPhysListEM : public G4VPhysicsConstructor
{
  public:
    MyPhysListEM(const G4String& name = "MyEMProc");
    ~MyPhysListEM();

    virtual void ConstructParticle();
    virtual void ConstructProcess();

  private:
    G4String fEmName;
};

#endif
