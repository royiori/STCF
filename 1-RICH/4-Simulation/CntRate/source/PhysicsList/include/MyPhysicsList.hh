#ifndef MyPhysicsList_h
#define MyPhysicsList_h 1

#include "G4VModularPhysicsList.hh"

class G4VPhysicsConstructor;
class MyDetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class MyPhysicsList : public G4VModularPhysicsList
{
public:
    MyPhysicsList(MyDetectorConstruction *det);
    ~MyPhysicsList();

    virtual void ConstructParticle();
    void SetCuts();

protected:
    void AddParameterisation();

private:
    MyDetectorConstruction *fDetector;
};

#endif
