#include "MyMagneticField.hh"

#include "G4GenericMessenger.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyMagneticField::MyMagneticField()
    : G4MagneticField(),
      fMessenger(nullptr), fBy(1.0 * tesla)
{
  // define commands for this class
  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyMagneticField::~MyMagneticField()
{
  delete fMessenger;
}

void MyMagneticField::GetFieldValue(const G4double[4], double *bField) const
{
  bField[0] = 0.;
  bField[1] = fBy;
  bField[2] = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyMagneticField::DefineCommands()
{
  // Define /B5/field command directory using generic messenger class
  fMessenger = new G4GenericMessenger(this,
                                      "/MySimulation/field/",
                                      "Field control");

  // fieldValue command
  auto &valueCmd = fMessenger->DeclareMethodWithUnit("value", "tesla",
                                                     &MyMagneticField::SetField,
                                                     "Set field strength.");
  valueCmd.SetParameterName("field", true);
  valueCmd.SetDefaultValue("1.");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
