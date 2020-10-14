//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file persistency/gdml/G04/src/G04SensitiveDetector.cc
/// \brief Implementation of the G04SensitiveDetector class
//
//
// $Id: G04SensitiveDetector.cc 68025 2013-03-13 13:43:46Z gcosmo $
//

#include "G04SensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4VProcess.hh"
#include "G4ios.hh"
#include "Verbose.hh"
//#include "g4root.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G04SensitiveDetector::G04SensitiveDetector(const G4String &name)
    : G4VSensitiveDetector(name)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G04SensitiveDetector::~G04SensitiveDetector()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G04SensitiveDetector::Initialize(G4HCofThisEvent *)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G04SensitiveDetector::ProcessHits(G4Step *, G4TouchableHistory *)
{
    return true;
    /*
    //G4double dE = aStep->GetTotalEnergyDeposit();
    //G4double stepl = aStep->GetStepLength();
    //G4bool isFirstDep = aStep->IsFirstStepInVolume();
    //G4double edelta = aStep->GetDeltaEnergy();
    //G4double eNonIoniz = aStep->GetNonIonizingEnergyDeposit();

    G4StepPoint *preStepPoint = aStep->GetPreStepPoint();
    G4StepPoint *postStepPoint = aStep->GetPostStepPoint();
    G4ThreeVector prePos = preStepPoint->GetPosition();
    G4ThreeVector postPos = postStepPoint->GetPosition();

    G4Track *aTrack = aStep->GetTrack();
    G4int trackID = aTrack->GetTrackID();
    G4int parentID = aTrack->GetParentID();
    G4int pdgID = aTrack->GetDefinition()->GetPDGEncoding();
    G4int StepNo = aTrack->GetCurrentStepNumber();
    G4ParticleDefinition *particleType = aTrack->GetDefinition();

    G4LogicalVolume *presentVolume = aStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

    //pro = track->GetCreatorProcess()->GetProcessName();
    //auto analysisManager = G4AnalysisManager::Instance();
    //analysisManager->FillNtupleDColumn(12, 0, fenergy);

    if (parentID == 0) //from incident particle
        return true;

    G4cout << "Processing hits ...." << presentVolume->GetName() << " " << presentVolume->GetName() << " " << prePos.z() << " " << aTrack->GetCreatorProcess()->GetProcessName() << G4endl;

    //1. position
    return true;
    */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G04SensitiveDetector::EndOfEvent(G4HCofThisEvent *)
{
}
