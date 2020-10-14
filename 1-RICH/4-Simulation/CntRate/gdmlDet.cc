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
// --------------------------------------------------------------
//      GEANT 4 - main
//
// --------------------------------------------------------------

#include <vector>

#include "MyPhysicsList.hh"
#include "MyActionInitialization.hh"
#include "MyDetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4GDMLParser.hh"
#include "G4VModularPhysicsList.hh"
#include "G4StepLimiterPhysics.hh"
#include "Randomize.hh"

#include "QBBC.hh"
int verbose = 0;
//#README, 打开输出verbose控制
//https://marketplace.visualstudio.com/items?itemName=Gruntfuggly.todo-tree
//用法： 1. cmd+shift+p, 
//      2. 选 TODO Tree, Add Tag
//      3. or 选settings，修改settings.json，改颜色或图标

int main(int argc, char **argv)
{
   // Detect interactive mode (if only one argument) and define UI session
   //
   G4UIExecutive *ui = 0;
   if (argc == 1)
      ui = new G4UIExecutive(argc, argv);

   // Choose the Random engine
   G4Random::setTheEngine(new CLHEP::RanecuEngine);

   G4RunManager *runManager = new G4RunManager;
   MyDetectorConstruction *detector = new MyDetectorConstruction();
   runManager->SetUserInitialization(detector);
   runManager->SetUserInitialization(new MyPhysicsList(detector));
   //runManager->SetUserInitialization(new QBBC);

   // User action initialization
   runManager->SetUserInitialization(new MyActionInitialization(detector));
   runManager->Initialize();

   // Initialize visualization
   G4VisManager *visManager = new G4VisExecutive;
   visManager->Initialize();

   // Get the pointer to the User Interface manager
   G4UImanager *UImanager = G4UImanager::GetUIpointer();

   // Process macro or start UI session
   if (!ui) // batch mode
   {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command + fileName);
   }
   else // interactive mode
   {
      UImanager->ApplyCommand("/control/execute vis.mac");
      ui->SessionStart();
      delete ui;
   }

   delete visManager;
   delete runManager;
}
