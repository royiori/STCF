#include "G4ProcessManager.hh"
#include "G4Decay.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

// particles

#include "G4ParticleDefinition.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

// EM physics
#include "G4EmDNAChemistry.hh"
#include "G4EmDNAPhysics.hh"
#include "G4EmDNAPhysicsActivator.hh"
#include "G4EmDNAPhysics_option1.hh"
#include "G4EmDNAPhysics_option2.hh"
#include "G4EmDNAPhysics_option3.hh"
#include "G4EmDNAPhysics_option4.hh"
#include "G4EmDNAPhysics_option5.hh"
#include "G4EmDNAPhysics_option7.hh"
#include "G4EmLEPTSPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmLivermorePolarizedPhysics.hh"
#include "G4EmLowEPPhysics.hh"
#include "G4EmModelActivator.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "G4EmStandardPhysicsGS.hh"
#include "G4EmStandardPhysicsSS.hh"
#include "G4EmStandardPhysicsWVI.hh"
#include "G4OpticalPhysics.hh"

// EM Extra
#include "G4BertiniElectroNuclearBuilder.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmMessenger.hh"
#include "G4PAIPhotModel.hh"

// Decay
#include "G4DecayPhysics.hh"
#include "G4MuonicAtomDecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4SpinDecayPhysics.hh"
#include "G4UnknownDecayPhysics.hh"

// Hadron physics
#include "G4IonElasticPhysics.hh"
#include "G4HadronHElasticPhysics.hh"
#include "G4ChargeExchangePhysics.hh"
#include "G4HadronDElasticPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4HadronElasticPhysicsXS.hh"
#include "G4HadronElasticPhysicsPHP.hh"
#include "G4HadronElasticPhysicsLEND.hh"

#include "G4VHadronPhysics.hh"
#include "G4HadronInelasticQBBC.hh"
#include "G4HadronPhysicsQGSP_BERT.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_FTFP_BERT.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4HadronPhysicsQGSP_BIC_AllHP.hh"
#include "G4HadronPhysicsQGS_BIC.hh"
#include "G4HadronPhysicsFTF_BIC.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4HadronPhysicsFTFP_BERT_HP.hh"
#include "G4HadronPhysicsFTFP_BERT_TRV.hh"
#include "G4HadronPhysicsFTFP_BERT_ATL.hh"
#include "G4HadronPhysicsINCLXX.hh"
#include "G4HadronPhysicsNuBeam.hh"
#include "G4HadronPhysicsShielding.hh"

#include "G4StoppingPhysics.hh"

#include "G4IonPhysics.hh"
#include "G4IonPhysicsPHP.hh"
#include "G4IonQMDPhysics.hh"
#include "G4IonINCLXXPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"

#include "G4FastSimulationPhysics.hh"
#include "G4ImportanceBiasing.hh"
#include "G4MinEkineCuts.hh"
#include "G4ParallelWorldPhysics.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4GenericBiasingPhysics.hh"
#include "G4MaxTimeCuts.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4SpecialCuts.hh"
#include "G4WeightWindowBiasing.hh"
#include "G4FastSimulationManagerProcess.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTable.hh"
#include "G4EmConfigurator.hh"
#include "G4RegionStore.hh"

// My physics process
#include "MyPhysListEM.hh"

#include "Verbose.hh"
#include "MyPhysicsList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyPhysicsList::MyPhysicsList() : G4VModularPhysicsList()
{
    if (verbose)
        G4cout << "====>MyPhysicsList::MyPhysicsList()" << G4endl;

    //-- EM Physics
    // options: (electromagnetic)
    // G4EmDNAChemistry.cc              G4EmDNAPhysics_option3.cc        G4EmLivermorePhysics.cc          G4EmStandardPhysics.cc           G4EmStandardPhysics_option2.cc
    // G4EmDNAPhysics.cc                G4EmDNAPhysics_option4.cc        G4EmLivermorePolarizedPhysics.cc G4EmStandardPhysicsGS.cc         G4EmStandardPhysics_option3.cc
    // G4EmDNAPhysicsActivator.cc       G4EmDNAPhysics_option5.cc        G4EmLowEPPhysics.cc              G4EmStandardPhysicsSS.cc         G4EmStandardPhysics_option4.cc
    // G4EmDNAPhysics_option1.cc        G4EmDNAPhysics_option7.cc        G4EmModelActivator.cc            G4EmStandardPhysicsWVI.cc        G4OpticalPhysics.cc
    // G4EmDNAPhysics_option2.cc        G4EmLEPTSPhysics.cc              G4EmPenelopePhysics.cc           G4EmStandardPhysics_option1.cc
    RegisterPhysics(new G4EmStandardPhysics(verbose));
    //RegisterPhysics(new MyPhysListEM("MyEMProc"));

    //-- Synchroton Radiation & Gamma+Nuclear(GN) Physics
    // options: (gamma_lepto_nuclear)
    // G4BertiniElectroNuclearBuilder.cc G4EmExtraPhysics.cc               G4EmMessenger.cc
    RegisterPhysics(new G4EmExtraPhysics(verbose));

    //-- EM Optics
    //fEMOptPhysicsList = new G4OpticalPhysics(verbose);
    //fEMOptPhysicsList->SetScintillation(scint);
    //fEMOptPhysicsList->SetFiniteRiseTime(true);
    G4OpticalPhysics *opticalPhysics = new G4OpticalPhysics();
    opticalPhysics->SetTrackSecondariesFirst(kCerenkov, true);
    RegisterPhysics(opticalPhysics);

    //-- Decays
    // options: (decay)
    // G4DecayPhysics.cc            G4MuonicAtomDecayPhysics.cc  G4RadioactiveDecayPhysics.cc G4SpinDecayPhysics.cc        G4UnknownDecayPhysics.cc
    RegisterPhysics(new G4DecayPhysics(verbose));

    //-- Hadron Physics
    // options: (hadron_elastic)
    // G4ChargeExchangePhysics.cc    G4HadronElasticPhysics.cc     G4HadronElasticPhysicsLEND.cc G4HadronElasticPhysicsXS.cc
    // G4HadronDElasticPhysics.cc    G4HadronElasticPhysicsHP.cc   G4HadronElasticPhysicsPHP.cc  G4HadronHElasticPhysics.cc
    RegisterPhysics(new G4HadronElasticPhysics(verbose));

    // options: (hadron_inelastic)
    // G4HadronInelasticQBBC.cc         G4HadronPhysicsFTFP_BERT_TRV.cc  G4HadronPhysicsQGSP_BERT.cc      G4HadronPhysicsQGSP_BIC_HP.cc    G4VHadronPhysics.cc
    // G4HadronPhysicsFTFP_BERT.cc      G4HadronPhysicsFTF_BIC.cc        G4HadronPhysicsQGSP_BERT_HP.cc   G4HadronPhysicsQGSP_FTFP_BERT.cc
    // G4HadronPhysicsFTFP_BERT_ATL.cc  G4HadronPhysicsINCLXX.cc         G4HadronPhysicsQGSP_BIC.cc       G4HadronPhysicsQGS_BIC.cc
    // G4HadronPhysicsFTFP_BERT_HP.cc   G4HadronPhysicsNuBeam.cc         G4HadronPhysicsQGSP_BIC_AllHP.cc G4HadronPhysicsShielding.cc
    //RegisterPhysics(new G4HadronInelasticQBBC(verbose));

    //-- Ion physics
    // options: (ion_elastic)
    //RegisterPhysics( new G4IonElasticPhysics(verbose));

    // options: (ions_inelastic)
    // G4IonBinaryCascadePhysics.cc G4IonINCLXXPhysics.cc        G4IonPhysics.cc              G4IonPhysicsPHP.cc           G4IonQMDPhysics.cc
    //RegisterPhysics(new G4IonPhysics(verbose));

    //-- Neutron tracking cut
    // options: (limiters)
    // G4FastSimulationPhysics.cc G4ImportanceBiasing.cc     G4MinEkineCuts.cc          G4ParallelWorldPhysics.cc  G4StepLimiterPhysics.cc
    // G4GenericBiasingPhysics.cc G4MaxTimeCuts.cc           G4NeutronTrackingCut.cc    G4SpecialCuts.cc           G4WeightWindowBiasing.cc
    //RegisterPhysics( new G4NeutronTrackingCut(verbose));

    // options: (stopping)
    // G4StoppingPhysics.cc
    RegisterPhysics( new G4StoppingPhysics(verbose) );

    G4StepLimiterPhysics *stepLimiterPhysics = new G4StepLimiterPhysics();
    RegisterPhysics(stepLimiterPhysics);

    //G4LossTableManager::Instance();
    //SetDefaultCutValue(1 * mm);
    //SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyPhysicsList::~MyPhysicsList()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyPhysicsList::ConstructParticle()
{
    if (verbose)
        G4cout << "====>MyPhysicsList::ConstructParticle()" << G4endl;

    G4BosonConstructor pBosonConstructor;
    pBosonConstructor.ConstructParticle();

    G4LeptonConstructor pLeptonConstructor;
    pLeptonConstructor.ConstructParticle();

    G4MesonConstructor pMesonConstructor;
    pMesonConstructor.ConstructParticle();

    G4BaryonConstructor pBaryonConstructor;
    pBaryonConstructor.ConstructParticle();

    G4IonConstructor pIonConstructor;
    pIonConstructor.ConstructParticle();

    G4ShortLivedConstructor pShortLivedConstructor;
    pShortLivedConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
void MyPhysicsList::ConstructProcess()
{
    if (verbose)
        G4cout << "====>MyPhysicsList::ConstructProcess()" << G4endl;

    AddTransportation();
    if (fEMPhysicsList != 0)
        fEMPhysicsList->ConstructProcess();
    if (fEMOptPhysicsList != 0)
        fEMOptPhysicsList->ConstructProcess();
    if (fEMExtraPhysicsList != 0)
        fEMExtraPhysicsList->ConstructProcess();
    if (fDecayPhysicsList != 0)
        fDecayPhysicsList->ConstructProcess();
    if (fHadronElasticPhysicsList != 0)
        fHadronElasticPhysicsList->ConstructProcess();
    if (fHadronInelasticPhysicsList != 0)
        fHadronInelasticPhysicsList->ConstructProcess();
    if (fStoppingPhysicsList != 0)
        fStoppingPhysicsList->ConstructProcess();
    if (fIonPhysicsList != 0)
        fIonPhysicsList->ConstructProcess();
    if (fNeutronPhysicsList != 0)
        fNeutronPhysicsList->ConstructProcess();

    //------------------------------
    // User defined processes
    AddParameterisation();
}
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyPhysicsList::SetCuts()
{
    /*
    G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(20. * eV, 100. * TeV);
    SetCutsWithDefault();

    G4Region *region = G4RegionStore::GetInstance()->GetRegion("Gas_Region");
    G4ProductionCuts *cuts = new G4ProductionCuts();
    cuts->SetProductionCut(1 * um, G4ProductionCuts::GetIndex("gamma"));
    cuts->SetProductionCut(1 * um, G4ProductionCuts::GetIndex("e-"));
    cuts->SetProductionCut(1 * um, G4ProductionCuts::GetIndex("e+"));
    if (region)
    {
        region->SetProductionCuts(cuts);
    }

    DumpCutValuesTable();
    */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyPhysicsList::AddParameterisation()
{
    if (verbose)
        G4cout << "====>MyPhysicsList::AddParameterisation()" << G4endl;

    /*
    //Create an instance of the G4FastSimulationManagerProcess
    G4FastSimulationManagerProcess *fastSimProcess_garfield = new G4FastSimulationManagerProcess("G4FSMP_garfield");
    theParticleTable = G4ParticleTable::GetParticleTable();
    G4ParticleTable::G4PTblDicIterator *theParticleIterator = theParticleTable->GetIterator();
    theParticleIterator->reset();

    while ((*theParticleIterator)())
    {
        G4ParticleDefinition *particle = theParticleIterator->value();
        G4ProcessManager *pmanager = particle->GetProcessManager();
        G4EmConfigurator *config = G4LossTableManager::Instance()->EmConfigurator();
        G4LossTableManager::Instance()->SetVerbose(1);

        //Select the particle for which the parametrisation applies
        G4PAIPhotModel *paiPhot = new G4PAIPhotModel(particle, "G4PAIModel");

        if (particle->GetParticleName() == "e-" || particle->GetParticleName() == "gamma")
            pmanager->AddDiscreteProcess(fastSimProcess_garfield);

        if (particle->GetParticleName() == "e-" || particle->GetParticleName() == "e+")
        {

            config->SetExtraEmModel(particle->GetParticleName(),
                                    "eIoni", paiPhot, "Gas_Region", 0 * MeV,
                                    1e+8 * MeV, paiPhot);
        }
        else if (particle->GetParticleName() == "mu-" || particle->GetParticleName() == "mu+")
        {

            config->SetExtraEmModel(particle->GetParticleName(),
                                    "muIoni", paiPhot, "Gas_Region", 0 * MeV,
                                    1e+8 * MeV, paiPhot);
        }
        else if (particle->GetParticleName() == "proton" || particle->GetParticleName() == "pi+" || particle->GetParticleName() == "pi-")
        {

            config->SetExtraEmModel(particle->GetParticleName(),
                                    "hIoni", paiPhot, "Gas_Region", 0 * MeV,
                                    1e+8 * MeV, paiPhot);
        }
        else if (particle->GetParticleName() == "alpha" || particle->GetParticleName() == "He3" || particle->GetParticleName() == "GenericIon")
        {

            config->SetExtraEmModel(particle->GetParticleName(),
                                    "ionIoni", paiPhot, "Gas_Region", 0 * MeV,
                                    1e+8 * MeV, paiPhot);
        }
    }
    */
}
