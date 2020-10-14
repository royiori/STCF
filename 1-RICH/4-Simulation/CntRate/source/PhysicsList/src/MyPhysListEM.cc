#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4PhysicsListHelper.hh"

// gamma processes
#include "G4PhotoElectricEffect.hh"
#include "G4PEEffectFluoModel.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4PenelopePhotoElectricModel.hh"
//#
#include "G4PolarizedPhotoElectricEffect.hh"
#include "G4PolarizedPEEffectModel.hh"
#include "G4LivermorePolarizedPhotoElectricModel.hh"

#include "G4ComptonScattering.hh"
#include "G4KleinNishinaCompton.hh"
#include "G4KleinNishinaModel.hh"
#include "G4HeatedKleinNishinaCompton.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4LivermoreComptonModifiedModel.hh"
#include "G4LowEPComptonModel.hh"
#include "G4PenelopeComptonModel.hh"
//#
#include "G4PolarizedCompton.hh"
#include "G4PolarizedComptonModel.hh"
#include "G4LivermorePolarizedComptonModel.hh"
#include "G4LowEPPolarizedComptonModel.hh"

#include "G4GammaConversion.hh"
#include "G4BetheHeitlerModel.hh"
#include "G4PairProductionRelModel.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4LivermoreGammaConversionModelRC.hh"
#include "G4PenelopeGammaConversionModel.hh"
//#
#include "G4PolarizedGammaConversion.hh"
#include "G4PolarizedGammaConversionModel.hh"
#include "G4LivermorePolarizedGammaConversionModel.hh"
//#
#include "G4BoldyshevTripletModel.hh"                //?
#include "G4LivermoreNuclearGammaConversionModel.hh" //?

#include "G4RayleighScattering.hh"
#include "G4LivermoreRayleighModel.hh"
#include "G4LivermorePolarizedRayleighModel.hh"
#include "G4PenelopeRayleighModel.hh"

// electron processes
#include "G4eMultipleScattering.hh"
#include "G4UrbanMscModel.hh"
#include "G4LowEWentzelVIModel.hh"
#include "G4WentzelVIModel.hh"
#include "G4GoudsmitSaundersonMscModel.hh"

#include "G4CoulombScattering.hh"
#include "G4eCoulombScatteringModel.hh"

#include "G4eIonisation.hh"
#include "G4MollerBhabhaModel.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4PenelopeIonisationModel.hh"
#include "G4PAIModel.hh"
#include "G4PAIPhotModel.hh"
//#
#include "G4ePolarizedIonisation.hh"
#include "G4PolarizedMollerBhabhaModel.hh"

#include "G4eBremsstrahlung.hh"
#include "G4SeltzerBergerModel.hh"
#include "G4eBremsstrahlungRelModel.hh"
#include "G4LivermoreBremsstrahlungModel.hh"
#include "G4PenelopeBremsstrahlungModel.hh"
//#
#include "G4ePolarizedBremsstrahlung.hh"
#include "G4ePolarizedBremsstrahlungModel.hh"

// positron processes
#include "G4eplusAnnihilation.hh"

// muon EM processes
#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

// hadron EM process
#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

// ion EM process
#include "G4ionIonisation.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4NuclearStopping.hh"

// distrubution
#include "G4UniversalFluctuation.hh"
//
#include "G4EmProcessOptions.hh"
#include "G4MscStepLimitType.hh"

#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"

#include "G4SystemOfUnits.hh"

#include "Verbose.hh"
#include "MyPhysListEM.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyPhysListEM::MyPhysListEM(const G4String &name)
    : G4VPhysicsConstructor(name)
{
    if (verbose)
        G4cout << "====>MyPhysListEM::MyPhysListEM()" << G4endl;

    SetVerboseLevel(verbose);
    fEmName = G4String(name);

    G4EmParameters *param = G4EmParameters::Instance();
    param->SetDefaults();
    param->SetVerbose(verbose);
    param->SetMinEnergy(100 * eV);
    param->SetMaxEnergy(1 * TeV);
    param->SetLowestElectronEnergy(100 * eV);
    param->SetNumberOfBinsPerDecade(20);
    param->ActivateAngularGeneratorForIonisation(true);
    param->SetMscRangeFactor(0.02);
    param->SetMuHadLateralDisplacement(true);
    param->SetMscStepLimitType(fUseDistanceToBoundary);
    param->SetFluo(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyPhysListEM::~MyPhysListEM()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyPhysListEM::ConstructParticle()
{

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MyPhysListEM::ConstructProcess()
{
    if (verbose)
        G4cout << "====>MyPhysListEM::ConstructProcess()" << G4endl;

    G4PhysicsListHelper *ph = G4PhysicsListHelper::GetPhysicsListHelper();

    // Add standard EM Processes
    //
    auto aParticleIterator = GetParticleIterator();

    aParticleIterator->reset();
    while ((*aParticleIterator)())
    {
        G4ParticleDefinition *particle = aParticleIterator->value();
        G4String particleName = particle->GetParticleName();

        //------------
        // 1. gamma
        if (particleName == "gamma")
        {

            //------------
            //
            //--> 1.1 photoelectric effect || polarizedPhotoElectric
            //G4PhotoElectricEffect *thePhotoElectricEffect = new G4PhotoElectricEffect();
            G4PolarizedPhotoElectricEffect* thePhotoElectricEffect = new G4PolarizedPhotoElectricEffect();

            //----> 1.1.0 default model: G4PEEffectFluoModel || G4PolarizedPEEffectModel

            //----> 1.1.1 G4LivermorePhotoElectricModel
            //
            //thePhotoElectricEffect->SetEmModel(new G4LivermorePhotoElectricModel(), 1);

            //----> 1.1.2 G4LivermorePolarizedPhotoElectricModel
            //
            G4double LivermoreHighEnergyLimit = 1.0*GeV;
            G4LivermorePolarizedPhotoElectricModel* theLivermorePhotoElectricModel = new G4LivermorePolarizedPhotoElectricModel();
            theLivermorePhotoElectricModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
            thePhotoElectricEffect->SetEmModel(theLivermorePhotoElectricModel);

            //----> 1.1.3 G4PenelopePhotoElectricModel
            //
            //G4double PenelopeHighEnergyLimit = 1.0*GeV;
            //G4PenelopePhotoElectricModel* thePEPenelopeModel = new G4PenelopePhotoElectricModel();
            //thePEPenelopeModel->SetHighEnergyLimit(PenelopeHighEnergyLimit);
            //thePhotoElectricEffect->SetEmModel(thePEPenelopeModel, 1);

            //----> other options:
            //
            // G4PEEffectFluoModel.cc			G4PolarizedPEEffectModel.cc
            // G4LivermorePhotoElectricModel.cc	G4LivermorePolarizedPhotoElectricGDModel.cc
            // G4PenelopePhotoElectricModel.cc
            // G4PhotoElectricAngularGeneratorPolarized.cc	G4PhotoElectricAngularGeneratorSauterGavrila.cc		G4PhotoElectricAngularGeneratorSimple.cc

            ph->RegisterProcess(thePhotoElectricEffect, particle);

            //------------
            //
            //--> 1.2 Compton scattering || polarized Compton scattering
            G4ComptonScattering *theComptonScattering = new G4ComptonScattering();
            //G4PolarizedCompton* theComptonScattering = new G4PolarizedCompton();

            //----> 1.2.0 default model: G4KleinNishinaCompton || G4PolarizedComptonModel

            //----> 1.2.1 G4LivermoreComptonModel
            //
            theComptonScattering->SetEmModel(new G4LivermoreComptonModel(), 1);

            //----> 1.2.2 G4LivermorePolarizedComptonModel
            //
            //G4LivermorePolarizedComptonModel* theLivermoreComptonModel = new G4LivermorePolarizedComptonModel();
            //theLivermoreComptonModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
            //theComptonScattering->AddEmModel(0, theLivermoreComptonModel);

            //----> 1.2.3 G4LowEPComptonModel
            //
            //G4LowEPComptonModel* theLowEPComptonModel = new G4LowEPComptonModel();
            //theLowEPComptonModel->SetHighEnergyLimit(20*MeV);
            //theComptonScattering->AddEmModel(0, theLowEPComptonModel);

            //----> 1.2.4 G4PenelopeComptonModel
            //
            //G4PenelopeComptonModel* theComptonPenelopeModel = new G4PenelopeComptonModel();
            //theComptonPenelopeModel->SetHighEnergyLimit(PenelopeHighEnergyLimit);
            //theComptonScattering->SetEmModel(theComptonPenelopeModel, 1);

            //----> other options:
            //
            // G4eInverseCompton.cc
            // G4KleinNishinaCompton.cc	G4HeatedKleinNishinaCompton.cc
            // G4LowEPComptonModel.cc		G4LowEPPolarizedComptonModel.cc
            // G4LivermoreComptonModel.cc	G4LivermoreComptonModifiedModel.cc	G4LivermorePolarizedComptonModel.cc
            // G4PenelopeComptonModel.cc

            ph->RegisterProcess(theComptonScattering, particle);

            //------------
            //
            //--> 1.3 gamma conversion || polarized gamma conversion
            G4GammaConversion *theGammaConversion = new G4GammaConversion();
            //G4PolarizedGammaConversionModel *theGammaConversion = new G4PolarizedGammaConversionModel();

            //----> 1.3.0 default model:   G4BetheHeitlerModel & G4PairProductionRelModel
            //                          || G4PolarizedGammaConversionModel

            //----> 1.3.1 G4LivermoreGammaConversionModel - Livermore model below 80 GeV
            //
            theGammaConversion->SetEmModel(new G4LivermoreGammaConversionModel(), 1);

            //----> 1.3.2 G4LivermorePolarizedGammaConversionModel
            //
            //G4LivermorePolarizedGammaConversionModel* theLivermoreGammaConversionModel = new G4LivermorePolarizedGammaConversionModel();
            //theLivermoreGammaConversionModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
            //theGammaConversion->AddEmModel(0, theLivermoreGammaConversionModel);

            //----> 1.3.4 G4PenelopeGammaConversionModel
            //
            //G4PenelopeGammaConversionModel* theGCPenelopeModel = new G4PenelopeGammaConversionModel();
            //theGammaConversion->SetEmModel(theGCPenelopeModel,1);

            //----> other options:
            //

            ph->RegisterProcess(theGammaConversion, particle);

            //------------
            //
            //--> 1.4 rayleigh - default Rayleigh scattering is Livermore
            G4RayleighScattering *theRayleigh = new G4RayleighScattering();

            //----> 1.4.1 default Rayleigh scattering is G4LivermoreRayleighModel

            //----> 1.4.2 G4LivermorePolarizedRayleighModel
            //
            //G4LivermorePolarizedRayleighModel* theRayleighModel = new G4LivermorePolarizedRayleighModel();
            //theRayleighModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
            //theRayleigh->AddEmModel(0, theRayleighModel);

            //----> 1.4.3 G4PenelopeRayleighModel
            //
            //G4PenelopeRayleighModel* theRayleighPenelopeModel = new G4PenelopeRayleighModel();
            ////theRayleighPenelopeModel->SetHighEnergyLimit(PenelopeHighEnergyLimit);
            //theRayleigh->SetEmModel(theRayleighPenelopeModel, 1);

            //----> other options:
            //

            ph->RegisterProcess(theRayleigh, particle);
        }
        //------------
        // 2. electron
        else if (particleName == "e-")
        {

            G4double highEnergyLimit = 100 * MeV;

            //------------
            //
            //--> 2.1 multiple scattering
            G4eMultipleScattering *msc = new G4eMultipleScattering;
            //msc->SetStepLimitType(fUseDistanceToBoundary);

            //----> 2.1.0 default model: G4UrbanMscModel

            //----> 2.1.1 use 1 model
            //
            // options:
            // G4UrbanMscModel, G4WentzelVIModel, G4LowEWentzelVIModel, G4GoudsmitSaundersonMscModel
            //G4LowEWentzelVIModel *msc1 = new G4LowEWentzelVIModel();
            //msc->SetEmModel(msc1, 1);

            //----> 2.1.2 use 2 models for different energy range
            //
            G4UrbanMscModel *msc1 = new G4UrbanMscModel();
            G4WentzelVIModel *msc2 = new G4WentzelVIModel();
            msc1->SetHighEnergyLimit(highEnergyLimit);
            msc2->SetLowEnergyLimit(highEnergyLimit);
            msc->AddEmModel(0, msc1);
            msc->AddEmModel(0, msc2);

            ph->RegisterProcess(msc, particle);

            //------------
            //
            //--> 2.2 coulomb scattering, default model: G4eCoulombScatteringModel
            G4CoulombScattering *ss = new G4CoulombScattering();

            G4eCoulombScatteringModel *ssm = new G4eCoulombScatteringModel();
            ssm->SetLowEnergyLimit(highEnergyLimit);
            ssm->SetActivationLowEnergyLimit(highEnergyLimit);
            ss->SetEmModel(ssm, 1);
            ss->SetMinKinEnergy(highEnergyLimit);

            ph->RegisterProcess(ss, particle);

            //------------
            //
            //--> 2.3 Ionisation and polarized ionisation
            G4eIonisation *eIoni = new G4eIonisation();
            //G4ePolarizedIonisation *eIoni = new G4ePolarizedIonisation();
            eIoni->SetStepFunction(0.2, 100 * um);

            //----> 2.3.1 default is G4MollerBhabhaModel || G4PolarizedMollerBhabhaModel

            //----> 2.3.2 Livermore should be used only for low energies
            //
            G4LivermoreIonisationModel *theIoniLivermore = new G4LivermoreIonisationModel();
            theIoniLivermore->SetHighEnergyLimit(0.1 * MeV);
            eIoni->AddEmModel(0, theIoniLivermore, new G4UniversalFluctuation());

            //----> 2.3.3 G4PenelopeIonisationModel
            //
            //G4PenelopeIonisationModel* theIoniPenelope = new G4PenelopeIonisationModel();
            //theIoniPenelope->SetHighEnergyLimit(PenelopeHighEnergyLimit);
            //eIoni->AddEmModel(0,theIoniPenelope,new G4UniversalFluctuation());

            //----> other options (PAI: PhotoAbsorption Ionization)
            // G4PAIModel.cc 	G4PAIPhotModel.cc
			//G4PAIPhotModel *thePAIphotModel = new G4PAIPhotModel();
            //thePAIphotModel->SetHighEnergyLimit(0.1 * MeV);
            //eIoni->AddEmModel(0, thePAIphotModel, new G4UniversalFluctuation());


            ph->RegisterProcess(eIoni, particle);

            //------------
            //
            //--> 2.4 Bremsstrahlung or polarized bremsstrahlung
            G4eBremsstrahlung *eBrem = new G4eBremsstrahlung();
            //G4ePolarizedBremsstrahlung* eBrem = new G4ePolarizedBremsstrahlung();

            //----> 2.4.1 default is G4SeltzerBergerModel (min ~ GeV) + G4eBremsstrahlungRelModel (GeV ~ max)
            //               || G4ePolarizedBremsstrahlungModel

            //----> 2.4.2 G4SeltzerBergerModel + Angular distribution
            //
            //G4VEmModel* theBrem = new G4SeltzerBergerModel();
            //theBrem->SetHighEnergyLimit(1*GeV);
            //theBrem->SetAngularDistribution(new G4Generator2BS());

            //----> 2.4.3 Livermore
            //
            //G4VEmModel* theBremLivermore = new G4LivermoreBremsstrahlungModel();
            //theBremLivermore->SetHighEnergyLimit(1*GeV);
            //theBremLivermore->SetAngularDistribution(new G4Generator2BS());
            //eBrem->SetEmModel(theBremLivermore,1);

            //----> 2.4.4 G4PenelopeBremsstrahlungModel
            //
            //G4PenelopeBremsstrahlungModel* theBremPenelope = new G4PenelopeBremsstrahlungModel();
            //theBremPenelope->SetHighEnergyLimit(PenelopeHighEnergyLimit);
            //eBrem->AddEmModel(0,theBremPenelope);

            ph->RegisterProcess(eBrem, particle);
        }
        //------------
        // 3. position
        else if (particleName == "e+")
        {

            G4double highEnergyLimit = 100 * MeV;

            //------------
            //
            //--> 3.1 multiple scattering
            G4eMultipleScattering *msc = new G4eMultipleScattering;
            //msc->SetStepLimitType(fUseDistanceToBoundary);

            //----> 3.1.0 default model: G4UrbanMscModel

            //----> 3.1.1 use 1 model
            //
            // options:
            // G4UrbanMscModel, G4WentzelVIModel, G4LowEWentzelVIModel, G4GoudsmitSaundersonMscModel
            //G4LowEWentzelVIModel *msc1 = new G4LowEWentzelVIModel();
            //msc->SetEmModel(msc1, 1);

            //----> 3.1.2 use 2 models for different energy range
            //
            G4UrbanMscModel *msc1 = new G4UrbanMscModel();
            G4WentzelVIModel *msc2 = new G4WentzelVIModel();
            msc1->SetHighEnergyLimit(highEnergyLimit);
            msc2->SetLowEnergyLimit(highEnergyLimit);
            msc->AddEmModel(0, msc1);
            msc->AddEmModel(0, msc2);

            ph->RegisterProcess(msc, particle);

            //------------
            //
            //--> 3.2 coulomb scattering, default model: G4eCoulombScatteringModel
            G4CoulombScattering *ss = new G4CoulombScattering();

            G4eCoulombScatteringModel *ssm = new G4eCoulombScatteringModel();
            ssm->SetLowEnergyLimit(highEnergyLimit);
            ssm->SetActivationLowEnergyLimit(highEnergyLimit);
            ss->SetEmModel(ssm, 1);
            ss->SetMinKinEnergy(highEnergyLimit);

            ph->RegisterProcess(ss, particle);

            //------------
            //
            //--> 3.3 Ionisation and polarized ionisation
            G4eIonisation *eIoni = new G4eIonisation();
            //G4ePolarizedIonisation *eIoni = new G4ePolarizedIonisation();
            eIoni->SetStepFunction(0.2, 100 * um);

            //----> 3.3.1 default is G4MollerBhabhaModel || G4PolarizedMollerBhabhaModel

            //----> 3.3.2 Livermore should be used only for low energies electron, not for positron

            //----> 3.3.3 G4PenelopeIonisationModel
            //
            //G4PenelopeIonisationModel* theIoniPenelope = new G4PenelopeIonisationModel();
            //theIoniPenelope->SetHighEnergyLimit(PenelopeHighEnergyLimit);
            //eIoni->AddEmModel(0,theIoniPenelope,new G4UniversalFluctuation());

            //----> other options (PAI: PhotoAbsorption Ionization)
            // G4PAIModel.cc  G4PAIPhotModel.cc

            ph->RegisterProcess(eIoni, particle);

            //------------
            //
            //--> 3.4 Bremsstrahlung or polarized bremsstrahlung
            G4eBremsstrahlung *eBrem = new G4eBremsstrahlung();
            //G4ePolarizedBremsstrahlung* eBrem = new G4ePolarizedBremsstrahlung();

            //----> 3.4.1 default is G4SeltzerBergerModel (min ~ GeV) + G4eBremsstrahlungRelModel (GeV ~ max)
            //               || G4ePolarizedBremsstrahlungModel

            //----> 3.4.2 G4SeltzerBergerModel + Angular distribution
            //
            //G4VEmModel* theBrem = new G4SeltzerBergerModel();
            //theBrem->SetHighEnergyLimit(1*GeV);
            //theBrem->SetAngularDistribution(new G4Generator2BS());

            //----> 3.4.3 Livermore
            //
            //G4VEmModel* theBremLivermore = new G4LivermoreBremsstrahlungModel();
            //theBremLivermore->SetHighEnergyLimit(1*GeV);
            //theBremLivermore->SetAngularDistribution(new G4Generator2BS());
            //eBrem->SetEmModel(theBremLivermore,1);

            //----> 3.4.4 G4PenelopeBremsstrahlungModel
            //
            //G4PenelopeBremsstrahlungModel* theBremPenelope = new G4PenelopeBremsstrahlungModel();
            //theBremPenelope->SetHighEnergyLimit(PenelopeHighEnergyLimit);
            //eBrem->AddEmModel(0,theBremPenelope);

            ph->RegisterProcess(eBrem, particle);

            //------------
            //
            //--> 3.5 Annihilation
            G4eplusAnnihilation *eAnn = new G4eplusAnnihilation();

            //----> 3.5.1 default is G4eeToTwoGammaModel

            //----> 3.5.2 G4PenelopeAnnihilationModel
            //G4PenelopeAnnihilationModel* theAnnPenelope = new G4PenelopeAnnihilationModel();
            //theAnnPenelope->SetHighEnergyLimit(PenelopeHighEnergyLimit);
            //eAnni->AddEmModel(0,theAnnPenelope);

            ph->RegisterProcess(eAnn, particle);
        }
        else if (particleName == "mu+" ||
                 particleName == "mu-")
        {

            G4MuIonisation *muIoni = new G4MuIonisation();
            muIoni->SetStepFunction(0.2, 50 * um);
            ph->RegisterProcess(muIoni, particle);

            G4MuMultipleScattering *mumsc = new G4MuMultipleScattering();
            mumsc->AddEmModel(0, new G4WentzelVIModel());
            ph->RegisterProcess(mumsc, particle);

            G4MuBremsstrahlung *mub = new G4MuBremsstrahlung();
            ph->RegisterProcess(mub, particle);

            G4MuPairProduction *mup = new G4MuPairProduction();
            ph->RegisterProcess(mup, particle);

            G4CoulombScattering *mucoul = new G4CoulombScattering();
            ph->RegisterProcess(mucoul, particle);
        }
        else if (particleName == "alpha" ||
                 particleName == "He3")
        {

            G4hMultipleScattering *msc = new G4hMultipleScattering();
            ph->RegisterProcess(msc, particle);

            G4ionIonisation *ionIoni = new G4ionIonisation();
            ionIoni->SetStepFunction(0.1, 10 * um);
            ph->RegisterProcess(ionIoni, particle);

            G4NuclearStopping *pnuc = new G4NuclearStopping();
            ph->RegisterProcess(pnuc, particle);
            pnuc->SetMaxKinEnergy(MeV);
        }
        else if (particleName == "pi-" ||
                 particleName == "pi+")
        {

            G4hMultipleScattering *pimsc = new G4hMultipleScattering();
            ph->RegisterProcess(pimsc, particle);

            G4hIonisation *hIoni = new G4hIonisation();
            hIoni->SetStepFunction(0.2, 50 * um);
            ph->RegisterProcess(hIoni, particle);

            G4hBremsstrahlung *pib = new G4hBremsstrahlung();
            ph->RegisterProcess(pib, particle);

            G4hPairProduction *pip = new G4hPairProduction();
            ph->RegisterProcess(pip, particle);
        }
        else if (particleName == "kaon+" ||
                 particleName == "kaon-")
        {

            G4hMultipleScattering *kmsc = new G4hMultipleScattering();
            ph->RegisterProcess(kmsc, particle);

            G4hIonisation *hIoni = new G4hIonisation();
            hIoni->SetStepFunction(0.2, 50 * um);
            ph->RegisterProcess(hIoni, particle);

            G4hBremsstrahlung *kb = new G4hBremsstrahlung();
            ph->RegisterProcess(kb, particle);

            G4hPairProduction *kp = new G4hPairProduction();
            ph->RegisterProcess(kp, particle);
        }
        else if (particleName == "proton" ||
                 particleName == "anti_proton")
        {

            G4hMultipleScattering *pmsc = new G4hMultipleScattering();
            ph->RegisterProcess(pmsc, particle);

            G4hIonisation *hIoni = new G4hIonisation();
            hIoni->SetStepFunction(0.2, 50 * um);
            ph->RegisterProcess(hIoni, particle);

            G4hBremsstrahlung *pb = new G4hBremsstrahlung();
            ph->RegisterProcess(pb, particle);

            G4hPairProduction *pp = new G4hPairProduction();
            ph->RegisterProcess(pp, particle);

            G4NuclearStopping *pnuc = new G4NuclearStopping();
            ph->RegisterProcess(pnuc, particle);
            pnuc->SetMaxKinEnergy(MeV);
        }
        else if (particleName == "GenericIon")
        {

            G4ionIonisation *ionIoni = new G4ionIonisation();
            ionIoni->SetEmModel(new G4IonParametrisedLossModel());
            ionIoni->SetStepFunction(0.1, 1 * um);
            ph->RegisterProcess(ionIoni, particle);

            G4hMultipleScattering *hmsc = new G4hMultipleScattering("ionmsc");
            ph->RegisterProcess(hmsc, particle);

            G4NuclearStopping *pnuc = new G4NuclearStopping();
            ph->RegisterProcess(pnuc, particle);
            pnuc->SetMaxKinEnergy(MeV);
        }
        else if (particleName == "B+" ||
                 particleName == "B-" ||
                 particleName == "D+" ||
                 particleName == "D-" ||
                 particleName == "Ds+" ||
                 particleName == "Ds-" ||
                 particleName == "anti_He3" ||
                 particleName == "anti_alpha" ||
                 particleName == "anti_deuteron" ||
                 particleName == "anti_lambda_c+" ||
                 particleName == "anti_omega-" ||
                 particleName == "anti_sigma_c+" ||
                 particleName == "anti_sigma_c++" ||
                 particleName == "anti_sigma+" ||
                 particleName == "anti_sigma-" ||
                 particleName == "anti_triton" ||
                 particleName == "anti_xi_c+" ||
                 particleName == "anti_xi-" ||
                 particleName == "deuteron" ||
                 particleName == "lambda_c+" ||
                 particleName == "omega-" ||
                 particleName == "sigma_c+" ||
                 particleName == "sigma_c++" ||
                 particleName == "sigma+" ||
                 particleName == "sigma-" ||
                 particleName == "tau+" ||
                 particleName == "tau-" ||
                 particleName == "triton" ||
                 particleName == "xi_c+" ||
                 particleName == "xi-")
        {

            G4hMultipleScattering *hmsc = new G4hMultipleScattering("ionmsc");
            ph->RegisterProcess(hmsc, particle);

            G4hIonisation *hIoni = new G4hIonisation();
            ph->RegisterProcess(hIoni, particle);

            G4NuclearStopping *pnuc = new G4NuclearStopping();
            ph->RegisterProcess(pnuc, particle);
            pnuc->SetMaxKinEnergy(MeV);
        }
        else if ((!particle->IsShortLived()) &&
                 (particle->GetPDGCharge() != 0.0) &&
                 (particle->GetParticleName() != "chargedgeantino"))
        {

            //all others charged particles except geantino
            ph->RegisterProcess(new G4hMultipleScattering(), particle);
            ph->RegisterProcess(new G4hIonisation(), particle);
        }
    }

    
    // Deexcitation
    //
    G4VAtomDeexcitation *de = new G4UAtomicDeexcitation();
    // default: all false
    de->SetFluo(true);   // active Deexcitation, default is false
    de->SetAuger(false); // Auger electron
    de->SetPIXE(false);  // Particle induced X-ray emission
    G4LossTableManager::Instance()->SetAtomDeexcitation(de);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
