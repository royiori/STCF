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
//
//

#include "Verbose.hh"
#include "MyDetectorConstruction.hh"
#include "G04SensitiveDetector.hh"
#include "G4SDManager.hh"
#include "G4GDMLParser.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MyDetectorConstruction::MyDetectorConstruction(const G4GDMLParser &parser)
    : G4VUserDetectorConstruction(),
      fParser(parser)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *MyDetectorConstruction::Construct()
{
    /*
    std::map<G4int, G4LogicalVolume*> fdEMap;;
    const G4LogicalVolumeStore* lvs = G4LogicalVolumeStore::GetInstance();
    std::vector<G4LogicalVolume*>::const_iterator lvciter;
    for( lvciter = lvs->begin(); lvciter != lvs->end(); lvciter++ )
    {
        G4cout<<"lv: "<<(*lvciter)->GetName()<<G4endl;

        G4GDMLAuxListType auxInfo = fParser.GetVolumeAuxiliaryInformation(*lvciter);
        //std::vector<G4GDMLAuxPairType>::const_iterator ipair = auxInfo.begin();
        auto ipair = auxInfo.begin();
        for( ipair = auxInfo.begin(); ipair != auxInfo.end(); ipair++ )
        {
            G4String str=ipair->type;
            G4String val=ipair->value;
            G4cout << "Auxiliary Information is found for Logical Volume :  "
            	<< (*lvciter)->GetName() << G4endl;
            G4cout << "Name of Auxiliary type is     :  " << str << G4endl;
            G4cout << "Associated Auxiliary value is :  " << val << G4endl;
            fdEMap.insert(std::make_pair(G4UIcommand::ConvertToInt(val), *lvciter));
        }
    }

    G4PhysicalVolumeStore* pvs = G4PhysicalVolumeStore::GetInstance();
    std::vector<G4VPhysicalVolume*>::iterator cite;

    for(cite=pvs->begin();cite!=pvs->end();cite++)
    {
        G4VPhysicalVolume* pv = *cite;
        G4cout<<"pv: "<<pv->GetName()<<G4endl;
        if(pv->GetName()=="RICH_0")
        {
            G4LogicalVolume * logvol = pv->GetLogicalVolume();
            while(logvol->GetNoDaughters()>0)
                {
                    G4cout<<logvol->GetName()<<" have ";
                    logvol = logvol->GetDaughter(0)->GetLogicalVolume();
                    G4cout<<logvol->GetName()<<G4endl;
                }
        }
    }
    */

    return fParser.GetWorldVolume();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// 根据gdml里的Aux标签找对应的值
//   例如在gdml文件里的 <auxiliary auxtype="nSDType" auxvalue="1"/>  返回1
//   如果没找到，返回空字符串

G4String MyDetectorConstruction::GetGDMLAuxValue(G4String label)
{
    const G4LogicalVolumeStore *lvs = G4LogicalVolumeStore::GetInstance();
    std::vector<G4LogicalVolume *>::const_iterator lvciter;
    for (lvciter = lvs->begin(); lvciter != lvs->end(); lvciter++)
    {
        G4String res = GetGDMLAuxValue(*lvciter, label);
        if (res != "")
            return res;
    }
    return G4String("");
}

G4String MyDetectorConstruction::GetGDMLAuxValue(G4LogicalVolume *lvciter, G4String label)
{
    G4GDMLAuxListType auxInfo = fParser.GetVolumeAuxiliaryInformation(lvciter);
    auto ipair = auxInfo.begin();
    for (ipair = auxInfo.begin(); ipair != auxInfo.end(); ipair++)
    {
        G4String str = ipair->type;
        G4String val = ipair->value;
        if (str == label)
            return val;
    }
    return G4String("");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// 先看ConstructSDandField的注释，里面介绍了gdml的构造
// 根据逻辑体的名字，找此逻辑体下（包括自身）定义为灵敏探测器SD的逻辑体
//   例如在gdml文件里，定义探测器PhysicalVolume为RICH_0，里面包含的读出逻辑体为ReadoutBoxVol，在gdml里对其定义：
//       <volume name="ReadoutBoxVol">
//			<materialref ref="G4_Cu" />
//			<solidref ref="ReadoutBox" />
//			<auxiliary auxtype="SensDet" auxvalue="RICH_readout"/>
//		</volume>
//   那么此函数将返回ReadoutBoxVol的指针。
//   没找到则返回NULL空指针。
//
// 多叉树循环遍历，将momlabel下的子逻辑体放入vector中
//
G4VPhysicalVolume *MyDetectorConstruction::FindSDComponent(G4String momlabel, G4String sdlabel)
{

    G4cout << "--> FindSDComponent: looking for physical volume " << momlabel << " for " << sdlabel << G4endl;

    /*
    G4VPhysicalVolume *momPhysical = NULL;

    //1. 根据momlabel找出momPhysicalVolume
    G4PhysicalVolumeStore *pvs = G4PhysicalVolumeStore::GetInstance();
    std::vector<G4VPhysicalVolume *>::iterator pvcite;
    for (pvcite = pvs->begin(); pvcite != pvs->end(); pvcite++)
    {
        if ((*pvcite)->GetName() == momlabel)
        {
            momPhysical = (*pvcite);
            break;
        }
    }

    if (momPhysical == NULL)
        return NULL;

    //2. 多叉树遍历，将momlabel下的子逻辑体放入vector中。 
    std::stack<G4VPhysicalVolume *> PVstack;
    G4VPhysicalVolume *PVnode;

    PVstack.push(momPhysical);

    while (!PVstack.empty())
    {
        PVnode = PVstack.top();
        PVstack.pop();

        G4String SDType = GetGDMLAuxValue(PVnode->GetLogicalVolume(), sdlabel);
        if (SDType != "")
        {
            G4cout << "  --> found " << PVnode->GetName() << " has aux_value " << SDType << " " << PVnode << " " << momPhysical << G4endl;
            return PVnode;
        }

        for (int i = 0; i < PVnode->GetLogicalVolume()->GetNoDaughters(); i++)
            PVstack.push(PVnode->GetLogicalVolume()->GetDaughter(i));
    }
    */
    return NULL;
}

void MyDetectorConstruction::ConstructSDandField()
{
    //------------------------------------------------
    // 定义SD
    //------------------------------------------------

    G4SDManager *SDman = G4SDManager::GetSDMpointer();
    G04SensitiveDetector *aRICH = new G04SensitiveDetector("RICH");
    SDman->AddNewDetector(aRICH);

    //------------------------------------------------
    // Sensitive detectors GDML 定义
    //------------------------------------------------
    //  在RICH.gdml里，定义如下关键词
    //        <volume name="ReadoutBoxVol">                     定义了读出体
    //      	    <materialref ref="G4_Cu" />                 定义了读出体的材质
    //      		<solidref ref="ReadoutBox" />               定义了读出体的结构尺寸
    //      		<auxiliary auxtype="SensDet" auxvalue="RICH"/>    SensDet表示这是个读出体，其值为RICH，这个名字必须和上面的SDMananger那里对应上
    //	      </volume>
    //------------------------------------------------
    // 注意，这里RICH的logicVolume只有一个，但是放置了好几次，因此会有多个PhysicsVolume
    //
    const G4GDMLAuxMapType *auxmap = fParser.GetAuxMap();
    for (G4GDMLAuxMapType::const_iterator iter = auxmap->begin(); iter != auxmap->end(); iter++)
    {
        for (G4GDMLAuxListType::const_iterator vit = (*iter).second.begin(); vit != (*iter).second.end(); vit++)
        {
            if ((*vit).type == "SensDet")
            {
                G4cout << "Attaching sensitive detector " << (*vit).value
                       << " to volume " << ((*iter).first)->GetName()
                       << G4endl << G4endl;

                G4VSensitiveDetector *mydet = SDman->FindSensitiveDetector((*vit).value);
                if (mydet)
                {
                    G4LogicalVolume *myvol = (*iter).first;
                    myvol->SetSensitiveDetector(mydet);
                }
            }
        }
    }

    //------------------------------------------------
    // Fields
    //------------------------------------------------
    // Create global magnetic field messenger.
    // Uniform magnetic field is then created automatically if
    // the field value is not zero.
    // active by mac file through:
    //       /globalField/setValue 1.0 0 0 tesla

    G4ThreeVector fieldValue;
    fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
    fMagFieldMessenger->SetVerboseLevel(1);

    // Register the field messenger for deleting
    G4AutoDelete::Register(fMagFieldMessenger);
}
