#ifndef MyDetectorSettings_H
#define MyDetectorSettings_H 1

#include <map>

#include "G4GDMLParser.hh"
class G4LogicalVolume;

#include "RadiatorDescription.hh"

class MyDetectorSettings
{

public:
    MyDetectorSettings();
    ~MyDetectorSettings();

    void ApplyAuxValue(G4GDMLParser *fParser);
    void ApplyAuxValue(const G4GDMLAuxListType *auxInfoList, G4LogicalVolume *vol = NULL);

    //#AuxXML 1. 定义Setting里的用户函数
    void setColor(G4LogicalVolume *vol, G4String value, G4String unit); //<auxiliary auxtype="setColor" auxvalue="kRed" auxunit="0"/> 参考ROOT的色盘
    void setAlpha(G4LogicalVolume *vol, G4String value); //<auxiliary auxtype="setAlpha" auxvalue="0.5"/>

    //for TRD settings
    RadiatorDescription *fRadiatorDescription = NULL;
    void setTRD_FoilThickness(G4LogicalVolume *vol, G4String value, G4String unit); //<auxiliary auxtype="setTRD_FoilThickness" auxvalue="0.006" auxunit="mm"/>
    void setTRD_GasThickness(G4LogicalVolume *vol, G4String value, G4String unit); //<auxiliary auxtype="setTRD_GasThickness" auxvalue="0.020" auxunit="mm"/>
    void setTRD_FoilMaterial(G4LogicalVolume *vol, G4String value); //<auxiliary auxtype="setTRD_FoilMaterial" auxvalue="G4_C"/>
    void setTRD_GasMaterial(G4LogicalVolume *vol, G4String value); //<auxiliary auxtype="setTRD_GasMaterial" auxvalue="G4_AIR"/>
    void setTRD_FoilNumber(G4LogicalVolume *vol, G4String value); //<auxiliary auxtype="setTRD_FoilNumber" auxvalue="2000"/>
    
    RadiatorDescription* GetRadiatorDescription() { return fRadiatorDescription; }
};

#endif
