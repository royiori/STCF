#ifndef MyDetectorReader_H
#define MyDetectorReader_H 1

#include <map>
#include "G4GDMLReadStructure.hh"

class G4VisAttributes;
class G4AssemblyVolume;

/// GDML reader for the color attributes

class MyDetectorReader : public G4GDMLReadStructure
{

public:
    MyDetectorReader();
    ~MyDetectorReader();

    void ExtensionRead(const xercesc::DOMElement *const element);
    void ColorRead(const xercesc::DOMElement *const element);

    G4VisAttributes *GetVisAttribute(const G4String &ref);

protected:
    virtual void VolumeRead(const xercesc::DOMElement *const);
    virtual void Volume_contentRead(const xercesc::DOMElement *const);

    G4LogicalVolume *MyFileRead(const xercesc::DOMElement *const fileElement);
    void MyPhysvolRead(const xercesc::DOMElement *const, G4AssemblyVolume *assembly = 0);

private:
    std::map<G4String, G4VisAttributes *> fAttribs;
};

#endif
