// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME libdIRICHDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "inc/MyGuiMainFrame.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void delete_MyGuiMainFrame(void *p);
   static void deleteArray_MyGuiMainFrame(void *p);
   static void destruct_MyGuiMainFrame(void *p);
   static void streamer_MyGuiMainFrame(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MyGuiMainFrame*)
   {
      ::MyGuiMainFrame *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MyGuiMainFrame >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MyGuiMainFrame", ::MyGuiMainFrame::Class_Version(), "inc/MyGuiMainFrame.h", 35,
                  typeid(::MyGuiMainFrame), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MyGuiMainFrame::Dictionary, isa_proxy, 16,
                  sizeof(::MyGuiMainFrame) );
      instance.SetDelete(&delete_MyGuiMainFrame);
      instance.SetDeleteArray(&deleteArray_MyGuiMainFrame);
      instance.SetDestructor(&destruct_MyGuiMainFrame);
      instance.SetStreamerFunc(&streamer_MyGuiMainFrame);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MyGuiMainFrame*)
   {
      return GenerateInitInstanceLocal((::MyGuiMainFrame*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::MyGuiMainFrame*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr MyGuiMainFrame::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MyGuiMainFrame::Class_Name()
{
   return "MyGuiMainFrame";
}

//______________________________________________________________________________
const char *MyGuiMainFrame::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MyGuiMainFrame*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MyGuiMainFrame::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MyGuiMainFrame*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MyGuiMainFrame::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MyGuiMainFrame*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MyGuiMainFrame::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MyGuiMainFrame*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void MyGuiMainFrame::Streamer(TBuffer &R__b)
{
   // Stream an object of class MyGuiMainFrame.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TGMainFrame::Streamer(R__b);
      R__b >> fSettingText1;
      R__b >> fCTab;
      int R__i;
      for (R__i = 0; R__i < 100; R__i++)
         R__b >> fCA[R__i];
      {
         map<int,TGTextButton*> &R__stl =  butList;
         R__stl.clear();
         TClass *R__tcl2 = TBuffer::GetClass(typeid(class TGTextButton *));
         if (R__tcl2==0) {
            Error("butList streamer","Missing the TClass object for class TGTextButton *!");
            return;
         }
         int R__i, R__n;
         R__b >> R__n;
         for (R__i = 0; R__i < R__n; R__i++) {
            int R__t;
            R__b >> R__t;
            TGTextButton* R__t2;
            R__t2 = (TGTextButton*)R__b.ReadObjectAny(R__tcl2);
            typedef int Value_t;
            std::pair<Value_t const, class TGTextButton * > R__t3(R__t,R__t2);
            R__stl.insert(R__t3);
         }
      }
      {
         map<int,TGSlider*> &R__stl =  sldList;
         R__stl.clear();
         TClass *R__tcl2 = TBuffer::GetClass(typeid(class TGSlider *));
         if (R__tcl2==0) {
            Error("sldList streamer","Missing the TClass object for class TGSlider *!");
            return;
         }
         int R__i, R__n;
         R__b >> R__n;
         for (R__i = 0; R__i < R__n; R__i++) {
            int R__t;
            R__b >> R__t;
            TGSlider* R__t2;
            R__t2 = (TGSlider*)R__b.ReadObjectAny(R__tcl2);
            typedef int Value_t;
            std::pair<Value_t const, class TGSlider * > R__t3(R__t,R__t2);
            R__stl.insert(R__t3);
         }
      }
      R__b >> fComboCmd;
      R__b >> fCommand;
      R__b >> fCommandBuf;
      R__b >> const_cast< int &>( NResPAGE );
      R__b.CheckByteCount(R__s, R__c, MyGuiMainFrame::IsA());
   } else {
      R__c = R__b.WriteVersion(MyGuiMainFrame::IsA(), kTRUE);
      TGMainFrame::Streamer(R__b);
      R__b << fSettingText1;
      R__b << fCTab;
      int R__i;
      for (R__i = 0; R__i < 100; R__i++)
         R__b << fCA[R__i];
      {
         map<int,TGTextButton*> &R__stl =  butList;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            map<int,TGTextButton*>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << ((*R__k).first );
            R__b << ((*R__k).second);
            }
         }
      }
      {
         map<int,TGSlider*> &R__stl =  sldList;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            map<int,TGSlider*>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b << ((*R__k).first );
            R__b << ((*R__k).second);
            }
         }
      }
      R__b << fComboCmd;
      R__b << fCommand;
      R__b << fCommandBuf;
      R__b << const_cast< int &>( NResPAGE );
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_MyGuiMainFrame(void *p) {
      delete ((::MyGuiMainFrame*)p);
   }
   static void deleteArray_MyGuiMainFrame(void *p) {
      delete [] ((::MyGuiMainFrame*)p);
   }
   static void destruct_MyGuiMainFrame(void *p) {
      typedef ::MyGuiMainFrame current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_MyGuiMainFrame(TBuffer &buf, void *obj) {
      ((::MyGuiMainFrame*)obj)->::MyGuiMainFrame::Streamer(buf);
   }
} // end of namespace ROOT for class ::MyGuiMainFrame

namespace ROOT {
   static TClass *maplEintcOTGTextButtonmUgR_Dictionary();
   static void maplEintcOTGTextButtonmUgR_TClassManip(TClass*);
   static void *new_maplEintcOTGTextButtonmUgR(void *p = 0);
   static void *newArray_maplEintcOTGTextButtonmUgR(Long_t size, void *p);
   static void delete_maplEintcOTGTextButtonmUgR(void *p);
   static void deleteArray_maplEintcOTGTextButtonmUgR(void *p);
   static void destruct_maplEintcOTGTextButtonmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,TGTextButton*>*)
   {
      map<int,TGTextButton*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,TGTextButton*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,TGTextButton*>", -2, "map", 898,
                  typeid(map<int,TGTextButton*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOTGTextButtonmUgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,TGTextButton*>) );
      instance.SetNew(&new_maplEintcOTGTextButtonmUgR);
      instance.SetNewArray(&newArray_maplEintcOTGTextButtonmUgR);
      instance.SetDelete(&delete_maplEintcOTGTextButtonmUgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOTGTextButtonmUgR);
      instance.SetDestructor(&destruct_maplEintcOTGTextButtonmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,TGTextButton*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,TGTextButton*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOTGTextButtonmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,TGTextButton*>*)0x0)->GetClass();
      maplEintcOTGTextButtonmUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOTGTextButtonmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOTGTextButtonmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,TGTextButton*> : new map<int,TGTextButton*>;
   }
   static void *newArray_maplEintcOTGTextButtonmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,TGTextButton*>[nElements] : new map<int,TGTextButton*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOTGTextButtonmUgR(void *p) {
      delete ((map<int,TGTextButton*>*)p);
   }
   static void deleteArray_maplEintcOTGTextButtonmUgR(void *p) {
      delete [] ((map<int,TGTextButton*>*)p);
   }
   static void destruct_maplEintcOTGTextButtonmUgR(void *p) {
      typedef map<int,TGTextButton*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,TGTextButton*>

namespace ROOT {
   static TClass *maplEintcOTGSlidermUgR_Dictionary();
   static void maplEintcOTGSlidermUgR_TClassManip(TClass*);
   static void *new_maplEintcOTGSlidermUgR(void *p = 0);
   static void *newArray_maplEintcOTGSlidermUgR(Long_t size, void *p);
   static void delete_maplEintcOTGSlidermUgR(void *p);
   static void deleteArray_maplEintcOTGSlidermUgR(void *p);
   static void destruct_maplEintcOTGSlidermUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const map<int,TGSlider*>*)
   {
      map<int,TGSlider*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(map<int,TGSlider*>));
      static ::ROOT::TGenericClassInfo 
         instance("map<int,TGSlider*>", -2, "map", 898,
                  typeid(map<int,TGSlider*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &maplEintcOTGSlidermUgR_Dictionary, isa_proxy, 0,
                  sizeof(map<int,TGSlider*>) );
      instance.SetNew(&new_maplEintcOTGSlidermUgR);
      instance.SetNewArray(&newArray_maplEintcOTGSlidermUgR);
      instance.SetDelete(&delete_maplEintcOTGSlidermUgR);
      instance.SetDeleteArray(&deleteArray_maplEintcOTGSlidermUgR);
      instance.SetDestructor(&destruct_maplEintcOTGSlidermUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::MapInsert< map<int,TGSlider*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const map<int,TGSlider*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *maplEintcOTGSlidermUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const map<int,TGSlider*>*)0x0)->GetClass();
      maplEintcOTGSlidermUgR_TClassManip(theClass);
   return theClass;
   }

   static void maplEintcOTGSlidermUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_maplEintcOTGSlidermUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,TGSlider*> : new map<int,TGSlider*>;
   }
   static void *newArray_maplEintcOTGSlidermUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) map<int,TGSlider*>[nElements] : new map<int,TGSlider*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_maplEintcOTGSlidermUgR(void *p) {
      delete ((map<int,TGSlider*>*)p);
   }
   static void deleteArray_maplEintcOTGSlidermUgR(void *p) {
      delete [] ((map<int,TGSlider*>*)p);
   }
   static void destruct_maplEintcOTGSlidermUgR(void *p) {
      typedef map<int,TGSlider*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class map<int,TGSlider*>

namespace {
  void TriggerDictionaryInitialization_RICHDict_Impl() {
    static const char* headers[] = {
"inc/MyGuiMainFrame.h",
0
    };
    static const char* includePaths[] = {
"/Users/chad/Work/src/programs/root/root-6.16.00/install/include",
"/Users/chad/Documents/GitHub/STCF/1-RICH/0-Decode/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "RICHDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$inc/MyGuiMainFrame.h")))  MyGuiMainFrame;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "RICHDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "inc/MyGuiMainFrame.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"MyGuiMainFrame", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("RICHDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_RICHDict_Impl, {}, classesHeaders, /*has no C++ module*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_RICHDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_RICHDict() {
  TriggerDictionaryInitialization_RICHDict_Impl();
}
