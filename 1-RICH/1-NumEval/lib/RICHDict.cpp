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
#include "inc/MyRICHDetector.h"
#include "inc/MyMainFrameGui.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_MyRICHDetector(void *p = 0);
   static void *newArray_MyRICHDetector(Long_t size, void *p);
   static void delete_MyRICHDetector(void *p);
   static void deleteArray_MyRICHDetector(void *p);
   static void destruct_MyRICHDetector(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MyRICHDetector*)
   {
      ::MyRICHDetector *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MyRICHDetector >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MyRICHDetector", ::MyRICHDetector::Class_Version(), "inc/MyRICHDetector.h", 14,
                  typeid(::MyRICHDetector), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MyRICHDetector::Dictionary, isa_proxy, 4,
                  sizeof(::MyRICHDetector) );
      instance.SetNew(&new_MyRICHDetector);
      instance.SetNewArray(&newArray_MyRICHDetector);
      instance.SetDelete(&delete_MyRICHDetector);
      instance.SetDeleteArray(&deleteArray_MyRICHDetector);
      instance.SetDestructor(&destruct_MyRICHDetector);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MyRICHDetector*)
   {
      return GenerateInitInstanceLocal((::MyRICHDetector*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::MyRICHDetector*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_MyMainFrameGui(void *p);
   static void deleteArray_MyMainFrameGui(void *p);
   static void destruct_MyMainFrameGui(void *p);
   static void streamer_MyMainFrameGui(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MyMainFrameGui*)
   {
      ::MyMainFrameGui *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MyMainFrameGui >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MyMainFrameGui", ::MyMainFrameGui::Class_Version(), "inc/MyMainFrameGui.h", 42,
                  typeid(::MyMainFrameGui), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MyMainFrameGui::Dictionary, isa_proxy, 16,
                  sizeof(::MyMainFrameGui) );
      instance.SetDelete(&delete_MyMainFrameGui);
      instance.SetDeleteArray(&deleteArray_MyMainFrameGui);
      instance.SetDestructor(&destruct_MyMainFrameGui);
      instance.SetStreamerFunc(&streamer_MyMainFrameGui);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MyMainFrameGui*)
   {
      return GenerateInitInstanceLocal((::MyMainFrameGui*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::MyMainFrameGui*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr MyRICHDetector::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MyRICHDetector::Class_Name()
{
   return "MyRICHDetector";
}

//______________________________________________________________________________
const char *MyRICHDetector::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MyRICHDetector*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MyRICHDetector::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MyRICHDetector*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MyRICHDetector::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MyRICHDetector*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MyRICHDetector::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MyRICHDetector*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr MyMainFrameGui::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MyMainFrameGui::Class_Name()
{
   return "MyMainFrameGui";
}

//______________________________________________________________________________
const char *MyMainFrameGui::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MyMainFrameGui*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MyMainFrameGui::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MyMainFrameGui*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MyMainFrameGui::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MyMainFrameGui*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MyMainFrameGui::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MyMainFrameGui*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void MyRICHDetector::Streamer(TBuffer &R__b)
{
   // Stream an object of class MyRICHDetector.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(MyRICHDetector::Class(),this);
   } else {
      R__b.WriteClassBuffer(MyRICHDetector::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_MyRICHDetector(void *p) {
      return  p ? new(p) ::MyRICHDetector : new ::MyRICHDetector;
   }
   static void *newArray_MyRICHDetector(Long_t nElements, void *p) {
      return p ? new(p) ::MyRICHDetector[nElements] : new ::MyRICHDetector[nElements];
   }
   // Wrapper around operator delete
   static void delete_MyRICHDetector(void *p) {
      delete ((::MyRICHDetector*)p);
   }
   static void deleteArray_MyRICHDetector(void *p) {
      delete [] ((::MyRICHDetector*)p);
   }
   static void destruct_MyRICHDetector(void *p) {
      typedef ::MyRICHDetector current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MyRICHDetector

//______________________________________________________________________________
void MyMainFrameGui::Streamer(TBuffer &R__b)
{
   // Stream an object of class MyMainFrameGui.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TGMainFrame::Streamer(R__b);
      R__b >> fSettingText;
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
      R__b >> fComboCmd;
      R__b >> fCommand;
      R__b >> fCommandBuf;
      R__b >> NPage;
      R__b.CheckByteCount(R__s, R__c, MyMainFrameGui::IsA());
   } else {
      R__c = R__b.WriteVersion(MyMainFrameGui::IsA(), kTRUE);
      TGMainFrame::Streamer(R__b);
      R__b << fSettingText;
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
      R__b << fComboCmd;
      R__b << fCommand;
      R__b << fCommandBuf;
      R__b << NPage;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_MyMainFrameGui(void *p) {
      delete ((::MyMainFrameGui*)p);
   }
   static void deleteArray_MyMainFrameGui(void *p) {
      delete [] ((::MyMainFrameGui*)p);
   }
   static void destruct_MyMainFrameGui(void *p) {
      typedef ::MyMainFrameGui current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_MyMainFrameGui(TBuffer &buf, void *obj) {
      ((::MyMainFrameGui*)obj)->::MyMainFrameGui::Streamer(buf);
   }
} // end of namespace ROOT for class ::MyMainFrameGui

namespace ROOT {
   static TClass *vectorlEstringgR_Dictionary();
   static void vectorlEstringgR_TClassManip(TClass*);
   static void *new_vectorlEstringgR(void *p = 0);
   static void *newArray_vectorlEstringgR(Long_t size, void *p);
   static void delete_vectorlEstringgR(void *p);
   static void deleteArray_vectorlEstringgR(void *p);
   static void destruct_vectorlEstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<string>*)
   {
      vector<string> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<string>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<string>", -2, "vector", 470,
                  typeid(vector<string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEstringgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<string>) );
      instance.SetNew(&new_vectorlEstringgR);
      instance.SetNewArray(&newArray_vectorlEstringgR);
      instance.SetDelete(&delete_vectorlEstringgR);
      instance.SetDeleteArray(&deleteArray_vectorlEstringgR);
      instance.SetDestructor(&destruct_vectorlEstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<string> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<string>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<string>*)0x0)->GetClass();
      vectorlEstringgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEstringgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string> : new vector<string>;
   }
   static void *newArray_vectorlEstringgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string>[nElements] : new vector<string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEstringgR(void *p) {
      delete ((vector<string>*)p);
   }
   static void deleteArray_vectorlEstringgR(void *p) {
      delete [] ((vector<string>*)p);
   }
   static void destruct_vectorlEstringgR(void *p) {
      typedef vector<string> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<string>

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = 0);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 470,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<double>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)0x0)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *vectorlETH2FmUgR_Dictionary();
   static void vectorlETH2FmUgR_TClassManip(TClass*);
   static void *new_vectorlETH2FmUgR(void *p = 0);
   static void *newArray_vectorlETH2FmUgR(Long_t size, void *p);
   static void delete_vectorlETH2FmUgR(void *p);
   static void deleteArray_vectorlETH2FmUgR(void *p);
   static void destruct_vectorlETH2FmUgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TH2F*>*)
   {
      vector<TH2F*> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TH2F*>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TH2F*>", -2, "vector", 470,
                  typeid(vector<TH2F*>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETH2FmUgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TH2F*>) );
      instance.SetNew(&new_vectorlETH2FmUgR);
      instance.SetNewArray(&newArray_vectorlETH2FmUgR);
      instance.SetDelete(&delete_vectorlETH2FmUgR);
      instance.SetDeleteArray(&deleteArray_vectorlETH2FmUgR);
      instance.SetDestructor(&destruct_vectorlETH2FmUgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TH2F*> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TH2F*>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETH2FmUgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TH2F*>*)0x0)->GetClass();
      vectorlETH2FmUgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETH2FmUgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETH2FmUgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TH2F*> : new vector<TH2F*>;
   }
   static void *newArray_vectorlETH2FmUgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TH2F*>[nElements] : new vector<TH2F*>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETH2FmUgR(void *p) {
      delete ((vector<TH2F*>*)p);
   }
   static void deleteArray_vectorlETH2FmUgR(void *p) {
      delete [] ((vector<TH2F*>*)p);
   }
   static void destruct_vectorlETH2FmUgR(void *p) {
      typedef vector<TH2F*> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TH2F*>

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
         instance("map<int,TGTextButton*>", -2, "map", 899,
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

namespace {
  void TriggerDictionaryInitialization_RICHDict_Impl() {
    static const char* headers[] = {
"inc/MyRICHDetector.h",
"inc/MyMainFrameGui.h",
0
    };
    static const char* includePaths[] = {
"/Users/chad/Work/src/programs/root/root-6.16.00/install/include",
"/Users/chad/Documents/GitHub/STCF/1-RICH/1-NumEval/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "RICHDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$inc/MyRICHDetector.h")))  MyRICHDetector;
class __attribute__((annotate("$clingAutoload$inc/MyMainFrameGui.h")))  MyMainFrameGui;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "RICHDict dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "inc/MyRICHDetector.h"
#include "inc/MyMainFrameGui.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"MyMainFrameGui", payloadCode, "@",
"MyRICHDetector", payloadCode, "@",
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
