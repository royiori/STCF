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
         instance("MyRICHDetector", ::MyRICHDetector::Class_Version(), "inc/MyRICHDetector.h", 13,
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
   static TClass *MyMainFrameGui_Dictionary();
   static void MyMainFrameGui_TClassManip(TClass*);
   static void delete_MyMainFrameGui(void *p);
   static void deleteArray_MyMainFrameGui(void *p);
   static void destruct_MyMainFrameGui(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MyMainFrameGui*)
   {
      ::MyMainFrameGui *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::MyMainFrameGui));
      static ::ROOT::TGenericClassInfo 
         instance("MyMainFrameGui", "inc/MyMainFrameGui.h", 24,
                  typeid(::MyMainFrameGui), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &MyMainFrameGui_Dictionary, isa_proxy, 0,
                  sizeof(::MyMainFrameGui) );
      instance.SetDelete(&delete_MyMainFrameGui);
      instance.SetDeleteArray(&deleteArray_MyMainFrameGui);
      instance.SetDestructor(&destruct_MyMainFrameGui);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MyMainFrameGui*)
   {
      return GenerateInitInstanceLocal((::MyMainFrameGui*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::MyMainFrameGui*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *MyMainFrameGui_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::MyMainFrameGui*)0x0)->GetClass();
      MyMainFrameGui_TClassManip(theClass);
   return theClass;
   }

   static void MyMainFrameGui_TClassManip(TClass* ){
   }

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
         instance("vector<string>", -2, "vector", 464,
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
         instance("vector<double>", -2, "vector", 464,
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
         instance("vector<TH2F*>", -2, "vector", 464,
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
class __attribute__((annotate(R"ATTRDUMP(Jet class)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$inc/MyRICHDetector.h")))  MyRICHDetector;
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
