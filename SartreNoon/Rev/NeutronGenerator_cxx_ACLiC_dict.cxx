// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dIhomedIbhaktadISartremINoondISartreNoondIRevdINeutronGenerator_cxx_ACLiC_dict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
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

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "/home/bhakta/Sartre-Noon/SartreNoon/Rev/NeutronGenerator.cxx"

// Header files passed via #pragma extra_include

   namespace ROOTDict {
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance();
      static TClass *ROOT_Dictionary();

      // Function generating the singleton type initializer
      inline ::ROOT::TGenericClassInfo *GenerateInitInstance()
      {
         static ::ROOT::TGenericClassInfo 
            instance("ROOT", 0 /*version*/, "Rtypes.h", 105,
                     ::ROOT::Internal::DefineBehavior((void*)nullptr,(void*)nullptr),
                     &ROOT_Dictionary, 0);
         return &instance;
      }
      // Insure that the inline function is _not_ optimized away by the compiler
      ::ROOT::TGenericClassInfo *(*_R__UNIQUE_DICT_(InitFunctionKeeper))() = &GenerateInitInstance;  
      // Static variable to force the class initialization
      static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstance(); R__UseDummy(_R__UNIQUE_DICT_(Init));

      // Dictionary for non-ClassDef classes
      static TClass *ROOT_Dictionary() {
         return GenerateInitInstance()->GetClass();
      }

   }

namespace ROOT {
   static void *new_NeutronGenerator(void *p = nullptr);
   static void *newArray_NeutronGenerator(Long_t size, void *p);
   static void delete_NeutronGenerator(void *p);
   static void deleteArray_NeutronGenerator(void *p);
   static void destruct_NeutronGenerator(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::NeutronGenerator*)
   {
      ::NeutronGenerator *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::NeutronGenerator >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("NeutronGenerator", ::NeutronGenerator::Class_Version(), "NeutronGenerator.h", 16,
                  typeid(::NeutronGenerator), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::NeutronGenerator::Dictionary, isa_proxy, 4,
                  sizeof(::NeutronGenerator) );
      instance.SetNew(&new_NeutronGenerator);
      instance.SetNewArray(&newArray_NeutronGenerator);
      instance.SetDelete(&delete_NeutronGenerator);
      instance.SetDeleteArray(&deleteArray_NeutronGenerator);
      instance.SetDestructor(&destruct_NeutronGenerator);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::NeutronGenerator*)
   {
      return GenerateInitInstanceLocal(static_cast<::NeutronGenerator*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::NeutronGenerator*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr NeutronGenerator::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *NeutronGenerator::Class_Name()
{
   return "NeutronGenerator";
}

//______________________________________________________________________________
const char *NeutronGenerator::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NeutronGenerator*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int NeutronGenerator::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::NeutronGenerator*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *NeutronGenerator::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NeutronGenerator*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *NeutronGenerator::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::NeutronGenerator*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void NeutronGenerator::Streamer(TBuffer &R__b)
{
   // Stream an object of class NeutronGenerator.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(NeutronGenerator::Class(),this);
   } else {
      R__b.WriteClassBuffer(NeutronGenerator::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_NeutronGenerator(void *p) {
      return  p ? new(p) ::NeutronGenerator : new ::NeutronGenerator;
   }
   static void *newArray_NeutronGenerator(Long_t nElements, void *p) {
      return p ? new(p) ::NeutronGenerator[nElements] : new ::NeutronGenerator[nElements];
   }
   // Wrapper around operator delete
   static void delete_NeutronGenerator(void *p) {
      delete (static_cast<::NeutronGenerator*>(p));
   }
   static void deleteArray_NeutronGenerator(void *p) {
      delete [] (static_cast<::NeutronGenerator*>(p));
   }
   static void destruct_NeutronGenerator(void *p) {
      typedef ::NeutronGenerator current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::NeutronGenerator

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = nullptr);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "vector", 389,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 0,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<double>","std::vector<double, std::allocator<double> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<double>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<double>*>(nullptr))->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete (static_cast<vector<double>*>(p));
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] (static_cast<vector<double>*>(p));
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace {
  void TriggerDictionaryInitialization_NeutronGenerator_cxx_ACLiC_dict_Impl() {
    static const char* headers[] = {
"NeutronGenerator.cxx",
nullptr
    };
    static const char* includePaths[] = {
"/usr/local/include",
"/usr/local/etc/",
"/usr/local/etc//cling",
"/usr/local/etc//cling/plugins/include",
"/usr/local/include/",
"/usr/local/include/",
"/home/bhakta/Sartre-Noon/SartreNoon/Rev/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "NeutronGenerator_cxx_ACLiC_dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$NeutronGenerator.h")))  __attribute__((annotate("$clingAutoload$NeutronGenerator.cxx")))  NeutronGenerator;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "NeutronGenerator_cxx_ACLiC_dict dictionary payload"

#ifndef __ACLIC__
  #define __ACLIC__ 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "NeutronGenerator.cxx"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"", payloadCode, "@",
"NeutronGenerator", payloadCode, "@",
"NeutronGenerator::HadronicIntModel_t", payloadCode, "@",
"NeutronGenerator::RunMode_t", payloadCode, "@",
"NeutronGenerator::fgIsA", payloadCode, "@",
"ROOT::GenerateInitInstance", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("NeutronGenerator_cxx_ACLiC_dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_NeutronGenerator_cxx_ACLiC_dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_NeutronGenerator_cxx_ACLiC_dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_NeutronGenerator_cxx_ACLiC_dict() {
  TriggerDictionaryInitialization_NeutronGenerator_cxx_ACLiC_dict_Impl();
}
