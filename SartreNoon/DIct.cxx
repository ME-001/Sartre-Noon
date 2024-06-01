// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME DIct
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

// Header files passed as explicit arguments
#include "NeutronGenerator.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static void *new_NeutronGenerator(void *p = nullptr);
   static void *newArray_NeutronGenerator(Long_t size, void *p);
   static void delete_NeutronGenerator(void *p);
   static void deleteArray_NeutronGenerator(void *p);
   static void destruct_NeutronGenerator(void *p);
   static void streamer_NeutronGenerator(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::NeutronGenerator*)
   {
      ::NeutronGenerator *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::NeutronGenerator >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("NeutronGenerator", ::NeutronGenerator::Class_Version(), "NeutronGenerator.h", 30,
                  typeid(::NeutronGenerator), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::NeutronGenerator::Dictionary, isa_proxy, 16,
                  sizeof(::NeutronGenerator) );
      instance.SetNew(&new_NeutronGenerator);
      instance.SetNewArray(&newArray_NeutronGenerator);
      instance.SetDelete(&delete_NeutronGenerator);
      instance.SetDeleteArray(&deleteArray_NeutronGenerator);
      instance.SetDestructor(&destruct_NeutronGenerator);
      instance.SetStreamerFunc(&streamer_NeutronGenerator);
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

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      void *ptr_fRunMode = (void*)&fRunMode;
      R__b >> *reinterpret_cast<Int_t*>(ptr_fRunMode);
      void *ptr_fHadronicIntModel = (void*)&fHadronicIntModel;
      R__b >> *reinterpret_cast<Int_t*>(ptr_fHadronicIntModel);
      R__b >> iEvent;
      R__b >> hInputRapidity;
      R__b >> hInputMass;
      R__b >> fRapMin;
      R__b >> fRapMax;
      R__b >> fMassMin;
      R__b >> fMassMax;
      lineString.Streamer(R__b);
      R__b >> const_cast< Double_t &>( neutron_M );
      R__b >> const_cast< Double_t &>( hbarc );
      R__b >> const_cast< Double_t &>( hbarcmev );
      R__b >> const_cast< Double_t &>( pi );
      R__b >> const_cast< Double_t &>( alpha );
      R__b >> const_cast< Double_t &>( nucleus_R );
      R__b >> const_cast< Int_t &>( nSteps_impactPar );
      R__b >> const_cast< Int_t &>( nSteps_GG );
      R__b >> const_cast< Int_t &>( nSteps_R );
      R__b >> const_cast< Int_t &>( nSteps_Phi );
      R__b >> const_cast< Int_t &>( nSteps_Energy );
      R__b >> const_cast< Int_t &>( maxNeutrons );
      R__b >> nFluxes;
      R__b >> nucleus_Z;
      R__b >> nucleus_A;
      R__b >> beamGamma;
      R__b >> gammaTarget;
      R__b >> neutronSepThr;
      R__b >> saturationEnergy;
      R__b.ReadStaticArray((double*)energyGamma_Xn);
      R__b.ReadStaticArray((double*)xSection_Xn);
      R__b >> gSection_Nn;
      R__b >> hSection_Nn;
      R__b >> gPhotonFluxTable;
      R__b >> gNucleusBreakupTable;
      R__b >> hTwoPhotonFluxModulationTable;
      R__b >> gNuclearThickness;
      R__b >> gHadronicProbabilityTable;
      R__b >> hXsectionXn;
      R__b >> gMeanNeutronN;
      R__b >> gWidthNeutronN;
      R__b >> fitMeanNeutronN;
      R__b >> fitWidthNeutronN;
      R__b >> hENDF_2D;
      R__b >> hENDF_1D;
      R__b >> fENDFFile;
      R__b >> hBranchingRatioMap;
      R__b >> hEventBreakupMap;
      R__b >> gUnitaryLeak;
      R__b >> fEventTree;
      fParticles->Streamer(R__b);
      R__b >> fOutputFile;
      R__b >> fQAhistList;
      R__b >> hNeutronMultiplicity;
      R__b >> hEnergyGen;
      R__b >> hKinEnergyGen;
      R__b >> hMomGen;
      R__b >> hPhiGen;
      R__b >> hThetaGen;
      R__b >> hEnergyBoosted;
      R__b >> hMomBoosted;
      R__b >> hPhiBoosted;
      R__b >> hThetaBoosted;
      R__b >> hNeutronRapidity;
      R__b >> hNeutronEta;
      R__b >> hEnergyBin;
      R__b >> hEnergyForNeutronMulti;
      R__b >> hProbabilityXn;
      R__b >> hPhotonK;
      R__b >> hRapidityVM;
      R__b >> hMassVM;
      R__b >> kStoreQA;
      R__b >> kStoreGen;
      R__b.CheckByteCount(R__s, R__c, NeutronGenerator::IsA());
   } else {
      R__c = R__b.WriteVersion(NeutronGenerator::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << (Int_t)fRunMode;
      R__b << (Int_t)fHadronicIntModel;
      R__b << iEvent;
      R__b << hInputRapidity;
      R__b << hInputMass;
      R__b << fRapMin;
      R__b << fRapMax;
      R__b << fMassMin;
      R__b << fMassMax;
      lineString.Streamer(R__b);
      R__b << const_cast< Double_t &>( neutron_M );
      R__b << const_cast< Double_t &>( hbarc );
      R__b << const_cast< Double_t &>( hbarcmev );
      R__b << const_cast< Double_t &>( pi );
      R__b << const_cast< Double_t &>( alpha );
      R__b << const_cast< Double_t &>( nucleus_R );
      R__b << const_cast< Int_t &>( nSteps_impactPar );
      R__b << const_cast< Int_t &>( nSteps_GG );
      R__b << const_cast< Int_t &>( nSteps_R );
      R__b << const_cast< Int_t &>( nSteps_Phi );
      R__b << const_cast< Int_t &>( nSteps_Energy );
      R__b << const_cast< Int_t &>( maxNeutrons );
      R__b << nFluxes;
      R__b << nucleus_Z;
      R__b << nucleus_A;
      R__b << beamGamma;
      R__b << gammaTarget;
      R__b << neutronSepThr;
      R__b << saturationEnergy;
      R__b.WriteArray(energyGamma_Xn, 625);
      R__b.WriteArray(xSection_Xn, 625);
      R__b << gSection_Nn;
      R__b << (TObject*)hSection_Nn;
      R__b << gPhotonFluxTable;
      R__b << gNucleusBreakupTable;
      R__b << (TObject*)hTwoPhotonFluxModulationTable;
      R__b << gNuclearThickness;
      R__b << gHadronicProbabilityTable;
      R__b << (TObject*)hXsectionXn;
      R__b << gMeanNeutronN;
      R__b << gWidthNeutronN;
      R__b << fitMeanNeutronN;
      R__b << fitWidthNeutronN;
      R__b << (TObject*)hENDF_2D;
      R__b << (TObject*)hENDF_1D;
      R__b << fENDFFile;
      R__b << (TObject*)hBranchingRatioMap;
      R__b << (TObject*)hEventBreakupMap;
      R__b << gUnitaryLeak;
      R__b << fEventTree;
      fParticles->Streamer(R__b);
      R__b << fOutputFile;
      R__b << fQAhistList;
      R__b << (TObject*)hNeutronMultiplicity;
      R__b << (TObject*)hEnergyGen;
      R__b << (TObject*)hKinEnergyGen;
      R__b << (TObject*)hMomGen;
      R__b << (TObject*)hPhiGen;
      R__b << (TObject*)hThetaGen;
      R__b << (TObject*)hEnergyBoosted;
      R__b << (TObject*)hMomBoosted;
      R__b << (TObject*)hPhiBoosted;
      R__b << (TObject*)hThetaBoosted;
      R__b << (TObject*)hNeutronRapidity;
      R__b << (TObject*)hNeutronEta;
      R__b << (TObject*)hEnergyBin;
      R__b << (TObject*)hEnergyForNeutronMulti;
      R__b << (TObject*)hProbabilityXn;
      R__b << (TObject*)hPhotonK;
      R__b << (TObject*)hRapidityVM;
      R__b << (TObject*)hMassVM;
      R__b << kStoreQA;
      R__b << kStoreGen;
      R__b.SetByteCount(R__c, kTRUE);
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
   // Wrapper around a custom streamer member function.
   static void streamer_NeutronGenerator(TBuffer &buf, void *obj) {
      ((::NeutronGenerator*)obj)->::NeutronGenerator::Streamer(buf);
   }
} // end of namespace ROOT for class ::NeutronGenerator

namespace {
  void TriggerDictionaryInitialization_DIct_Impl() {
    static const char* headers[] = {
"NeutronGenerator.h",
nullptr
    };
    static const char* includePaths[] = {
"/usr/local/include/",
"/home/bhakta/ModSartre/sartre/src/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "DIct dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
class __attribute__((annotate("$clingAutoload$NeutronGenerator.h")))  NeutronGenerator;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "DIct dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "NeutronGenerator.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"NeutronGenerator", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("DIct",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_DIct_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_DIct_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_DIct() {
  TriggerDictionaryInitialization_DIct_Impl();
}
