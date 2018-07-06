// E/p analysis for run 2
// Joakim Olsson (joakim.olsson@cern.ch)

// c++ include(s):
#include <iostream>
#include <typeinfo>
#include <sstream>

// EL include(s):
#include <EventLoop/Job.h>
#include <EventLoop/Worker.h>

// EDM include(s):
#include "AthContainers/ConstDataVector.h"
#include "InDetTrackSelectionTool/InDetTrackSelectionTool.h"

// package include(s):
#include "xAODEventInfo/EventInfo.h"
#include "xAODAnaHelpers/HelperFunctions.h"
#include "EoverPAnalysis/TrackVertexSelection.h"
#include "AsgTools/MessageCheck.h"

// ROOT include(s):
#include "TFile.h"
#include "TObjArray.h"
#include "TObjString.h"

// this is needed to distribute the algorithm to the workers
ClassImp(TrackVertexSelection)


TrackVertexSelection :: TrackVertexSelection (std::string className) :
    Algorithm(className),
    m_cutflowHist(nullptr),
    m_cutflowHistW(nullptr),
    m_trk_cutflowHist_1(nullptr),
    m_trkSelection(nullptr)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
  Info("TrackVertexSelection()", "Calling constructor");

  // read debug flag from .config file
  m_debug         = false;
  m_useCutFlow    = true;

  // input container to be read from TEvent or TStore
  m_inContainerName  = "";

  // decorate selected objects that pass the cuts
  m_decorateSelectedObjects = true;
  // additional functionality : create output container of selected objects
  //                            using the SG::VIEW_ELEMENTS option
  //                            decorating and output container should not be mutually exclusive
  m_createSelectedContainer = false;
  // if requested, a new container is made using the SG::VIEW_ELEMENTS option
  m_outContainerName        = "";
  // if only want to look at a subset of object
  m_nToProcess              = -1;

  // cuts
  m_pass_max                = -1;
  m_pass_min                = -1;

  m_cutLevel                = "LoosePrimary"; 
  m_minPt                   = -1.;
  m_maxAbsEta               = 1e8;
  m_maxZ0SinTheta           = 1e8;
  m_maxZ0                   = 1e8;
  m_maxD0                   = 1e8;
  m_minNPixelHits           = -1;
  m_minNSctHits             = -1;
  m_minNSiHits              = -1;
  m_minNTrtHits             = -1;

  m_passAuxDecorKeys        = "";
  m_failAuxDecorKeys        = "";

}

EL::StatusCode TrackVertexSelection :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.

  Info("setupJob()", "Calling setupJob");

  job.useXAOD ();
  xAOD::Init( "TrackVertexSelection" ).ignore(); // call before opening first file

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TrackVertexSelection :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.

  Info("histInitialize()", "Calling histInitialize");
  ANA_CHECK( xAH::Algorithm::algInitialize());
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TrackVertexSelection :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed

  Info("fileExecute()", "Calling fileExecute");

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TrackVertexSelection :: changeInput (bool /*firstFile*/)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.

  Info("changeInput()", "Calling changeInput");

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TrackVertexSelection :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  if(m_useCutFlow) {
    TFile *file = wk()->getOutputFile ("cutflow");
    m_cutflowHist  = (TH1D*)file->Get("cutflow");
    m_cutflowHistW = (TH1D*)file->Get("cutflow_weighted");
    m_cutflow_bin  = m_cutflowHist->GetXaxis()->FindBin("track selection");
    m_cutflowHistW->GetXaxis()->FindBin("track selection");

    // retrieve the object cutflow
    m_trk_cutflowHist_1  = (TH1D*)file->Get("cutflow_trks_1");
    m_trk_cutflow_all    = m_trk_cutflowHist_1->GetXaxis()->FindBin("all");
    m_trk_cutflow_accept = m_trk_cutflowHist_1->GetXaxis()->FindBin("pass trk selection");
  }

  // parse and split by comma
  std::string token;
  std::istringstream ss(m_passAuxDecorKeys);
  while(std::getline(ss, token, ',')){
    m_passKeys.push_back(token);
  }
  ss.clear();
  ss.str(m_failAuxDecorKeys);
  while(std::getline(ss, token, ',')){
    m_failKeys.push_back(token);
  }

  // initialize and configure the track selection tool
  //------------------------------------------------------
  m_trkSelection = new InDet::InDetTrackSelectionTool("TrackSelection");
  ANA_CHECK( m_trkSelection->setProperty("CutLevel", m_cutLevel.c_str()) );
  ANA_CHECK( m_trkSelection->setProperty("minPt", static_cast<double>(m_minPt)*1000) ); 
  ANA_CHECK( m_trkSelection->setProperty("maxAbsEta", static_cast<double>(m_maxAbsEta)) ); 
  // ANA_CHECK( m_trkSelection->setProperty("minPt", 400), "failed to set minPt property"); 
  // ANA_CHECK( m_trkSelection->setProperty("maxAbsEta", 2.5), "failed to set maxAbsEta property"); 
  ANA_CHECK( m_trkSelection->setProperty("maxZ0SinTheta", static_cast<double>(m_maxZ0SinTheta)) ); 
  ANA_CHECK( m_trkSelection->setProperty("maxD0", static_cast<double>(m_maxD0)) );
  ANA_CHECK( m_trkSelection->setProperty("maxZ0", static_cast<double>(m_maxZ0)) );
  ANA_CHECK( m_trkSelection->setProperty("minNPixelHits", static_cast<int>(m_minNPixelHits)) ); 
  ANA_CHECK( m_trkSelection->setProperty("minNSctHits", static_cast<int>(m_minNSctHits)) ); 
  ANA_CHECK( m_trkSelection->setProperty("minNSiHits", static_cast<int>(m_minNSiHits)) ); 
  ANA_CHECK( m_trkSelection->setProperty("maxTrtEtaAcceptance", 0.0) ); 
  ANA_CHECK( m_trkSelection->setProperty("maxEtaForTrtHitCuts", 2.0) ); 
  if (m_minNTrtHits =! -1) ANA_CHECK( m_trkSelection->setProperty("minNTrtHits", static_cast<int>(m_minNTrtHits)) ); 
  ANA_CHECK( m_trkSelection->initialize()); 

  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();

  Info("initialize()", "Number of events in file: %lld ", m_event->getEntries() );

  m_numEvent      = 0;
  m_numObject     = 0;
  m_numEventPass  = 0;
  m_weightNumEventPass  = 0;
  m_numObjectPass = 0;

  Info("initialize()", "TrackVertexSelection Interface succesfully initialized!" );

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TrackVertexSelection :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  if(m_debug) Info("execute()", "Applying Track Selection... ");

  // retrieve event
  const xAOD::EventInfo* eventInfo(nullptr);
  ANA_CHECK( HelperFunctions::retrieve(eventInfo, m_eventInfoContainerName, m_event, m_store));

  // MC event weight
  float mcEvtWeight(1.0);
  static SG::AuxElement::Accessor< float > mcEvtWeightAcc("mcEventWeight");
  if ( ! mcEvtWeightAcc.isAvailable( *eventInfo ) ) {
    Error("execute()  ", "mcEventWeight is not available as decoration! Aborting" );
    return EL::StatusCode::FAILURE;
  }
  mcEvtWeight = mcEvtWeightAcc( *eventInfo );

  m_numEvent++;

  // get the collection from TEvent or TStore
  const xAOD::TrackParticleContainer* inTracks(nullptr);
  ANA_CHECK( HelperFunctions::retrieve(inTracks, m_inContainerName, m_event, m_store));

  // get primary vertex
  const xAOD::VertexContainer *vertices(nullptr);
  ANA_CHECK( HelperFunctions::retrieve(vertices, "PrimaryVertices", m_event, m_store));
  const xAOD::Vertex *pvx = HelperFunctions::getPrimaryVertex(vertices);


  // create output container (if requested) - deep copy

  ConstDataVector<xAOD::TrackParticleContainer>* selectedTracks = 0;
  if(m_createSelectedContainer) {
    selectedTracks = new ConstDataVector<xAOD::TrackParticleContainer>(SG::VIEW_ELEMENTS);
  }

  int nPass(0); int nObj(0);
  for(const auto trk : *inTracks){

    // if only looking at a subset of tracks make sure all are decorrated
    if( m_nToProcess > 0 && nObj >= m_nToProcess ) {
      if(m_decorateSelectedObjects) {
        trk->auxdecor< char >( "passSel" ) = -1;
      } else {
        break;
      }
      continue;
    }

    nObj++;
    int passSel = m_trkSelection->accept(*trk, pvx);  
    if(m_decorateSelectedObjects) {
      trk->auxdecor< char >( "passSel" ) = passSel;
    }

    if(passSel) {
      nPass++;
      if(m_createSelectedContainer) {
        selectedTracks->push_back( trk );
      }
    }
  }

  m_numObject     += nObj;
  m_numObjectPass += nPass;

  // apply event selection based on minimal/maximal requirements on the number of objects per event passing cuts
  if( m_pass_min > 0 && nPass < m_pass_min ) {
    wk()->skipEvent();
    return EL::StatusCode::SUCCESS;
  }
  if( m_pass_max > 0 && nPass > m_pass_max ) {
    wk()->skipEvent();
    return EL::StatusCode::SUCCESS;
  }

  // add output container to TStore
  if( m_createSelectedContainer ) {
    ANA_CHECK( m_store->record( selectedTracks, m_outContainerName ) );
  }

  m_numEventPass++;
  m_weightNumEventPass += mcEvtWeight;

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode TrackVertexSelection :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.

  if(m_debug) Info("postExecute()", "Calling postExecute");

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TrackVertexSelection :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.
  //
  if ( m_useCutFlow ) {
    Info("histFinalize()", "Filling cutflow");
    m_cutflowHist ->SetBinContent( m_cutflow_bin, m_numEventPass        );
    m_cutflowHistW->SetBinContent( m_cutflow_bin, m_weightNumEventPass  );

    // fill cutflow bin 'all' before any cut
    if(m_useCutFlow) m_trk_cutflowHist_1->Fill(m_trk_cutflow_all, m_numObject);
    if(m_useCutFlow) m_trk_cutflowHist_1->Fill(m_trk_cutflow_accept, m_numObjectPass );
  }

  Info("finalize()", "Deleting tool instances...");

  if ( m_trkSelection ) {
    delete m_trkSelection; m_trkSelection = nullptr;
  }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode TrackVertexSelection :: histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.
  ANA_CHECK( xAH::Algorithm::algFinalize() );
  return EL::StatusCode::SUCCESS;
}
