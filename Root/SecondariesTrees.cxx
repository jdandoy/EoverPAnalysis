/*
  Root/SecondariesTrees.cxx

  M. LeBlanc (Arizona) <matt.leblanc@cern.ch>
  for Run 2 E/p studies on identified resonances.
*/

#include <AsgTools/MessageCheck.h>
#include <EoverPAnalysis/SecondariesTrees.h>

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <EventLoop/OutputStream.h>

#include "TFile.h"
#include "TSystem.h"
#include "TLorentzVector.h"

// EDM include(s):
#include "xAODCaloEvent/CaloClusterContainer.h"
#include "xAODCaloEvent/CaloCluster.h"
#include "xAODCore/ShallowCopy.h"
#include "AthContainers/ConstDataVector.h"
#include "xAODTracking/Vertex.h"

#include "xAODTracking/VertexContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODJet/JetAuxContainer.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODEventInfo/EventInfo.h"

SecondariesTrees :: SecondariesTrees (const std::string& name,
					  ISvcLocator *pSvcLocator)
  : EL::AnaAlgorithm (name, pSvcLocator)
{
  declareProperty( "label", m_label = "nominal",
                   "A label for the histograms" );

  declareProperty( "isData", m_isData = false,
                   "Is this data?" );

  declareProperty( "MessageFrequency", m_MessageFrequency = 1000,
		   "Frequency of debug messages" );

  declareProperty( "VertexContainer", m_VertexContainer = "PrimaryVertices",
                   "Input vertex container?" );
}

SecondariesTrees :: ~SecondariesTrees()
{
  // vector output branches need to be cleaned up
  //delete m_jetPt;
  //delete m_jetEta;
  //delete m_jetPhi;
  //delete m_jetM;
}

StatusCode SecondariesTrees :: initialize ()
{
  m_event = wk()->xaodEvent();

  // Histogram booking

  ANA_CHECK (book (TH1F ("h_npv_"+m_label, "h_npv_"+m_label, 60, 0, 60))); // NPV
  ANA_CHECK (book (TH1F ("h_mu_"+m_label, "h_mu_"+m_label, 80, 0, 80))); // mu

  ANA_CHECK (book (TH1F ("h_nConstituents_"+m_label, "h_nConstituents_"+m_label, 60, 0, 60))); // nConstituents

  ANA_CHECK (book (TH2F ("h2_DTRatioVsClusterE_"+m_label, "h2_DTRatioVsClusterE_"+m_label, 40, 0, 20, 120, 0, 1.2))); // DigiTruth bread-n-butter
  
  // Output tree
  // This left mostly as an example

  ANA_CHECK (book (TTree(TString("Tree_"+m_label).Data(), "Output Tree")) );
  TTree* t = tree (TString("Tree_"+m_label).Data());

  t->Branch("DSID",&m_DSID);
  t->Branch("weight",&m_weight);
  t->Branch("mcWeight",&m_mcWeight);
  t->Branch("mu",&m_mu);
  t->Branch("npv",&m_npv);
  
  t->Branch("vertex_N",&m_vertex_N);
  t->Branch("vertex_pt",&m_vertex_pt);
  t->Branch("vertex_eta",&m_vertex_eta);
  t->Branch("vertex_phi",&m_vertex_phi);
  t->Branch("vertex_mass",&m_vertex_mass);
  t->Branch("vertex_massErr",&m_vertex_massErr);
  t->Branch("vertex_Rxy",&m_vertex_Rxy);
  t->Branch("vertex_chiSquared",&m_vertex_chiSquared);
  t->Branch("vertex_numberDoF",&m_vertex_numberDoF);
  t->Branch("vertex_dr",&m_vertex_dr);
  
  t->Branch("track1_pt",&m_track1_pt);
  t->Branch("track2_pt",&m_track2_pt);
  t->Branch("track1_etaID",&m_track1_etaID);
  t->Branch("track1_phiID",&m_track1_phiID);
  t->Branch("track2_etaID",&m_track2_etaID);
  t->Branch("track2_phiID",&m_track2_phiID);
    
  return StatusCode::SUCCESS;
}

StatusCode SecondariesTrees :: execute ()
{
  static int count = 0;
  bool debug=false;
  if(count==0)
    {
      debug = true;
      std::cout << "SecondariesTrees :: execute() BEGIN" << std::endl;
    }
  count++;
  if(count % m_MessageFrequency == 0) debug = true; // Defaults to print every 1000 events

  if(debug) std::cout << "SecondariesTrees :: execute()\tProcessing event " << count << ". "  << std::endl;

  // NPV and mu
  int npv = 0;
  const xAOD::VertexContainer* PrimaryVertices = 0;
  m_event->retrieve(PrimaryVertices,m_VertexContainer.Data());
  for(const auto vertex : *PrimaryVertices)
    if(vertex->trackParticleLinks().size() >= 2)
      npv++;  
  m_npv = npv;
  m_vertex_N = npv;

  if(debug) std::cout << "SecondariesTrees :: execute()\t" << npv << " primary vertices!" << std::endl; 

  // Truth particles
  const xAOD::TruthParticleContainer* TruthParticles = 0;
  if(!m_isData) m_event->retrieve(TruthParticles, "TruthParticles");
  
  // Event weights
  const xAOD::EventInfo* eventInfo = 0;
  if( !m_event->retrieve( eventInfo, "EventInfo" ).isSuccess()) return EL::StatusCode::FAILURE;
  
  Float_t weight = 1.0;

  if(!m_isData)
    {
      Float_t mc_weight = eventInfo->mcEventWeight();
      Int_t dsid = eventInfo->mcChannelNumber();
      
      weight *= mc_weight;
      
      float xs = 1.;
      
      // R21
      if(dsid==361020) xs=(78420000000000./15527000.); // xSec*FE/evnts in fb !!!!!
      if(dsid==361021) xs=(78420000000000./15997000.);
      if(dsid==361022) xs=(2433200000000.*3.3423E-04/15983000.);
      if(dsid==361023) xs=(2433200000000.*3.3423E-04/15983000.);
      if(dsid==361024) xs=(254630000.*5.3138E-04/15978500.);
      if(dsid==361025) xs=(4553500.*9.2409E-04/15991500.);
      if(dsid==361026) xs=(257530.000*9.4092E-04/17834000.);
      if(dsid==361027) xs=(16215.*3.9280E-04/15118500.);
      if(dsid==361028) xs=(625.030*6.2503E-04/15997000.);
      if(dsid==361029) xs=(19.639*1.2076E-02/14552500.);
      if(dsid==361030) xs=(1.196*5.9087E-03/15998000.);
      if(dsid==361031) xs=(0.042*2.6761E-03/15973000.);
      if(dsid==361032) xs=(0.0010367*4.2592E-04/15070000.);
      
      weight *= xs;

      m_mcWeight = mc_weight;
      m_weight = weight;
    }
  
  for(const auto vertex : *PrimaryVertices)
    {
      //std::cout << vertex->trackParticleLinks().size() << " tracks on vertex" << std::endl;         
      TLorentzVector tlv_vertex;
      if(vertex->trackParticleLinks().size()!=2) continue;
      else
	{
	  //This is the track pT, and coordinates of the track in the ID	  	  	  
	  m_track1_pt = vertex->trackParticle(0)->pt()/1.e3;
	  m_track1_etaID = vertex->trackParticle(0)->eta();
	  m_track1_phiID = vertex->trackParticle(0)->phi();
	  //m_track1_m = vertex->trackParticle(0)->m()/1.e3;
	  
	  m_track2_pt = vertex->trackParticle(1)->pt()/1.e3;
	  m_track2_etaID = vertex->trackParticle(1)->eta();
          m_track2_phiID = vertex->trackParticle(1)->phi();
	  //m_track2_m = vertex->trackParticle(1)->m()/1.e3;

	  if(m_label=="Lambda")
	    {
	      if(m_track1_pt > m_track2_pt) m_track1_m = 0.938272;
	      else if (m_track1_pt < m_track2_pt) m_track2_m = 0.938272;
	    }
	  if(m_label=="Ks")
            {
              m_track1_m = 0.13957;
              m_track2_m = 0.13957;
            }
	  if(m_label=="Phi")
            {
              m_track1_m = 0.493677;
              m_track2_m = 0.493677;
            }
	}
      
      TLorentzVector tlv_track1;
      tlv_track1.SetPtEtaPhiM(m_track1_pt, 
			      m_track1_etaID,
			      m_track1_phiID,
			      m_track1_m);
      
      TLorentzVector tlv_track2;
      tlv_track2.SetPtEtaPhiM(m_track2_pt,
                              m_track2_etaID,
                              m_track2_phiID,
                              m_track2_m);

      tlv_vertex.Clear();
      tlv_vertex += tlv_track1;
      tlv_vertex += tlv_track2;

      if(debug) std::cout << "Track1 4vec:\tpt: " << tlv_track1.Pt() <<"\teta: "
			  << tlv_track1.Eta() << "\tphi: " << tlv_track1.Phi()
			  << "\tm: " << tlv_track1.M()
			  << std::endl;
      
      if(debug) std::cout << "Track2 4vec:\tpt: " << tlv_track2.Pt() <<"\teta: "
                          << tlv_track2.Eta() << "\tphi: " << tlv_track2.Phi()
                          << "\tm: " << tlv_track2.M()
                          << std::endl;
    
      if(debug) std::cout << "Vertex 4vec:\tpt: " << tlv_vertex.Pt() <<"\teta: "
			  << tlv_vertex.Eta() << "\tphi: " << tlv_vertex.Phi()
			  << "\tm: " << tlv_vertex.M()
			  << std::endl;

      static SG::AuxElement::ConstAccessor< float > acc_chiSquared("chiSquared");
      if(acc_chiSquared.isAvailable(*vertex) ) m_vertex_chiSquared = acc_chiSquared(*vertex);
      
      static SG::AuxElement::ConstAccessor< float > acc_numberDoF("numberDoF");
      if(acc_numberDoF.isAvailable(*vertex) ) m_vertex_numberDoF = acc_numberDoF(*vertex);

      if(m_label=="Lambda")
	{
	  static SG::AuxElement::ConstAccessor< float > acc_Lambda_mass("Lambda_mass");
	  static SG::AuxElement::ConstAccessor< float > acc_Lambda_massErr("Lambda_massErr");
	  if( acc_Lambda_mass.isAvailable(*vertex) ) m_vertex_mass = acc_Lambda_mass(*vertex);
	  if( acc_Lambda_massErr.isAvailable(*vertex) ) m_vertex_massErr = acc_Lambda_massErr(*vertex);
	}

      if(m_label=="Ks")
	{
          static SG::AuxElement::ConstAccessor< float > acc_Ks_mass("Ks_mass");
          static SG::AuxElement::ConstAccessor< float > acc_Ks_massErr("Ks_massErr");
          if( acc_Ks_mass.isAvailable(*vertex) ) m_vertex_mass = acc_Ks_mass(*vertex);
          if( acc_Ks_massErr.isAvailable(*vertex) ) m_vertex_massErr = acc_Ks_massErr(*vertex);
	}
      
      if(m_label=="Phi")
	{
          static SG::AuxElement::ConstAccessor< float > acc_Phi_mass("Phi_mass");
          static SG::AuxElement::ConstAccessor< float > acc_Phi_massErr("Phi_massErr");
          if( acc_Phi_mass.isAvailable(*vertex) ) m_vertex_mass = acc_Phi_mass(*vertex);
          if( acc_Phi_massErr.isAvailable(*vertex) ) m_vertex_massErr = acc_Phi_massErr(*vertex);
	}

      m_vertex_pt = tlv_vertex.Pt();
      m_vertex_eta = tlv_vertex.Eta();
      m_vertex_phi = tlv_vertex.Phi();

      m_vertex_dr = tlv_track1.DeltaR(tlv_track2);

      m_vertex_Rxy = vertex->position().perp();

      // Fill for every vertex??
      tree (TString("Tree_"+m_label).Data())->Fill ();
    }  
    
  // fill the branches  
  //tree (TString("Tree_"+m_label).Data())->Fill (); 

  return StatusCode::SUCCESS;
}


StatusCode SecondariesTrees :: finalize ()
{
  return StatusCode::SUCCESS;
}
