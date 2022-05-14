// rootupler for AOD (the information on the associate tracks of the PV) is store different in other formats

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "TLorentzVector.h"
#include "TVector3.h"
#include "TTree.h"
#include <vector>
#include <sstream>

class HIchicRootupler:public edm::EDAnalyzer {
      public:
	explicit HIchicRootupler(const edm::ParameterSet &);
	~HIchicRootupler() override {};
	static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

      private:
	void analyze(const edm::Event &, const edm::EventSetup &) override;
	bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);
	UInt_t getTriggerBits(const edm::Event &);

	std::string file_name;
        edm::EDGetTokenT<pat::CompositeCandidateCollection> chic_;
        edm::EDGetTokenT<pat::CompositeCandidateCollection> psi_;
        edm::EDGetTokenT<pat::CompositeCandidateCollection> refit1_;
        edm::EDGetTokenT<reco::VertexCollection>            primaryVertices_;
        edm::EDGetTokenT<edm::TriggerResults>               triggerResults_;
        edm::EDGetTokenT<reco::BeamSpot> thebeamspot_;

	bool isMC_,bestCandidateOnly_;
        std::vector<std::string> FilterNames_;

	UInt_t    run;
        ULong64_t event;
        UInt_t    lumiblock;

	TLorentzVector chic_p4;
	TLorentzVector dimuon_p4;
	TLorentzVector muonP_p4;
	TLorentzVector muonM_p4;
	TLorentzVector photon_p4;

	TLorentzVector rf1S_chic_p4;
	Double_t invm1S;
        Double_t rf1S_vProb,vProb;
	Double_t y1S_nsigma;

	Double_t ele_lowerPt_pt;
	Double_t ele_higherPt_pt;
        Double_t deltapi0, distMinA, vtxProb;

	Double_t ctpv, ctpv_error, ctbs, ctbs_error;
	Double_t rf1S_ctpv, rf1S_ctpv_error, rf1S_ctbs, rf1S_ctbs_error;
	Double_t conv_vertex;
	Double_t dz;

        Double_t z_pv, dz_pv_11, dz_pv_12, dz_pv_21, dz_pv_22;
        Double_t sigma_tk1_vtx, sigma_tk2_vtx;
        Double_t dPhiTvx, vtxChi2;
        Int_t    vtxNDoF, nTracks;
        Double_t conv_vertex_rho;
        Int_t    conv_tkvtx_comp, conv_comp_ihits, conv_high_purity, conv_reject_pi0;
        Int_t    conv_algo, conv_qual_high_purity, conv_qual_generaltks;

        Double_t tk1_chi2, tk2_chi2, tk1_ndof, tk2_ndof;

	UInt_t photon_flags;
	UInt_t numPrimaryVertices,numDiMuons,numChis,numKVFChis;
	UInt_t trigger;
	UInt_t rf1S_rank;
        UInt_t nTrk_Vtx_Q, nTrk_Vtx;

        TVector3       pvtx_s;
        TVector3       pvtx_b;
        TVector3       psivtx;
        TVector3       chivtx;

	TTree *chic_tree;

	Int_t chic_pdgId;
	Int_t psi_pdgId;
	TLorentzVector gen_chic_p4;
	TLorentzVector gen_psi_p4;
        TLorentzVector gen_dimuon_p4;
	TLorentzVector gen_photon_p4;
	TLorentzVector gen_muonP_p4;
	TLorentzVector gen_muonM_p4;
        TVector3       gen_pvtx;
        TVector3       gen_psivtx;
        TVector3       gen_chivtx;

        edm::EDGetTokenT<reco::GenParticleCollection> genCands_;

        TTree *psi_tree;
        TLorentzVector mumu_p4, muP_p4, muM_p4;
        UInt_t mumu_rank;

  int evt_total;
  int evt_pass;

};

static const double pi0_mass =  0.134977;
static const double y1SMass  =  3.0969;

/*
// 2011 par
static const double Y_sig_par_A = 0.058;
static const double Y_sig_par_B = 0.047;
static const double Y_sig_par_C = 0.22;
*/

// 2012 par
static const double Y_sig_par_A = 62.62;
static const double Y_sig_par_B = 56.3;
static const double Y_sig_par_C = -20.77;

HIchicRootupler::HIchicRootupler(const edm::ParameterSet & iConfig): 
chic_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("chic_cand"))),
psi_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("psi_cand"))),
refit1_(consumes<pat::CompositeCandidateCollection>(iConfig.getParameter < edm::InputTag > ("refit1S"))),
primaryVertices_(consumes<reco::VertexCollection>(iConfig.getParameter < edm::InputTag > ("primaryVertices"))),
triggerResults_(consumes<edm::TriggerResults>(iConfig.getParameter < edm::InputTag > ("TriggerResults"))), 
thebeamspot_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpotTag"))),
isMC_(iConfig.getParameter < bool > ("isMC")),
bestCandidateOnly_(iConfig.getParameter < bool > ("only_best")),
FilterNames_(iConfig.getParameter<std::vector<std::string>>("FilterNames"))
{
    edm::Service < TFileService > fs;
    chic_tree = fs->make < TTree > ("chicTree", "Tree of mumugamma");

    chic_tree->Branch("run",      &run,      "run/i");
    chic_tree->Branch("event",    &event,    "event/l");
    chic_tree->Branch("lumiblock",&lumiblock,"lumiblock/i");


    chic_tree->Branch("chic_p4",   "TLorentzVector", &chic_p4);
    chic_tree->Branch("dimuon_p4", "TLorentzVector", &dimuon_p4);
    chic_tree->Branch("muonP_p4",  "TLorentzVector", &muonP_p4);
    chic_tree->Branch("muonM_p4",  "TLorentzVector", &muonM_p4);
    chic_tree->Branch("photon_p4", "TLorentzVector", &photon_p4);

    chic_tree->Branch("rf1S_chic_p4", "TLorentzVector", &rf1S_chic_p4);
    chic_tree->Branch("pvtx_s",       "TVector3",       &pvtx_s);
    chic_tree->Branch("pvtx_b",       "TVector3",       &pvtx_b);
    chic_tree->Branch("psivtx",       "TVector3",       &psivtx);
    chic_tree->Branch("chivtx",       "TVector3",       &chivtx);

    chic_tree->Branch("invm1S",       &invm1S,          "invm1S/D");
    chic_tree->Branch("rf1S_vProb",   &rf1S_vProb,      "rf1S_vProb/D");
    chic_tree->Branch("y1S_nsigma",   &y1S_nsigma,      "y1S_nsigma/D");

    chic_tree->Branch("ele_lowerPt_pt",  &ele_lowerPt_pt,  "ele_lowerPt_pt/D");
    chic_tree->Branch("ele_higherPt_pt", &ele_higherPt_pt, "ele_higherPt_pt/D");

    chic_tree->Branch("deltapi0",        &deltapi0,        "deltapi0/D");
    chic_tree->Branch("distMinA",        &distMinA,        "distMinA/D");
    chic_tree->Branch("vtxProb",         &vtxProb,         "vtxProb/D");

    chic_tree->Branch("dPhiTvx",         &dPhiTvx,         "dPhiTvx/D");
    chic_tree->Branch("vtxChi2",         &vtxChi2,         "vtxChi2/D");
    chic_tree->Branch("vtxNDoF",         &vtxNDoF,         "vtxNDoF/i");
    chic_tree->Branch("nTracks",         &nTracks,         "nTracks/i");

    chic_tree->Branch("z_pv",            &z_pv,            "z_pv/D");
    chic_tree->Branch("dz_pv_11",        &dz_pv_11,        "dz_pv_11/D");
    chic_tree->Branch("dz_pv_12",        &dz_pv_12,        "dz_pv_12/D");
    chic_tree->Branch("dz_pv_21",        &dz_pv_21,        "dz_pv_21/D");
    chic_tree->Branch("dz_pv_22",        &dz_pv_22,        "dz_pv_22/D");

    chic_tree->Branch("sigma_tk1_vtx",   &sigma_tk1_vtx,   "sigma_tk1_vtx/D");
    chic_tree->Branch("sigma_tk2_vtx",   &sigma_tk2_vtx,   "sigma_tk2_vtx/D");
    chic_tree->Branch("tk1_chi2",        &tk1_chi2,        "tk1_chi2/D");
    chic_tree->Branch("tk2_chi2",        &tk2_chi2,        "tk2_chi2/D");
    chic_tree->Branch("tk1_ndof",        &tk1_ndof,        "tk1_ndof/D");
    chic_tree->Branch("tk2_ndof",        &tk2_ndof,        "tk2_ndof/D");

    chic_tree->Branch("conv_vertex_rho",       &conv_vertex_rho,       "conv_vertex_rho/D");
    chic_tree->Branch("conv_tkvtx_comp",       &conv_tkvtx_comp,       "conv_tkvtx_comp/i");
    chic_tree->Branch("conv_comp_ihits",       &conv_comp_ihits,       "conv_comp_ihits/i");
    chic_tree->Branch("conv_high_purity",      &conv_high_purity,      "conv_high_purity/i");
    chic_tree->Branch("conv_reject_pi0",       &conv_reject_pi0,       "conv_reject_pi0/i");
    chic_tree->Branch("conv_algo",             &conv_algo,             "conv_algo/i");
    chic_tree->Branch("conv_qual_high_purity", &conv_qual_high_purity, "conv_qual_high_purity/i");
    chic_tree->Branch("conv_qual_generaltks",  &conv_qual_generaltks,  "conv_qual_generaltks/i");

    chic_tree->Branch("vProb",          &vProb,          "vProb/D");
    chic_tree->Branch("ctpv",           &ctpv,           "ctpv/D");
    chic_tree->Branch("ctpv_error",     &ctpv_error,     "ctpv_error/D");
    chic_tree->Branch("ctbs",           &ctbs,           "ctbs/D");
    chic_tree->Branch("ctbs_error",     &ctbs_error,     "ctbs_error/D");
    chic_tree->Branch("rf1S_ctpv",      &rf1S_ctpv,      "rf1S_ctpv/D");
    chic_tree->Branch("rf1S_ctpv_error",&rf1S_ctpv_error,"rf1S_ctpv_error/D");
    chic_tree->Branch("rf1S_ctbs",      &rf1S_ctbs,      "rf1S_ctbs/D");
    chic_tree->Branch("rf1S_ctbs_error",&rf1S_ctbs_error,"rf1S_ctbs_error/D");
    chic_tree->Branch("conv_vertex",    &conv_vertex,    "conv_vertex/D");
    chic_tree->Branch("dz",             &dz,             "dz/D");

    chic_tree->Branch("photon_flags",   &photon_flags,   "photon_flags/i");

    chic_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/i");
    chic_tree->Branch("numDiMuons",         &numDiMuons,         "numDiMuons/i");
    chic_tree->Branch("numChis",            &numChis,            "numChis/i");
    chic_tree->Branch("numKVFChis",         &numKVFChis,         "numKVFChis/i");

    chic_tree->Branch("trigger",            &trigger,            "trigger/i");
    chic_tree->Branch("nTrk_Vtx_Q",         &nTrk_Vtx_Q,         "nTrk_Vtx_Q/i");
    chic_tree->Branch("nTrk_Vtx",           &nTrk_Vtx,           "nTrk_Vtx/i");
    chic_tree->Branch("rf1S_rank",          &rf1S_rank,          "rf1S_rank/i");

    if (isMC_) {
       chic_tree->Branch("chic_pdgId",    &chic_pdgId,       "chic_pdgId/I");
       chic_tree->Branch("psi_pdgId",     &psi_pdgId,        "psi_pdgId/I");
       chic_tree->Branch("gen_chic_p4",   "TLorentzVector",  &gen_chic_p4);
       chic_tree->Branch("gen_psi_p4",    "TLorentzVector",  &gen_psi_p4);
       chic_tree->Branch("gen_dimuon_p4", "TLorentzVector",  &gen_dimuon_p4);
       chic_tree->Branch("gen_photon_p4", "TLorentzVector",  &gen_photon_p4);
       chic_tree->Branch("gen_muonP_p4",  "TLorentzVector",  &gen_muonP_p4);
       chic_tree->Branch("gen_muonM_p4",  "TLorentzVector",  &gen_muonM_p4);
       chic_tree->Branch("gen_pvtx",      "TVector3",        &gen_pvtx);
       chic_tree->Branch("gen_psivtx",    "TVector3",        &gen_psivtx);
       chic_tree->Branch("gen_chivtx",    "TVector3",        &gen_chivtx);
    }

    psi_tree = fs->make<TTree>("psiTree","Tree of dimuons");
    psi_tree->Branch("run",       &run,             "run/i");
    psi_tree->Branch("event",     &event,           "event/l");
    psi_tree->Branch("lumiblock", &lumiblock,       "lumiblock/i");
    psi_tree->Branch("numDiMuons",&numDiMuons,         "numDiMuons/i");
    psi_tree->Branch("mumu_p4",   "TLorentzVector", &mumu_p4);
    psi_tree->Branch("muP_p4",    "TLorentzVector", &muP_p4);
    psi_tree->Branch("muM_p4",    "TLorentzVector", &muM_p4);
    psi_tree->Branch("vProb",     &vProb,           "vProb/D");
    psi_tree->Branch("ctpv",      &ctpv,            "ctpv/D");
    psi_tree->Branch("ctpv_error",&ctpv_error,      "ctpv_error/D");
    psi_tree->Branch("ctbs",      &ctbs,            "ctbs/D");
    psi_tree->Branch("ctbs_error",&ctbs_error,      "ctbs_error/D");
    psi_tree->Branch("pvtx_s",    "TVector3",       &pvtx_s);
    psi_tree->Branch("pvtx_b",    "TVector3",       &pvtx_b);
    psi_tree->Branch("psivtx",    "TVector3",       &psivtx);
    psi_tree->Branch("trigger",   &trigger,         "trigger/i");
    psi_tree->Branch("nTrk_Vtx_Q",         &nTrk_Vtx_Q,         "nTrk_Vtx_Q/i");
    psi_tree->Branch("nTrk_Vtx",           &nTrk_Vtx,           "nTrk_Vtx/i");
    psi_tree->Branch("numPrimaryVertices", &numPrimaryVertices, "numPrimaryVertices/i");
    psi_tree->Branch("mumu_rank", &mumu_rank,       "mumu_rank/i"); 
    if (isMC_) {
       psi_tree->Branch("chic_pdgId",    &chic_pdgId,       "chic_pdgId/I");
       psi_tree->Branch("psi_pdgId",     &psi_pdgId,        "psi_pdgId/I");
       psi_tree->Branch("gen_chic_p4",   "TLorentzVector",  &gen_chic_p4);
       psi_tree->Branch("gen_psi_p4",    "TLorentzVector",  &gen_psi_p4);
       psi_tree->Branch("gen_dimuon_p4", "TLorentzVector",  &gen_dimuon_p4);
       psi_tree->Branch("gen_photon_p4", "TLorentzVector",  &gen_photon_p4);
       psi_tree->Branch("gen_muonP_p4",  "TLorentzVector",  &gen_muonP_p4);
       psi_tree->Branch("gen_muonM_p4",  "TLorentzVector",  &gen_muonM_p4);
       psi_tree->Branch("gen_pvtx",      "TVector3",        &gen_pvtx);
       psi_tree->Branch("gen_psivtx",    "TVector3",        &gen_psivtx);
       psi_tree->Branch("gen_chivtx",    "TVector3",        &gen_chivtx);
    }
    genCands_ = consumes<reco::GenParticleCollection>((edm::InputTag)"genParticles");
    evt_total = 0;
    evt_pass  = 0;

}

//Check recursively if any ancestor of particle is the given one
bool HIchicRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   if (particle->numberOfMothers() && isAncestor(ancestor,particle->mother(0))) return true;
   return false;
}

UInt_t HIchicRootupler::getTriggerBits(const edm::Event& iEvent ) {
   UInt_t trigger = 0;
   edm::Handle<edm::TriggerResults> triggerresults;
   iEvent.getByToken(triggerResults_, triggerresults);
   if (triggerresults.isValid()) {
      const edm::TriggerNames & TheTriggerNames = iEvent.triggerNames(*triggerresults);
      for (unsigned int i = 0; i < FilterNames_.size(); i++) {
         for (int version = 1; version < 99; version++) {
            std::stringstream ss;
            ss << FilterNames_[i] << "_v" << version;
            unsigned int bit = TheTriggerNames.triggerIndex(edm::InputTag(ss.str()).label());
            if (bit < triggerresults->size() && triggerresults->accept(bit) && !triggerresults->error(bit)) {
               trigger += (1<<i);
               break;
            }
         }
      }
   } else std::cout << "MMGrootupler::getTriggerBits: *** NO triggerResults found *** " << iEvent.id().run() << "," << iEvent.id().event() << std::endl;
   return trigger;
}

// ------------ method called for each event  ------------
void HIchicRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  edm::Handle < pat::CompositeCandidateCollection >chic_cand_handle;
  iEvent.getByToken(chic_, chic_cand_handle);

  edm::Handle < pat::CompositeCandidateCollection >psi_hand;
  iEvent.getByToken(psi_, psi_hand);

  edm::Handle < pat::CompositeCandidateCollection >refit1S_handle;
  iEvent.getByToken(refit1_, refit1S_handle);

  edm::Handle < reco::VertexCollection  >primaryVertices_handle;
  iEvent.getByToken(primaryVertices_, primaryVertices_handle);

  edm::Handle < edm::TriggerResults > triggerResults_handle;
  iEvent.getByToken(triggerResults_, triggerResults_handle);

  reco::Vertex theBeamSpotV;
  edm::Handle<reco::BeamSpot> theBeamSpot;
  iEvent.getByToken(thebeamspot_, theBeamSpot);
  reco::BeamSpot bs = *theBeamSpot;
  theBeamSpotV = reco::Vertex(bs.position(), bs.covariance3D());

  numPrimaryVertices = primaryVertices_handle->size();
  numDiMuons         = (psi_hand.isValid() && !psi_hand->empty()) ? psi_hand->size() : 0;
  numChis            = (chic_cand_handle.isValid() && !chic_cand_handle->empty()) ? chic_cand_handle->size() : 0;
  numKVFChis         = (refit1S_handle.isValid() && !refit1S_handle->empty()) ? refit1S_handle->size() : 0 ;

  run       = iEvent.id().run();
  event     = iEvent.id().event();
  lumiblock = iEvent.id().luminosityBlock();

  pat::CompositeCandidate chic_cand;
  pat::CompositeCandidate refit1S;

  edm::Handle<reco::GenParticleCollection> pruned;
  iEvent.getByToken(genCands_,pruned);

  evt_total++;
  if (isMC_ && pruned.isValid()) {
   gen_chic_p4.SetPtEtaPhiM(0, 0, 0, 0);
   gen_psi_p4.SetPtEtaPhiM(0, 0, 0, 0);
   gen_dimuon_p4.SetPtEtaPhiM(0, 0, 0, 0);
   chic_pdgId = 0;
   for (size_t i=0; i<pruned->size(); i++) {
      int p_id = abs((*pruned)[i].pdgId());
      int p_status = (*pruned)[i].status();
      psi_pdgId = 0;
      int foundit = 0;
      if ( ( p_id == 20443 || p_id == 445 || p_id == 10441) && p_status == 2)  psi_pdgId = 443;
      if (psi_pdgId > 0) {
         chic_pdgId = p_id;
         foundit++;
         const reco::Candidate * pwave = &(*pruned)[i];
         gen_chic_p4.SetPtEtaPhiM(pwave->pt(),pwave->eta(),pwave->phi(),pwave->mass());
         gen_pvtx.SetXYZ(pwave->vx(),pwave->vy(),pwave->vz());
         for (size_t j=0; j<pwave->numberOfDaughters(); j++) {
            const reco::Candidate *dau = pwave->daughter(j);
            if (dau->pdgId() == psi_pdgId && dau->status() == 2) {
               gen_psi_p4.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
               gen_chivtx.SetXYZ(dau->vx(),dau->vy(),dau->vz());
               uint nmuons = 0;
               for (size_t k=0; k<dau->numberOfDaughters(); k++) {
                  const reco::Candidate *gdau = dau->daughter(k);
                  if (gdau->pdgId() == 13 && gdau->status()==1) {
                     nmuons++;
                     gen_muonM_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
                     gen_psivtx.SetXYZ(gdau->vx(),gdau->vy(),gdau->vz());
                  } else {
                     if (gdau->pdgId() == -13 && gdau->status()==1) {
                        nmuons++;
                        gen_muonP_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
                     } 
                  }
               }
               if (nmuons == 2 ) {
                  foundit += 3;                                  // found complete dimuon decay
                  gen_dimuon_p4 = gen_muonM_p4 + gen_muonP_p4;   // will account fsr
               }
            } else {
               if (dau->pdgId() == 22 && dau->status() ==1) { 
                  foundit++;
                  gen_photon_p4.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
               }  else std::cout << "Rootupler: unexpected pdg_id " << dau->pdgId() << " (" << run << "," << event << ")" << std::endl; 
            }
            if (foundit == 5 ) break;                             // decay found !
         } 
      }
      if (chic_pdgId && psi_pdgId && foundit==5) break;        // just one decay of this kind is expected
      else chic_pdgId = 0;
   } 
   if (!chic_pdgId)  std::cout << "Rootupler does not found the given decay " << run << "," << event << std::endl;
   else  evt_pass++;
  }

    trigger = getTriggerBits(iEvent);   //grab Trigger informations

    //bool bestCandidateOnly_ = false;
    rf1S_rank = 0;
    photon_flags = 0;
    nTrk_Vtx_Q = 0;
    nTrk_Vtx = 0;
    // grabbing chi inforamtion
    if (chic_cand_handle.isValid() && !chic_cand_handle->empty()) {

       unsigned int csize = chic_cand_handle->size();
       if (bestCandidateOnly_) csize = 1;

       for (unsigned int i = 0; i < csize; i++) {
	   chic_cand = chic_cand_handle->at(i);

	   chic_p4.SetPtEtaPhiM(chic_cand.pt(), chic_cand.eta(), chic_cand.phi(), chic_cand.mass());
	   dimuon_p4.SetPtEtaPhiM(chic_cand.daughter("dimuon")->pt(), chic_cand.daughter("dimuon")->eta(), 
                                  chic_cand.daughter("dimuon")->phi(), chic_cand.daughter("dimuon")->mass());

	   photon_p4.SetPtEtaPhiM(chic_cand.daughter("photon")->pt(), chic_cand.daughter("photon")->eta(), 
                                  chic_cand.daughter("photon")->phi(), chic_cand.daughter("photon")->mass());

	   reco::Candidate::LorentzVector vP = chic_cand.daughter("dimuon")->daughter("muon1")->p4();
	   reco::Candidate::LorentzVector vM = chic_cand.daughter("dimuon")->daughter("muon2")->p4();

           const pat::CompositeCandidate *ThePhoton = (dynamic_cast<const pat::CompositeCandidate *>(chic_cand.daughter("photon")));

           conv_vertex_rho = ThePhoton->userFloat("conv_vertex_rho");
           conv_tkvtx_comp = ThePhoton->userInt("conv_tkvtx_comp");
           conv_comp_ihits = ThePhoton->userInt("conv_comp_ihits");
           conv_high_purity = ThePhoton->userInt("conv_high_purity");
           conv_reject_pi0 = ThePhoton->userInt("conv_reject_pi0");
           conv_algo = ThePhoton->userInt("conv_algo");
           conv_qual_high_purity = ThePhoton->userInt("conv_qual_high_purity");
           conv_qual_generaltks = ThePhoton->userInt("conv_qual_generaltks");

           deltapi0 = ThePhoton->userFloat("deltapi0");
           distMinA = ThePhoton->userFloat("distMinA");
           vtxProb  = ThePhoton->userFloat("vtxProb");

	   if (chic_cand.daughter("dimuon")->daughter("muon1")->charge() < 0) {
	      vP = chic_cand.daughter("dimuon")->daughter("muon2")->p4();
	      vM = chic_cand.daughter("dimuon")->daughter("muon1")->p4();
	   }

	   muonP_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
	   muonM_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());

	   Double_t ele1_pt = (ThePhoton->userData<reco::Track>("track0"))->pt();
	   Double_t ele2_pt = (ThePhoton->userData<reco::Track>("track1"))->pt();

	   if (ele1_pt > ele2_pt) {
	      ele_higherPt_pt = ele1_pt;
	      ele_lowerPt_pt = ele2_pt;
	   } else {
	      ele_higherPt_pt = ele2_pt;
	      ele_lowerPt_pt = ele1_pt;
	   }

           vProb = (dynamic_cast < pat::CompositeCandidate * >(chic_cand.daughter("dimuon")))->userFloat("vProb");
	   ctpv = (dynamic_cast < pat::CompositeCandidate * >(chic_cand.daughter("dimuon")))->userFloat("ppdlPV");
	   ctpv_error = (dynamic_cast < pat::CompositeCandidate * >(chic_cand.daughter("dimuon")))->userFloat("ppdlErrPV");
           ctbs = (dynamic_cast < pat::CompositeCandidate * >(chic_cand.daughter("dimuon")))->userFloat("ppdlBS");
           ctbs_error = (dynamic_cast < pat::CompositeCandidate * >(chic_cand.daughter("dimuon")))->userFloat("ppdlErrBS");
	   photon_flags = (UInt_t) ThePhoton->userInt("flags");
          
           //const reco::Vertex *thePV_closest = ThePhoton->userData<reco::Vertex>("closet_pv_z");
           const reco::Vertex *thePV_dz1 = ThePhoton->userData<reco::Vertex>("closet_pv_d1");
           const reco::Vertex *thePV_dz2 = ThePhoton->userData<reco::Vertex>("closet_pv_d2");

           z_pv = ThePhoton->userFloat("pv_z"); //zOfPrimaryVertexFromTracks(thePV_closest->position());
           const reco::Track *theTk0 = ThePhoton->userData<reco::Track>("track0");
           const reco::Track *theTk1 = ThePhoton->userData<reco::Track>("track1");

           dz_pv_11 = theTk0->dz(thePV_dz1->position());
           dz_pv_21 = theTk0->dz(thePV_dz2->position());
           dz_pv_12 = theTk1->dz(thePV_dz1->position());
           dz_pv_22 = theTk1->dz(thePV_dz2->position());

           double dzError0_ = theTk0->dzError();
           dzError0_ = sqrt(dzError0_*dzError0_+ thePV_dz1->covariance(2,2));
           sigma_tk1_vtx = fabs(dz_pv_11)/dzError0_;

           double dzError1_ = theTk1->dzError();
           dzError1_ = sqrt(dzError1_*dzError1_+ thePV_dz1->covariance(2,2));
           sigma_tk2_vtx = fabs(dz_pv_12)/dzError1_;
           
           tk1_chi2 = theTk0->normalizedChi2();
           tk2_chi2 = theTk1->normalizedChi2();

           tk1_ndof = theTk0->ndof();
           tk2_ndof = theTk1->ndof();
 
           dPhiTvx  = ThePhoton->userFloat("dPhiTvx");
           vtxChi2  = ThePhoton->userFloat("vtxChi2");
           vtxNDoF  = ThePhoton->userInt("vtxNDoF");
           nTracks  = ThePhoton->userInt("nTracks");

	   conv_vertex = ThePhoton->vertex().rho();
	   dz = chic_cand.userFloat("dz");

           UInt_t nTrk_Vtx_Q_tmp = 0;
           UInt_t nTrk_Vtx_tmp = 0;
           const reco::Vertex *thePrimaryV =  (dynamic_cast < pat::CompositeCandidate * >(chic_cand.daughter("dimuon")))->userData<reco::Vertex>("PVwithmuons");
           pvtx_s.SetXYZ(thePrimaryV->position().x(), thePrimaryV->position().y(), thePrimaryV->position().z());
           pvtx_b.SetXYZ(theBeamSpotV.position().x(), theBeamSpotV.position().y(), theBeamSpotV.position().z());
           const reco::Vertex *thePsiV =  (dynamic_cast < pat::CompositeCandidate * >(chic_cand.daughter("dimuon")))->userData<reco::Vertex>("commonVertex");
           psivtx.SetXYZ(thePsiV->position().x(), thePsiV->position().y(), thePsiV->position().z());

           std::vector<reco::TrackBaseRef>::const_iterator itPVtrack = thePrimaryV->tracks_begin();
           for (; itPVtrack != thePrimaryV->tracks_end(); ++itPVtrack) {
               const reco::Track &track = **itPVtrack;
               if (track.pt() < 0.4) continue;
               if (fabs(track.eta())>2.4) continue;
               nTrk_Vtx_Q_tmp++;
               if (track.quality(reco::TrackBase::highPurity)) nTrk_Vtx_tmp++;
           }
           nTrk_Vtx_Q = nTrk_Vtx_Q_tmp;
           nTrk_Vtx   = nTrk_Vtx_tmp;

	   // 2012 parameterization
	   double sigma = Y_sig_par_A + Y_sig_par_B * pow(fabs(dimuon_p4.Rapidity()), 2) + 
                          Y_sig_par_C * pow(fabs(dimuon_p4.Rapidity()), 3);

	   y1S_nsigma = fabs(dimuon_p4.M() - y1SMass) / sigma;
           double QValue = chic_p4.M() - dimuon_p4.M();

           invm1S = QValue + y1SMass;

           unsigned int j = csize + 99;
           if (refit1S_handle.isValid()  && !refit1S_handle->empty()) { 
              for (unsigned int k = 0; k < numKVFChis; k++) {
                  if ((refit1S_handle->at(k)).userInt("Index") == int(i)) {
                     j = k;
                     break;
                  }
              }
           }
	   if ( j <= numKVFChis ) {
	      refit1S = refit1S_handle->at(j);
	      rf1S_chic_p4.SetPtEtaPhiM(refit1S.pt(), refit1S.eta(), refit1S.phi(), refit1S.mass());
	      rf1S_vProb = refit1S.userFloat("vProb");
              rf1S_ctpv = refit1S.userFloat("ctauPV");
              rf1S_ctpv_error = refit1S.userFloat("ctauErrPV");
              rf1S_ctbs = refit1S.userFloat("ctauBS");
              rf1S_ctbs_error = refit1S.userFloat("ctauErrBS");
              const reco::Vertex *theChiV = refit1S.userData<reco::Vertex>("commonVertex");
              chivtx.SetXYZ(theChiV->position().x(), theChiV->position().y(), theChiV->position().z());
	   } else {
	      rf1S_chic_p4.SetPtEtaPhiM(chic_cand.pt(), chic_cand.eta(), chic_cand.phi(), invm1S);
	      rf1S_vProb = -2.;
	   }	// if rf1S is valid
	   chic_tree->Fill();
           rf1S_rank++;
	}		// for i on chic_cand_handle
    } //else std::cout << "no valid chi handle" << std::endl;
    
    mumu_rank = 0;
    if (psi_hand.isValid() && !psi_hand->empty()) {
      for (unsigned int i=0; i< psi_hand->size(); i++) {
        pat::CompositeCandidate psi_ = psi_hand->at(i);
        mumu_p4.SetPtEtaPhiM(psi_.pt(), psi_.eta(), psi_.phi(), psi_.mass());
	ctpv = psi_.userFloat("ppdlPV");
	ctpv_error = psi_.userFloat("ppdlErrPV");
        ctbs = psi_.userFloat("ppdlBS");
        ctbs_error = psi_.userFloat("ppdlErrBS");
        vProb = psi_.userFloat("vProb");

        reco::Candidate::LorentzVector vP = psi_.daughter("muon1")->p4();
        reco::Candidate::LorentzVector vM = psi_.daughter("muon2")->p4();
        if (psi_.daughter("muon1")->charge() < 0) {
           vP = psi_.daughter("muon2")->p4();
           vM = psi_.daughter("muon1")->p4();
        }

        UInt_t nTrk_Vtx_Q_tmp = 0;
        UInt_t nTrk_Vtx_tmp = 0;
        const reco::Vertex *thePrimaryV =  psi_.userData<reco::Vertex>("PVwithmuons");
        pvtx_s.SetXYZ(thePrimaryV->position().x(), thePrimaryV->position().y(), thePrimaryV->position().z());
        pvtx_b.SetXYZ(theBeamSpotV.position().x(), theBeamSpotV.position().y(), theBeamSpotV.position().z());
        const reco::Vertex *thePsiV =  psi_.userData<reco::Vertex>("commonVertex");
        psivtx.SetXYZ(thePsiV->position().x(), thePsiV->position().y(), thePsiV->position().z());

        std::vector<reco::TrackBaseRef>::const_iterator itPVtrack = thePrimaryV->tracks_begin();
        for (; itPVtrack != thePrimaryV->tracks_end(); ++itPVtrack) {
            const reco::Track &track = **itPVtrack;
            if (track.pt() < 0.4) continue;
            if (fabs(track.eta())>2.4) continue;
            nTrk_Vtx_Q_tmp++;
            if (track.quality(reco::TrackBase::highPurity)) nTrk_Vtx_tmp++;
        }
        nTrk_Vtx_Q = nTrk_Vtx_Q_tmp;
        nTrk_Vtx = nTrk_Vtx_tmp;

        muP_p4.SetPtEtaPhiM(vP.pt(), vP.eta(), vP.phi(), vP.mass());
        muM_p4.SetPtEtaPhiM(vM.pt(), vM.eta(), vM.phi(), vM.mass());
        psi_tree->Fill();
        mumu_rank++;
	break;  // just one combination per event, the bestonly
      }
    } 
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void HIchicRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HIchicRootupler);
