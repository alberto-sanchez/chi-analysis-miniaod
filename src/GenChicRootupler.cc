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

class GenChicRootupler:public edm::EDAnalyzer {
      public:
	explicit GenChicRootupler(const edm::ParameterSet &);
	~GenChicRootupler() override {};
	static void fillDescriptions(edm::ConfigurationDescriptions & descriptions);

      private:
	void analyze(const edm::Event &, const edm::EventSetup &) override;
	bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);

	std::string file_name;

	UInt_t    run;
        ULong64_t event;

	TTree *chic_tree;

	Int_t chic_pdgId;
	Int_t psi_pdgId;
	TLorentzVector gen_chic_p4;
	TLorentzVector gen_psi_p4;
        TLorentzVector gen_dimuon_p4;
	TLorentzVector gen_photon_p4;
	TLorentzVector gen_muonP_p4;
	TLorentzVector gen_muonM_p4;

        edm::EDGetTokenT<reco::GenParticleCollection> genCands_;

};

GenChicRootupler::GenChicRootupler(const edm::ParameterSet & iConfig)
{
    edm::Service < TFileService > fs;
    chic_tree = fs->make < TTree > ("genTree", "Tree of mumugamma");

    chic_tree->Branch("run",      &run,      "run/i");
    chic_tree->Branch("event",    &event,    "event/l");

    chic_tree->Branch("chic_pdgId",    &chic_pdgId,       "chic_pdgId/I");
    chic_tree->Branch("psi_pdgId",     &psi_pdgId,        "psi_pdgId/I");
    chic_tree->Branch("gen_chic_p4",   "TLorentzVector",  &gen_chic_p4);
    chic_tree->Branch("gen_psi_p4",    "TLorentzVector",  &gen_psi_p4);
    chic_tree->Branch("gen_dimuon_p4", "TLorentzVector",  &gen_dimuon_p4);
    chic_tree->Branch("gen_photon_p4", "TLorentzVector",  &gen_photon_p4);
    chic_tree->Branch("gen_muonP_p4",  "TLorentzVector",  &gen_muonP_p4);
    chic_tree->Branch("gen_muonM_p4",  "TLorentzVector",  &gen_muonM_p4);

    genCands_ = consumes<reco::GenParticleCollection>((edm::InputTag)"genParticles");

}

//Check recursively if any ancestor of particle is the given one
bool GenChicRootupler::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle ) return true;
   if (particle->numberOfMothers() && isAncestor(ancestor,particle->mother(0))) return true;
   return false;
}

// ------------ method called for each event  ------------
void GenChicRootupler::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup) {

  run       = iEvent.id().run();
  event     = iEvent.id().event();

  edm::Handle<reco::GenParticleCollection> pruned;
  iEvent.getByToken(genCands_,pruned);

  if (pruned.isValid()) {
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
         for (size_t j=0; j<pwave->numberOfDaughters(); j++) {
            const reco::Candidate *dau = pwave->daughter(j);
            if (dau->pdgId() == psi_pdgId && dau->status() == 2) {
               gen_psi_p4.SetPtEtaPhiM(dau->pt(),dau->eta(),dau->phi(),dau->mass());
               uint nmuons = 0;
               for (size_t k=0; k<dau->numberOfDaughters(); k++) {
                  const reco::Candidate *gdau = dau->daughter(k);
                  if (gdau->pdgId() == 13 && gdau->status()==1) {
                     nmuons++;
                     gen_muonM_p4.SetPtEtaPhiM(gdau->pt(),gdau->eta(),gdau->phi(),gdau->mass());
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
      else {
         chic_pdgId = 0;
         psi_pdgId = 0;
      }
   } 
   if (!chic_pdgId)  std::cout << "Rootupler: does not found the given decay " << run << "," << event << std::endl;
   else  chic_tree->Fill();
  } else std::cout << "Rootupler: GenParticleCollection  not valid " << run << "," << event << std::endl;

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void GenChicRootupler::fillDescriptions(edm::ConfigurationDescriptions & descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenChicRootupler);
