// system include files
#include <memory>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>

#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TString.h>
#include <TObjString.h>
#include <TCanvas.h>
#include <TVector3.h>
#include <TMatrixD.h>
#include <TRandom.h>
#include <TMath.h>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"

#include "DataFormats/Candidate/interface/Candidate.h"

#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include <Math/Functions.h>
#include <Math/SVector.h>
#include <Math/SMatrix.h>


//
// class decleration
//

#define PI 3.1416
using namespace std;

class D0TreeMatchGen : public edm::EDAnalyzer {
public:
  explicit D0TreeMatchGen(const edm::ParameterSet&);
  ~D0TreeMatchGen();


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  // ----------member data ---------------------------
    
    TTree* D0para;
    TTree* D0para_GEN;
    
    TH1D* hDR_d1;
    TH1D* hDpt_d1;
    TH1D* hDR_d2;
    TH1D* hDpt_d2;
    
    TH1D* hNtrk;
    TH1D* hZvtx;
    TH1D* hXvtx;
    TH1D* hYvtx;
    TH1D* hNd0;
    
    int nMult_ass_good;
    float bestvx;
    float bestvy;
    float bestvz;
    
    double multMax_;
    double multMin_;
    double deltaR_;
    
    float pt;
    float eta;
    float y;
    float mass;
    float dzos1;
    float dzos2;
    float dxyos1;
    float dxyos2;
    float VtxProb;
    int nhit1;
    int nhit2;
    float dlos;
    float dl;
    float dlerror;
    float agl;
    bool trkquality1;
    bool trkquality2;
    float pt1;
    float pt2;
    float ptErr1;
    float ptErr2;
    float p1;
    float p2;
    float eta1;
    float eta2;
    float charge1;
    float charge2;
    float vtxChi2;
    float ndf;
    float agl_abs;
    float agl2D;
    float agl2D_abs;
    float dlos2D;
    float H2dedx1;
    float H2dedx2;
    float T4dedx1;
    float T4dedx2;
    float trackdca;
    float d0dca;
    float trkChi1;
    float trkChi2;
    bool isSwap;
    bool isPrompt;
    int idmom_reco;
    bool matchGEN;    

    float pt_gen;
    float eta_gen;
    int status_gen;
    int nDau_gen;
    int idmom;
    
    vector< vector<double> > *pVect;
    vector<double> *Dvector1;
    vector<double> *Dvector2;
    vector<int> *pVectIDmom;
    
    edm::EDGetTokenT<reco::VertexCollection> tok_offlinePV_;
    edm::EDGetTokenT<reco::TrackCollection> tok_generalTrk_;
    edm::EDGetTokenT<reco::VertexCompositeCandidateCollection> recoVertexCompositeCandidateCollection_Token_;
    edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token1_;
    edm::EDGetTokenT<edm::ValueMap<reco::DeDxData> > Dedx_Token2_;
    
    edm::EDGetTokenT<reco::GenParticleCollection> tok_genParticle_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

D0TreeMatchGen::D0TreeMatchGen(const edm::ParameterSet& iConfig)
{

  //now do what ever initialization is needed
    multMax_ = iConfig.getUntrackedParameter<double>("multMax", 0.0);
    multMin_ = iConfig.getUntrackedParameter<double>("multMin", 999.9);
    deltaR_ = iConfig.getUntrackedParameter<double>("deltaR", 0.5);

    tok_offlinePV_ = consumes<reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices"));
    tok_generalTrk_ = consumes<reco::TrackCollection>(edm::InputTag("generalTracks"));
    recoVertexCompositeCandidateCollection_Token_ = consumes<reco::VertexCompositeCandidateCollection>(edm::InputTag("generalD0CandidatesNew:D0"));
    
    Dedx_Token1_ = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxHarmonic2"));
    Dedx_Token2_ = consumes<edm::ValueMap<reco::DeDxData> >(edm::InputTag("dedxTruncated40"));
    
    tok_genParticle_ = consumes<reco::GenParticleCollection>(edm::InputTag("genParticles"));

}


D0TreeMatchGen::~D0TreeMatchGen()
{
 
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
D0TreeMatchGen::analyze(const edm::Event& iEvent, const edm::EventSetup&
iSetup)
{
    using std::vector;
    using namespace edm;
    
    // select on requirement of valid vertex
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(tok_offlinePV_,vertices);
    bestvz=-999.9; bestvx=-999.9; bestvy=-999.9;
    double bestvzError=-999.9, bestvxError=-999.9, bestvyError=-999.9;
    const reco::Vertex & vtx = (*vertices)[0];
    bestvz = vtx.z(); bestvx = vtx.x(); bestvy = vtx.y();
    bestvzError = vtx.zError(); bestvxError = vtx.xError(); bestvyError = vtx.yError();
    
    //if(bestvz < -15.0 || bestvz>15.0) return;
    
    hZvtx->Fill(bestvz);
    hXvtx->Fill(bestvx);
    hYvtx->Fill(bestvy);

    edm::Handle<reco::VertexCompositeCandidateCollection> v0candidates;
    iEvent.getByToken(recoVertexCompositeCandidateCollection_Token_,v0candidates);
    if(!v0candidates.isValid()) return;
    
    const reco::VertexCompositeCandidateCollection * v0candidates_ks = v0candidates.product();
    
    edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle1;
    iEvent.getByToken(Dedx_Token1_, dEdxHandle1);

    edm::Handle<edm::ValueMap<reco::DeDxData> > dEdxHandle2;
    iEvent.getByToken(Dedx_Token2_, dEdxHandle2);
    
    //gen par
    edm::Handle<reco::GenParticleCollection> genpars;
    iEvent.getByToken(tok_genParticle_,genpars);
    
    pVect = new vector< vector<double>>;
    pVectIDmom = new vector<int>;
    
    int nD=0;
    
    for(unsigned it=0; it<genpars->size(); ++it){
        
        const reco::GenParticle & trk = (*genpars)[it];
        
        int id = trk.pdgId();
        if(fabs(id)!=421) continue; //check is D0
        nD = nD+1;
        
        pt_gen = trk.pt();
        eta_gen = trk.eta();
        status_gen = trk.status();
        nDau_gen = trk.numberOfDaughters();
        idmom = -77;
        
        if(trk.numberOfMothers()!=0)
        {
            const reco::Candidate * mom = trk.mother();
            idmom = mom->pdgId();
            if(mom->numberOfMothers()!=0)
            {
                const reco::Candidate * mom1 = mom->mother();
                if(mom1->pdgId()>100) idmom = mom1->pdgId();
                if(mom1->numberOfMothers()!=0)
                {
                    const reco::Candidate * mom2 = mom1->mother();
                    if(mom2->pdgId()>100) idmom = mom2->pdgId();
                }
            }
        }
        
        D0para_GEN->Fill();
        
        if(trk.numberOfDaughters()!=2) continue; //check 2-pron decay
 
        const reco::Candidate * Dd1 = trk.daughter(0);
        const reco::Candidate * Dd2 = trk.daughter(1);
        
        if(Dd1->charge()==Dd2->charge()) continue; //check opposite charge daughter
        if(Dd1->status()!=1 || Dd2->status()!=1) continue; //check stable daughter
        if(!(fabs(Dd1->pdgId())==321 && fabs(Dd2->pdgId())==211) && !(fabs(Dd2->pdgId())==321 && fabs(Dd1->pdgId())==211)) continue; //check pi-K daughter
        
        Dvector1 = new vector<double>;
        Dvector2 = new vector<double>;
        
        Dvector1->push_back(Dd1->pt());
        Dvector1->push_back(Dd1->eta());
        Dvector1->push_back(Dd1->phi());
        Dvector1->push_back(Dd1->charge());
        Dvector1->push_back(Dd1->mass());
        
        Dvector2->push_back(Dd2->pt());
        Dvector2->push_back(Dd2->eta());
        Dvector2->push_back(Dd2->phi());
        Dvector2->push_back(Dd2->charge());
        Dvector2->push_back(Dd2->mass());
        
        pVect->push_back(*Dvector1);
        pVect->push_back(*Dvector2);
        
        pVectIDmom->push_back(idmom);
        
        delete Dvector1;
        delete Dvector2;
    }

    hNd0->Fill(nD);

    //track selection, Ntrkoffline couting
    edm::Handle<reco::TrackCollection> tracks;
    iEvent.getByToken(tok_generalTrk_, tracks);
    
    nMult_ass_good = 0;
    for(unsigned it=0; it<tracks->size(); ++it){
        
        const reco::Track & trk = (*tracks)[it];
        
        math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzvtx = trk.dz(bestvtx);
        double dxyvtx = trk.dxy(bestvtx);
        double dzerror = sqrt(trk.dzError()*trk.dzError()+bestvzError*bestvzError);
        double dxyerror = sqrt(trk.d0Error()*trk.d0Error()+bestvxError*bestvyError);
        
        if(!trk.quality(reco::TrackBase::highPurity)) continue;
        if(fabs(trk.ptError())/trk.pt()>0.10) continue;
        if(fabs(dzvtx/dzerror) > 3) continue;
        if(fabs(dxyvtx/dxyerror) > 3) continue;
        
        double eta = trk.eta();
        double pt  = trk.pt();
        
        if(fabs(eta)>2.4) continue;
        if(pt<=0.4) continue;
        nMult_ass_good++;
    }
    
    hNtrk->Fill(nMult_ass_good);
    
    for(unsigned it=0; it<v0candidates_ks->size(); ++it){
        
        const reco::VertexCompositeCandidate & trk = (*v0candidates_ks)[it];

        double secvz=-999.9, secvx=-999.9, secvy=-999.9;
        //double secvzError=-999.9, secvxError=-999.9, secvyError=-999.9;
        //pt,mass
        eta = trk.eta();
        //if(eta>2.4 || eta<-2.4) continue;
        y = trk.rapidity();
        pt = trk.pt();
        double px = trk.px();
        double py = trk.py();
        double pz = trk.pz();
        mass = trk.mass();
        
        const reco::Candidate * d1 = trk.daughter(0);
        const reco::Candidate * d2 = trk.daughter(1);
        
        //Gen match
        matchGEN = false;
        int nGenDau = (int)pVect->size();
        isSwap = false;
        isPrompt = true;
        idmom_reco = -77;
        
        for(int i=0;i<nGenDau;i++)
        {
            vector<double> Dvector1 = (*pVect)[i]; //get GEN daugther vector
            if(d1->charge()!=Dvector1.at(3)) continue; //check match charge
            double deltaR = sqrt(pow(d1->eta()-Dvector1.at(1),2)+pow(d1->phi()-Dvector1.at(2),2));
            hDR_d1->Fill(deltaR);
            hDpt_d1->Fill(fabs((d1->pt()-Dvector1.at(0))/d1->pt()));
            
            if(deltaR > deltaR_) continue; //check deltaR matching
            if(fabs((d1->pt()-Dvector1.at(0))/d1->pt()) > 0.5) continue; //check deltaPt matching
            double d1massGEN = Dvector1.at(4);
            double d1mass = d1->mass();
            double d2massGEN=0, d2mass=0;
            //check dau2
            if(i%2==0)
            {
                vector<double> Dvector2 = (*pVect)[i+1]; //get GEN daugther vector for track2
                if(d2->charge()!=Dvector2.at(3)) continue; //check match charge
                double deltaR = sqrt(pow(d2->eta()-Dvector2.at(1),2)+pow(d2->phi()-Dvector2.at(2),2));
                hDR_d2->Fill(deltaR);
                hDpt_d2->Fill(fabs((d2->pt()-Dvector2.at(0))/d2->pt()));
                
                if(deltaR > deltaR_) continue; //check deltaR matching
                if(fabs((d2->pt()-Dvector2.at(0))/d2->pt()) > 0.5) continue; //check deltaPt matching
                d2massGEN = Dvector2.at(4);
                d2mass = d2->mass();
                
                matchGEN = true; //matched gen
            }
            if(i%2==1)
            {
                vector<double> Dvector2 = (*pVect)[i-1]; //get GEN daugther vector for track2
                if(d2->charge()!=Dvector2.at(3)) continue; //check match charge
                double deltaR = sqrt(pow(d2->eta()-Dvector2.at(1),2)+pow(d2->phi()-Dvector2.at(2),2));
                hDR_d2->Fill(deltaR);
                hDpt_d2->Fill(fabs((d2->pt()-Dvector2.at(0))/d2->pt()));

                if(deltaR > deltaR_) continue; //check deltaR matching
                if(fabs((d2->pt()-Dvector2.at(0))/d2->pt()) > 0.5) continue; //check deltaPt matching
                d2massGEN = Dvector2.at(4);
                d2mass = d2->mass();
                
                matchGEN = true; //matched gen
            }
            
            //check swap
            if(abs(d1massGEN - d1mass)>0.01 || abs(d2massGEN - d2mass)>0.01) isSwap = true;
            
            //check prompt & record mom id
            idmom_reco = pVectIDmom->at(i/2);
            if(fabs(idmom_reco)==521 || fabs(idmom_reco)==511 || fabs(idmom_reco)==531 || fabs(idmom_reco)==541) isPrompt = false;
            
        }
        
        //if(!matchGEN) continue;
        
        double pxd1 = d1->px();
        double pyd1 = d1->py();
        double pzd1 = d1->pz();
        //double pd1 = d1->p();
        double pxd2 = d2->px();
        double pyd2 = d2->py();
        double pzd2 = d2->pz();
        //double pd2 = d2->p();
        
        TVector3 dauvec1(pxd1,pyd1,pzd1);
        TVector3 dauvec2(pxd2,pyd2,pzd2);
        
        auto dau1 = d1->get<reco::TrackRef>();
        auto dau2 = d2->get<reco::TrackRef>();
        
        //DCA between track
        edm::ESHandle<TransientTrackBuilder> theB;
        iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB);

        const TransientTrackBuilder *builder_;
        builder_ = theB.product();
        const reco::TransientTrack dt1 = builder_->build(dau1);
        const reco::TransientTrack dt2 = builder_->build(dau2);
        
        if (!dt1.impactPointTSCP().isValid() || !dt2.impactPointTSCP().isValid()) continue;
        FreeTrajectoryState const & posState = dt1.impactPointTSCP().theState();
        FreeTrajectoryState const & negState = dt2.impactPointTSCP().theState();
        ClosestApproachInRPhi cApp;
        cApp.calculate(posState, negState);
        if (!cApp.status()) continue;
        trackdca = std::abs(cApp.distance());
        
        //trk quality
        
        trkquality1 = dau1->quality(reco::TrackBase::highPurity);
        trkquality2 = dau2->quality(reco::TrackBase::highPurity);
        
        if(!dau1->quality(reco::TrackBase::highPurity)) continue;
        if(!dau2->quality(reco::TrackBase::highPurity)) continue;
        
        //trk dEdx
        H2dedx1 = -999.9;
        H2dedx2 = -999.9;

        if(dEdxHandle1.isValid()){
            const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle1.product();
            H2dedx1 = dEdxTrack[dau1].dEdx();
            H2dedx2 = dEdxTrack[dau2].dEdx();
        }
        
        T4dedx1 = -999.9;
        T4dedx2 = -999.9;
        
        if(dEdxHandle2.isValid()){
            const edm::ValueMap<reco::DeDxData> dEdxTrack = *dEdxHandle2.product();
            T4dedx1 = dEdxTrack[dau1].dEdx();
            T4dedx2 = dEdxTrack[dau2].dEdx();
        }

        //track pt
        pt1 = d1->pt();
        pt2 = d2->pt();
        
        //track momentum
        p1 = d1->p();
        p2 = d2->p();
        
        //track eta
        eta1 = d1->eta();
        eta2 = d2->eta();
        
        //track charge
        charge1 = d1->charge();
        charge2 = d2->charge();

        //track Chi2
        trkChi1 = dau1->normalizedChi2();
        trkChi2 = dau2->normalizedChi2();
        
        //track pT error
        ptErr1 = dau1->ptError();
        ptErr2 = dau2->ptError();

        
        //vertexCovariance 00-xError 11-y 22-z
        secvz = trk.vz(); secvx = trk.vx(); secvy = trk.vy();
        //secvzError = sqrt(trk.vertexCovariance(2,2)); secvxError = sqrt(trk.vertexCovariance(0,0)); secvyError = sqrt(trk.vertexCovariance(1,1));
        
        //trkNHits
        nhit1 = dau1->numberOfValidHits();
        nhit2 = dau2->numberOfValidHits();

        //DCA
        math::XYZPoint bestvtx(bestvx,bestvy,bestvz);
        
        double dzbest1 = dau1->dz(bestvtx);
        double dxybest1 = dau1->dxy(bestvtx);
        double dzerror1 = sqrt(dau1->dzError()*dau1->dzError()+bestvzError*bestvzError);
        double dxyerror1 = sqrt(dau1->d0Error()*dau1->d0Error()+bestvxError*bestvyError);
    
        dzos1 = dzbest1/dzerror1;
        dxyos1 = dxybest1/dxyerror1;
        
        double dzbest2 = dau2->dz(bestvtx);
        double dxybest2 = dau2->dxy(bestvtx);
        double dzerror2 = sqrt(dau2->dzError()*dau2->dzError()+bestvzError*bestvzError);
        double dxyerror2 = sqrt(dau2->d0Error()*dau2->d0Error()+bestvxError*bestvyError);
        
        dzos2 = dzbest2/dzerror2;
        dxyos2 = dxybest2/dxyerror2;

        //vtxChi2
        vtxChi2 = trk.vertexChi2();
        ndf = trk.vertexNdof();
        VtxProb = TMath::Prob(vtxChi2,ndf);

        //PAngle
        TVector3 ptosvec(secvx-bestvx,secvy-bestvy,secvz-bestvz);
        TVector3 secvec(px,py,pz);
        
        TVector3 ptosvec2D(secvx-bestvx,secvy-bestvy,0);
        TVector3 secvec2D(px,py,0);
            
        agl = cos(secvec.Angle(ptosvec));
        agl_abs = secvec.Angle(ptosvec);

        agl2D = cos(secvec2D.Angle(ptosvec2D));
        agl2D_abs = secvec2D.Angle(ptosvec2D);

        //Decay length 3D
        typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
        typedef ROOT::Math::SVector<double, 3> SVector3;
        typedef ROOT::Math::SVector<double, 6> SVector6;
        
        SMatrixSym3D totalCov = vtx.covariance() + trk.vertexCovariance();
        SVector3 distanceVector(secvx-bestvx,secvy-bestvy,secvz-bestvz);
        
        dl = ROOT::Math::Mag(distanceVector);
        dlerror = sqrt(ROOT::Math::Similarity(totalCov, distanceVector))/dl;
        
        dlos = dl/dlerror;
        
        //Decay length 2D
        SVector6 v1(vtx.covariance(0,0), vtx.covariance(0,1),vtx.covariance(1,1),0,0,0);
        SVector6 v2(trk.vertexCovariance(0,0), trk.vertexCovariance(0,1),trk.vertexCovariance(1,1),0,0,0);
        
        SMatrixSym3D sv1(v1);
        SMatrixSym3D sv2(v2);
        
        SMatrixSym3D totalCov2D = sv1 + sv2;
        SVector3 distanceVector2D(secvx-bestvx,secvy-bestvy,0);

        double dl2D = ROOT::Math::Mag(distanceVector2D);
        double dl2Derror = sqrt(ROOT::Math::Similarity(totalCov2D, distanceVector2D))/dl2D;
        
        dlos2D = dl2D/dl2Derror;
        
        //D0 DCA
        d0dca = dl*sin(agl_abs);
        
        //Fill
        D0para->Fill();
        
    }
}


// ------------ method called once each job just before starting event
//loop  ------------
void 
D0TreeMatchGen::beginJob()
{
    edm::Service<TFileService> fs;
        
    TH1D::SetDefaultSumw2();
    
    hDR_d1 = fs->make< TH1D>("hDR_d1","hDR_d1",1000,0,10);
    hDpt_d1 = fs->make< TH1D>("hDpt_d1","hDpt_d1",500,0,5);
    
    hDR_d2 = fs->make< TH1D>("hDR_d2","hDR_d2",1000,0,10);
    hDpt_d2 = fs->make< TH1D>("hDpt_d2","hDpt_d2",500,0,5);

    hNtrk = fs->make< TH1D>("Ntrkoffline","Ntrkoffline",500,0,500);
    hZvtx = fs->make< TH1D>("Zvtx","Zvtx",6000,-30,30);
    hXvtx = fs->make< TH1D>("Xvtx","Xvtx",6000,-30,30);
    hYvtx = fs->make< TH1D>("Yvtx","Yvtx",6000,-30,30);
    hNd0 = fs->make< TH1D>("Nd0","Nd0",100,0,100);
    
    D0para = fs->make< TTree>("D0para","D0para");
    
    //Event info
    //D0para->Branch("Ntrkoffline",&nMult_ass_good);
    //D0para->Branch("vtxX",&bestvx);
    //D0para->Branch("vtxY",&bestvy);
    //D0para->Branch("vtxZ",&bestvz);
    
    //SV info
    D0para->Branch("pT",&pt);
    D0para->Branch("eta",&eta);
    D0para->Branch("y",&y);
    D0para->Branch("mass",&mass);
    D0para->Branch("VtxProb",&VtxProb);
    //D0para->Branch("VtxChi2",&vtxChi2);
    //D0para->Branch("VtxNDoF",&ndf);
    //D0para->Branch("3DCosPointingAngle",&agl);
    D0para->Branch("3DPointingAngle",&agl_abs);
    //D0para->Branch("2DCosPointingAngle",&agl2D);
    //D0para->Branch("2DPointingAngle",&agl2D_abs);
    D0para->Branch("3DDecayLengthSignificance",&dlos);
    D0para->Branch("3DDecayLength",&dl);
    D0para->Branch("3DDecayLengthError",&dlerror);
    //D0para->Branch("2DDecayLengthSignificance",&dlos2D);
    //D0para->Branch("TrackDCA",&trackdca);
    D0para->Branch("D0DCA",&d0dca);
    

    //daugther info
    //D0para->Branch("zDCASignificanceDaugther1",&dzos1);
    //D0para->Branch("zDCASignificanceDaugther2",&dzos2);
    //D0para->Branch("xyDCASignificanceDaugther1",&dxyos1);
    //D0para->Branch("xyDCASignificanceDaugther2",&dxyos2);
    D0para->Branch("NHitD1",&nhit1);
    D0para->Branch("NHitD2",&nhit2);
    //D0para->Branch("HighPurityDaugther1",&trkquality1);
    //D0para->Branch("HighPurityDaugther2",&trkquality2);
    D0para->Branch("pTD1",&pt1);
    D0para->Branch("pTD2",&pt2);
    D0para->Branch("pTerrD1",&ptErr1);
    D0para->Branch("pTerrD2",&ptErr2);
    D0para->Branch("pD1",&p1);
    D0para->Branch("pD2",&p2);
    D0para->Branch("EtaD1",&eta1);
    D0para->Branch("EtaD2",&eta2);
    //D0para->Branch("chargeD1",&charge1);
    //D0para->Branch("chargeD2",&charge2);
    //D0para->Branch("dedxHarmonic2D1",&H2dedx1);
    //D0para->Branch("dedxHarmonic2D2",&H2dedx2);
    //D0para->Branch("dedxTruncated40Daugther1",&T4dedx1);
    //D0para->Branch("dedxTruncated40Daugther2",&T4dedx2);
    //D0para->Branch("normalizedChi2Daugther1",&trkChi1);
    //D0para->Branch("normalizedChi2Daugther2",&trkChi2);
    D0para->Branch("isSwap",&isSwap);
    D0para->Branch("isPrompt",&isPrompt);
    D0para->Branch("idmom",&idmom_reco);
    D0para->Branch("matchGEN",&matchGEN);

    D0para_GEN = fs->make< TTree>("D0para_GEN","D0para_GEN");
    
    //GEN info
    D0para_GEN->Branch("pT",&pt_gen);
    D0para_GEN->Branch("eta",&eta_gen);
    D0para_GEN->Branch("status",&status_gen);
    D0para_GEN->Branch("nDau",&nDau_gen);
    D0para_GEN->Branch("MotherID",&idmom);
}

// ------------ method called once each job just after ending the event
//loop  ------------
void 
D0TreeMatchGen::endJob() {
    
}

//define this as a plug-in
DEFINE_FWK_MODULE(D0TreeMatchGen);






