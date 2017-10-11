#include "marlin/Processor.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Track.h"
#include "EVENT/MCParticle.h"
#include "lcio.h"
#include <vector>
#include "IMPL/LCCollectionVec.h"
#include "TFile.h"
#include "TTree.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "TLorentzVector.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"

using namespace lcio ;
	/** PhotonChecks:<br>
 *
 * (modelled after GammaGammaCandidateTruthFilter processor)
 * 
 * @author Justin Anguiano, University of Kansas
 * 
 */


 class PhotonChecks : public marlin::Processor {

 public:

 virtual marlin::Processor*  newProcessor() { return new PhotonChecks ; }

  PhotonChecks() ;

  /** Called at the beginning of the job before anything is read.
   *  Use to initialize the proscessor, e.g. book histograms.
   */
  virtual void init() ;
  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;

  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ;


  /** Called after data processing for clean up.
   */
  virtual void end() ;

  bool FindMCParticles( LCEvent* evt);
  bool FindPfos( LCEvent* evt );
  int getpfoindex();
  double getPhi(TLorentzVector v);
  private:
  int nEvt;
  std::vector<MCParticle*> _mcpartvec;
  std::vector<ReconstructedParticle*> _pfovec;
  int   _printing;

  double PhotonEnergy;
  double xlow;
  double xup;
  double pullup;
  double pulldown;
  int nbins;
  
  TFile* f;
  TH1D* hE;
  TH1D* hPx;
  TH1D* hPy;
  TH1D* hPz;
  TH1D* htheta;
  TH1D* hphi;
  TH1D* hEmc;
  TH1D* hPxmc;
  TH1D* hPymc;
  TH1D* hPzmc;
  TH1D* hthetamc;
  TH1D* hphimc;
  TH1D* hEpull;
  TH1D* hthetapull;
  TH1D* hphipull;
  TH1D* hnPfos;
  TH1D* hnPfosmc;
  TH2D*  h_theta_pull;
  TH2D* h_thetaerr_pull;  
  TProfile* h_deltacostheta_costheta;
  TProfile* h_deltatheta;
  TProfile* h_deltaphi;
  TProfile* h_costheta_costhetamc;
  TProfile* h_deltatheta_deltaphi;
  std::string _inputParticleCollectionName;
  std::string _mcParticleCollectionName;
  

  std::string m_rootFile;
};
