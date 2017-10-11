#include "PhotonChecks.h"

PhotonChecks aPhotonChecks ;


PhotonChecks::PhotonChecks() : Processor("PhotonChecks") {

  // modify processor description
  _description = "prints photons from reco particles and " ;


  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "Printing" ,
                              "Print certain messages"  ,
                              _printing,
                               (int)1 ) ;
  registerProcessorParameter("PhotonEnergy" ,
				"The photon energy" ,
				PhotonEnergy,
				(double)10.0);
  
   registerProcessorParameter("Xlow" ,
				"histogram lower bound",
				xlow,
				(double)5.0);

  registerProcessorParameter("Xup" ,
				"histogram upper bound",
				xup,
				(double)5.0);
  registerProcessorParameter("Pullup",
				"pull upper bound",
				pullup,
				(double)3.0);
  registerProcessorParameter("Pulldown",
				"pull lower bound",
				pulldown,
				(double)-3.0);

  registerProcessorParameter("Nbins" ,
				"histogram bins",
				nbins,
				(int)100);
  registerProcessorParameter("rootfilename",
				"the output root file name",
				m_rootFile,
				std::string("photons.root") );				

	//MCParticlesofDecay
  registerInputCollection( LCIO::MCPARTICLE,
                           "MCParticleCollection" ,
                           "Name of the MCParticle input collection"  ,
                           _mcParticleCollectionName ,
                           std::string("MCDecayParticles") ) ;
 
 std::string inputParticleCollectionName = "PandoraPFOs";
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     "InputParticleCollectionName" , 
			     "Input Particle Collection Name "  ,
			     _inputParticleCollectionName,
			      inputParticleCollectionName);



}


void PhotonChecks::init() {

  streamlog_out(DEBUG) << "   init called  "
                       << std::endl ;

  // usually a good idea to
  printParameters() ;
  nEvt = 0;



//set up a TH1Ds for plots
	 f = new TFile(m_rootFile.c_str(),"RECREATE");


      //widen photon energy because its only used for bins
      PhotonEnergy = PhotonEnergy + 3.0;	

	double pi = 3.14159;

	//reco
	 hE = new TH1D("hE","Photon Energy;GeV;EventsPerBin",nbins,xlow,xup);
	 hPx = new TH1D("hPx","Photon Px;GeV;EventsPerBin",nbins,-PhotonEnergy,PhotonEnergy);
	 hPy = new TH1D("hPy","Photon Py;GeV;EventsPerBin",nbins,-PhotonEnergy,PhotonEnergy);
	 hPz = new TH1D("hPz","Photon Pz;GeV;EventsPerBin",nbins,-PhotonEnergy,PhotonEnergy);
	 htheta = new TH1D("htheta","Photon theta;rad;EventsPerBin",nbins,0,3.14159);
	 hphi = new TH1D("hphi","Photon phi;rad;EventsPerBin",nbins,-pi,3.14159);

	//mc
	 hEmc = new TH1D("hEmc","True Photon Energy;GeV;EventsPerBin",nbins,xlow,xup);
	 hPxmc = new TH1D("hPxmc","True Photon Px;GeV;EventsPerBin",nbins,-PhotonEnergy,PhotonEnergy);
	 hPymc = new TH1D("hPymc","True Photon Py;GeV;EventsPerBin",nbins,-PhotonEnergy,PhotonEnergy);
	 hPzmc = new TH1D("hPzmc","True Photon Pz;GeV;EventsPerBin",nbins,-PhotonEnergy,PhotonEnergy);
	 hthetamc = new TH1D("hthetamc","True Photon theta;rad;EventsPerBin",nbins,0,3.14159);
	 hphimc = new TH1D("hphimc","True Photon phi;rad;EventsPerBin",nbins,-pi,3.14159);

	//pulls
	 hEpull = new TH1D("hEpull","Photon Energy Pull",nbins,pulldown,pullup);
	// hPxpull = new TH1D("hPxpull","Photon Px Pull",nbins,-3,3);
	// hPypull = new TH1D("hPypull","Photon Py Pull",nbins,-3,3);
	// hPzpull = new TH1D("hPzpull","Photon Pz Pull",nbins,-3,3);
	 hthetapull = new TH1D("hthetapull","Photon theta pull #frac{#theta-#theta_{mc}}{#sigma_{#theta}}",nbins,pulldown,pullup);
	 hphipull = new TH1D("hphipull","Photon phi pull #frac{#phi-#phi_{mc}}{#sigma_{#phi}}",nbins,pulldown,pullup);

	 hnPfos = new TH1D("hnPfos","N Pfos Reconstructed",7,0,7);
	 hnPfosmc = new TH1D("hnPfosmc","N mc pfos",7,0,7);

//         TH2D (const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup)
	h_theta_pull= new TH2D("h_theta_pull", "measured theta and pull value #frac{#theta-#theta_{mc}}{#sigma_{#theta}}",10 , 0,3.141, 10, -6, 6);
        h_thetaerr_pull = new TH2D("h_thetaerr_pull","measured theta err and pull value",10, 0,0.001, 10,-6,6);

	h_deltacostheta_costheta = new TProfile("h_deltacostheta_costheta","cos#theta - cos#theta_{mc} as a function of cos #theta_{mc};cos #theta_{mc}; #Delta cos #theta",30,-1.0,1.0,0.0,1);
	h_deltatheta = new TProfile("h_deltatheta","#theta - #theta_{mc}  as a function of  cos#theta_{mc};cos #theta_{mc};#Delta #theta",30,-1.0,1.0,-0.005,0.005);

	h_deltaphi = new TProfile("h_deltaphi","#phi - #phi_{mc} as a function of #phi_{mc};#phi_{mc};#Delta #phi",30,-pi,pi,-0.005,0.005);

	 h_costheta_costhetamc = new TProfile("h_costheta_costhetamc","cos #theta as a function cos #theta_{mc}",10,-1,1,-1,1);

	h_deltatheta_deltaphi = new TProfile("h_deltatheta_deltaphi","#Delta #phi as a function of #Delta #theta;#Delta #theta;#Delta #phi",30,-0.005,0.005,-0.005,0.005);
}

void PhotonChecks::processRunHeader( LCRunHeader* run) {

}
bool PhotonChecks::FindPfos( LCEvent* evt ) {

  bool tf = false;

  // clear old vector
  _pfovec.clear();
  typedef const std::vector<std::string> StringVec ;
  StringVec* strVec = evt->getCollectionNames() ;
  for(StringVec::const_iterator name=strVec->begin(); name!=strVec->end(); name++){
 // std::cout<<"collection item: "<< *name << " with " << _inputParticleCollectionName<<std::endl;    
    if(*name==_inputParticleCollectionName){ 
//	std::cout<< "found matching collection name for photon" <<std::endl;
      LCCollection* col = evt->getCollection(*name);
      unsigned int nelem = col->getNumberOfElements();
//	std::cout<<" number of elements "<<nelem<<std::endl;
      tf = true;
      for(unsigned int i=0;i<nelem;i++){
	ReconstructedParticle* recoPart = dynamic_cast<ReconstructedParticle*>(col->getElementAt(i));
	_pfovec.push_back(recoPart);
      }
    }
  }

  if(_printing>1)std::cout << "Find PFOs : " << tf << std::endl; 

  return tf;
}

bool PhotonChecks::FindMCParticles( LCEvent* evt ){
          bool tf = false;

  // clear old vector
  _mcpartvec.clear();
  typedef const std::vector<std::string> StringVec ;
  StringVec* strVec = evt->getCollectionNames() ;
  for(StringVec::const_iterator name=strVec->begin(); name!=strVec->end(); name++){
    if(*name==_mcParticleCollectionName){
      LCCollection* col = evt->getCollection(*name);
      unsigned int nelem = col->getNumberOfElements();
      tf = true;
      for(unsigned int i=0;i<nelem;i++){
        MCParticle* mcPart = dynamic_cast<MCParticle*>(col->getElementAt(i));
        _mcpartvec.push_back(mcPart);
      }
    }
  }

  if(_printing>1)std::cout << "Find MCParticles : " << tf << std::endl;

  return tf;
}
double PhotonChecks::getPhi(TLorentzVector v){
//	if(v.Phi() < 0) return (v.Phi() + 2*3.14159);
//	else return (v.Phi());
 return v.Phi();
}

//return the index in pfovec of the particle that matches the mc particle
int PhotonChecks::getpfoindex(){
	if(_mcpartvec.size() > 1) return -1 ;
	if(_mcpartvec.size() == 0) return -1;	
	int closest_match_index=-1;
	double theta_residual=-1;
	double phi_residual=-1;
	double e_residual=-1;
	double tempresidual1=-1;
	double tempresidual2=-1;
        double tempresidual3=-1;
	TLorentzVector mc,rec ;
	mc.SetPxPyPzE(_mcpartvec[0]->getMomentum()[0],_mcpartvec[0]->getMomentum()[1],_mcpartvec[0]->getMomentum()[2],_mcpartvec[0]->getEnergy()); 
	for(int i=0; i<_pfovec.size(); i++){
		//compare angles
		rec.SetPxPyPzE(_pfovec[i]->getMomentum()[0],_pfovec[i]->getMomentum()[1],_pfovec[i]->getMomentum()[2],_pfovec[i]->getEnergy());
		if(closest_match_index==-1){
			closest_match_index = 0;
			theta_residual = abs(rec.Theta() - mc.Theta());
			phi_residual = abs(getPhi(rec) - getPhi(mc));
 			e_residual = abs(rec.E() - mc.E());
		}
		tempresidual1 = abs(rec.Theta() - mc.Theta());
		tempresidual2 = abs(getPhi(rec) - getPhi(mc));
		tempresidual3 = abs(rec.E() - mc.E());
		if((tempresidual1+tempresidual2+tempresidual3) < (theta_residual+phi_residual+e_residual)){
			closest_match_index=i;
			theta_residual = tempresidual1;
			phi_residual = tempresidual2;
			e_residual = tempresidual3;
		}
		
	}
	return closest_match_index;
}
void PhotonChecks::processEvent( LCEvent * evt ) {


  FindMCParticles(evt);
  FindPfos(evt);

  streamlog_out(MESSAGE) << " start processing event " << std::endl;
 TLorentzVector gamma, gammamc;

	if(_pfovec.size() != 1) return;	
	 hnPfos->Fill(_pfovec.size());
	std::cout<<"Pandora Photons:"<<std::endl;
	for(int i=0; i<_pfovec.size(); i++){
		gamma.SetPxPyPzE(_pfovec[i]->getMomentum()[0],_pfovec[i]->getMomentum()[1],_pfovec[i]->getMomentum()[2],_pfovec[i]->getEnergy());
		hE->Fill( gamma.E());
		hPx->Fill(gamma.Px());
		hPy->Fill(gamma.Py());
		hPz->Fill(gamma.Pz());
		htheta->Fill(gamma.Theta());
		hphi->Fill(getPhi(gamma));
		if(_printing>3){
			std::cout<<"Reconstructed Photon "<<i<<std::endl;
			std::cout<<"(px,py,pz,E) "<<gamma.Px()<<" "<<gamma.Py()<<" "<<gamma.Pz()<<" "<<gamma.E()<<std::endl;
			std::cout<<"(E,theta,phi)"<<gamma.E()<<" "<<gamma.Theta()<<" "<<getPhi(gamma)<<std::endl;
			std::cout<<"Mass "<<gamma.M()<<std::endl;
		}
		
	}
	hnPfosmc->Fill(_mcpartvec.size());
	std::cout<<"Filtered MC Photons:"<<std::endl;
	for(int i=0; i<_mcpartvec.size(); i++){
		gammamc.SetPxPyPzE(_mcpartvec[i]->getMomentum()[0],_mcpartvec[i]->getMomentum()[1],_mcpartvec[i]->getMomentum()[2],_mcpartvec[i]->getEnergy());
		hEmc->Fill( gammamc.E());
		hPxmc->Fill(gammamc.Px());
		hPymc->Fill(gammamc.Py());
		hPzmc->Fill(gammamc.Pz());
		hthetamc->Fill(gammamc.Theta());
		hphimc->Fill(getPhi(gammamc));
		if(_printing>3){
			std::cout<<"MC Photon "<<i<<std::endl;
			std::cout<<"(px,py,pz,E) "<<gammamc.Px()<<" "<<gammamc.Py()<<" "<<gammamc.Pz()<<" "<<gammamc.E()<<std::endl;
			std::cout<<"(E,theta,phi)"<<gammamc.E()<<" "<<gammamc.Theta()<<" "<<getPhi(gammamc)<<std::endl;
			std::cout<<"MC Mass "<<gammamc.M()<<std::endl;
		}

	}

	//pulls
	int pfoindex = getpfoindex();
	//errors
	double dE,dtheta,dphi,dEmc,dthetamc,dphimc;
	if(pfoindex != -1){
		gamma.SetPxPyPzE(_pfovec[pfoindex]->getMomentum()[0],_pfovec[pfoindex]->getMomentum()[1],_pfovec[pfoindex]->getMomentum()[2],_pfovec[pfoindex]->getEnergy());
		gammamc.SetPxPyPzE(_mcpartvec[0]->getMomentum()[0],_mcpartvec[0]->getMomentum()[1],_mcpartvec[0]->getMomentum()[2],_mcpartvec[0]->getEnergy());
	//lets assume 1 mm resolution	
	//copy paste stuff dTheta dPhi (this is actually spatial resolution)
//	double _dTheta = 1.0;//1mm res
//	double _dPhi = 1.0;//1mm res dx and dy

		dE = 0.20*sqrt(gamma.E());
		dtheta = 0.001/sqrt(gamma.E());
		dphi = 0.001/sqrt(gamma.E());
								
				

		dEmc = 0.20*sqrt(gammamc.E());
		dthetamc = 0.001/sqrt(gammamc.E());
		dphimc = 0.001/sqrt(gammamc.E());  

	/*	        double phi = gamma.Phi();
        double theta = gamma.Theta();
        double dx = _dTheta/std::sqrt(gamma.E());//temp stuff to see if this mode actually works
        double dy = _dPhi/std::sqrt(gamma.E());
        double r = 2000; // assume resolution in mm (based on scale of detector thingies) say r is approx 2 meters out
        //trying new spherical model based on resoultion from x,y 
        //eventually rename dtheta dphi to dx and dy
        double a1 = cos(phi)*cos(phi)*dx*dx - sin(phi)*sin(phi)*dy*dy;
        double b1 = ( pow(cos(phi),4) - pow(sin(phi),4))*r*r*cos(theta)*cos(theta) ;
        double sigma1 = std::sqrt(a1/b1); // this is for dtheta

        double a2 = sin(phi)*sin(phi)*dx*dx - cos(phi)*cos(phi)*dy*dy;
        double b2 = ( pow(sin(phi),4) - pow(cos(phi),4))*r*r*sin(theta)*sin(theta);
        double sigma2 = std::sqrt(a2/b2) ; //this is for dphi
	
		double dtheta = sigma1;
		double dphi = sigma2;
*/
		
		hEpull->Fill((gamma.E() - gammamc.E())/dE);
		hthetapull->Fill((gamma.Theta() - gammamc.Theta())/dtheta);
		std::cout<<"theta pull "<<(gamma.Theta()- gammamc.Theta())/ dtheta<<std::endl;
		hphipull->Fill((getPhi(gamma) - getPhi(gammamc))/dphi);

		h_theta_pull->Fill(gamma.Theta(),  (gamma.Theta()- gammamc.Theta())/dtheta);
		h_thetaerr_pull->Fill(dtheta, (gamma.Theta()- gammamc.Theta())/dtheta);
		h_deltacostheta_costheta->Fill(gammamc.CosTheta(),gamma.CosTheta()-gammamc.CosTheta() );
		h_deltatheta->Fill(gammamc.CosTheta(), gamma.Theta()-gammamc.Theta());
		h_deltaphi->Fill(gammamc.Phi(), gamma.Phi() - gammamc.Phi());
		std::cout<<"residual "<<gamma.CosTheta()-gammamc.CosTheta()<<std::endl;
		std::cout<<"costhetamc"<<gammamc.CosTheta();

		h_costheta_costhetamc->Fill(gammamc.CosTheta(), gamma.CosTheta());

		h_deltatheta_deltaphi->Fill(gamma.Theta()-gammamc.Theta(), gamma.Phi()-gammamc.Phi());
	}

  nEvt++;

  

 std::cout << "======================================== event " << nEvt << std::endl ;

}


//void PhotonChecks::check( LCEvent * evt ) {

//}


void PhotonChecks::end(){
 f->Write();
}






