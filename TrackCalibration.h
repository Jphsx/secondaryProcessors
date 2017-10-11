#include "marlin/Processor.h"

#include "EVENT/Track.h"

#include "lcio.h"
#include <vector>
#include "IMPL/LCCollectionVec.h"
#include "TFile.h"

#include "IMPL/ParticleIDImpl.h"

#include "IMPL/TrackImpl.h"
using namespace lcio ;
	/** TrackCalibration:<br>
 *
 * (modelled after GammaGammaCandidateTruthFilter processor)
 * 
 * @author Justin Anguiano, University of Kansas
 * 
 */


 class TrackCalibration : public marlin::Processor {

 public:

 virtual marlin::Processor*  newProcessor() { return new TrackCalibration ; }

  TrackCalibration() ;

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

 
  bool FindTracks(LCEvent* evt);
 
  double safeAcos(double x);
 
  private:
  int nEvt;
  
   std::vector<Track*> _trackvec;
  int   _printing;

   double  _D0Errcalibration;
   double  _Z0Errcalibration;
  double _OmeErrcalibration;
  double _PhiErrcalibration;
  double _TanLErrcalibration;
  
// _inputTrackCollectionName 
  std::string _outputTrackCollectionName;
  std::string _inputTrackCollectionName;
  std::string m_rootFile;
};
