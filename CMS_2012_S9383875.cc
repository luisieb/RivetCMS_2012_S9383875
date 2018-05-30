// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Tools/Logging.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/HeavyHadrons.hh"
#include "Rivet/Projections/HadronicFinalState.hh"

/// @todo Include more projections as required, e.g. ChargedFinalState, FastJets, ZFinder...

namespace Rivet {
   
      
   class CMS_2012_S9383875 : public Analysis {
      public:
      /// @name Constructors etc.
                   //@{
                   
       /// Constructor
        CMS_2012_S9383875() : Analysis("CMS_2012_S9383875"){}
               
        double Bjet18 ;
        double Bjet32 ;
        double Bjetmuon ;
                
                   //@}
                   
                   
      public:
                
       /// @name Analysis methods
       //@{
                  
       /// Book histograms and initialise projections before the run
      
       void init() {
                      
        /// @todo Initialise and register projections here
            const FinalState cnfs(-4, 4);
            addProjection(cnfs, "FS");
            addProjection(FastJets(cnfs, FastJets::ANTIKT, 0.5), "Jets");

            IdentifiedFinalState mufs(cnfs);
            mufs.acceptId(PID::MUON);
            mufs.acceptId(PID::ANTIMUON);
            addProjection(mufs, "Muons");
                      
            //b-Hadrons
            addProjection(HeavyHadrons(Cuts::abseta < 2.5 && Cuts::pT > 5*GeV), "BHadrons");

            /// @todo Book histograms here, e.g.:
            _h_dsigdpty05 = bookHisto1D(4, 1, 1);
            _h_dsigdpty10 = bookHisto1D(5, 1, 1);
            _h_dsigdpty15 = bookHisto1D(6, 1, 1);
            _h_dsigdpty20 = bookHisto1D(7, 1, 1);
            _h_dsigdpty22 = bookHisto1D(8, 1, 1);
            _h_dsigdpt    = bookHisto1D(9, 1, 1);
            _h_dsigdy     = bookHisto1D(11, 1, 1);
		      		
            _h_p1d        = bookProfile1D("NBHad_vs_b-jet-pT",60,15,2000);	   

            vector<double> Ptbining {15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507, 548, 592, 638, 686, 737, 790, 846, 905, 967,1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000 }; 

//,2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832,4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};

            vector<double> NBHadBinning;
                     
            for (int i = 0; i < 11; ++i)
                NBHadBinning.push_back(i-0.5);

            pt_N = bookHisto2D("BHad_mult_vs_b-jet-pT", Ptbining, NBHadBinning,"BHad_mult_vs_b-jet-pT","pt","N");
            ptbj_ptBHad = bookHisto2D("pT_b-jet_vs_pT-BHad",Ptbinig, Ptbning,"pT_b-jet_vs_pT-BHad","b-jet_pT","BHad_pT");          
            _h_dsigdphi = bookHisto1D("dsigmadphi_jets",50,0,M_PI,"dsigma/dphi for jets", "phi","dsigma/dphi");
		      _h_dsigdphi_bjets = bookHisto1D("dsigmadphi_bjets",40,0,M_PI,"dsigma/dphi for b-jets", "phi","dsigma/dphi");                        
                      
            Bjet18 = 0.;
            Bjet32 = 0.;
            Bjetmuon = 0.;
                     
            }
                
                 
      /// Perform the per-event analysis
      void analyze(const Event& event) {

            const double weight = event.weight();
                      
             /// @todo Do the event by event analysis here
            
             const FastJets& fastjets = applyProjection<FastJets>(event, "Jets");
             const Jets jets = fastjets.jetsByPt(10.);
                      
             const FinalState& muons = applyProjection<FinalState>(event, "Muons");
             const Particles& bHadrons = applyProjection<HeavyHadrons>(event, "BHadrons").bHadrons();
                     

             bool muongoodevent=false;
                      
             bool onebtag=false;
            
             //For b-jets 
             double deltaphi = 0.; 
             double leadj_phi = 0.;
	          double subleadj_phi = 0.;
             
             //For jets
             double deltaphi_jets = 0.;
             double leadj_phi_jet = 0.;
             double subleadj_phi_jet = 0.;

             unsigned int num_jets = 0;            
             unsigned int num_b_jets = 0;                       
           
             foreach (const Jet& j, jets) {
                               
                  bool btag=false;
                  int BHad =0;
			    	   ++num_jets;
                  
                  if (num_jets < 3){

                        if (num_jets == 1)
                        leadj_phi_jet = j.momentum().phi();
                        else if(num_jets == 2)
                        subleadj_phi_jet = j.momentum().phi();

                        deltaphi_jets = deltaPhi(leadj_phi_jet,subleadj_phi_jet);
                       // cout<<"********************* num_jets ="<<num_jets<<" ******************************"<<endl;
                       // cout<<"deltaphi_jets= "<<deltaphi_jets<<endl;
                  }

                  const double pT_BHad; 
                  bool BHad_lead = true;
                 
                  foreach(const Particle& b, bHadrons) { 

                         double difference = deltaR(j.momentum(), b.momentum());
                                 
                          if(  difference <0.4) {
                               btag=true;
                               onebtag=true;  
				                   ++BHad;
                               if(BHad_lead){
                                  pT_BHad = b.momentum().pT();
                                  BHad_lead = false; 
			          cout << "*****************  PT_BHad = "<< pT_BHad <<"***********************"<<endl;		
                               }
                           }
                  }

				
                           
                  if(btag){

                        num_b_jets += 1;
                               
			               if (num_b_jets < 3){  	
                           
                              if (num_b_jets == 1) 
				                  leadj_phi = j.momentum().phi();
                              else if(num_b_jets == 2)
				                  subleadj_phi = j.momentum().phi();
                                  
			                     deltaphi = deltaPhi(leadj_phi,subleadj_phi);
                             // cout<<"*************** "<<"num_b_jets =" <<num_b_jets <<" *****************************"<<endl;
                             // cout <<"deltaphi_bjets ="<< deltaphi<<endl;   
                        }

                        const double ptB = j.momentum().pT();
                        const double yB = j.momentum().rapidity();
                        
                         if( fabs(yB) < 0.5) _h_dsigdpty05->fill( ptB, weight );
                         else if( fabs(yB) > 0.5 && fabs(yB) < 1.0) _h_dsigdpty10->fill( ptB, weight );
                         else if( fabs(yB) > 1.0 && fabs(yB) < 1.5) _h_dsigdpty15->fill( ptB, weight );
                         else if( fabs(yB) > 1.5 && fabs(yB) < 2.0) _h_dsigdpty20->fill( ptB, weight );
                         else if( fabs(yB) > 2.0 && fabs(yB) < 2.2) _h_dsigdpty22->fill( ptB, weight );
                         if ( fabs(yB) < 2.2 && ptB > 18. ) Bjet18 = Bjet18 + weight;
                         if ( fabs(yB) < 2.2 && ptB > 32. ) Bjet32 = Bjet32 + weight;
                         
                        // cout<<"num_b-jets = "<< num_b_jets<<endl;
                        // cout<<"num_b-hadrons ="<< BHad << endl;

                         pt_N->fill(ptB,BHad,weight);
                        _h_p1d -> fill(ptB,BHad, weight);
                         ptbj_ptBHad -> fill(ptB, pT_BHad, weight);            

                   }


                   if(j.momentum().pT()>30 && fabs(j.momentum().rapidity())<2.4){

                        foreach (const Particle& muon, muons.particles()) {
                                  if ( fabs(muon.momentum().eta()) < 2.4 && muon.momentum().pT() > 9 ) muongoodevent=true;
                         }

                        const double ptB= j.momentum().pT();
                        const double yB= j.momentum().rapidity();

                        if(muongoodevent && onebtag){
                              Bjetmuon+=weight;
                              _h_dsigdpt->fill( ptB, weight );
                              _h_dsigdy->fill( fabs(yB), weight );
                         }
                   }
            }
           
            if( num_b_jets > 2) 
            _h_dsigdphi_bjets -> fill(deltaphi,weight);
            if( num_jets > 2)
            _h_dsigdphi -> fill(deltaphi_jets,weight);
      }
                   
       /// Normalise histograms etc., after the run
      void finalize() {

      /// @todo Normalise, scale and otherwise manipulate histograms here

      // scale(_h_YYYY, crossSection()/sumOfWeights());# norm to cross section

            cout << " CMS_2012_S9383875 [pb] gen xsec = " << crossSection() << endl;
            cout << " CMS_2012_S9383875  Bjet > 18 GeV  xsec = " << Bjet18/sumOfWeights()*crossSection() << endl;
            cout << " CMS_2012_S9383875  Bjet > 32 GeV  xsec = " << Bjet32/sumOfWeights()*crossSection() << endl;

            double invlumi = crossSection()/picobarn/sumOfWeights();
                      
            scale(_h_dsigdpty05, invlumi);
            scale(_h_dsigdpty10, invlumi);
            scale(_h_dsigdpty15, invlumi);
            scale(_h_dsigdpty20, invlumi);
            scale(_h_dsigdpty22, invlumi/0.4);

            double invlumiNano = crossSection()/nanobarn/sumOfWeights();

            cout << " CMS_2012_S9383875  Bjet > 30 GeV with muon  xsec = " << Bjetmuon/sumOfWeights()*crossSection()/nanobarn << endl;
           
            scale(_h_dsigdpt, invlumiNano);
            //scale(_h_dsigdy, invlumiNano/0.5);
            scale(_h_dsigdy, invlumiNano);

      }

      //@}

   private:
                // Data members like post-cuts event weight counters go here

   private:

      /// @name Histograms
      //@{

      Histo1DPtr _h_dsigdpty05;
      Histo1DPtr _h_dsigdpty10;
      Histo1DPtr _h_dsigdpty15;
      Histo1DPtr _h_dsigdpty20;
      Histo1DPtr _h_dsigdpty22;
      Histo1DPtr _h_dsigdpt;
      Histo1DPtr _h_dsigdy;
      
      Profile1DPtr _h_p1d;
      Histo2DPtr pt_N;							      
      Histo2DPtr ptbj_ptBHad;
      
      Histo1DPtr _h_dsigdphi_bjets;
      Histo1DPtr _h_dsigdphi;

                //@}
 };

   // The hook for the plugin system
   DECLARE_RIVET_PLUGIN(CMS_2012_S9383875);
   
}
