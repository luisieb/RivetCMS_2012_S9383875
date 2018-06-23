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

//double x = xmin;
//vector<double> steps;
//enum Stepmode {LINEAR, LOGARITHMIC};
//switch(logx)
//{
//   case LINEAR:
//      {
//         const double steplength = (xmax-xmin)/double(nsteps);
//         for(int i = 0; i<=nsteps; i++)
//         {
//            steps.push_back(x);
//            x+=steplength;
//         }
//         cout << program_name << "\tDoing " << nsteps << " linear nsteps (step length = " << steplength << ")" << endl;
//      }
//      break;
//   case LOGARITHMIC:
//      {
//         const double steplength = (log10(xmax/xmin))/double(nsteps);
//         for(int i = 0; i<=nsteps; i++)
//         {
//            steps.push_back(x);
//            x*=pow(10,steplength);
//         }
//         cout << program_name << "\tDoing " << nsteps << " logarithmic nsteps (logstep = " << steplength << ")" << endl;
//      }
//      break;
//}
//return steps;

namespace Rivet {
   
      
   class CMS_2012_S9383875 : public Analysis {
      public:
      /// @name Constructors etc.
                   //@{
                   
       /// Constructor
        CMS_2012_S9383875() : Analysis("CMS_2012_S9383875"){}
               
      /*  double Bjet18 ;
        double Bjet32 ;
        double Bjetmuon ;*/
                
                   //@}
                   
                   
      public:
                
       /// @name Analysis methods
       //@{
                  
       /// Book histograms and initialise projections before the run
      
       void init() {
                      
        /// @todo Initialise and register projections here
            const FinalState cnfs(-4, 4);
            addProjection(cnfs, "FS");
            addProjection(FastJets(cnfs, FastJets::ANTIKT, 0.4), "Jets");

            IdentifiedFinalState mufs(cnfs);
            mufs.acceptId(PID::MUON);
            mufs.acceptId(PID::ANTIMUON);
           // addProjection(mufs, "Muons");
                      
            //b-Hadrons
            addProjection(HeavyHadrons(Cuts::abseta < 2.5 && Cuts::pT > 5*GeV), "BHadrons");

            /// @todo Book histograms here, e.g.:
       /*     _h_dsigdpty05 = bookHisto1D(4, 1, 1);
            _h_dsigdpty10 = bookHisto1D(5, 1, 1);
            _h_dsigdpty15 = bookHisto1D(6, 1, 1);
            _h_dsigdpty20 = bookHisto1D(7, 1, 1);
            _h_dsigdpty22 = bookHisto1D(8, 1, 1);
            _h_dsigdpt    = bookHisto1D(9, 1, 1);
            _h_dsigdy     = bookHisto1D(11, 1, 1);
		      		
            _h_p1d        = bookProfile1D("NBHad_vs_b-jet-pT",60,15,2000);*/	   


       /****     vector<double> Ptbining {0,3,6,9,12,15, 18, 21, 24, 28, 32, 37, 43, 49, 56, 64, 74, 84, 97, 114, 133, 153, 174, 196, 220, 245, 272, 300, 330, 362, 395, 430, 468,507, 548, 592, 638, 686, 737, 790, 846, 905, 967,1032, 1101, 1172, 1248, 1327, 1410, 1497, 1588, 1684, 1784, 1890, 2000 }; 

//,2116, 2238, 2366, 2500, 2640, 2787, 2941, 3103, 3273, 3450, 3637, 3832,4037, 4252, 4477, 4713, 4961, 5220, 5492, 5777, 6076, 6389, 6717, 7000};

            vector<double> NBHadBinning;
                     
            for (int i = 0; i < 11; ++i)
                NBHadBinning.push_back(i-0.5);

            pt_N = bookHisto2D("BHad_mult_vs_b-jet-pT", Ptbining, NBHadBinning,"BHad_mult_vs_b-jet-pT","pt","N");
          
            ptbj_ptBHad = bookHisto2D("pT_b-jet_vs_pT-BHad",Ptbining, Ptbining,"pT_b-jet_vs_pT-BHad","b-jet_pT","BHad_pT");          
            ptbj_ptallBHad =  bookHisto2D("pT_b-jet_vs_pT-allBHad",Ptbining, Ptbining,"pT_b-jet_vs_pT-BHad","b-jet_pT","BHad_pT"); */
            const int binning = 20;           

            _h_dsigdphi = bookHisto1D("dsigmadphi_jets",binning,0,M_PI,"dsigma/dphi for jets", "phi","dsigma/dphi");
            _h_dsigdphi_bjets = bookHisto1D("dsigmadphi_bjets",binning,0,M_PI,"dsigma/dphi for b-jets", "phi","dsigma/dphi");                        
           
            // Histograms for deltaphi slices in pT
            //10 < pT < 20             
            _h_dsigdphi_10pT20 = bookHisto1D("dsigmadphi_jets_10pT20",binning,0,M_PI,"dsigma/dphi for jets (10<pT<20)", "phi","dsigma/dphi");
            _h_dsigdphi_bjets_10pT20 = bookHisto1D("dsigmadphi_bjets_10pT20",binning,0,M_PI,"dsigma/dphi for b-jets (10<pT<20)", "phi","dsigma/dphi");                         
            // 20 < pT < 50             
            _h_dsigdphi_20pT50 = bookHisto1D("dsigmadphi_jets_20pT50",binning,0,M_PI,"dsigma/dphi for jets (20<pT<50)", "phi","dsigma/dphi");
	         _h_dsigdphi_bjets_20pT50 = bookHisto1D("dsigmadphi_bjets_20pT50",binning,0,M_PI,"dsigma/dphi for b-jets (20<pT<50)", "phi","dsigma/dphi");                        
            // 50 < pT < 100             
            _h_dsigdphi_50pT100 = bookHisto1D("dsigmadphi_jets_50pT100",binning,0,M_PI,"dsigma/dphi for jets (50<pT<100)", "phi","dsigma/dphi");
            _h_dsigdphi_bjets_50pT100= bookHisto1D("dsigmadphi_bjets_50pT100",binning,0,M_PI,"dsigma/dphi for b-jets (50<pT<100)", "phi","dsigma/dphi");                        
            // 100 < pT < 200            
            _h_dsigdphi_100pT200 = bookHisto1D("dsigmadphi_jets_100pT200",binning,0,M_PI,"dsigma/dphi for jets (100<pT<200)", "phi","dsigma/dphi");
		      _h_dsigdphi_bjets_100pT200= bookHisto1D("dsigmadphi_bjets_100pT200",binning,0,M_PI,"dsigma/dphi for b-jets (100<pT<200)", "phi","dsigma/dphi");
             // 200 < pT < 400            
            _h_dsigdphi_200pT400 = bookHisto1D("dsigmadphi_jets_200pT400",binning,0,M_PI,"dsigma/dphi for jets (200<pT<400)", "phi","dsigma/dphi");
		      _h_dsigdphi_bjets_200pT400= bookHisto1D("dsigmadphi_bjets_200pT400",binning,0,M_PI,"dsigma/dphi for b-jets (200<pT<400)", "phi","dsigma/dphi");
             // 400 < pT < 800            
            _h_dsigdphi_400pT800 = bookHisto1D("dsigmadphi_jets_400pT800",binning,0,M_PI,"dsigma/dphi for jets (400<pT<800)", "phi","dsigma/dphi");
		      _h_dsigdphi_bjets_400pT800= bookHisto1D("dsigmadphi_bjets_400pT800",binning,0,M_PI,"dsigma/dphi for b-jets (400<pT<800)", "phi","dsigma/dphi");
             // 800 < pT < 1200            
            _h_dsigdphi_800pT1200 = bookHisto1D("dsigmadphi_jets_800pT1200",binning,0,M_PI,"dsigma/dphi for jets (800<pT<1200)", "phi","dsigma/dphi");
		      _h_dsigdphi_bjets_800pT1200= bookHisto1D("dsigmadphi_bjets_800pT1200",binning,0,M_PI,"dsigma/dphi for b-jets (800<pT<1200)", "phi","dsigma/dphi");
             // 1200 < pT < 1600            
            _h_dsigdphi_1200pT1600 = bookHisto1D("dsigmadphi_jets_1200pT1600",binning,0,M_PI,"dsigma/dphi for jets (1200<pT<1600)", "phi","dsigma/dphi");
		      _h_dsigdphi_bjets_1200pT1600= bookHisto1D("dsigmadphi_bjets_1200pT1600",binning,0,M_PI,"dsigma/dphi for b-jets (1200<pT<1600)", "phi","dsigma/dphi");
             // 1600 < pT < 2000            
            _h_dsigdphi_1600pT2000 = bookHisto1D("dsigmadphi_jets_1600pT2000",binning,0,M_PI,"dsigma/dphi for jets (1600<pT<2000)", "phi","dsigma/dphi");
	         _h_dsigdphi_bjets_1600pT2000= bookHisto1D("dsigmadphi_bjets_1600pT2000",binning,0,M_PI,"dsigma/dphi for b-jets (1600<pT<2000)", "phi","dsigma/dphi");
            // 1600 < pT < 2000            
            _h_dsigdphi_2000pT = bookHisto1D("dsigmadphi_jets_2000pT",binning,0,M_PI,"dsigma/dphi for jets (pT>2000)", "phi","dsigma/dphi");
	         _h_dsigdphi_bjets_2000pT= bookHisto1D("dsigmadphi_bjets_2000pT",binning,0,M_PI,"dsigma/dphi for b-jets (pT>2000)", "phi","dsigma/dphi");

            
            const int nsteps = 20;
            const double xmin = 0.0001;
            const double xmax = 10.;
            double x = xmin;  
            vector<double> steps;
            const double steplength = (log10(xmax/xmin))/nsteps;
            for(int i = 0; i<=nsteps; i++){
            	steps.push_back(x);
                x*=pow(10,steplength);
            }
            
            //phi* histograms
	       _h_dsigmadphistar_jets = bookHisto1D("phistar_jets",steps,"phistar for jets", "phi*","dsigma/dphi*");
	       _h_dsigmadphistar_bjets = bookHisto1D("phistar_bjets",steps,"phistar for b-jets", "phi*","dsigma/dphi*");	     

          /*  Bjet18 = 0.;
            Bjet32 = 0.;
            Bjetmuon = 0.;*/
                     
            }
                


           //function which gives nB and maximal pT of bHadrons close to jet j
           void getBmaxPtAndNBs(const Jet &j, const Particles& bHadrons, int &nB, double &pT_Max) {
              pT_Max = 0;
              nB = 0;
              foreach(const Particle& b, bHadrons) { 
                 double difference = deltaR(j.momentum(), b.momentum());
                 if(difference < 0.4) {
                    pT_Max = max(pT_Max, b.momentum().pT());
                    ++nB;
                 }
              }
           }

           double getBmaxPt(const Jet &j, const Particles& bHadrons) {
              int nB;
              double pT_Max;
              getBmaxPtAndNBs(j, bHadrons, nB,pT_Max);
              return pT_Max;
           }


           int nBs(const Jet &j, const Particles& bHadrons) {
              int nB;
              double pT_Max;
              getBmaxPtAndNBs(j, bHadrons, nB,pT_Max);
              return nB;
           }




                 
      /// Perform the per-event analysis
      void analyze(const Event& event) {

            const double weight = event.weight();
                      
             /// @todo Do the event by event analysis here
            
             const FastJets& fastjets = applyProjection<FastJets>(event, "Jets");
             const Jets jets = fastjets.jetsByPt(10.);
                      
           //  const FinalState& muons = applyProjection<FinalState>(event, "Muons");
             const Particles& bHadrons = applyProjection<HeavyHadrons>(event, "BHadrons").bHadrons();
                     

            // bool muongoodevent=false;
            // bool onebtag      =false;
            
//             foreach (const Jet& j, jets) {

                  //number of Bhadrons close to jet j
           //       int BHad = nBs(j, bHadrons); 
                  
                  //filling the BHad_pT of all BHad vs b-jet pT
          /*       foreach(const Particle& b, bHadrons) { 
                     double difference = deltaR(j.momentum(), b.momentum());
                     if(difference < 0.4) {
                        ptbj_ptallBHad -> fill(j.momentum().pT(),b.momentum().pT(),weight);       
                        
                       //cout <<"pT_jet = " << j.momentum().pT() <<endl;
                       //cout <<"pT_BHad = "<< b.momentum().pT() << endl;

                     }
                 }*/

            /*      if(BHad > 0) { //is bJet
                        onebtag = true;
                        const double ptB = j.momentum().pT();
                        const double yBa = abs(j.momentum().rapidity());
                        
                         if( yBa < 0.5)                   _h_dsigdpty05->fill( ptB, weight );
                         else if( yBa > 0.5 && yBa < 1.0) _h_dsigdpty10->fill( ptB, weight );
                         else if( yBa > 1.0 && yBa < 1.5) _h_dsigdpty15->fill( ptB, weight );
                         else if( yBa > 1.5 && yBa < 2.0) _h_dsigdpty20->fill( ptB, weight );
                         else if( yBa > 2.0 && yBa < 2.2) _h_dsigdpty22->fill( ptB, weight );

                         if ( yBa < 2.2 && ptB > 18. ) Bjet18 +=  weight;
                         if ( yBa < 2.2 && ptB > 32. ) Bjet32 +=  weight;

                         double pT_BHad_Max = getBmaxPt(j, bHadrons);
                        //cout << "pT_BHad_Max = "<< pT_BHad_Max <<endl;  

                         pt_N ->fill(ptB, BHad,weight);
                        _h_p1d->fill(ptB, BHad, weight);
                         ptbj_ptBHad -> fill(ptB, pT_BHad_Max, weight);            
                        
                   }*/


                  //Analysis with muons
               /*    if(j.momentum().pT()>30 && fabs(j.momentum().rapidity())<2.4){

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
                   }*/
//            }
           
            //Delta Phi for jets above 100 GeV   
            if(jets.size() >= 2 && jets[1].pt() > 10*GeV) {
               const double deltaphi = deltaPhi(jets[0].phi(),jets[1].phi());
	            const double pt = jets[0].pt();

      	       // Compute phi*
               const double phi_acop = M_PI - deltaphi;
               const double costhetastar = tanh( 0.5 * (jets[0].eta() - jets[1].eta()) );
               const double sin2thetastar = (costhetastar > 1) ? 0.0 : (1.0 - sqr(costhetastar));
               const double phistar = tan(0.5 * phi_acop) * sqrt(sin2thetastar);

               _h_dsigdphi -> fill(deltaphi, weight);
               _h_dsigmadphistar_jets -> fill(phistar,weight);

               if( pt > 10 && pt < 20 )
                  _h_dsigdphi_10pT20 -> fill(deltaphi, weight);      
               else if( pt > 20 && pt < 50 )
                  _h_dsigdphi_20pT50 -> fill(deltaphi, weight);
               else if( pt > 50 && pt < 100 )
                  _h_dsigdphi_50pT100 -> fill(deltaphi, weight); 
               else if( pt > 100 && pt < 200 )
                  _h_dsigdphi_100pT200 -> fill(deltaphi, weight);      
               else if( pt > 200 && pt < 400 )
                  _h_dsigdphi_200pT400 -> fill(deltaphi, weight);     
               else if( pt > 400 && pt < 800 )
                  _h_dsigdphi_400pT800 -> fill(deltaphi, weight);     
               else if( pt > 800 && pt < 1200 )
                  _h_dsigdphi_800pT1200 -> fill(deltaphi, weight);     
               else if( pt > 1200 && pt < 1600 )
                  _h_dsigdphi_1200pT1600 -> fill(deltaphi, weight);
               else if( pt > 1600 && pt < 2000 )
                  _h_dsigdphi_1600pT2000 -> fill(deltaphi, weight);
	            else if( pt > 2000 )
                  _h_dsigdphi_2000pT -> fill(deltaphi, weight);	
   
              
           
               if(nBs(jets[0], bHadrons) && nBs(jets[1], bHadrons)) {
                  
		             _h_dsigdphi_bjets -> fill(deltaphi, weight);
		             _h_dsigmadphistar_bjets -> fill(phistar,weight);	
              
                  if( pt > 10 && pt < 20 )
                      _h_dsigdphi_bjets_10pT20 -> fill(deltaphi, weight);      
                  else if( pt > 20 && pt < 50 )
                      _h_dsigdphi_bjets_20pT50 -> fill(deltaphi, weight);
                  else if( pt > 50 && pt < 100 )
                      _h_dsigdphi_bjets_50pT100 -> fill(deltaphi, weight); 
                  else if( pt > 100 && pt < 200 )
                      _h_dsigdphi_bjets_100pT200 -> fill(deltaphi, weight);      
                  else if( pt > 200 && pt < 400 )
                      _h_dsigdphi_bjets_200pT400 -> fill(deltaphi, weight);     
                  else if( pt > 400 && pt < 800 )
                      _h_dsigdphi_bjets_400pT800 -> fill(deltaphi, weight);     
                  else if( pt > 800 && pt < 1200 )
                      _h_dsigdphi_bjets_800pT1200 -> fill(deltaphi, weight);     
                  else if( pt > 1200 && pt < 1600 )
                      _h_dsigdphi_bjets_1200pT1600 -> fill(deltaphi, weight);
                  else if( pt > 1600 && pt < 2000 )
                      _h_dsigdphi_bjets_1600pT2000 -> fill(deltaphi, weight);
                  else if( pt > 2000 )
                      _h_dsigdphi_bjets_2000pT -> fill(deltaphi, weight);	
                 
               }
            }

      }
                   
       /// Normalise histograms etc., after the run
      void finalize() {

      /// @todo Normalise, scale and otherwise manipulate histograms here

      // scale(_h_YYYY, crossSection()/sumOfWeights());# norm to cross section

       /*     cout << " CMS_2012_S9383875 [pb] gen xsec = " << crossSection() << endl;
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
            scale(_h_dsigdy, invlumiNano);*/
          
      }

      //@}

   private:
                // Data members like post-cuts event weight counters go here

   private:

      /// @name Histograms
      //@{

     /* Histo1DPtr _h_dsigdpty05;
      Histo1DPtr _h_dsigdpty10;
      Histo1DPtr _h_dsigdpty15;
      Histo1DPtr _h_dsigdpty20;
      Histo1DPtr _h_dsigdpty22;
      Histo1DPtr _h_dsigdpt;
      Histo1DPtr _h_dsigdy;
      
      Profile1DPtr _h_p1d;
      Histo2DPtr pt_N;							      
      Histo2DPtr ptbj_ptBHad;
      Histo2DPtr ptbj_ptallBHad; */

       Histo1DPtr _h_dsigdphi_bjets;
       Histo1DPtr _h_dsigdphi;

       // 10 < pT < 20         
       Histo1DPtr _h_dsigdphi_10pT20;
       Histo1DPtr _h_dsigdphi_bjets_10pT20;
       // 20 < pT < 50         
       Histo1DPtr _h_dsigdphi_20pT50;
       Histo1DPtr _h_dsigdphi_bjets_20pT50; 
       // 50 < pT < 100         
       Histo1DPtr _h_dsigdphi_50pT100;
		 Histo1DPtr _h_dsigdphi_bjets_50pT100; 
       // 100 < pT < 200         
       Histo1DPtr _h_dsigdphi_100pT200;
		 Histo1DPtr _h_dsigdphi_bjets_100pT200;
       // 200 < pT < 400            
       Histo1DPtr _h_dsigdphi_200pT400;
		 Histo1DPtr _h_dsigdphi_bjets_200pT400;
       // 400 < pT < 800            
       Histo1DPtr _h_dsigdphi_400pT800;
		 Histo1DPtr _h_dsigdphi_bjets_400pT800;
       // 800 < pT < 1200            
       Histo1DPtr _h_dsigdphi_800pT1200;
	    Histo1DPtr _h_dsigdphi_bjets_800pT1200;
       // 1200 < pT < 1600            
       Histo1DPtr _h_dsigdphi_1200pT1600;
       Histo1DPtr _h_dsigdphi_bjets_1200pT1600;
       // 1600 < pT < 2000            
       Histo1DPtr _h_dsigdphi_1600pT2000;
       Histo1DPtr _h_dsigdphi_bjets_1600pT2000;
       // pT > 2000            
       Histo1DPtr _h_dsigdphi_2000pT;
       Histo1DPtr _h_dsigdphi_bjets_2000pT;
       
       //phistar histograms
       Histo1DPtr _h_dsigmadphistar_jets;
       Histo1DPtr _h_dsigmadphistar_bjets;
             

                //@}
 };

   // The hook for the plugin system
   DECLARE_RIVET_PLUGIN(CMS_2012_S9383875);
   
}
