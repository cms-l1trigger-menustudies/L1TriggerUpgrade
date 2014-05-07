//
//  Script depends on the following libraries being loaded:
//          L1Ntupe_C.so
//          L1Menu2015_C.so
//          EvaluateL1Menu_C.so
//  (The libraries are loaded with init.C)
//
//  example:
//
//   linux> root init.C
//   root [1] .x runL1Menu2015.C

// Version 6

std::string InputRootFileName="L1UpgradeNtuples_v24_Tauola_TTbar_PU140_list.txt";
  //"L1UpgradeNtuples_v24_SingleElectron_PU140_list.txt";
  //"L1UpgradeNtuples_v24_Tauola_TTbar_PU140_list.txt";
//L1UpgradeNtuples_v24_SingleElectron_PU140_list.txt";

L1UpgradeNtuple *myNtpl;

test(){
  
   gSystem->Exec("date"); //beginning time stamp
	
// Instantiate the L1 Menu Evaluation Code and open the ROOT file
//   myL1Menu = new L1Menu2015(0,targetLumi,NumberOfLumiSections,InstLumi,InputRootFileName,0.,ZeroBiasPrescale,false);
   myNtpl = new L1UpgradeNtuple;
   myNtpl->OpenWithList(InputRootFileName);
	
// Turn off Dumping of Events
//   myL1Menu->SetDumpEvents(0);	
   myNtpl->testPerformance();

   return;

// If threshold plots need to be made, Run the job and create them
   totalRate = 0.;
   if(makeThrPlots>0) {

// Dump Some Events to screen	
//	    myL1Menu->SetDumpEvents(25);
		 
       printf("\n ---> Making Rate vs Threshold Plots \n");
       totalRate = myL1Menu->RunL1Ana(L1MenuFile,lsFileName,jobTag,numEvts,makeThrPlots,selectDataInput);

// Set the file name to correspond to the rates just calculated
       thresholdPlotsFile = "L1RateHist_";
       thresholdPlotsFile += jobTag;
       thresholdPlotsFile += "Thr";
       thresholdPlotsFile += makeThrPlots;
       thresholdPlotsFile += "_rates.root";
       
// Turn off Dumping of Events for Iteration Phase
//       myL1Menu->SetDumpEvents(0);

   }
 
// Instantiate the code that will adjust thresholds
    myEval = new EvaluateL1Menu();

// Load the rate vs threshold plots
    myEval->LoadThresholdPlots(thresholdPlotsFile);
	
// Load the L1 Menu	
    myEval->LoadL1Menu(L1MenuFile);

// Determine the thresholds
    myEval->DetermineThresholds();
		
// Write out a temporary L1 Menu file with new thresholds 
    TString tmpMenu = "L1Menu_Tmp_";
	 tmpMenu += jobTag;
	 tmpMenu += ".txt";  
    myEval->WriteL1Menu(tmpMenu);	

// Turn off making the rate vs threshold plots...already have them	
    makeThrPlots = 0;

// Get total rate with new thresholds
    totalRate = myL1Menu->RunL1Ana(tmpMenu,lsFileName,jobTag,numEvts,makeThrPlots,selectDataInput);

// Save a copy of the first pass of the total rates (before any bandwidth scaling)
    TString fName = "L1Rates_";
	 fName += jobTag;
	 TString fNameSave = fName;
	 fName += ".txt";
	 fNameSave  += "_Default.txt";
	 TString cmd = "mv "; cmd += fName; cmd += " "; cmd += fNameSave;
	 gSystem->Exec(cmd.Data());
	 
// Now iterate to fill the bandwidth if necessary
    int numIter = 0;
//	printf("Starting Interation %f %f %f \n",totalRate,targetTotalRate,fabs(totalRate-targetTotalRate));
	
   while(fabs(totalRate-targetTotalRate)>targetTolerance && numIter<maxIter) {

// Scale the bandwidths to get the total rate near the target
      double scaleFactor = targetTotalRate/totalRate;
      myEval->ScaleBandwidth(scaleFactor);

// Determine new thresholds for the new bandwidth
      myEval->DetermineThresholds();
	
// Write out new menu 
      myEval->WriteL1Menu(tmpMenu);

// Get total rate with new menu
      totalRate = myL1Menu->RunL1Ana(tmpMenu,lsFileName,jobTag,numEvts,makeThrPlots,selectDataInput);
//	   printf("Total Rate Value %f after iteration %i \n",totalRate,numIter);

// Count interations
      numIter++;
   }

// Write out final L1 Menu
   TString finalMenu = "L1Menu_";
   finalMenu += jobTag;
   finalMenu += "_Final.txt";
   myEval->WriteL1Menu(finalMenu);	
		
   printf("Job %s ending...",jobTag.Data());	
   gSystem->Exec("date"); //ending time stamp

   
}
