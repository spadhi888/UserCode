//-----------------------------------------------
// Takes the text files and dumps the contents 
// in an ntuple.
// It's OK if a text file does not exist.
// However, the text files need to have the 
// same order of gluino-mass LSP-mass etc...
// Meaning that the lines must correspond
// If they do not, the program will barf
//
// It is advisable to "sort" the textfiles ahead 
// of time, eg,
// sort T1tttt_landsOUT_SR1.txt > SR1_sorted.txt
// sort T1tttt_landsOUT_SR2.txt > SR2_sorted.txt
// sort T1tttt_landsOUT_SR3.txt > SR3_sorted.txt
// sort T1tttt_landsOUT_SR4.txt > SR4_sorted.txt
// sort T1tttt_landsOUT_SR5.txt > SR5_sorted.txt
// sort T1tttt_landsOUT_SR6.txt > SR6_sorted.txt
// sort T1tttt_landsOUT_SR7.txt > SR7_sorted.txt
//
// Also reads in the gluino xsection txt file.  
// Assumes that the gluino mass increases monotonycally
// from oe line to the next
//
//
//  Usage: 
//  root> .L makeNtuple()
//  root> makeNtuple()
//------------------------------------------- 

//
//

void makeNtuple(){

  //-----------------------------------------------
  // The name of the file with the output ntuple
  //-----------------------------------------------
  char* ntFile = "ntuple.root";


  //-----------------------------------------------
  // Raed in the gluino pair production xsection  
  //-----------------------------------------------
  float xsec_mass[2000];
  float xsec[2000];
  float xsecup[2000];
  float xsecdwn[2000];
  int nxsec=0;
  ifstream fromfile("../../xsec/gluino_xsec.txt");
  if ( !fromfile.good() ) {
    cout << "gluino xsec file does not exist" << endl ;
    return;
  }
  cout << "Reading gluino xsection" << endl;
  int xs=0;
  while (!fromfile.eof()) {
    TString d1, d2, d3, d4;
    fromfile >> d1 >> xsec_mass[xs] >> d2 >> xsec[xs] >> d3 
	     >> xsecup[xs] >> d4 >> xsecdwn[xs];
    xs++;
    nxsec++;
  }
  fromfile.close();
  nxsec = nxsec-1;  //last line is junk

  for (int i=0; i<nxsec; i++) {
    cout<<xsec_mass[i]<<" "<<xsec[i]<<" "<<xsecup[i]<<" "<< xsecdwn[i]<< endl;
  }

  //--------------------------------------------------------
  // Check that the gluino mass in the xsection file increases monotonically
  //---------------------------------------------------------
  float previous = xsec_mass[0];
  for (int i=1; i<nxsec; i++) {
    if (xsec_mass[i] < previous) {
      cout << "Mass in xsection file does not increase monotonically..Abort"<<endl;
      return;
    }
    previous = xsec_mass[i];
  }

  //-----------------------------------------------
  // Define the variables in the txt files
  //-----------------------------------------------
  float glmass[7][1000];
  float lspmass[7][1000];
  float eff[7][1000];
  float err[7][1000];
  float obslim[7][1000];
  float explim[7][1000];
  bool  usedflag[7][1000];  // a flag to be used later
  int nentries[7];          // number of entries in the file
  
  //-----------------------------------------------
  // The file names, ordered by signal region...
  //-----------------------------------------------
  vector<TString> filenames;
  filenames.push_back("../txt/t2ttww_sr1.txt");
  filenames.push_back("../txt/t2ttww_sr2.txt");
  filenames.push_back("../txt/t2ttww_sr3.txt");
  filenames.push_back("../txt/t2ttww_sr4.txt");
  filenames.push_back("../txt/t2ttww_sr5.txt");
  filenames.push_back("../txt/t2ttww_sr6.txt");
  filenames.push_back("../txt/t2ttww_sr7.txt");

  //-----------------------------------------------
  // Loop over signal regions and load the arrays from the txt files
  //-----------------------------------------------
  for (int ireg=0; ireg<7; ireg++) {
    nentries[ireg] = 0;
    ifstream fromfile(filenames.at(ireg));

    // if the file does not exist, go to the next one
    if ( !fromfile.good() ) {
      cout << "file " << filenames.at(ireg) << " does not exist" << endl ;
      continue;
    }
    // the file exist, read from it
    cout << "Reading from " << filenames.at(ireg) << endl;
    int j=0;
    while (!fromfile.eof()) {
      fromfile >> glmass[ireg][j]
	       >> lspmass[ireg][j]
	       >> eff[ireg][j]
	       >> err[ireg][j]
	       >> obslim[ireg][j]
	       >> explim[ireg][j];      
      j++;
      nentries[ireg]++;
    }
    fromfile.close();
    // the last entry is junk
    nentries[ireg] = nentries[ireg]-1;
  }
  //-------------------------------------------
  // Now check that the line order is OK
  // And that the number of lines is consistent
  //--------------------------------------------
  bool bad=false;
  int numbLines = 0;
  for (int ii=0; ii<7; ii++) {
    for (int jj=ii+1; jj<7; jj++) {
      if (nentries[ii]*nentries[jj] !=0 ) {
	numbLines = nentries[ii];
	if (nentries[ii] != nentries[jj]) {
	  cout << "Mismatched lines: file " << filenames.at(ii);
          cout << " has " << nentries[ii] << " lines" << endl;
	  cout << "                  file " << filenames.at(jj);
          cout << " has " << nentries[jj] << " lines" << endl;
	  bad=true;
	}
	for (int line=0; line<nentries[ii]; line++) {
	  if (glmass[ii][line] != glmass[jj][line]) {
	    cout << "Gluino mass on line " << line  
	          << " of files " << filenames.at(ii) << " and " 
	         << filenames.at(jj) << " do not match" << endl;
	    bad=true;
	  }
	  if (lspmass[ii][line] != lspmass[jj][line]) {
	    cout << "LSP mass on line " << line 
	         << " of files " << filenames.at(ii) << " and " 
	         << filenames.at(jj) << " do not match" << endl;
	    bad=true;
	  }
	}
      }
    }
  }
  if (bad) return;
  cout << " The file formats match...Now fill ntuple" <<endl;

  //-----------------------------------------------------
  // Now the content of the files is loaded in arrays
  // We are going to stick everything in an ntuple
  //-----------------------------------------------------
  TFile* babyFile_ = TFile::Open(Form("%s", ntFile), "RECREATE");
  babyFile_->cd();
  babyTree_ = new TTree("tree", "A Baby Ntuple");
  babyTree_->SetDirectory(0);

  // define the variables
  float glmass_;
  float lspmass_;
  float xsec_;
  float xsecup_;
  float xsecdwn_;
  int   bestsr_ = 0;   // best signal region

  float effsr1_ = -1.;
  float effsr2_ = -1.;
  float effsr3_ = -1.;
  float effsr4_ = -1.;
  float effsr5_ = -1.;
  float effsr6_ = -1.;
  float effsr7_ = -1.;
  float effsrb_ = -1.;   // b = best SR

  float errsr1_ = -1.;
  float errsr2_ = -1.;
  float errsr3_ = -1.;
  float errsr4_ = -1.;
  float errsr5_ = -1.;
  float errsr6_ = -1.;
  float errsr7_ = -1.;
  float errsrb_ = -1.;   // b = best SR

  float obslimsr1_ = -1.;
  float obslimsr2_ = -1.;
  float obslimsr3_ = -1.;
  float obslimsr4_ = -1.;
  float obslimsr5_ = -1.;
  float obslimsr6_ = -1.;
  float obslimsr7_ = -1.;
  float obslimsrb_ = -1.;   // b = best SR

  float explimsr1_ = -1.;
  float explimsr2_ = -1.;
  float explimsr3_ = -1.;
  float explimsr4_ = -1.;
  float explimsr5_ = -1.;
  float explimsr6_ = -1.;
  float explimsr7_ = -1.;
  float explimsrb_ = -1.;   // b = best SR


  // put the branches in the tree
  babyTree_->Branch("glmass",  &glmass_);
  babyTree_->Branch("lspmass", &lspmass_);
  babyTree_->Branch("xsec",    &xsec_);
  babyTree_->Branch("xsecup",  &xsecup_);
  babyTree_->Branch("xsecdwn", &xsecdwn_);
  babyTree_->Branch("bestsr",  &bestsr_);


  babyTree_->Branch("effsr1", &effsr1_);
  babyTree_->Branch("effsr2", &effsr2_);
  babyTree_->Branch("effsr3", &effsr3_);
  babyTree_->Branch("effsr4", &effsr4_);
  babyTree_->Branch("effsr5", &effsr5_);
  babyTree_->Branch("effsr6", &effsr6_);
  babyTree_->Branch("effsr7", &effsr7_);
  babyTree_->Branch("effsrb", &effsrb_);   // b = best SR

  babyTree_->Branch("errsr1", &errsr1_);
  babyTree_->Branch("errsr2", &errsr2_);
  babyTree_->Branch("errsr3", &errsr3_);
  babyTree_->Branch("errsr4", &errsr4_);
  babyTree_->Branch("errsr5", &errsr5_);
  babyTree_->Branch("errsr6", &errsr6_);
  babyTree_->Branch("errsr7", &errsr7_);
  babyTree_->Branch("errsrb", &errsrb_);   // b = best SR

  babyTree_->Branch("obslimsr1", &obslimsr1_);
  babyTree_->Branch("obslimsr2", &obslimsr2_);
  babyTree_->Branch("obslimsr3", &obslimsr3_);
  babyTree_->Branch("obslimsr4", &obslimsr4_);
  babyTree_->Branch("obslimsr5", &obslimsr5_);
  babyTree_->Branch("obslimsr6", &obslimsr6_);
  babyTree_->Branch("obslimsr7", &obslimsr7_);
  babyTree_->Branch("obslimsrb", &obslimsrb_);   // b = best SR

  babyTree_->Branch("explimsr1", &explimsr1_);
  babyTree_->Branch("explimsr2", &explimsr2_);
  babyTree_->Branch("explimsr3", &explimsr3_);
  babyTree_->Branch("explimsr4", &explimsr4_);
  babyTree_->Branch("explimsr5", &explimsr5_);
  babyTree_->Branch("explimsr6", &explimsr6_);
  babyTree_->Branch("explimsr7", &explimsr7_);
  babyTree_->Branch("explimsrb", &explimsrb_);  // b = best SR

  // fill the variables.  
  // Note: the best region is the one with the smallest
  // ratio of expected limit and efficiency
  cout << "Filling an ntuple with " << numbLines << " entries" << endl;
  for (int il=0; il<numbLines; il++) {
    float best = 999999.;
    int ibest = -1;
    if (nentries[0] != 0) {
      glmass_  = glmass[0][il];
      lspmass_ = lspmass[0][il];
      effsr1_ = eff[0][il];
      errsr1_ = eff[0][il];
      obslimsr1_ = obslim[0][il];
      explimsr1_ = explim[0][il];
      float temp = explimsr1_/effsr1_;
      if (temp < best) {
	best = temp;
	ibest = 1;
      }
    }
    if (nentries[1] != 0) {
      glmass_  = glmass[1][il];
      lspmass_ = lspmass[1][il];
      effsr2_ = eff[1][il];
      errsr2_ = eff[1][il];
      obslimsr2_ = obslim[1][il];
      explimsr2_ = explim[1][il];
      float temp = explimsr2_/effsr2_;
      if (temp < best) {
	best = temp;
	ibest = 2;
      }
    }
    if (nentries[2] != 0) {
      glmass_  = glmass[2][il];
      lspmass_ = lspmass[2][il];
      effsr3_ = eff[2][il];
      errsr3_ = eff[2][il];
      obslimsr3_ = obslim[2][il];
      explimsr3_ = explim[2][il];
      float temp = explimsr3_/effsr3_;
      if (temp < best) {
	best = temp;
	ibest = 3;
      }
    }
    if (nentries[3] != 0) {
      glmass_  = glmass[3][il];
      lspmass_ = lspmass[3][il];
      effsr4_ = eff[3][il];
      errsr4_ = eff[3][il];
      obslimsr4_ = obslim[3][il];
      explimsr4_ = explim[3][il];
      float temp = explimsr4_/effsr4_;
      if (temp < best) {
	best = temp;
	ibest = 4;
      }
    }
    if (nentries[4] != 0) {
      glmass_  = glmass[4][il];
      lspmass_ = lspmass[4][il];
      effsr5_ = eff[4][il];
      errsr5_ = eff[4][il];
      obslimsr5_ = obslim[4][il];
      explimsr5_ = explim[4][il];
      float temp = explimsr5_/effsr5_;
      if (temp < best) {
	best = temp;
	ibest = 5;
      }
    }
    if (nentries[5] != 0) {
      glmass_  = glmass[5][il];
      lspmass_ = lspmass[5][il];
      effsr6_ = eff[5][il];
      errsr6_ = eff[5][il];
      obslimsr6_ = obslim[5][il];
      explimsr6_ = explim[5][il];
      float temp = explimsr6_/effsr6_;
      if (temp < best) {
	best = temp;
	ibest = 6;
      }
    }
    if (nentries[6] != 0) {
      glmass_  = glmass[6][il];
      lspmass_ = lspmass[6][il];
      effsr7_ = eff[6][il];
      errsr7_ = eff[6][il];
      obslimsr7_ = obslim[6][il];
      explimsr7_ = explim[6][il];
      float temp = explimsr7_/effsr7_;
      if (temp < best) {
	best = temp;
	ibest = 7;
      }
    }
    bestsr_ = ibest;
    effsrb_ = eff[ibest-1][il];
    errsrb_ = eff[ibest-1][il];
    obslimsrb_ = obslim[ibest-1][il];
    explimsrb_ = explim[ibest-1][il];

    // Fill the cross-section
    bool done = false;
    for (int xs=0; xs<nxsec-1; xs++) {
      if ( (glmass_ >= xsec_mass[xs] && glmass_ < xsec_mass[xs+1]) ||
           (glmass_ <  xsec_mass[xs] && xs==0) ) {

	float dd1 = glmass_ - xsec_mass[xs];
	float dd2 = xsec_mass[xs+1] - xsec_mass[xs];
	xsec_     = xsec[xs]    + (xsec[xs+1]    - xsec[xs])   *(dd1/dd2);
        xsecup_   = xsecup[xs]  + (xsecup[xs+1]  - xsecup[xs]) *(dd1/dd2);
        xsecdwn_  = xsecdwn[xs] + (xsecdwn[xs+1] - xsecdwn[xs])*(dd1/dd2);
	done = true;
	break;
      }
    }
    if (!done) {
      int xs = nxsec - 1; 
      float dd1 = glmass_ - xsec_mass[xs];
      float dd2 = xsec_mass[xs] - xsec_mass[xs-1];
      xsec_     = xsec[xs]    + (xsec[xs]    - xsec[xs-1])   *(dd1/dd2);
      xsecup_   = xsecup[xs]  + (xsecup[xs]  - xsecup[xs-1]) *(dd1/dd2);
      xsecdwn_  = xsecdwn[xs] + (xsecdwn[xs] - xsecdwn[xs-1])*(dd1/dd2);
    }


    // Now fill the ntuple
    babyTree_->Fill(); 

  }  
  // Now close the file

    babyFile_->cd();
    babyTree_->Write();
    babyFile_->Close();


}

