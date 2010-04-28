#ifndef SuPageSetup_cxx
#define SuPageSetup_cxx
// Conference Note 5931-CONF - Preliminary results - to be submitted soon 
// Please contact the D0 NP conveners: duperrin@cppm.in2p3.fr, gusbroo@nevis.columbia.edu
#include <iostream>
#include <sstream>
#include "TRandom.h"
#include "TStyle.h"
#include "TPave.h"
#include "TPaveLabel.h"
#include "TText.h"
#include "TPostScript.h"
// Is 25 pad's enough ?
#define NX_MAX 5
#define NY_MAX 5
TPad* Pad[NX_MAX][NY_MAX];
//Margins setup
#define TVSPACE  0.05 // space for title
#define H_MARGIN 0.03
#define V_MARGIN 0.05
#define HSPACE   0.02
#define VSPACE   0.02
#define PWW     (1.-2.*H_MARGIN)
#define PWH     (1.-2.*V_MARGIN-TVSPACE)
#define PW_XLOW (0.+H_MARGIN)
#define PW_YLOW (0.+V_MARGIN)
#define PW_XHG  (PW_XLOW+PWW)
#define PW_YHG  (PW_YLOW+PWH)

void  drawH1(TH1F* hist,char* hist_title, char* X_title, Int_t color, Int_t same);
void  drawH2(TH2F* hist,char* hist_title, char* X_title, char* Y_title, Int_t color);
void  paintH1(TH1F* hist,char* hist_title, char* X_title, Int_t color, Int_t same);
void  paintH2(TH2F* hist,char* hist_title, char* X_title, char* Y_title, Int_t color);
void  testLayout();
//void  NewPage(TPostScript* ps,Int_t nx,Int_t ny,Int_t cpad,char* page_title,char* meeting_id);
void lq3style(Int_t nx=1, Int_t ny=1);
Int_t SuPageSetup(Int_t nx,Int_t ny,Int_t cpad,const char* page_title,char* meeting_id,Int_t pstyle=1);
Int_t SuPageSetup(Int_t nx,Int_t ny, Int_t cpad,
		  const char* page_title, char* meeting_id,Int_t pstyle)
{
   void draw_page_title(const char*);
   void draw_page_basement(char*);
   void SetGlobal(Int_t nx,Int_t ny);
   void SetGpad(TPad* pad,Int_t iy,Int_t ix);

    Double_t  page_width=1.;
    Double_t  page_height=1.;
    Double_t  page_xlow=0.;
    Double_t  page_ylow=0.;
   

// reject unresonoble division attempts
   if(!(nx>0&&ny>0&&nx<=NX_MAX&&ny<=NY_MAX)){
    std::cout << nx << " X " << ny << 
     " ??? Are you crazy or joking? Check your eyes..." << std::endl;
    return -1;
   }
 // pads parameters
   if (pstyle>1){
    page_width=PWW;
    page_height=PWH;
    page_xlow=PW_XLOW;
    page_ylow=PW_YLOW;
    if(cpad) cpad=10; //  color
// page title
    draw_page_title(page_title);
// page basement
    draw_page_basement(meeting_id);
   } 
   Double_t  p_width =page_width  -(nx-1.)*HSPACE;
   Double_t  p_height=page_height -(ny-1.)*VSPACE;
   p_width  = p_width/nx;
   p_height = p_height/ny;

//    if(cpad) cpad=10; //  color
// // page title
//    draw_page_title(page_title);
// // page basement
//    draw_page_basement(meeting_id);
// //   
//   SetGlobal(nx,ny);
   lq3style(nx,ny);
// subpads layout
    for(Int_t irow=0;irow!=ny;irow++){
      for(Int_t iclmn=0;iclmn!=nx;iclmn++){
        Double_t px_low=iclmn       *(HSPACE+p_width) +page_xlow;
        Double_t py_low=(ny-irow-1) *(VSPACE+p_height)+page_ylow;
        Double_t px_hgh=px_low + p_width;
        Double_t py_hgh=py_low + p_height;
        char *name = new char [10];
        sprintf(name,"p%d_%d\0",iclmn+1,irow+1);
        TPad* pad = (TPad*)gROOT->FindObject(name); 
        if(pad) delete pad;
        pad=new TPad(name,name,px_low,py_low,px_hgh,py_hgh,cpad);
        Pad[irow][iclmn] = pad;

        SetGpad(Pad[irow][iclmn],nx,ny);
        Pad[irow][iclmn]->Draw();
        delete [] name;
       }
    }
 return 0;
//
}
//========================
void SetGpad(TPad* pad,Int_t iy,Int_t ix){

Double_t Lmrg=0.18-(0.12/NY_MAX)*iy;
Double_t Bmrg=0.15-0.02/NX_MAX*ix;
Double_t Tmrg=0.14-0.01/NX_MAX*ix;
     pad->Draw();
     pad->SetGridx(1);
     pad->SetGridy(1);
     pad->SetLeftMargin(Lmrg);
     pad->SetRightMargin(0.04);
     pad->SetTopMargin(Tmrg);
     pad->SetBottomMargin(Bmrg);
}
//---------------------------------------------------------------------------//
void SetGlobal(Int_t nx,Int_t ny){
  //===========
Double_t KLblSize  =0.070- (nx+ny)*0.005;
Double_t KTtlSize  =0.070- (ny+ny)*0.005;
Double_t KTtlW     = 0.80  - (nx+ny)*0.05;
Double_t KTtlH     = 0.070 - (nx+ny)*0.005;
//Double_t KStW      = 0.03  + (nx+ny)*0.1;
//Double_t KStH      = 0.03  + (nx+ny)*0.06;
Double_t Lbl_ofs_x = 0.004 - 0.003/NY_MAX*(nx+ny);
Double_t Lbl_ofs_y = 0.015 + (nx+ny)*0.006;
//Double_t Lbl_ofs_y = 0.030 + (nx+ny)*0.006;
Double_t Ttl_ofs_x = 1.2*Lbl_ofs_y;
Double_t Ttl_ofs_y = 2.4*Lbl_ofs_x;


gStyle->SetLineWidth(1);
gStyle->SetFrameLineWidth(1);
gStyle->SetFillStyle(3001);

gStyle->SetHistLineWidth(2);
//gStyle->SetHistLineColor(57);    
gStyle->SetHistLineColor(221);    //dark blue
gStyle->SetLabelSize(KLblSize,"XYZ");
gStyle->SetLabelOffset(Lbl_ofs_x,"X");
gStyle->SetLabelOffset(Lbl_ofs_y,"Y");

gStyle->SetLabelColor(kBlack,"XYZ");
gStyle->SetTitleSize(KTtlSize,"XYZ");
gStyle->SetTitleOffset(Ttl_ofs_x ,"X");
gStyle->SetTitleOffset(Ttl_ofs_y ,"Y");

gStyle->SetTitleW(KTtlW);           // w of title bar
gStyle->SetTitleH(KTtlH);           // h of title bar 

gStyle->SetTitleColor(60);          // axis title fill color 

gStyle->SetTitleTextColor(50);      // red  hist title fill color 
gStyle->SetTitleY(0.97);            // title position 
//gStyle->SetPadColor(10);
gStyle->SetPadColor(0);
gStyle->SetCanvasColor(0);
gStyle->SetPadTickX(1);
gStyle->SetPadTickY(1);
gStyle->SetHistFillColor(226);      // dark ultramarin
gStyle->SetHistFillColor(19);       // light grey
gStyle->SetNdivisions(106,"XY");
gStyle->SetStatFormat("6.2g");
gStyle->SetFitFormat("6.3g");
gStyle->SetOptStat(111);            // entries, mean value, rms
gStyle->SetOptDate(1);  // date in bottom left corner
// cout << "Ttl_ofs_x= " << Ttl_ofs_x << ", Ttl_ofs_y= " <<  Ttl_ofs_y << endl; 
}

// Page header
//---------------------------------------------------------------------------//
void draw_page_title( const char* title )
{
// Draw title
   TPaveLabel* page_title = (TPaveLabel*)gROOT->FindObject("SUpageTitle"); 
   if(page_title) delete page_title;
   page_title = new TPaveLabel(0.1,0.92,0.9,0.96,title);
   page_title->SetName("SUpageTitle");
   page_title->SetFillColor(19);
   // title->SetTextFont(32);
   // title->SetTextColor(226); //ultramarine
   // title->SetTextColor(106);  //magenta
   page_title->SetTextColor(225);  //magenta
   page_title->Draw();
}
// bottom bar
//---------------------------------------------------------------------------//
void draw_page_basement(char* meeting_id)
{
//  Date
    Double_t sign_x=0.03;
    Double_t sign_y=0.025;
    Float_t  sign_tsize=0.028;
    Int_t    sign_font=32;
    Int_t    sign_color=1;
    Int_t    sign_align=11;
 
    gStyle->SetDateX(sign_x);
    gStyle->SetDateY(sign_y);
    gStyle->GetAttDate()->SetTextFont(sign_font);
    gStyle->GetAttDate()->SetTextSize(sign_tsize);
    gStyle->GetAttDate()->SetTextColor(sign_color);
    gStyle->GetAttDate()->SetTextAlign(sign_align);
    gStyle->SetOptDate(0);  // date in bottom left corner
    
//  Affilation
    sign_x=0.81;
    TText*   aff       = (TText*)gROOT->FindObject("SU_AFF"); 
    if(aff)       delete aff;
    aff= new TText(sign_x,sign_y,"S.Uzunyan, NIU");
    aff->SetName("SU_AFF");
    aff->SetTextSize(sign_tsize);
    aff->SetTextColor(sign_color);
    aff->SetTextFont(sign_font);
    aff->SetTextAlign(12);
    aff->Draw();
//  Meet_name
   Double_t meet_name_x=0.50;
   Double_t meet_name_y=0.023;
   TText*   meet_name = (TText*)gROOT->FindObject("SUmeetName"); 
   if(meet_name) delete meet_name; 
   meet_name= new TText(meet_name_x,meet_name_y,meeting_id);
   meet_name->SetName("SUmeetName");
   meet_name->SetTextSize(sign_tsize*1.3);
   meet_name->SetTextColor(50);
   meet_name->SetTextFont(22);
   meet_name->SetTextAlign(22);
   meet_name->Draw();
}
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
void paintH1(TH1F* hist, char* H_title, char* X_title, Int_t color, Int_t  same ){
  //=======
  if(!hist){
   std::cout << " drawHistogram has been called with zero pointer to" << std::endl;
   std::cout << " histogram %s : skipping..."                         << std::endl;
   return;
  }

  hist->UseCurrentStyle();
  hist->SetXTitle(X_title);
  hist->SetTitle(H_title);

  if(color>=0) hist->SetFillColor(color);
     if(same>=0) { hist->Draw("same"); return ; }
     hist->Draw();
  return;
}
//---------------------------------------------------------------------------//
void drawH1(TH1F* hist, char* H_title, char* X_title, Int_t color, Int_t  same ){
  //=======
  if(!hist){
   std::cout << " drawHistogram has been called with zero pointer to" << std::endl;
   std::cout << " histogram %s : skipping..."                         << std::endl;
   return;
  }

  hist->UseCurrentStyle();
  hist->SetXTitle(X_title);
  hist->SetTitle(H_title);
  // hist->SetFillStyle(3001);
  if(color>=0) hist->SetFillColor(color);
  if(same>=0) { hist->Draw("hist same"); return ; }
  hist->Draw("hist");
  return;
}
//---------------------------------------------------------------------------//
void drawH2(TH2F * hist,  char* H_title, char* X_title, char* Y_title, Int_t color){
  //=======
  if(!hist){
   std::cout << " drawHistogram has been called with zero pointer to" << std::endl;
   std::cout << " histogram %s : skipping..."                         << std::endl;
   return;
  }
   hist->UseCurrentStyle();

   hist->SetFillColor(226);  
   hist->SetTitle(H_title);
   hist->SetXTitle(X_title);
   hist->SetYTitle(Y_title);
    //Size and position of the axis numbers
   hist->SetLabelSize(  1.2*(0.065 - 4*0.005),       "X");
   hist->SetLabelSize(  1.2*(0.065 - 4*0.005),       "Y");
   hist->SetLabelOffset(    (0.004- 0.003/(NY_MAX+NX_MAX)*4.0),"X");
   hist->SetLabelOffset(5.2*(0.004- 0.003/(NY_MAX+NX_MAX)*4.0),"Y");
    //X,Y Title 
   hist->SetTitleOffset(1.2 *(0.004- 0.003/(NY_MAX+NX_MAX)*4.0),"X");
   hist->SetTitleOffset(16.2*(0.004- 0.003/(NY_MAX+NX_MAX)*4.0),"Y");
   hist->SetTitleSize( (0.075 - 4*0.005), "X");
   hist->SetTitleSize(      (0.075 - 4*0.005)                  ,"Y");


  if(color>=0) hist->SetFillColor(color);
  hist->Draw("box");  
  
  return;
}
void paintH2(TH2F * hist,  char* H_title, char* X_title, char* Y_title, Int_t color){
  //=======
  if(!hist){
   std::cout << " drawHistogram has been called with zero pointer to" << std::endl;
   std::cout << " histogram %s : skipping..."                         << std::endl;
   return;
  }
   hist->UseCurrentStyle();

   hist->SetFillColor(226);  
   hist->SetTitle(H_title);
   hist->SetXTitle(X_title);
   hist->SetYTitle(Y_title);
//     //Size and position of the axis numbers
//    hist->SetLabelSize(  1.2*(0.065 - 4*0.005),       "X");
//    hist->SetLabelSize(  1.2*(0.065 - 4*0.005),       "Y");
//    hist->SetLabelOffset(    (0.004- 0.003/(NY_MAX+NX_MAX)*4.0),"X");
//    hist->SetLabelOffset(5.2*(0.004- 0.003/(NY_MAX+NX_MAX)*4.0),"Y");
//     //X,Y Title 
//    hist->SetTitleOffset(1.2 *(0.004- 0.003/(NY_MAX+NX_MAX)*4.0),"X");
//    hist->SetTitleOffset(16.2*(0.004- 0.003/(NY_MAX+NX_MAX)*4.0),"Y");
//    hist->SetTitleSize( (0.075 - 4*0.005), "X");
//    hist->SetTitleSize( (0.075 - 4*0.005) ,"Y");


  if(color>=0) hist->SetFillColor(color);
  //  hist->Draw("box");  
  
  return;
}

//---------------------------------------------------------------------------//
void testLayout(){
  //Draw test histogram in 1x1, 1x2, 2x1, 2x2, 1x3, 2x3, 3x3 layouts
  //Create test hists
   const Int_t   Nentries=10000;
   const Int_t   N_CHANNELS = 100;
   const Float_t N_CH_LO = -3.;
   const Float_t N_CH_HI =  3.;

    TH1F* test_1D = (TH1F*)gROOT->FindObject("test_1D"); 
    if(test_1D)       delete test_1D;
    TH2F* test_2D = (TH2F*)gROOT->FindObject("test_2D"); 
    if(test_2D)       delete test_2D;

    test_1D = new TH1F("1D hist test", "1D hist test",N_CHANNELS,N_CH_LO,N_CH_HI);
    test_2D = new TH2F("2D hist test", "2D hist test",
                               N_CHANNELS,N_CH_LO,N_CH_HI, N_CHANNELS,N_CH_LO,N_CH_HI);
   //fill histos
   for(Int_t i=0;i!=Nentries;++i){
     test_1D->Fill(gRandom->Landau(-1,0.5 ));  
     test_2D->Fill(gRandom->Landau(-1,0.5 ),
                   gRandom->Landau(-1,0.5 ));  
   }
  //Define Meeting name
  char meeting_id[]="*** Test  ***\0";
  //Create Canvas 
  TCanvas* c1=new TCanvas("c1","L2muEffAnalyse", 1,0,675,925);
  // create PS file
  TPostScript* ps = (TPostScript*)gROOT->FindObject("PageExamples.ps"); 
  if(ps) ps->Close();
  ps = new TPostScript("PageExamples.ps",111);
    gStyle->SetNdivisions(106,"XY");
    gStyle->SetStatFormat("6.2g");
    gStyle->SetFitFormat("6.3g");
    gStyle->SetOptStat(11);            //  entries, mean value, rms
  // Plot test pages
  // 1x1 =============================
    ps->NewPage();
    c1->Clear();
    SuPageSetup(1,1,10,"Test layout 1x1",meeting_id);
    Pad[0][0]->cd();    
    drawH2(test_2D ,"#Test 2-D hist", "#X, tugriks", "Y, tugriks ", -1); 
    c1->Update();
  // 1x2 =============================
    ps->NewPage();
    c1->Clear();
    SuPageSetup(1,2,10,"Test layout 1x2",meeting_id);
    Pad[0][0]->cd();    
    drawH1(test_1D ,"#Test 1-D hist","#X, tugriks", -1, -1); 
    Pad[1][0]->cd();    
    drawH2(test_2D ,"#Test 2-D hist","#X, tugriks", "Y, tugriks ", -1); 
    c1->Update();
  // 2x1 =============================
    ps->NewPage();
    c1->Clear();
    SuPageSetup(2,1,10,"Test layout 2x1",meeting_id);
    Pad[0][0]->cd();    
    drawH1(test_1D ,"#Test 1-D hist","#X, tugriks", -1, -1); 
    Pad[0][1]->cd();    
    drawH2(test_2D ,"#Test 2-D hist","#X, tugriks", "Y, tugriks ", -1); 
    c1->Update();
  // 1x3 =============================
    ps->NewPage();
    c1->Clear();
    SuPageSetup(1,3,10,"Test layout 1x3",meeting_id);
    Pad[0][0]->cd();    
    drawH1(test_1D ,"#Test 1-D hist","#X, tugriks", -1, -1); 

    Pad[1][0]->cd();    
    drawH1(test_1D ,"#Test 1-D hist","#X, tugriks", -1, -1); 

    Pad[2][0]->cd();    
    drawH2(test_2D ,"#Test 2-D hist","#X, tugriks", "Y, tugriks ", -1); 
    c1->Update();
  // 2x3 =============================
    ps->NewPage();
    c1->Clear();
    SuPageSetup(2,3,10,"Test layout 2x3",meeting_id);
    //row 1
    Pad[0][0]->cd();    
    drawH1(test_1D ,"#Test 1-D hist","#X, tugriks", -1, -1); 

    Pad[1][0]->cd();    
    drawH1(test_1D ,"#Test 1-D hist","#X, tugriks", -1, -1); 

    Pad[2][0]->cd();    
    drawH2(test_2D ,"#Test 1-D hist","#X, tugriks", "Y, tugriks ", -1); 
    //row 2
    Pad[0][1]->cd();    
    drawH1(test_1D ,"#Test 1-D hist","#X, tugriks", -1, -1); 

    Pad[1][1]->cd();    
    drawH1(test_1D ,"#Test 1-D hist","#X, tugriks", -1, -1); 

    Pad[2][1]->cd();    
    drawH2(test_2D ,"#Test 2-D hist","#X, tugriks", "Y, tugriks ", -1); 
    c1->Update();
  // 3x2 =============================
    ps->NewPage();
    c1->Clear();
    SuPageSetup(3,2,10,"Test layout 3x2",meeting_id);
    //row 1
    Pad[0][0]->cd();    
    drawH1(test_1D ,"#Test 1-D hist","#X, tugriks", -1, -1); 

    Pad[1][0]->cd();    
    drawH2(test_2D ,"#Test 1-D hist","#X, tugriks", "Y, tugriks ", -1); 
    //row 2
    Pad[0][1]->cd();    
    drawH1(test_1D ,"#Test 1-D hist","#X, tugriks", -1, -1); 

    Pad[1][1]->cd();    
    drawH2(test_2D ,"#Test 1-D hist","#X, tugriks", "Y, tugriks ", -1); 
    //row 3
    Pad[0][2]->cd();    
    drawH1(test_1D ,"#Test 1-D hist","#X, tugriks", -1, -1); 

    Pad[1][2]->cd();    
    drawH2(test_2D ,"#Test 2-D hist","#X, tugriks", "Y, tugriks ", -1); 
    c1->Update();
  // 3x3 =============================
    ps->NewPage();
    c1->Clear();
    SuPageSetup(3,3,10,"Test layout 3x3",meeting_id);
    //row 1
    Pad[0][0]->cd();    
    drawH1(test_1D ,"#Test 1-D hist","#X, tugriks", -1, -1); 

    Pad[1][0]->cd();    
    drawH2(test_2D ,"#Test 2-D hist","#X, tugriks", "Y, tugriks ", -1); 

    Pad[2][0]->cd();    
    drawH2(test_2D ,"#Test 2-D hist","#X, tugriks", "Y, tugriks ", 33); 
    //row 2
    Pad[0][1]->cd();    
    drawH1(test_1D ,"#Test 1-D hist","#X, tugriks", -1, -1); 

    Pad[1][1]->cd();    
    drawH2(test_2D ,"#Test 2-D hist","#X, tugriks", "Y, tugriks ", -1); 

    Pad[2][1]->cd();    
    drawH2(test_2D ,"#Test 2-D hist","#X, tugriks", "Y, tugriks ", 33); 
    //row 3
    Pad[0][2]->cd();    
    drawH1(test_1D ,"#Test 1-D hist","#X, tugriks", -1, -1); 

    Pad[1][2]->cd();    
    drawH2(test_2D ,"#Test 2-D hist","#X, tugriks", "Y, tugriks ", -1); 
    c1->Update();

    Pad[2][2]->cd();    
    drawH2(test_2D ,"#Test 2-D hist","#X, tugriks", "Y, tugriks ", 33); 
    c1->Update();
    //==================================================
   ps->Close(); 
   delete test_1D;
   delete test_2D;
   delete ps;
   delete c1;
   return;
}
void lq3style(Int_t nx,   Int_t ny) {
  //===========
Double_t KLblSize  = 0.055 - (nx+ny)*0.0001;
//Double_t KTtlSize  = 0.070 - (ny+ny)*0.005;
Double_t KTtlSize  = KLblSize;
Double_t KTtlW     = 0.75  - (nx+ny)*0.0001;
Double_t KTtlH     = 0.070 - (nx+ny)*0.0001;
//Double_t KStW      = 0.03  + (nx+ny)*0.1;
//Double_t KStH      = 0.03  + (nx+ny)*0.06;
Double_t Lbl_ofs_x = 0.015 - 0.003/NX_MAX*(nx+ny);
Double_t Lbl_ofs_y = 0.015 - 0.003/NY_MAX*(nx+ny);
Double_t Ttl_ofs_x = 1.1+KTtlSize + Lbl_ofs_x;
//Double_t Ttl_ofs_y = Lbl_ofs_y+2.4*KLblSize;
Double_t Ttl_ofs_y = 1.2+KTtlSize + 10.0*Lbl_ofs_y;;

gROOT->SetStyle();
//
gStyle->SetCanvasDefH(460);
gStyle->SetCanvasDefW(460);
//

Double_t Lmrg=0.18-(0.12/NY_MAX)*ny;
//Double_t Bmrg=2*Ttl_ofs_x+3*Lbl_ofs_x;
Double_t Bmrg=0.13-0.02/NX_MAX*nx;
Double_t Tmrg=0.14-0.01/NX_MAX*nx;
//cout << "Lmargun : " << gStyle->GetPadLeftMargin()  << " Tmargun : " << Tmrg  << " Bmargun : "<<  Bmrg << endl;;

gStyle->SetPadGridY(1);
gStyle->SetPadGridX(1);
gStyle->SetPadBottomMargin(Bmrg);
gStyle->SetPadTopMargin(0.8*Bmrg);
gStyle->SetPadLeftMargin(0.9*Bmrg);

gStyle->SetPadRightMargin(0.9*Bmrg);

//gStyle->SetPadRightMargin(0.9.*Bmrg);

//cout << "Lmargun : " << gStyle->GetPadLeftMargin()   << " Tmargun : " << gStyle->GetPadTopMargin() 
//     << " Bmargun : "<< gStyle->GetPadBottomMargin() << " Rmargun : " << gStyle->GetPadRightMargin()  << endl;
// cout << "Lbl_ofs_x " << Lbl_ofs_x << "  Ttl_ofs_x " << Ttl_ofs_x << " KTtlSize " << KTtlSize << endl;
//gStyle->SetLineWidth(1);
//gStyle->SetFrameLineWidth(1);
gStyle->SetFillStyle(3001);

gStyle->SetHistLineWidth(2);
//gStyle->SetHistLineColor(57);    
gStyle->SetHistLineColor(221);    //dark blue
gStyle->SetLabelSize(KLblSize,"XYZ");
gStyle->SetLabelOffset(Lbl_ofs_x,"X");
gStyle->SetLabelOffset(Lbl_ofs_y,"Y");

//gStyle->SetLabelColor(kBlack,"XYZ");
gStyle->SetTitleSize(KTtlSize,"XYZ");
gStyle->SetTitleOffset(Ttl_ofs_x ,"X");
gStyle->SetTitleOffset(Ttl_ofs_y ,"Y");

gStyle->SetTitleW(KTtlW);           // w of title bar
gStyle->SetTitleH(KTtlH);           // h of title bar 

//gStyle->SetTitleColor(102);          // axis title fill color 

gStyle->SetTitleTextColor(102);     // red  hist title fill color 
gStyle->SetTitleY(0.99);            // title position 
gStyle->SetPadColor(0);
gStyle->SetCanvasColor(0);
gStyle->SetPadTickX(1);
//gStyle->SetPadTickY(1);
//gStyle->SetHistFillColor(226);      // dark ultramarin
//gStyle->SetHistFillColor(19);       // light grey
gStyle->SetHistFillColor(0);       // light grey
gStyle->SetNdivisions(206,"XY");
//gStyle->SetNdivisions(106,"XY");
gStyle->SetStatFormat("6.2g");
gStyle->SetFitFormat("6.3g");
gStyle->SetOptStat(111);            // entries, mean value, rms
//gStyle->SetOptDate(0);  // date in bottom left corner
//cout << "Ttl_ofs_x= " << Ttl_ofs_x << ", Ttl_ofs_y= " <<  Ttl_ofs_y << endl; 
gStyle->SetOptDate(1);  // date in bottom left corner
}
#endif  // SuPageSetup_cxx

