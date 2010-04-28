#include "SuPageSetup.cpp"
Sbottom2d_conf(Int_t istyle=1)
{
// Conference Note 5931-CONF - Preliminary results - to be submitted soon 
// Please contact the D0 NP conveners: duperrin@cppm.in2p3.fr, gusbroo@nevis.columbia.edu
  /////////////////////////////////////////////////////////////////////
   TCanvas *c1 = new TCanvas("c1", " ",0,0,775,675);
   SuPageSetup(1,1,10,  "MUJET/MHT exclusion plots ", "NP", istyle);
   Pad[0][0]->cd();
   gPad->SetGridy(0);
   gPad->SetGridx(0);
   gPad->SetTopMargin(0.05);

   gStyle->SetCanvasColor(0);
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetPadColor(0);
   gStyle->SetPadBorderMode(0);
   gStyle->SetFrameBorderMode(0);


   float left_margin=0.15;
   float right_margin=0.05;
   float top_margin=0.05;
   float bottom_margin=0.15;
   int canvas_color=0;
   int frame_color=0;

   gPad->SetLeftMargin(left_margin);
   gPad->SetRightMargin(right_margin);
   gPad->SetTopMargin(top_margin);
   gPad->SetBottomMargin(bottom_margin);
   gPad->SetFillColor(canvas_color);
   gPad->SetFrameFillColor(frame_color);
   
   gPad->SetBorderMode(0);
   gPad->SetBorderSize(1);
   gPad->SetFrameBorderMode(0);

   gPad->SetFillColor(0);
   gPad->SetFillStyle(4000);
   gPad->SetTickx();
   gPad->SetTicky();
   gPad->SetFrameFillStyle(0);
   gPad->SetFrameFillStyle(0);
   

   ////////////////////////////////
   // Histo - axe etc.
   ////////////////////////////////

   TH1 *H11 = new TH1F("H11","",230,0,270);

   int line_width=3;
   int line_color=1;
   int line_style=1; 
   int fill_style=1001;
   int fill_color=50;
   float y_min=-1111.;
   float y_max=-1111.;
   int ndivx=508;
   int ndivy=508;
   int marker_style=20;
   int marker_color=1;
   float marker_size=1.2;
   int optstat=0;

   H11->SetLineWidth(line_width);
   H11->SetLineColor(line_color);
   H11->SetLineStyle(line_style);
   H11->SetFillColor(fill_color);
   H11->SetFillStyle(fill_style);
   H11->SetMaximum(y_max);
   H11->SetMinimum(y_min);
   H11->GetXaxis()->SetNdivisions(ndivx);
   H11->GetYaxis()->SetNdivisions(ndivy);

   H11->SetMarkerStyle(marker_style);
   H11->SetMarkerColor(marker_color);
   H11->SetMarkerSize(marker_size);
   H11->SetStats(optstat);

   H11->SetMinimum(0);
   H11->SetMaximum(122);
   H11->SetDirectory(0);
   H11->SetStats(0);

   H11->GetXaxis()->SetTitle("Sbottom Mass (GeV)");
   H11->GetYaxis()->SetTitle("Neutralino Mass (GeV)");
  ////////////////////////////////
   // D0 RUN II - 4fb-1 Green contour
   ////////////////////////////////
   
   ////////////////////////////////////////////////////////////
   //Observed
   ////////////////////////////////////////////////////////////

   TGraph *graph_o4 = new TGraph(17);
   graph_o4->SetFillColor(3);
   graph_o4->SetLineWidth(4);
   graph_o4->SetLineColor(3);
   graph_o4->SetPoint(0,88, 0);
   graph_o4->SetPoint(1,88,60);
   graph_o4->SetPoint(2,90,65);
   graph_o4->SetPoint(3,100,70); 
   graph_o4->SetPoint(4,105,75);
   graph_o4->SetPoint(5,110,80);
   graph_o4->SetPoint(6,122,85);
   graph_o4->SetPoint(7,137,90);
   graph_o4->SetPoint(8,170,95);
   graph_o4->SetPoint(9,207,95); 
   graph_o4->SetPoint(10,215,90); 
   graph_o4->SetPoint(11,219,85);
   graph_o4->SetPoint(12,239,80);
   graph_o4->SetPoint(13,243,75);
   graph_o4->SetPoint(14,247,60); 
   graph_o4->SetPoint(15,251,35); 
   graph_o4->SetPoint(16,255,0.); 

   
   graph_o4->SetHistogram(H11);
   graph_o4->Draw("f");


   // (4fb observed)
   TGraph *graphL_o4 = new TGraph(17);
   graphL_o4->SetLineWidth(4);
   graphL_o4->SetLineColor(46);
   graphL_o4->SetPoint(0,88, 0);
   graphL_o4->SetPoint(1,88,60);
   graphL_o4->SetPoint(2,90,65);
   graphL_o4->SetPoint(3,100,70); 
   graphL_o4->SetPoint(4,105,75);
   graphL_o4->SetPoint(5,110,80);
   graphL_o4->SetPoint(6,122,85);
   graphL_o4->SetPoint(7,137,90);
   graphL_o4->SetPoint(8,170,95);
   graphL_o4->SetPoint(9,207,95); 
   graphL_o4->SetPoint(10,215,90); 
   graphL_o4->SetPoint(11,219,85); 
   graphL_o4->SetPoint(12,239,80);
   graphL_o4->SetPoint(13,243,75);
   graphL_o4->SetPoint(14,247,60); 
   graphL_o4->SetPoint(15,251,35); 
   graphL_o4->SetPoint(16,255,0.); 

   graphL_o4->SetHistogram(H11);
   graphL_o4->Draw("l");

   // (4fb expected, magenta
   TGraph *graphL_e4 = new TGraph(17);
   graphL_e4->SetLineWidth(4);
   graphL_e4->SetLineColor(46);
   graphL_e4->SetLineStyle(5);
   graphL_e4->SetPoint(0,88,0);
   graphL_e4->SetPoint(1,88,60);
   graphL_e4->SetPoint(2,90,65);
   graphL_e4->SetPoint(3,100,70); 
   graphL_e4->SetPoint(4,105,75);
   graphL_e4->SetPoint(5,110,80);
   graphL_e4->SetPoint(6,122,85);
   graphL_e4->SetPoint(7,137,90);
   graphL_e4->SetPoint(8,170,95);
   graphL_e4->SetPoint(9,214,95); 
   graphL_e4->SetPoint(10,232,90); 
   graphL_e4->SetPoint(11,228,85); 
   graphL_e4->SetPoint(12,239,80);
   graphL_e4->SetPoint(13,243,75);
   graphL_e4->SetPoint(14,247,60); 
   graphL_e4->SetPoint(15,251,35); 
   graphL_e4->SetPoint(16,255,0.); 

   graphL_e4->SetHistogram(H11);
   graphL_e4->Draw("l");

   ////////////////////////////////
   // D0 RUN II - Blue contour
   ////////////////////////////////
   
   ////////////////////////////////////////////////////////////
   //Single b-tag observed
   ////////////////////////////////////////////////////////////

   TGraph *graph_so = new TGraph(19);
   graph_so->SetFillColor(7);
   graph_so->SetLineWidth(4);
   graph_so->SetLineStyle(4);
   graph_so->SetLineColor(4);
   graph_so->SetPoint(0,93,0);
   graph_so->SetPoint(1,93,60);
   graph_so->SetPoint(2,100,64);
   graph_so->SetPoint(3,107,70);
   graph_so->SetPoint(4,120,73);
   graph_so->SetPoint(5,122,75);
   graph_so->SetPoint(6,130,77);
   graph_so->SetPoint(7,140,79);
   graph_so->SetPoint(8,150,80); 
   graph_so->SetPoint(9,155,79); 
   graph_so->SetPoint(10,160,83); 
   graph_so->SetPoint(11,170,85);
   graph_so->SetPoint(12,175,86);
   graph_so->SetPoint(13,180,87);
   graph_so->SetPoint(14,185,88);
   graph_so->SetPoint(15,205,87);
   graph_so->SetPoint(16,215,75); 
   graph_so->SetPoint(17,222,60); 
   graph_so->SetPoint(18,222,0.); 
   
   graph_so->SetHistogram(H11);
    graph_so->Draw("f");


   // (Single b-tag observed)
   TGraph *graphL_so = new TGraph(19);
   graphL_so->SetLineWidth(4);
   graphL_so->SetLineColor(4);
   graphL_so->SetPoint(0,93,0);
   graphL_so->SetPoint(1,93,60);
   graphL_so->SetPoint(2,100,64);
   graphL_so->SetPoint(3,107,70);
   graphL_so->SetPoint(4,120,73);
   graphL_so->SetPoint(5,122,75);
   graphL_so->SetPoint(6,130,77);
   graphL_so->SetPoint(7,140,79);
   graphL_so->SetPoint(8,150,80); 
   graphL_so->SetPoint(9,155,79); 
   graphL_so->SetPoint(10,160,83); 
   graphL_so->SetPoint(11,170,85);
   graphL_so->SetPoint(12,175,86);
   graphL_so->SetPoint(13,180,87);
   graphL_so->SetPoint(14,185,88);
   graphL_so->SetPoint(15,205,87);
   graphL_so->SetPoint(16,215,75); 
   graphL_so->SetPoint(17,222,60); 
   graphL_so->SetPoint(18,222,0.); 

   graphL_so->SetHistogram(H11);
    graphL_so->Draw("l");

   // single b-tagging expected
   TGraph *graphL_se = new TGraph(16);
   graphL_se->SetLineWidth(4);
   graphL_se->SetLineColor(4);
   graphL_se->SetLineStyle(2);
   graphL_se->SetPoint(0,100,0);
   graphL_se->SetPoint(1,100,64);
   graphL_se->SetPoint(2,107,70); 
   graphL_se->SetPoint(3,120,74);
   graphL_se->SetPoint(4,122,75);
   graphL_se->SetPoint(5,130,81);
   graphL_se->SetPoint(6,140,82); 
   graphL_se->SetPoint(7,155,82); 
   graphL_se->SetPoint(8,160,84);
   graphL_se->SetPoint(9,170,79); 
   graphL_se->SetPoint(10,175,75); 
   graphL_se->SetPoint(11,183,65); 
   graphL_se->SetPoint(12,185,61); 
   graphL_se->SetPoint(13,191,50);
   graphL_se->SetPoint(14,194,35);
   graphL_se->SetPoint(15,194,0);

   graphL_se->SetHistogram(H11);
   //  graphL_se->Draw("l");

   ////////////////////////////////
   // CDF RUN II - Yellow contour
   ////////////////////////////////

   TGraph *graph0 = new TGraph(9);
   graph0->SetName("");
   graph0->SetTitle("");
   graph0->SetFillColor(5);
   graph0->SetLineWidth(0);
   graph0->SetPoint(0,90,0);
   graph0->SetPoint(1,90,60);
   graph0->SetPoint(2,110,80);
   graph0->SetPoint(3,140,80);
   graph0->SetPoint(4,160,82.5);
   graph0->SetPoint(5,180,70);
   graph0->SetPoint(6,190,60);
   graph0->SetPoint(7,193,40);
   graph0->SetPoint(8,195,0);
 
   graph0->SetHistogram(H11);
   graph0->Draw("f");

   // red line
   TGraph *graph1 = new TGraph(16);
   graph1->SetName("");
   graph1->SetTitle("");
   graph1->SetFillColor(1);
   graph1->SetLineWidth(3);
   graph1->SetLineStyle(9);
   graph1->SetLineColor(2);
   graph1->SetPoint(0,90,0);
   graph1->SetPoint(1,90,60);
   graph1->SetPoint(2,110,80);
   graph1->SetPoint(3,140,80);
   graph1->SetPoint(4,160,82.5);
   graph1->SetPoint(5,180,70);
   graph1->SetPoint(6,190,60);
   graph1->SetPoint(7,193,40);
   graph1->SetPoint(8,195,0);


   graph1->SetHistogram(H11);
   graph1->Draw("al");



   ////////////////////////////////
   // CDF RUN I - Yellow contour
   ////////////////////////////////

   TGraph *graph0 = new TGraph(16);
   graph0->SetName("");
   graph0->SetTitle("");
   //   graph0->SetFillColor(5);
   graph0->SetFillColor(41);
   //graph0->SetFillStyle(1);
   graph0->SetLineWidth(0);
   //graph0->SetLineColor(2);
   graph0->SetPoint(0,34,0);
   graph0->SetPoint(1,35,9);
   graph0->SetPoint(2,51,29);
   graph0->SetPoint(3,59.5,35.5);
   graph0->SetPoint(4,70,42);
   graph0->SetPoint(5,85,56);
   graph0->SetPoint(6,93.5,60.5);
   graph0->SetPoint(7,102,66);
   graph0->SetPoint(8,111,70);
   graph0->SetPoint(9,118,72.5);
   graph0->SetPoint(10,130,75);
   graph0->SetPoint(11,140,68);
   graph0->SetPoint(12,143,55.5);
   graph0->SetPoint(13,146,50);
   graph0->SetPoint(14,146,9);
   graph0->SetPoint(15,148,0);
   
   graph0->SetHistogram(H11);
   graph0->Draw("f");

   // red line
   TGraph *graph1 = new TGraph(16);
   graph1->SetName("");
   graph1->SetTitle("");
   graph1->SetFillColor(1);
   graph1->SetLineWidth(3);
   graph1->SetLineColor(2);
   graph1->SetLineStyle(9);
   graph1->SetPoint(0,34,0);
   graph1->SetPoint(1,35,9);
   graph1->SetPoint(2,51,29);
   graph1->SetPoint(3,59.5,35.5);
   graph1->SetPoint(4,70,42);
   graph1->SetPoint(5,85,56);
   graph1->SetPoint(6,93.5,60.5);
   graph1->SetPoint(7,102,66);
   graph1->SetPoint(8,111,70);
   graph1->SetPoint(9,118,72.5);
   graph1->SetPoint(10,130,75);
   graph1->SetPoint(11,140,68);
   graph1->SetPoint(12,143,55.5);
   graph1->SetPoint(13,146,50);
   graph1->SetPoint(14,146,9);
   graph1->SetPoint(15,148,0);
   

   graph1->SetHistogram(H11);
   graph1->Draw("al");

   ////////////////////////////////
   // D0 RUN I - Magenta contour
   ////////////////////////////////

   TGraph *graph0DR1 = new TGraph(8);
   graph0DR1->SetName("");
   graph0DR1->SetTitle("");
   //   graph0DR1->SetFillColor(5);
   graph0DR1->SetFillColor(69);
   //graph0DR1->SetFillStyle(1);
   graph0DR1->SetLineWidth(0);
   //graph0->SetLineColor(2);
   graph0DR1->SetPoint(0,45,0);
   graph0DR1->SetPoint(1,45,17);
   graph0DR1->SetPoint(2,55,23);
   graph0DR1->SetPoint(3,70,39);
   graph0DR1->SetPoint(4,85,47);
   graph0DR1->SetPoint(5,100,40);
   graph0DR1->SetPoint(6,115,20);
   graph0DR1->SetPoint(7,115,0);
   
   graph0DR1->SetHistogram(H11);
   graph0DR1->Draw("f");

   // red line
   TGraph *graph1DR1 = new TGraph(8);
   graph1DR1->SetName("");
   graph1DR1->SetTitle("");
   graph1DR1->SetLineWidth(3);
   graph1DR1->SetLineColor(63);
   graph1DR1->SetPoint(0,45,0);
   graph1DR1->SetPoint(1,45,17);
   graph1DR1->SetPoint(2,55,23);
   graph1DR1->SetPoint(3,70,39);
   graph1DR1->SetPoint(4,85,47);
   graph1DR1->SetPoint(5,100,40);
   graph1DR1->SetPoint(6,115,20);
   graph1DR1->SetPoint(7,115,0);

   graph1DR1->SetHistogram(H11);
   graph1DR1->Draw("al");


   ////////////////////////////////
   // ADLO
   ////////////////////////////////
   // http://lepsusy.web.cern.ch/lepsusy/www/squarks_summer04/stop_combi_208_final.html

   TGraph *graph2 = new TGraph(5);
   graph2->SetLineWidth(3);
   graph2->SetLineColor(1);
   graph2->SetPoint(0,100,0);
   graph2->SetPoint(1,100,87.);
   graph2->SetPoint(2,99,88.);
   graph2->SetPoint(3,97,88.);
   graph2->SetPoint(4,37,28);

   graph2->SetHistogram(H11);
   
   graph2->Draw("");
   //redraw RunII
   //   graphL_o4->Draw("l");
   //
   ////////////////////////////////
   // TLegend
   ////////////////////////////////
 
   TLegend *leg = new TLegend(.74,.82,.89,.92,NULL,"brNDC");
   leg->SetTextSize(0.04);
   leg->SetFillStyle(0);
   leg->SetFillColor(0);
   leg->SetTextFont(22);
   leg->SetTextAlign(32);
   leg->SetBorderSize(0);
   leg->AddEntry(graphL_o4,"Observed","l");
   leg->AddEntry(graphL_e4,"Expected","l");
    leg->Draw("same");

   ////////////////////////////////
   // Text
   ////////////////////////////////
    TLatex* D0RunII = new TLatex(99,115,"D0 Run II Preliminary (4 fb^{-1})");
    //    D0RunII->SetNDC();
    D0RunII->SetTextSize(0.055);
    D0RunII->SetTextColor(1);
    D0RunII->SetTextFont(22);
    D0RunII->SetTextAlign(22);
    D0RunII->Draw();
// D0 Run I
   TLatex *   tex = new TLatex(66,25,"D0");
   tex->SetTextColor(63);
   tex->SetTextSize(0.032);
   tex->SetLineWidth(2);
   tex->Draw();

   TLatex *   tex = new TLatex(66,20,"Run I");
   tex->SetTextColor(63);
   tex->SetTextSize(0.032);
   tex->SetLineWidth(2);
   tex->Draw();

   TLatex *   tex = new TLatex(66,15,"92 pb^{-1}");
   tex->SetTextColor(63);
   tex->SetTextSize(0.032);
   tex->SetLineWidth(2);
   tex->Draw();

   // CDF Run I
   TLatex *   tex = new TLatex(115.1,40,"CDF");
   tex->SetTextColor(2);
   tex->SetTextSize(0.032);
   tex->SetLineWidth(2);
   tex->Draw();

   TLatex *   tex = new TLatex(115.1,35,"Run I");
   tex->SetTextColor(2);
   tex->SetTextSize(0.032);
   tex->SetLineWidth(2);
   tex->Draw();

   TLatex *   tex = new TLatex(115.1,30,"88 pb^{-1}");
   tex->SetTextColor(2);
   tex->SetTextSize(0.032);
   tex->SetLineWidth(2);
   tex->Draw();

  // CDF Run II
   TLatex *   tex = new TLatex(150.1,55,"CDF");
   tex->SetTextColor(2);
   tex->SetTextSize(0.032);
   tex->SetLineWidth(2);
   tex->Draw();

   TLatex *   tex = new TLatex(150.1,50,"Run IIa");
   tex->SetTextColor(2);
   tex->SetTextSize(0.032);
   tex->SetLineWidth(2);
   tex->Draw();

   TLatex *   tex = new TLatex(150.1,45,"295 pb^{-1}");
   tex->SetTextColor(2);
   tex->SetTextSize(0.032);
   tex->SetLineWidth(2);
   tex->Draw();




  // D0 Run II
   TLatex *   tex = new TLatex(179,81,"D0");
   tex->SetTextColor(4);
   //   tex->SetTextColor(9);
   tex->SetTextSize(0.032);
   tex->SetLineWidth(2);
   tex->Draw();
   TLatex *   tex = new TLatex(179,76,"Run IIa");
   tex->SetTextColor(4);
   //   tex->SetTextColor(9);
   tex->SetTextSize(0.032);
   tex->SetLineWidth(2);
   tex->Draw();

   TLatex *   tex = new TLatex(179,71,"310 pb^{-1}");
   tex->SetTextColor(4);
   tex->SetTextSize(0.032);
   tex->SetLineWidth(2);
   tex->Draw();

   TLatex *   tex = new TLatex(69,55,"LEP");
   tex->SetTextColor(1);
   tex->SetTextSize(0.032);
   tex->SetLineWidth(2);
   tex->Draw();

   TLatex *   tex = new TLatex(62,50,"#sqrt{s}=208 GeV");
   tex->SetTextColor(1);
   tex->SetTextSize(0.025);
   tex->SetLineWidth(2);
   tex->Draw();

   //
   TLine* line_b1 = new TLine(5.,0.,127.,122.);
   line_b1->SetLineWidth(2);
   line_b1->SetLineColor(1);
   line_b1->SetLineStyle(2);
   line_b1->Draw("same");


   TLatex *   tex = new TLatex(40,40,"M(#tilde{b}) = M(b) + M(#chi^{0}_{1})");
   tex->SetTextColor(1);
   tex->SetTextSize(0.03);
   tex->SetLineWidth(2);
   tex->SetTextAngle(63.);
   tex->Draw();


   TPad *name = new TPad("name", "title",0.174017,0.88563,0.60,0.932551);
   //name->Draw();
   name->cd();
   name->Range(0,0,1,1);
   name->SetFillColor(0);
   name->SetBorderMode(0);
   name->SetBorderSize(2);
   name->SetFrameBorderMode(0);
   name->Modified();


   Pad[0][0]->cd();
   //   gPad->cd();
   gPad->Modified();
   gPad->Update();

//    TText *text = new TText(0.,0.2,"D0 Run II Preliminary");
//    text->SetNDC();
//    text->SetTextSize(1.);
//    text->Draw();


//     char char_sigma_def[500];
//    //    sprintf(char_sigma_def,"p#bar{p}#rightarrowLQ+#bar{LQ}+X, NLO, #sqrt{s}=1.96 TeV");
//     sprintf(char_sigma_def,"D0 Run II Preliminary");
//   TLatex* sigma_def = new TLatex(0.65,.81,char_sigma_def);
 
 
   gPad->RedrawAxis();
   gPad->Update();   
 

   gPad->SaveAs("SB_contour.eps");
   gPad->SaveAs("SB_contour.png");
}
