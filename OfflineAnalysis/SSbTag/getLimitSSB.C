//----------------------------------------------
// Simple function to return the limit for  
// the SSB analysis as an interpolation between 
// different values of the acceptance uncertainty
//
// inputs
// region = a number between 1 and 7
// uncer  = uncertainty on sig acceptance (units are such that 0.1 means 10%)
// observed + true if obs limit, otherwise expected limit
//
float getLimitSSb(int region, float uncert, bool observed=true){

  // check that the region makes sens
  if (region < 1 || region > 7) {
    cout << "Bad region = " << region << endl;
    return -1.;
  }

  // The uncertainties for which we calculated the limits
  float defunc[6];
  defunc[0] = 0.10;
  defunc[1] = 0.15;
  defunc[2] = 0.20;
  defunc[3] = 0.25;
  defunc[4] = 0.30;
  defunc[5] = 0.35;

  // The precalculated limits
  // First index is the region
  // second index is the uncertainty
  // LAST index is observed (index=0) or expected (index=1)
  float deflimit[7][6][2];
  // Region 1
  deflimit[0][0][0] =  7.36065; 
  deflimit[0][1][0] =  7.47573;
  deflimit[0][2][0] =  7.63392;
  deflimit[0][3][0] =  7.84519;
  deflimit[0][4][0] =  8.06697;
  deflimit[0][5][0] =  8.32292;
  deflimit[0][0][1] =  7.66041;
  deflimit[0][1][1] =  7.78795;
  deflimit[0][2][1] =  7.96758;
  deflimit[0][3][1] =  8.18392;
  deflimit[0][4][1] =  8.41051;
  deflimit[0][5][1] =  8.69804;


  // Region 2
  deflimit[1][0][0] = 6.86233;   
  deflimit[1][1][0] = 7.00144;
  deflimit[1][2][0] = 7.17087;
  deflimit[1][3][0] = 7.38516;
  deflimit[1][4][0] = 7.61948;
  deflimit[1][5][0] = 7.89678;
  deflimit[1][0][1] = 6.13321;
  deflimit[1][1][1] = 6.26641;
  deflimit[1][2][1] = 6.37852;
  deflimit[1][3][1] = 6.53499;
  deflimit[1][4][1] = 6.70809;
  deflimit[1][5][1] = 6.97797;

  // Region 3 
  deflimit[2][0][0] = 5.11835;
  deflimit[2][1][0] = 5.23902;
  deflimit[2][2][0] = 5.38619;
  deflimit[2][3][0] = 5.56763;
  deflimit[2][4][0] = 5.76217;
  deflimit[2][5][0] = 5.97349;
  deflimit[2][0][1] = 3.72018;
  deflimit[2][1][1] = 3.8084;
  deflimit[2][2][1] = 3.87847;
  deflimit[2][3][1] = 3.9698;
  deflimit[2][4][1] = 4.08307;
  deflimit[2][5][1] = 4.20746;

  // Region 4
  deflimit[3][0][0] = 7.25207;
  deflimit[3][1][0] = 7.41479;
  deflimit[3][2][0] = 7.61226;
  deflimit[3][3][0] = 7.85479;
  deflimit[3][4][0] = 8.14328;
  deflimit[3][5][0] = 8.44737;
  deflimit[3][0][1] = 5.89297;
  deflimit[3][1][1] = 5.99519;
  deflimit[3][2][1] = 6.12458;
  deflimit[3][3][1] = 6.27872;
  deflimit[3][4][1] = 6.472;
  deflimit[3][5][1] = 6.69073;

  // Region 5
  deflimit[4][0][0] = 4.66266;
  deflimit[4][1][0] = 4.74443;
  deflimit[4][2][0] = 4.86318;
  deflimit[4][3][0] = 4.99548;
  deflimit[4][4][0] = 5.14248;
  deflimit[4][5][0] = 5.30601;
  deflimit[4][0][1] = 4.4103;
  deflimit[4][1][1] = 4.4857;
  deflimit[4][2][1] = 4.58763;
  deflimit[4][3][1] = 4.69626;
  deflimit[4][4][1] = 4.82989;
  deflimit[4][5][1] = 4.973;

  // Region 6
  deflimit[5][0][0] = 2.82262;
  deflimit[5][1][0] = 2.7939;
  deflimit[5][2][0] = 2.7806;
  deflimit[5][3][0] = 2.78247;
  deflimit[5][4][0] = 2.80842;
  deflimit[5][5][0] = 2.84401;
  deflimit[5][0][1] = 2.95897;
  deflimit[5][1][1] = 2.98792;
  deflimit[5][2][1] = 3.03174;
  deflimit[5][3][1] = 3.0862;
  deflimit[5][4][1] = 3.16545;
  deflimit[5][5][1] = 3.24383;

  // Region 7
  deflimit[6][0][0] = 2.83979;
  deflimit[6][1][0] = 2.80244;
  deflimit[6][2][0] = 2.78913;
  deflimit[6][3][0] = 2.77535;
  deflimit[6][4][0] = 2.79348;
  deflimit[6][5][0] = 2.81533;
  deflimit[6][0][1] = 2.88701;
  deflimit[6][1][1] = 2.87543;
  deflimit[6][2][1] = 2.88644;
  deflimit[6][3][1] = 2.89446;
  deflimit[6][4][1] = 2.93053;
  deflimit[6][5][1] = 2.9846;


  float answer;
  int j = 1;
  if (observed) j=0;
  for (int i=0; i<5; i++) {
    float uncstep = defunc[i+1]-defunc[i]; 
    float dlimit  = deflimit[region-1][i+1][j] - deflimit[region-1][i][j];
    float dunc    = uncert - defunc[i];
    
    if (uncert >= defunc[i] && uncert < defunc[i+1]) {
      answer = deflimit[region-1][i][j] + dlimit*dunc/uncstep;
      return answer;
    } else if (i==0 && uncert < defunc[i]) {
      answer = deflimit[region-1][i][j] + dlimit*dunc/uncstep;
      return answer;
    } else if (i==4 && uncert >= defunc[i+1]) {
      dunc = uncert - defunc[i+1];      
      answer = deflimit[region-1][i+1][j] + dlimit*dunc/uncstep;
      return answer;
    }
  }
  cout << "Should never get here" << endl;
  return -1.;
}
