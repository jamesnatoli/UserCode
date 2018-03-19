// This is a header file for the slimanalysis research
// Define everything here
#if !defined SLIMANALYSIS_2017_H
#define SLIMANALYSIS_2017_H
#include <string>
#include <vector>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <TH1F.h>
#include <TH2F.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TLegend.h>
#include <TLine.h>
#include <TLatex.h>
#include <TFitResult.h>
#include <TMath.h>
#include <initializer_list>
#include "edges.h"

// Correction to wire-chamber positions
const float wc_xab = +0.642158;
const float wc_xbc = +0.209902;
const float wc_xca = -0.830777;
const float wc_yab = -8.84928 ;
const float wc_ybc = +18.0467 ;
const float wc_yca = -9.20219 ;

// (Possibly) using sigma of delta-position to clean up more
const float wc_xab_res = 1.16192 ;
const float wc_xbc_res = 1.05491 ;
const float wc_xca_res = 1.54024 ;
const float wc_yab_res = 1.2041  ;
const float wc_ybc_res = 0.892824;
const float wc_yca_res = 1.33447 ;

const int NUMCHAN = 8;

const vector<int> TIMESLICES {5, 6, 7, 8, 9};

// Fiducial region: for each tile, 4 points are needed (8 numbers)
// TB: top, bottom; LR: left, right
// x,y BL; x,y BR; x,y TL; x,y TR
const float fiducialX[NUMCHAN][4] = {
  {-33,  30, -46,  17}, // EJ-260
  {-25,  35, -40,  28}, // EJ-260 2P
  {-30,  35, -37,  30}, // EJ-200 
  // {-10,10,-10,10},   // Scint-X #signma
  // {  7,18,-6,  5},   // Scint-X finger
  {-32, -20, -42, -30}, // SCSN-81 finger 1
  {-16,  -4, -27, -15}, // SCSN-81 finger 2
  {  0,  13, -11,   2}, // SCSN-81 finger 3
  { 17,  30,   5,  18}, // SCSN-81 finger 4
  {-38,  32, -45,  22}  // SCSN-81 #sigma

  // old numbers
  /*{-42,35,-37,38},
  {-37,40,-37,40},
  {-50,29,-49,30},
  {-50,28,-50,28},
  {-52,27,-51,28} */
};
const float fiducialY[NUMCHAN][4] = {
  {-52,-38,15,29}, // EJ_260
  {-52,-40,25,41}, // EJ-260 2P
  {-47,-39,34,42}, // EJ-200
  // {-10,-10,10,10}, // Scint-X #signma
  // {-42,-41,20,22}, // Scint-X finger
  {-50,-48,23,25}, // SCSN-81 finger 1
  {-50,-48,25,27}, // SCSN-81 finger 2
  {-50,-48,28,30}, // SCSN-81 finger 3
  {-46,-44,30,32}, // SCSN-81 finger 4
  {-52,-45,20,33}  // SCSN-81 #sigma

  // old numbers
  /*{-48,-50,25,24},
  {-45,-45,33,33},
  {-45,-47,34,32},
  {-49,-49,32,32},
  {-41,-42,36,35}*/
};

// Identify channels we need to use                                                                                                                         
struct channel {
  int chan;
  int ieta;
  int idepth;
  string name;
};

vector<channel> channels = 
  {{6,23,5,"EJ-260"},
   {9,23,6,"EJ-260_2P"},
   {12,23,7,"EJ-200"},
   // {3,23,4,"Scint-XS"},
   // {1,23,2,"Scint-XF"},
   {4,24,4,"SCSN-81F1"},
   {7,24,5,"SCSN-81F2"},
   {10,24,6,"SCSN-81F3"},
   {13,24,7,"SCSN-81F4"},
   {5,25,4,"SCSN-81S"}
  };

const string entry[NUMCHAN] = {
  "EJ-260",
  "EJ-260 2P",
  "EJ-200",
  // "Scint-X #sigma",
  // "Scint-X finger",
  "SCSN-81 F1",
  "SCSN-81 F2",
  "SCSN-81 F3",
  "SCSN-81 F4",
  "SCSN-81 #sigma"
};

const int color[NUMCHAN] = {
  kBlack,
  kGreen,
  kBlue,
  // kBlue+1,
  // kGreen+1,
  kRed,
  kRed+1,
  kRed+2,
  kRed+3,
  kViolet
};

const int style[NUMCHAN] = {
  kDotted,
  kSolid,
  kSolid,
  // kSolid,                                                                                                                                   
  // kDashed,                                                                                                                                  
  kDotted,
  kSolid,
  kSolid,
  kSolid,
  kDashed
};

const int NCH = 16;
const int NTS = 10;

float rot_fiducialX[NUMCHAN][4];
float rot_fiducialY[NUMCHAN][4];
double thetas[NUMCHAN];

// distance from wire-chamber A to plastic tiles, in [mm]
const double z_ex = 7300;

// Debug this?
bool isFiducial(int tile, float x, float y);
void fill_Rot_Array();
double calc_theta(int channel_num);
double rotate_Point(double point_X, double point_Y, int channel_num, char xy);


#endif
