#ifndef GLOBALS_H
#define GLOBALS_H
// Some enums to be used by plotting macros

#include "cross_sections-WAA.h"


// MC Samples definition ////////////////////////////////////////////////////
enum {qqWAAEI = 0, // qq' -> WAA -> e nu ISR                            //  0
      qqWAAMI,     // qq' -> WAA -> mu nu ISR                           //  1
      qqWAATI,     // qq' -> WAA -> tau nu ISR                          //  2
      qqWAAEF,     // qq' -> WAA -> e nu ISR                            //  3
      qqWAAMF,     // qq' -> WAA -> mu nu ISR                           //  4
      qqWAATF,     // qq' -> WAA -> tau nu ISR                          //  5

      // ZAA
      qqZAAE,                                                           //  6 
      qqZAAM,							        //  7 
      qqZAAT,							        //  8

      // QCD + top background
      TTBAR,							        //  9
      TTAA,                                                             // 10
      TtE,        // t-channel single-top -> e,mu,tau 		        // 11
      TtM,							        // 12
      TtT,							        // 13
      TsE,        // s-channel single-top 			        // 14
      TsM,							        // 15
      TsT,							        // 16
      WT,         // inclusive W + single-top			        // 17

      // AA + jets
      qqAA_10TO25,						        // 18
      qqAA_25TO250,						        // 19
      qqAA_250TOINF,						        // 20
      DIPHOTONJETS,						        // 21

      // VV + jets
      qqWZ2L2Q,							        // 22
      qqWZ3L1N,							        // 23
      qqZZ4E,							        // 24
      qqZZ4M,							        // 25
      qqZZ4T,							        // 26
      qqZZ2E2M,							        // 27
      qqZZ2E2T,							        // 28
      qqZZ2M2T,							        // 29
      qqWW2L2N,							        // 30

      // WWA
      qqWWA,							        // 31

      // WA + jets
      qqWA,							        // 32

      // ZA + jets
      qqZA,							        // 33

      // Z + jets
      qqZ,							        // 34

      // W + jets
      qqW,							        // 35

      // A + jets
      qqA,							        // 36

      // VVV
      qqZZZ,							        // 37
      qqWZZ,							        // 38
      qqWWZ,							        // 39
      qqWWW                                                             // 40
};

// Make arrays of enums; MCID will loop on the number of collective samples,
// and all the individual samples contained in each collective-sample array
// will be added together
// We need the individual samples only when opening the root files; once
// we have the histograms scaled and added together, we can just use
// the collective samples

const int MCID = 40+1;
const char* mcsample[MCID]={
  "qqwaa_enu_isr",
  "qqwaa_munu_isr",
  "qqwaa_taunu_isr",
  "qqwaa_enu_fsr",
  "qqwaa_munu_fsr",
  "qqwaa_taunu_fsr",

  "qqzza_ee",
  "qqzza_mumu",
  "qqzza_tautau",

  "tt",
  "ttaa",
  "st_t_enu",
  "st_t_munu",
  "st_t_taunu",
  "st_s_enu",
  "st_s_munu",
  "st_s_taunu",
  "wt",

  "qqaa_10to25",
  "qqaa_25to250",
  "qqaa_250toinf",
  "aaj",

  "wz_lnuqq",
  "wz_lnull",
  //"wz_qqll",
  //"wzp_taunull",
  //"wzp_lnutautau",
  //"wzp_taunutautau",
  //"wzp_qqtautau",
  "zz_eeee",
  "zz_mumumumu",
  "zz_tautautautau",
  "zz_eemumu",
  "zz_eetautau",
  "zz_mumutautau",
  "ww_lnulnu",
  //"zz_all",
  //"zz_llqq",
  //"zz_llll",
  //"zz_lnulnu",
  //"zz_lltautau",

  "wwa",

  "gw",
  //"gw_e",
  //"gw_mu",
  //"gw_tau",
  "gz",
  //"zgamma_ee",
  //"zgamma_mumu",
  //"zgamma_tautau",

  "w",
  //"w_enu_np0",
  //"w_enu_np1",
  //"w_enu_np2",
  //"w_enu_np3",
  //"w_enu_np4",
  //"w_enu_np5",
  //"w_munu_np0",
  //"w_munu_np1",
  //"w_munu_np2",
  //"w_munu_np3",
  //"w_munu_np4",
  //"w_munu_np5",
  //"w_taunu_np0",
  //"w_taunu_np1",
  //"w_taunu_np2",
  //"w_taunu_np3",
  //"w_taunu_np4",
  //"w_taunu_np5",

  "z",
  //"z_ee_np0",
  //"z_ee_np1",
  //"z_ee_np2",
  //"z_ee_np3",
  //"z_ee_np4",
  //"z_ee_np5",
  //"z_mumu_np0",
  //"z_mumu_np1",
  //"z_mumu_np2",
  //"z_mumu_np3",
  //"z_mumu_np4",
  //"z_mumu_np5",
  //"z_tautau_np0",
  //"z_tautau_np1",
  //"z_tautau_np2",
  //"z_tautau_np3",
  //"z_tautau_np4",
  //"z_tautau_np5",
  //"drellyan_tt",
  //"drellyan_ee",
  //"drellyan_mm",

  "zzz",
  "wzz",
  "wwz",
  "www"
};

// Cross-section * filter eff * BR * trig eff * rec eff * ....

const double scaleFactor[MCID] = {
  qqwaa_enu_isr_xsec,
  qqwaa_munu_isr_xsec,
  qqwaa_taunu_isr_xsec,
  qqwaa_enu_fsr_xsec,
  qqwaa_munu_fsr_xsec,
  qqwaa_taunu_fsr_xsec,

  qqzaa_ee_xsec,
  qqzaa_mumu_xsec,
  qqzaa_tautau_xsec,

  tt_xsec*filter_tt_eff,
  ttaa_xsec,
  st_t_enu_xsec,
  st_t_munu_xsec,
  st_t_taunu_xsec,
  st_s_enu_xsec,
  st_s_munu_xsec,
  st_s_taunu_xsec,
  wt_xsec,

  qqaa_10to25_xsec  *filter_qqaa_10to25_eff,
  qqaa_25to250_xsec *filter_qqaa_25to250_eff,
  qqaa_250toinf_xsec*filter_qqaa_250toinf_eff,
  aaj_xsec,

  wz_lnuqq_xsec,
  wz_lnull_xsec,

  zz_eeee_xsec,
  zz_mumumumu_xsec,
  zz_tautautautau_xsec,
  zz_eemumu_xsec,
  zz_eetautau_xsec,
  zz_mumutautau_xsec,
  ww_lnulnu_xsec,

  wwa_xsec,

  gw_xsec,
  gz_xsec,
  
  w_xsec,

  z_xsec,
  
  zzz_xsec,
  wzz_xsec,
  wwz_xsec,
  www_xsec

};

#endif // GLOBALS_H
