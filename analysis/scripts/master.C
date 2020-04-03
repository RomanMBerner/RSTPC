void master() {
  gROOT->ProcessLine(".x LoadLibs.C");
  gROOT->ProcessLine(".L RSTPC_T2wrapper.cc");
  gROOT->ProcessLine(".L TestPulses.C");
  gROOT->ProcessLine(".L PCA.C");
  gROOT->ProcessLine(".L Residuals.C");
  gROOT->ProcessLine(".L E_field_correction.C");
  gROOT->ProcessLine(".x StraightLines.C");
}
