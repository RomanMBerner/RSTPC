{
//gInterpreter->AddIncludePath("/home/francescop/ArCube/analysis/code");
gInterpreter->AddIncludePath("/home/rberner/analysis/RSTPC_analysis");
//gInterpreter->AddIncludePath("/home/francescop/analysis/code");
gInterpreter->AddIncludePath("/home/rberner/software/eigen/eigen-eigen-5a0156e40feb"); // should be before AddDynamicPath

gSystem->AddDynamicPath("/home/rberner/ACLiC_build"); // looks for .so files

gSystem->Load("DigitalFilters.so");
gSystem->Load("HistoManipulators.so");
gSystem->Load("RSTPC_Globals.so");
gSystem->Load("RSTPC_Analyser.so");
gSystem->Load("MppcTreeWrapper.so");
gSystem->Load("RSTPC_T1wrapper.so");
gSystem->Load("RSTPC_Hits.so");
gSystem->Load("RSTPC_RunProcessor.so");
}
