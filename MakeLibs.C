void MakeLibs()
{
	gSystem->SetBuildDir(gSystem->GetBuildDir(), true);
	
	gInterpreter->AddIncludePath("/home/francescop/ArCube/analysis/code");
	gInterpreter->AddIncludePath("/home/francescop/analysis/include");
	
	//gSystem->Unload("RSTPC_Globals");
	//gSystem->Unload("RSTPC_Analyser");
	//gSystem->Unload("MppcTreeWrapper");
	//gSystem->Unload("RSTPC_RunProcessor");
	//gSystem->Unload("RSTPC_T1wrapper");
	
	gSystem->CompileMacro("src/DigitalFilters.cc", "gkf", "DigitalFilters");
	gSystem->CompileMacro("src/HistoManipulators.cc", "gkf", "HistoManipulators");
}
