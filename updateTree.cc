using namespace std;

// Put into the same folder as the trees to be updated and the histograms to be used one folder above this

void upd() {


	vector<string> sampleNameList;
	sampleNameList.push_back("sw8026v1001.processed.4T_aMCatNLO_FXFX_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.AStar_aMCatNLO_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.AStar_Madgraph_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.MergedData.root");
	sampleNameList.push_back("sw8026v1001.processed.ST_antitop_tChannel_4f_Powheg_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.ST_antitop_tWChannel_5f_Powheg_NoFullyHadronicDecays_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.ST_sChannel_4f_aMCatNLO_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.ST_top_tChannel_4f_Powheg_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.ST_top_tWllChannel_5f_Powheg_MadGraph_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.ST_top_tWChannel_5f_Powheg_NoFullyHadronicDecays_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.TT_Dilepton_aMCatNLO_FXFX_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.TT_Dilepton_Powheg_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.TTG_aMCatNLO_FXFX_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.TTHTobb_Powheg_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.TTHToNonbb_Powheg_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.TTJets_Dilepton_Madgraph_MLM_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.TT_Powheg_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.TT_Semileptonic_Powheg_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.TTWToLNu_aMCatNLO_FXFX_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.TTWToQQ_aMCatNLO_FXFX_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.TTZToLLNuNu_aMCatNLO_FXFX_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.TTZToQQ_aMCatNLO_FXFX_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.TZQ_LL_aMCatNLO_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.VH_ToNonbb_aMCatNLO_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.WJetsToLNu_aMCatNLO_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.WWTo1L1Nu2Q_aMCatNLO_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.WWTo2L2Nu_Powheg_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.WWToLNuQQ_Powheg_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.WWZ_aMCatNLO_FXFX_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.WZTo1L1Nu2Q_aMCatNLO_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.WZTo1L3Nu_aMCatNLO_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.WZTo2L2Q_aMCatNLO_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.WZTo3LNu_aMCatNLO_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.WZTo3LNu_Powheg_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.WZZ_aMCatNLO_FXFX_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.ZJets_Madgraph_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.ZZTo2L2Nu_Powheg_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.ZZTo2L2Q_aMCatNLO_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.ZZTo2Q2Nu_aMCatNLO_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.ZZTo4L_Powheg_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.ZZTo4Q_aMCatNLO_Summer16_25ns.root");
	sampleNameList.push_back("sw8026v1001.processed.ZZZ_aMCatNLO_FXFX_Summer16_25ns.root");
	//~ sampleNameList.push_back("sw8026v1001.processed.T6bbllslepton_msbottom_900_mneutralino_150.root");
	//~ sampleNameList.push_back("sw8026v1001.processed.T6bbllslepton_msbottom_900_mneutralino_500.root");
	//~ sampleNameList.push_back("sw8026v1001.processed.T6bbllslepton_msbottom_1000_mneutralino_300.root");
	//~ sampleNameList.push_back("sw8026v1001.processed.T6bbllslepton_msbottom_1200_mneutralino_200.root");
	//~ sampleNameList.push_back("sw8026v1001.processed.T6bbllslepton_msbottom_1200_mneutralino_1000.root");

	
	TFile *electronScaleFactorFile = new TFile("../ScaleFactorsElectrons.root");							
	electronScaleFactorFile->ls();
	TH2F * electronIDScaleFactorHisto = (TH2F*)electronScaleFactorFile->Get("GsfElectronToMVATightTightIP2DSIP3D4");
	TH2F * electronIsoScaleFactorHisto = (TH2F*)electronScaleFactorFile->Get("MVAVLooseElectronToMini");
	TH2F * electronConvMissHitScaleFactorHisto = (TH2F*)electronScaleFactorFile->Get("MVATightElectronToConvVetoIHit0");
	
	TFile *electronTrackScaleFactorFile = new TFile("../TrackScaleFactorsElectrons.root");							
	electronTrackScaleFactorFile->ls();
	TH2F * electronTrackScaleFactorHisto = (TH2F*)electronTrackScaleFactorFile->Get("EGamma_SF2D");
	
	TFile *muonIDScaleFactorFile = new TFile("../ScaleFactorMuonID.root");							
	muonIDScaleFactorFile->ls();
	TH2F * muonIDScaleFactorHisto = (TH2F*)muonIDScaleFactorFile->Get("SF");
	
	TFile *muonIsoScaleFactorFile = new TFile("../ScaleFactorMuonMiniIso.root");							
	muonIsoScaleFactorFile->ls();
	TH2F * muonIsoScaleFactorHisto = (TH2F*)muonIsoScaleFactorFile->Get("SF");
	
	TFile *muonIP2DScaleFactorFile = new TFile("../ScaleFactorMuonIP2D.root");							
	muonIP2DScaleFactorFile->ls();
	TH2F * muonIP2DScaleFactorHisto = (TH2F*)muonIP2DScaleFactorFile->Get("SF");
	
	TFile *muonSIP3DScaleFactorFile = new TFile("../ScaleFactorMuonSIP3D.root");							
	muonSIP3DScaleFactorFile->ls();
	TH2F * muonSIP3DcaleFactorHisto = (TH2F*)muonSIP3DScaleFactorFile->Get("SF");
	
	//~ TFile *muonTrackScaleFactorFile = new TFile("../TrackScaleFactorsMuons.root");							
	//~ muonTrackScaleFactorFile->ls();
	//~ TH1F * muonTrackScaleFactorEtaHistoInput = (TH1F*)muonTrackScaleFactorEtaGraph->GetHistogram();
	//~ TH1F * muonTrackScaleFactorVtxHisto = (TH1F*)muonTrackScaleFactorVtxGraph->GetHistogram();
	
	TH1F * muonTrackScaleFactorEtaHisto = new TH1F("muonTrackScaleFactorEtaHisto","muonTrackScaleFactorEta",12,0.,2.4);
	TH1F * muonTrackScaleFactorVtxHisto = new TH1F("muonTrackScaleFactorVtxHisto","muonTrackScaleFactorVtx",23,0.,46);
	
	muonTrackScaleFactorEtaHisto->SetBinContent(1,0.9970);
	muonTrackScaleFactorEtaHisto->SetBinContent(2,0.9977);
	muonTrackScaleFactorEtaHisto->SetBinContent(3,0.9981);
	muonTrackScaleFactorEtaHisto->SetBinContent(4,0.9978);
	muonTrackScaleFactorEtaHisto->SetBinContent(5,0.9980);
	muonTrackScaleFactorEtaHisto->SetBinContent(6,0.9972);
	muonTrackScaleFactorEtaHisto->SetBinContent(7,0.9962);
	muonTrackScaleFactorEtaHisto->SetBinContent(8,0.9955);
	muonTrackScaleFactorEtaHisto->SetBinContent(9,0.9958);
	muonTrackScaleFactorEtaHisto->SetBinContent(10,0.9939);
	muonTrackScaleFactorEtaHisto->SetBinContent(11,0.9929);
	muonTrackScaleFactorEtaHisto->SetBinContent(12,0.9873);
	
	muonTrackScaleFactorVtxHisto->SetBinContent(1,0.9995);
	muonTrackScaleFactorVtxHisto->SetBinContent(2,0.9993);
	muonTrackScaleFactorVtxHisto->SetBinContent(3,0.9994);
	muonTrackScaleFactorVtxHisto->SetBinContent(4,0.9991);
	muonTrackScaleFactorVtxHisto->SetBinContent(5,0.9991);
	muonTrackScaleFactorVtxHisto->SetBinContent(6,0.9985);
	muonTrackScaleFactorVtxHisto->SetBinContent(7,0.9978);
	muonTrackScaleFactorVtxHisto->SetBinContent(8,0.9972);
	muonTrackScaleFactorVtxHisto->SetBinContent(9,0.9963);
	muonTrackScaleFactorVtxHisto->SetBinContent(10,0.9959);
	muonTrackScaleFactorVtxHisto->SetBinContent(11,0.9949);
	muonTrackScaleFactorVtxHisto->SetBinContent(12,0.9941);
	muonTrackScaleFactorVtxHisto->SetBinContent(13,0.9932);
	muonTrackScaleFactorVtxHisto->SetBinContent(14,0.9929);
	muonTrackScaleFactorVtxHisto->SetBinContent(15,0.9913);
	muonTrackScaleFactorVtxHisto->SetBinContent(16,0.9902);
	muonTrackScaleFactorVtxHisto->SetBinContent(17,0.9878);
	muonTrackScaleFactorVtxHisto->SetBinContent(18,0.9877);
	muonTrackScaleFactorVtxHisto->SetBinContent(19,0.9854);
	muonTrackScaleFactorVtxHisto->SetBinContent(20,0.9882);
	muonTrackScaleFactorVtxHisto->SetBinContent(21,0.9745);
	muonTrackScaleFactorVtxHisto->SetBinContent(22,0.9676);
	muonTrackScaleFactorVtxHisto->SetBinContent(23,0.9466);
	
	TFile *FastSimElectronIDScaleFactorFile = new TFile("../FastSimScaleFactorElectronID.root");
	TH2D * FastSimElectronIDScaleFactorHisto = (TH2D*)FastSimElectronIDScaleFactorFile->Get("histo2D");
	
	TFile *FastSimElectronIsoScaleFactorFile = new TFile("../FastSimScaleFactorElectronIso.root");
	TH2D * FastSimElectronIsoScaleFactorHisto = (TH2D*)FastSimElectronIsoScaleFactorFile->Get("histo2D");
	
	TFile *FastSimElectronConvVetoScaleFactorFile = new TFile("../FastSimScaleFactorElectronConvVeto.root");
	TH2D * FastSimElectronConvVetoScaleFactorHisto = (TH2D*)FastSimElectronConvVetoScaleFactorFile->Get("histo2D");
	
	
	TFile *FastSimMuonIDScaleFactorFile = new TFile("../FastSimScaleFactorMuonID.root");
	TH2D * FastSimMuonIDScaleFactorHisto = (TH2D*)FastSimMuonIDScaleFactorFile->Get("histo2D");
	
	TFile *FastSimMuonIsoScaleFactorFile = new TFile("../FastSimScaleFactorMuonIso.root");
	TH2D * FastSimMuonIsoScaleFactorHisto = (TH2D*)FastSimMuonIsoScaleFactorFile->Get("histo2D");
	
	TFile *FastSimMuonIP2DScaleFactorFile = new TFile("../FastSimScaleFactorMuonIP2D.root");
	TH2D * FastSimMuonIP2DScaleFactorHisto = (TH2D*)FastSimMuonIP2DScaleFactorFile->Get("histo2D");
	
	TFile *FastSimMuonSIP3DScaleFactorFile = new TFile("../FastSimScaleFactorMuonSIP3D.root");
	TH2D * FastSimMuonSIP3DScaleFactorHisto = (TH2D*)FastSimMuonSIP3DScaleFactorFile->Get("histo2D");
			
	float leptonFullSimScaleFactor1,leptonFullSimScaleFactor2;
	float leptonFullSimScaleFactorErr1,leptonFullSimScaleFactorErr2;
	//~ float leptonFastSimScaleFactor1,leptonFastSimScaleFactor2;
	//~ float leptonFastSimScaleFactorErr1,leptonFastSimScaleFactorErr2;
	float tempPt,tempPtTrack,tempEta;
	float diMuPt1,diMuPt2,diMuEta1,diMuEta2;
	float diElePt1,diElePt2,diEleEta1,diEleEta2;
	float EMuPt1,EMuPt2,EMuEta1,EMuEta2;
	float diEleIdEff1,diEleIdEff2;
	//~ float diMuIdEff1,diMuIdEff2;
	//~ float EMuIdEff1,EMuIdEff2;
	
	int diEleHBHENoiseFilter,diEleHBHENoiseIsoFilter,diEleCSCTightHaloFilter,diEleGoodVertices,diEleEcalDeadCellTriggerPrimitiveFilter,diEleEEBadScFilter;
	int diMuHBHENoiseFilter,diMuHBHENoiseIsoFilter,diMuCSCTightHaloFilter,diMuGoodVertices,diMuEcalDeadCellTriggerPrimitiveFilter,diMuEEBadScFilter;
	int EMuHBHENoiseFilter,EMuHBHENoiseIsoFilter,EMuCSCTightHaloFilter,EMuGoodVertices,EMuEcalDeadCellTriggerPrimitiveFilter,EMuEEBadScFilter;
	int diMuNVertices,EMuNVertices,tempNVertices;
	//~ Long64_t diEleHBHENoiseFilter,diEleHBHENoiseIsoFilter,diEleCSCTightHaloFilter,diEleGoodVertices,diEleEcalDeadCellTriggerPrimitiveFilter,diEleEEBadScFilter;
	//~ Long64_t diMuHBHENoiseFilter,diMuHBHENoiseIsoFilter,diMuCSCTightHaloFilter,diMuGoodVertices,diMuEcalDeadCellTriggerPrimitiveFilter,diMuEEBadScFilter;
	//~ Long64_t EMuHBHENoiseFilter,EMuHBHENoiseIsoFilter,EMuCSCTightHaloFilter,EMuGoodVertices,EMuEcalDeadCellTriggerPrimitiveFilter,EMuEEBadScFilter;
	//~ Long64_t diMuNVertices,EMuNVertices,tempNVertices;
	
	bool HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v, HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v,HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v,HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v;
	bool HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v, HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v, HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v,HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v,HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v, HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v, HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v, HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v, HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v, HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v;
	bool HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v,HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v,HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v,HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v,HLT_Mu27_TkMu8_v,HLT_Mu30_TkMu11_v;
		
	int metFilterSummary,triggerSummary;
	
	Long64_t nentries;
	string sampleName;
	
	for (vector<std::string>::const_iterator fname = sampleNameList.begin(); fname != sampleNameList.end(); ++fname){						
		// cout << i.first() << endl;
		cout << *fname << endl;
		sampleName = (*fname).c_str();
		TFile *f = new TFile((*fname).c_str(),"update");
		f->cd("cutsV34DileptonFinalTrees");
				
		TTree *DiElectronTree = (TTree*)gDirectory->Get("EEDileptonTree");
		
		TBranch *diElectronFullSimScaleFactor1Branch = DiElectronTree->Branch("leptonFullSimScaleFactor1",&leptonFullSimScaleFactor1,"leptonFullSimScaleFactor1/F");
		TBranch *diElectronFullSimScaleFactor2Branch = DiElectronTree->Branch("leptonFullSimScaleFactor2",&leptonFullSimScaleFactor2,"leptonFullSimScaleFactor2/F");
		TBranch *diElectronFullSimScaleFactorErr1Branch = DiElectronTree->Branch("leptonFullSimScaleFactorErr1",&leptonFullSimScaleFactorErr1,"leptonFullSimScaleFactorErr1/F");
		TBranch *diElectronFullSimScaleFactorErr2Branch = DiElectronTree->Branch("leptonFullSimScaleFactorErr2",&leptonFullSimScaleFactorErr2,"leptonFullSimScaleFactorErr2/F");
		//~ TBranch *diElectronFastSimScaleFactor1Branch = DiElectronTree->Branch("leptonFastSimScaleFactor1",&leptonFastSimScaleFactor1,"leptonFastSimScaleFactor1/F");
		//~ TBranch *diElectronFastSimScaleFactor2Branch = DiElectronTree->Branch("leptonFastSimScaleFactor2",&leptonFastSimScaleFactor2,"leptonFastSimScaleFactor2/F");
		//~ TBranch *diElectronFastSimScaleFactorErr1Branch = DiElectronTree->Branch("leptonFastSimScaleFactorErr1",&leptonFastSimScaleFactorErr1,"leptonFastSimScaleFactorErr1/F");
		//~ TBranch *diElectronFastSimScaleFactorErr2Branch = DiElectronTree->Branch("leptonFastSimScaleFactorErr2",&leptonFastSimScaleFactorErr2,"leptonFastSimScaleFactorErr2/F");
		TBranch *diElectronMetFilterBranch = DiElectronTree->Branch("metFilterSummary",&metFilterSummary,"metFilterSummary/I");
		TBranch *diElectronTriggerBranch = DiElectronTree->Branch("triggerSummary",&triggerSummary,"triggerSummary/I");
		
		DiElectronTree->SetBranchAddress("pt1",&diElePt1);
		DiElectronTree->SetBranchAddress("pt2",&diElePt2);		
		DiElectronTree->SetBranchAddress("eta1",&diEleEta1);
		DiElectronTree->SetBranchAddress("eta2",&diEleEta2);
	
		DiElectronTree->SetBranchAddress("Flag_HBHENoiseFilter",&diEleHBHENoiseFilter);
		DiElectronTree->SetBranchAddress("Flag_HBHENoiseIsoFilter",&diEleHBHENoiseIsoFilter);
		DiElectronTree->SetBranchAddress("Flag_globalTightHalo2016Filter",&diEleCSCTightHaloFilter);
		DiElectronTree->SetBranchAddress("Flag_goodVertices",&diEleGoodVertices);
		DiElectronTree->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter",&diEleEcalDeadCellTriggerPrimitiveFilter);
		DiElectronTree->SetBranchAddress("Flag_eeBadScFilter",&diEleEEBadScFilter);
		
		DiElectronTree->SetBranchAddress("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",&HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v);
		DiElectronTree->SetBranchAddress("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",&HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v);
		DiElectronTree->SetBranchAddress("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v",&HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v);
		DiElectronTree->SetBranchAddress("HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v",&HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v);	
					
		
		nentries = DiElectronTree->GetEntries();
		if (sampleName.find("MergedData") == string::npos){
			for (Long64_t i = 0; i < nentries; i++){
				DiElectronTree->GetEntry(i);
				leptonFullSimScaleFactor1 = 1.;
				leptonFullSimScaleFactor2 = 1.;
				leptonFullSimScaleFactorErr1 = 0.;
				leptonFullSimScaleFactorErr2 = 0.;
				
				//~ leptonFastSimScaleFactor1 = 1.;
				//~ leptonFastSimScaleFactor2 = 1.;
				//~ leptonFastSimScaleFactorErr1 = 0.;
				//~ leptonFastSimScaleFactorErr2 = 0.;
				
				if (diElePt1 > 200.) tempPt = 199.;
				else tempPt = diElePt1;
				if (diElePt1 > 500.)tempPtTrack = 499.;
				else if (diElePt1 < 25.)tempPtTrack = 26.;
				else tempPtTrack = diElePt1;
				tempEta = fabs(diEleEta1);
					
				leptonFullSimScaleFactor1 = leptonFullSimScaleFactor1 * electronIDScaleFactorHisto->GetBinContent(electronIDScaleFactorHisto->GetXaxis()->FindBin(tempPt),electronIDScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				leptonFullSimScaleFactor1 = leptonFullSimScaleFactor1 * electronIsoScaleFactorHisto->GetBinContent(electronIsoScaleFactorHisto->GetXaxis()->FindBin(tempPt),electronIsoScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				leptonFullSimScaleFactor1 = leptonFullSimScaleFactor1 * electronConvMissHitScaleFactorHisto->GetBinContent(electronConvMissHitScaleFactorHisto->GetXaxis()->FindBin(tempPt),electronConvMissHitScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				leptonFullSimScaleFactor1 = leptonFullSimScaleFactor1 * electronTrackScaleFactorHisto->GetBinContent(electronTrackScaleFactorHisto->GetXaxis()->FindBin(diEleEta1),electronTrackScaleFactorHisto->GetYaxis()->FindBin(tempPtTrack));
				
				leptonFullSimScaleFactorErr1 = electronIDScaleFactorHisto->GetBinError(electronIDScaleFactorHisto->GetXaxis()->FindBin(tempPt),electronIDScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				leptonFullSimScaleFactorErr1 = (leptonFullSimScaleFactorErr1**2 + electronIsoScaleFactorHisto->GetBinError(electronIsoScaleFactorHisto->GetXaxis()->FindBin(tempPt),electronIsoScaleFactorHisto->GetYaxis()->FindBin(tempEta))**2)**0.5;
				leptonFullSimScaleFactorErr1 = (leptonFullSimScaleFactorErr1**2 + electronConvMissHitScaleFactorHisto->GetBinError(electronConvMissHitScaleFactorHisto->GetXaxis()->FindBin(tempPt),electronConvMissHitScaleFactorHisto->GetYaxis()->FindBin(tempEta))**2)**0.5;
				
				leptonFullSimScaleFactorErr1 = (leptonFullSimScaleFactorErr1**2 + 0.01**2)**0.5;
				
				//~ leptonFastSimScaleFactor1 = leptonFastSimScaleFactor1 * FastSimElectronIDScaleFactorHisto->GetBinContent(FastSimElectronIDScaleFactorHisto->GetXaxis()->FindBin(tempPt),FastSimElectronIDScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				//~ leptonFastSimScaleFactor1 = leptonFastSimScaleFactor1 * FastSimElectronIsoScaleFactorHisto->GetBinContent(FastSimElectronIsoScaleFactorHisto->GetXaxis()->FindBin(tempPt),FastSimElectronIsoScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				//~ leptonFastSimScaleFactor1 = leptonFastSimScaleFactor1 * FastSimElectronConvVetoScaleFactorHisto->GetBinContent(FastSimElectronConvVetoScaleFactorHisto->GetXaxis()->FindBin(tempPt),FastSimElectronConvVetoScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				//~ 
				//~ leptonFastSimScaleFactorErr1 = 0.02;
				
				if (diElePt2 > 200.) tempPt = 199.;
				else tempPt = diElePt2;
				if (diElePt2 > 500.)tempPtTrack = 499.;
				else if (diElePt2 < 25.)tempPtTrack = 26.;
				else tempPtTrack = diElePt2;
				tempEta = fabs(diEleEta2);
					
				
				leptonFullSimScaleFactor2 = leptonFullSimScaleFactor2 * electronIDScaleFactorHisto->GetBinContent(electronIDScaleFactorHisto->GetXaxis()->FindBin(tempPt),electronIDScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				leptonFullSimScaleFactor2 = leptonFullSimScaleFactor2 * electronIsoScaleFactorHisto->GetBinContent(electronIsoScaleFactorHisto->GetXaxis()->FindBin(tempPt),electronIsoScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				leptonFullSimScaleFactor2 = leptonFullSimScaleFactor2 * electronConvMissHitScaleFactorHisto->GetBinContent(electronConvMissHitScaleFactorHisto->GetXaxis()->FindBin(tempPt),electronConvMissHitScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				leptonFullSimScaleFactor2 = leptonFullSimScaleFactor2 * electronTrackScaleFactorHisto->GetBinContent(electronTrackScaleFactorHisto->GetXaxis()->FindBin(diEleEta2),electronTrackScaleFactorHisto->GetYaxis()->FindBin(tempPtTrack));
	
				leptonFullSimScaleFactorErr2 = electronIDScaleFactorHisto->GetBinError(electronIDScaleFactorHisto->GetXaxis()->FindBin(tempPt),electronIDScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				leptonFullSimScaleFactorErr2 = (leptonFullSimScaleFactorErr2**2 + electronIsoScaleFactorHisto->GetBinError(electronIsoScaleFactorHisto->GetXaxis()->FindBin(tempPt),electronIsoScaleFactorHisto->GetYaxis()->FindBin(tempEta))**2)**0.5;
				leptonFullSimScaleFactorErr2 = (leptonFullSimScaleFactorErr2**2 + electronConvMissHitScaleFactorHisto->GetBinError(electronConvMissHitScaleFactorHisto->GetXaxis()->FindBin(tempPt),electronConvMissHitScaleFactorHisto->GetYaxis()->FindBin(tempEta))**2)**0.5;
				
				leptonFullSimScaleFactorErr2 = (leptonFullSimScaleFactorErr2**2 + 0.01**2)**0.5;
				
								
				//~ leptonFastSimScaleFactor2 = leptonFastSimScaleFactor2 * FastSimElectronIDScaleFactorHisto->GetBinContent(FastSimElectronIDScaleFactorHisto->GetXaxis()->FindBin(tempPt),FastSimElectronIDScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				//~ leptonFastSimScaleFactor2 = leptonFastSimScaleFactor2 * FastSimElectronIsoScaleFactorHisto->GetBinContent(FastSimElectronIsoScaleFactorHisto->GetXaxis()->FindBin(tempPt),FastSimElectronIsoScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				//~ leptonFastSimScaleFactor2 = leptonFastSimScaleFactor2 * FastSimElectronConvVetoScaleFactorHisto->GetBinContent(FastSimElectronConvVetoScaleFactorHisto->GetXaxis()->FindBin(tempPt),FastSimElectronConvVetoScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				//~ 
				//~ leptonFastSimScaleFactorErr2 = 0.02;
								
				if (diEleHBHENoiseFilter > 0 && diEleHBHENoiseIsoFilter > 0 && diEleCSCTightHaloFilter > 0 && diEleGoodVertices > 0 && diEleEcalDeadCellTriggerPrimitiveFilter > 0) metFilterSummary = 1;
				//~ if (diEleHBHENoiseFilter > 0 && diEleHBHENoiseIsoFilter > 0 && diEleGoodVertices > 0 && diEleEcalDeadCellTriggerPrimitiveFilter > 0) metFilterSummary = 1;
				else metFilterSummary = 0;
				
				if (HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v > 0 || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v > 0 || HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v > 0 || HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v > 0 ) triggerSummary = 1;
				else triggerSummary = 0;
	
				
				diElectronFullSimScaleFactor1Branch->Fill();
				diElectronFullSimScaleFactor2Branch->Fill();
				diElectronFullSimScaleFactorErr1Branch->Fill();
				diElectronFullSimScaleFactorErr2Branch->Fill();
				//~ diElectronFastSimScaleFactor1Branch->Fill();
				//~ diElectronFastSimScaleFactor2Branch->Fill();
				//~ diElectronFastSimScaleFactorErr1Branch->Fill();
				//~ diElectronFastSimScaleFactorErr2Branch->Fill();
				diElectronMetFilterBranch->Fill();
				diElectronTriggerBranch->Fill();
			}
			
		}
		
		else{
			for (Long64_t i = 0; i < nentries; i++){
				DiElectronTree->GetEntry(i);
				leptonFullSimScaleFactor1 = 1.;
				leptonFullSimScaleFactor2 = 1.;
				leptonFullSimScaleFactorErr1 = 0.;
				leptonFullSimScaleFactorErr2 = 0.;
				
				//~ leptonFastSimScaleFactor1 = 1.;
				//~ leptonFastSimScaleFactor2 = 1.;
				//~ leptonFastSimScaleFactorErr1 = 0.;
				//~ leptonFastSimScaleFactorErr2 = 0.;				
				
				if (diEleHBHENoiseFilter > 0 && diEleHBHENoiseIsoFilter > 0 && diEleCSCTightHaloFilter > 0 && diEleGoodVertices > 0 && diEleEcalDeadCellTriggerPrimitiveFilter > 0 && diEleEEBadScFilter > 0) metFilterSummary = 1;
				else metFilterSummary = 0;	
				
				if (HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v > 0 || HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v > 0 || HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v > 0 || HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v > 0 ) triggerSummary = 1;
				else triggerSummary = 0;
	
				
				diElectronFullSimScaleFactor1Branch->Fill();
				diElectronFullSimScaleFactor2Branch->Fill();
				diElectronFullSimScaleFactorErr1Branch->Fill();
				diElectronFullSimScaleFactorErr2Branch->Fill();
				//~ diElectronFastSimScaleFactor1Branch->Fill();
				//~ diElectronFastSimScaleFactor2Branch->Fill();
				//~ diElectronFastSimScaleFactorErr1Branch->Fill();
				//~ diElectronFastSimScaleFactorErr2Branch->Fill();
				diElectronMetFilterBranch->Fill();	
				diElectronTriggerBranch->Fill();		
				
			}
			
		}
			
	
		DiElectronTree->Write(0,TObject::kOverwrite);
		
		TTree *EMuTree = (TTree*)gDirectory->Get("EMuDileptonTree");

		
		TBranch *EMuFullSimScaleFactor1Branch = EMuTree->Branch("leptonFullSimScaleFactor1",&leptonFullSimScaleFactor1,"leptonFullSimScaleFactor1/F");
		TBranch *EMuFullSimScaleFactor2Branch = EMuTree->Branch("leptonFullSimScaleFactor2",&leptonFullSimScaleFactor2,"leptonFullSimScaleFactor2/F");
		TBranch *EMuFullSimScaleFactorErr1Branch = EMuTree->Branch("leptonFullSimScaleFactorErr1",&leptonFullSimScaleFactorErr1,"leptonFullSimScaleFactorErr1/F");
		TBranch *EMuFullSimScaleFactorErr2Branch = EMuTree->Branch("leptonFullSimScaleFactorErr2",&leptonFullSimScaleFactorErr2,"leptonFullSimScaleFactorErr2/F");
		TBranch *EMuMetFilterBranch = EMuTree->Branch("metFilterSummary",&metFilterSummary,"metFilterSummary/I");
		//~ TBranch *EMuFastSimScaleFactor1Branch = EMuTree->Branch("leptonFastSimScaleFactor1",&leptonFastSimScaleFactor1,"leptonFastSimScaleFactor1/F");
		//~ TBranch *EMuFastSimScaleFactor2Branch = EMuTree->Branch("leptonFastSimScaleFactor2",&leptonFastSimScaleFactor2,"leptonFastSimScaleFactor2/F");
		//~ TBranch *EMuFastSimScaleFactorErr1Branch = EMuTree->Branch("leptonFastSimScaleFactorErr1",&leptonFastSimScaleFactorErr1,"leptonFastSimScaleFactorErr1/F");
		//~ TBranch *EMuFastSimScaleFactorErr2Branch = EMuTree->Branch("leptonFastSimScaleFactorErr2",&leptonFastSimScaleFactorErr2,"leptonFastSimScaleFactorErr2/F");
		TBranch *EMuMetFilterBranch = EMuTree->Branch("metFilterSummary",&metFilterSummary,"metFilterSummary/I");
		TBranch *EMuMetTriggerBranch = EMuTree->Branch("triggerSummary",&triggerSummary,"triggerSummary/I");
		
		
		
		EMuTree->SetBranchAddress("pt1",&EMuPt1);
		EMuTree->SetBranchAddress("pt2",&EMuPt2);		
		EMuTree->SetBranchAddress("eta1",&EMuEta1);
		EMuTree->SetBranchAddress("eta2",&EMuEta2);
		EMuTree->SetBranchAddress("nVertices",&EMuNVertices);
		
		EMuTree->SetBranchAddress("Flag_HBHENoiseFilter",&EMuHBHENoiseFilter);
		EMuTree->SetBranchAddress("Flag_HBHENoiseIsoFilter",&EMuHBHENoiseIsoFilter);
		EMuTree->SetBranchAddress("Flag_globalTightHalo2016Filter",&EMuCSCTightHaloFilter);
		EMuTree->SetBranchAddress("Flag_goodVertices",&EMuGoodVertices);
		EMuTree->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter",&EMuEcalDeadCellTriggerPrimitiveFilter);
		EMuTree->SetBranchAddress("Flag_eeBadScFilter",&EMuEEBadScFilter);	
		
		EMuTree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",&HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v);	
		EMuTree->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v",&HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v);	
		EMuTree->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v",&HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v);	
		EMuTree->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",&HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v);	
		EMuTree->SetBranchAddress("HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v",&HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v);	
		EMuTree->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v",&HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v);	
		EMuTree->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v",&HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v);	
		EMuTree->SetBranchAddress("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v",&HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v);	
		EMuTree->SetBranchAddress("HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v",&HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v);	
		EMuTree->SetBranchAddress("HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v",&HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v);
		
		
		nentries = EMuTree->GetEntries();
		
		if (sampleName.find("MergedData") == string::npos){
		
			for (Long64_t i = 0; i < nentries; i++){
				EMuTree->GetEntry(i);
				leptonFullSimScaleFactor1 = 1.;
				leptonFullSimScaleFactor2 = 1.;
				leptonFullSimScaleFactorErr1 = 0.;
				leptonFullSimScaleFactorErr2 = 0.03;
				
				//~ leptonFastSimScaleFactor1 = 1.;
				//~ leptonFastSimScaleFactor2 = 1.;
				//~ leptonFastSimScaleFactorErr1 = 0.02;
				//~ leptonFastSimScaleFactorErr2 = 0.02;
				
				
					
				if (EMuPt1 > 200.) tempPt = 199.;
				else tempPt = EMuPt1;
				if (EMuPt1 > 500.)tempPtTrack = 499.;
				else if (EMuPt1 < 25.)tempPtTrack = 26.;
				else tempPtTrack = EMuPt1;
				tempEta = fabs(EMuEta1);
					
				leptonFullSimScaleFactor1 = leptonFullSimScaleFactor1 * electronIDScaleFactorHisto->GetBinContent(electronIDScaleFactorHisto->GetXaxis()->FindBin(tempPt),electronIDScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				leptonFullSimScaleFactor1 = leptonFullSimScaleFactor1 * electronIsoScaleFactorHisto->GetBinContent(electronIsoScaleFactorHisto->GetXaxis()->FindBin(tempPt),electronIsoScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				leptonFullSimScaleFactor1 = leptonFullSimScaleFactor1 * electronConvMissHitScaleFactorHisto->GetBinContent(electronConvMissHitScaleFactorHisto->GetXaxis()->FindBin(tempPt),electronConvMissHitScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				leptonFullSimScaleFactor1 = leptonFullSimScaleFactor1 * electronTrackScaleFactorHisto->GetBinContent(electronTrackScaleFactorHisto->GetXaxis()->FindBin(EMuEta1),electronTrackScaleFactorHisto->GetYaxis()->FindBin(tempPtTrack));
				
				leptonFullSimScaleFactorErr1 = electronIDScaleFactorHisto->GetBinError(electronIDScaleFactorHisto->GetXaxis()->FindBin(tempPt),electronIDScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				leptonFullSimScaleFactorErr1 = (leptonFullSimScaleFactorErr1**2 + electronIsoScaleFactorHisto->GetBinError(electronIsoScaleFactorHisto->GetXaxis()->FindBin(tempPt),electronIsoScaleFactorHisto->GetYaxis()->FindBin(tempEta))**2)**0.5;
				leptonFullSimScaleFactorErr1 = (leptonFullSimScaleFactorErr1**2 + electronConvMissHitScaleFactorHisto->GetBinError(electronConvMissHitScaleFactorHisto->GetXaxis()->FindBin(tempPt),electronConvMissHitScaleFactorHisto->GetYaxis()->FindBin(tempEta))**2)**0.5;
				
				leptonFullSimScaleFactorErr1 = (leptonFullSimScaleFactorErr1**2 + 0.01**2)**0.5;
				
				//~ leptonFastSimScaleFactor1 = leptonFastSimScaleFactor1 * FastSimElectronIDScaleFactorHisto->GetBinContent(FastSimElectronIDScaleFactorHisto->GetXaxis()->FindBin(tempPt),FastSimElectronIDScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				//~ leptonFastSimScaleFactor1 = leptonFastSimScaleFactor1 * FastSimElectronIsoScaleFactorHisto->GetBinContent(FastSimElectronIsoScaleFactorHisto->GetXaxis()->FindBin(tempPt),FastSimElectronIsoScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				//~ leptonFastSimScaleFactor1 = leptonFastSimScaleFactor1 * FastSimElectronConvVetoScaleFactorHisto->GetBinContent(FastSimElectronConvVetoScaleFactorHisto->GetXaxis()->FindBin(tempPt),FastSimElectronConvVetoScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				
				
				if (EMuPt2 > 120.)tempPt = 119.;
				else tempPt = EMuPt2;
				tempEta = fabs(EMuEta2);
					
				if (EMuNVertices > 45)tempNVertices = 45;
				else tempNVertices = EMuNVertices;
				
					
				leptonFullSimScaleFactor2 = leptonFullSimScaleFactor2 * muonIDScaleFactorHisto->GetBinContent(muonIDScaleFactorHisto->GetXaxis()->FindBin(tempPt),muonIDScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				leptonFullSimScaleFactor2 = leptonFullSimScaleFactor2 * muonIsoScaleFactorHisto->GetBinContent(muonIsoScaleFactorHisto->GetXaxis()->FindBin(tempPt),muonIsoScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				leptonFullSimScaleFactor2 = leptonFullSimScaleFactor2 * muonIP2DScaleFactorHisto->GetBinContent(muonIP2DScaleFactorHisto->GetXaxis()->FindBin(tempPt),muonIP2DScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				leptonFullSimScaleFactor2 = leptonFullSimScaleFactor2 * muonSIP3DcaleFactorHisto->GetBinContent(muonSIP3DcaleFactorHisto->GetXaxis()->FindBin(tempPt),muonSIP3DcaleFactorHisto->GetYaxis()->FindBin(tempEta));
				
				leptonFullSimScaleFactor2 = leptonFullSimScaleFactor2 * muonTrackScaleFactorEtaHisto->GetBinContent(muonTrackScaleFactorEtaHisto->GetXaxis()->FindBin(tempEta));
				leptonFullSimScaleFactor2 = leptonFullSimScaleFactor2 * muonTrackScaleFactorVtxHisto->GetBinContent(muonTrackScaleFactorVtxHisto->GetXaxis()->FindBin(tempNVertices));
				
				//~ leptonFastSimScaleFactor2 = leptonFastSimScaleFactor2 * FastSimMuonIDScaleFactorHisto->GetBinContent(FastSimMuonIDScaleFactorHisto->GetXaxis()->FindBin(tempPt),FastSimMuonIDScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				//~ leptonFastSimScaleFactor2 = leptonFastSimScaleFactor2 * FastSimMuonIsoScaleFactorHisto->GetBinContent(FastSimMuonIsoScaleFactorHisto->GetXaxis()->FindBin(tempPt),FastSimMuonIsoScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				//~ leptonFastSimScaleFactor2 = leptonFastSimScaleFactor2 * FastSimMuonIP2DScaleFactorHisto->GetBinContent(FastSimMuonIP2DScaleFactorHisto->GetXaxis()->FindBin(tempPt),FastSimMuonIP2DScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				//~ leptonFastSimScaleFactor2 = leptonFastSimScaleFactor2 * FastSimMuonSIP3DScaleFactorHisto->GetBinContent(FastSimMuonSIP3DScaleFactorHisto->GetXaxis()->FindBin(tempPt),FastSimMuonSIP3DScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				
				
									
				if (EMuHBHENoiseFilter > 0 && EMuHBHENoiseIsoFilter > 0 && EMuCSCTightHaloFilter > 0 && EMuGoodVertices > 0 && EMuEcalDeadCellTriggerPrimitiveFilter > 0) metFilterSummary = 1;
				//~ if (EMuHBHENoiseFilter > 0 && EMuHBHENoiseIsoFilter > 0 && EMuGoodVertices > 0 && EMuEcalDeadCellTriggerPrimitiveFilter > 0) metFilterSummary = 1;
				else metFilterSummary = 0;
				
				if (HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v > 0 || HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v > 0 || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v > 0 || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v > 0 || HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v > 0 || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v > 0 || HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v > 0 || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v > 0 || HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v > 0 || HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v > 0) triggerSummary = 1;
				else triggerSummary = 0;
	
				
				EMuFullSimScaleFactor1Branch->Fill();
				EMuFullSimScaleFactor2Branch->Fill();
				EMuFullSimScaleFactorErr1Branch->Fill();
				EMuFullSimScaleFactorErr2Branch->Fill();
				//~ EMuFastSimScaleFactor1Branch->Fill();
				//~ EMuFastSimScaleFactor2Branch->Fill();
				//~ EMuFastSimScaleFactorErr1Branch->Fill();
				//~ EMuFastSimScaleFactorErr2Branch->Fill();				
				EMuMetFilterBranch->Fill();
				EMuMetTriggerBranch->Fill();
			}
		}
		else{
			leptonFullSimScaleFactor1 = 1.;
			leptonFullSimScaleFactor2 = 1.;
			leptonFullSimScaleFactorErr1 = 0.;
			leptonFullSimScaleFactorErr2 = 0.;
			//~ leptonFastSimScaleFactor1 = 1.;
			//~ leptonFastSimScaleFactor2 = 1.;
			//~ leptonFastSimScaleFactorErr1 = 0.;
			//~ leptonFastSimScaleFactorErr2 = 0.;
			
			for (Long64_t i = 0; i < nentries; i++){
				EMuTree->GetEntry(i);			
				
				if (EMuHBHENoiseFilter > 0 && EMuHBHENoiseIsoFilter > 0 && EMuCSCTightHaloFilter > 0 && EMuGoodVertices > 0 && EMuEcalDeadCellTriggerPrimitiveFilter > 0 && EMuEEBadScFilter > 0) metFilterSummary = 1;
				else metFilterSummary = 0;
				
				if (HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v > 0 || HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v > 0 || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v > 0 || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v > 0 || HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v > 0 || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v > 0 || HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v > 0 || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v > 0 || HLT_Mu30_Ele30_CaloIdL_GsfTrkIdVL_v > 0 || HLT_Mu33_Ele33_CaloIdL_GsfTrkIdVL_v > 0) triggerSummary = 1;
				else triggerSummary = 0;
	
				
				EMuFullSimScaleFactor1Branch->Fill();
				EMuFullSimScaleFactor2Branch->Fill();
				EMuFullSimScaleFactorErr1Branch->Fill();
				EMuFullSimScaleFactorErr2Branch->Fill();
				//~ EMuFastSimScaleFactor1Branch->Fill();
				//~ EMuFastSimScaleFactor2Branch->Fill();
				//~ EMuFastSimScaleFactorErr1Branch->Fill();
				//~ EMuFastSimScaleFactorErr2Branch->Fill();			
				EMuMetFilterBranch->Fill();
				EMuMetTriggerBranch->Fill();
				
			}
			
		}
		
		EMuTree->Write(0,TObject::kOverwrite);
		
		TTree *DiMuonTree = (TTree*)gDirectory->Get("MuMuDileptonTree");
		
		TBranch *diMuonFullSimScaleFactor1Branch = DiMuonTree->Branch("leptonFullSimScaleFactor1",&leptonFullSimScaleFactor1,"leptonFullSimScaleFactor1/F");
		TBranch *diMuonFullSimScaleFactor2Branch = DiMuonTree->Branch("leptonFullSimScaleFactor2",&leptonFullSimScaleFactor2,"leptonFullSimScaleFactor2/F");
		TBranch *diMuonFullSimScaleFactorErr1Branch = DiMuonTree->Branch("leptonFullSimScaleFactorErr1",&leptonFullSimScaleFactorErr1,"leptonFullSimScaleFactorErr1/F");
		TBranch *diMuonFullSimScaleFactorErr2Branch = DiMuonTree->Branch("leptonFullSimScaleFactorErr2",&leptonFullSimScaleFactorErr2,"leptonFullSimScaleFactorErr2/F");
		TBranch *diMuonMetFilterBranch = DiMuonTree->Branch("metFilterSummary",&metFilterSummary,"metFilterSummary/I");
		//~ TBranch *diMuonFastSimScaleFactor1Branch = DiMuonTree->Branch("leptonFastSimScaleFactor1",&leptonFastSimScaleFactor1,"leptonFastSimScaleFactor1/F");
		//~ TBranch *diMuonFastSimScaleFactor2Branch = DiMuonTree->Branch("leptonFastSimScaleFactor2",&leptonFastSimScaleFactor2,"leptonFastSimScaleFactor2/F");
		//~ TBranch *diMuonFastSimScaleFactorErr1Branch = DiMuonTree->Branch("leptonFastSimScaleFactorErr1",&leptonFastSimScaleFactorErr1,"leptonFastSimScaleFactorErr1/F");
		//~ TBranch *diMuonFastSimScaleFactorErr2Branch = DiMuonTree->Branch("leptonFastSimScaleFactorErr2",&leptonFastSimScaleFactorErr2,"leptonFastSimScaleFactorErr2/F");
		TBranch *diMuonMetFilterBranch = DiMuonTree->Branch("metFilterSummary",&metFilterSummary,"metFilterSummary/I");
		TBranch *diMuonTriggerBranch = DiMuonTree->Branch("triggerSummary",&triggerSummary,"triggerSummary/I");	
		
		DiMuonTree->SetBranchAddress("pt1",&diMuPt1);
		DiMuonTree->SetBranchAddress("pt2",&diMuPt2);		
		DiMuonTree->SetBranchAddress("eta1",&diMuEta1);
		DiMuonTree->SetBranchAddress("eta2",&diMuEta2);		
		DiMuonTree->SetBranchAddress("nVertices",&diMuNVertices);		
		
		DiMuonTree->SetBranchAddress("Flag_HBHENoiseFilter",&diMuHBHENoiseFilter);
		DiMuonTree->SetBranchAddress("Flag_HBHENoiseIsoFilter",&diMuHBHENoiseIsoFilter);
		DiMuonTree->SetBranchAddress("Flag_globalTightHalo2016Filter",&diMuCSCTightHaloFilter);
		DiMuonTree->SetBranchAddress("Flag_goodVertices",&diMuGoodVertices);
		DiMuonTree->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter",&diMuEcalDeadCellTriggerPrimitiveFilter);
		DiMuonTree->SetBranchAddress("Flag_eeBadScFilter",&diMuEEBadScFilter);
		
		DiMuonTree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v",&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v);
		DiMuonTree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v",&HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v);
		DiMuonTree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",&HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v);
		DiMuonTree->SetBranchAddress("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v",&HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v);
		DiMuonTree->SetBranchAddress("HLT_Mu27_TkMu8_v",&HLT_Mu27_TkMu8_v);
		DiMuonTree->SetBranchAddress("HLT_Mu30_TkMu11_v",&HLT_Mu30_TkMu11_v);
			
				
		nentries = DiMuonTree->GetEntries();
		if (sampleName.find("MergedData") == string::npos){
			for (Long64_t i = 0; i < nentries; i++){
				DiMuonTree->GetEntry(i);
				leptonFullSimScaleFactor1 = 1.;
				leptonFullSimScaleFactor2 = 1.;
				leptonFullSimScaleFactorErr1 = 0.03;
				leptonFullSimScaleFactorErr2 = 0.03;
				
				//~ leptonFastSimScaleFactor1 = 1.;
				//~ leptonFastSimScaleFactor2 = 1.;
				//~ leptonFastSimScaleFactorErr1 = 0.02;
				//~ leptonFastSimScaleFactorErr2 = 0.02;
				
				
				
				if (diMuNVertices > 45)tempNVertices = 45;
				else tempNVertices = diMuNVertices;
				
				if (diMuPt1 > 120.)tempPt = 119.;
				else tempPt = diMuPt1;
				tempEta = fabs(diMuEta1);
					
				leptonFullSimScaleFactor1 = leptonFullSimScaleFactor1 * muonIDScaleFactorHisto->GetBinContent(muonIDScaleFactorHisto->GetXaxis()->FindBin(tempPt),muonIDScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				leptonFullSimScaleFactor1 = leptonFullSimScaleFactor1 * muonIsoScaleFactorHisto->GetBinContent(muonIsoScaleFactorHisto->GetXaxis()->FindBin(tempPt),muonIsoScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				leptonFullSimScaleFactor1 = leptonFullSimScaleFactor1 * muonIP2DScaleFactorHisto->GetBinContent(muonIP2DScaleFactorHisto->GetXaxis()->FindBin(tempPt),muonIP2DScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				leptonFullSimScaleFactor1 = leptonFullSimScaleFactor1 * muonSIP3DcaleFactorHisto->GetBinContent(muonSIP3DcaleFactorHisto->GetXaxis()->FindBin(tempPt),muonSIP3DcaleFactorHisto->GetYaxis()->FindBin(tempEta));
				
				leptonFullSimScaleFactor1 = leptonFullSimScaleFactor1 * muonTrackScaleFactorEtaHisto->GetBinContent(muonTrackScaleFactorEtaHisto->GetXaxis()->FindBin(tempEta));
				leptonFullSimScaleFactor1 = leptonFullSimScaleFactor1 * muonTrackScaleFactorVtxHisto->GetBinContent(muonTrackScaleFactorVtxHisto->GetXaxis()->FindBin(tempNVertices));
				
				//~ leptonFastSimScaleFactor1 = leptonFastSimScaleFactor1 * FastSimMuonIDScaleFactorHisto->GetBinContent(FastSimMuonIDScaleFactorHisto->GetXaxis()->FindBin(tempPt),FastSimMuonIDScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				//~ leptonFastSimScaleFactor1 = leptonFastSimScaleFactor1 * FastSimMuonIsoScaleFactorHisto->GetBinContent(FastSimMuonIsoScaleFactorHisto->GetXaxis()->FindBin(tempPt),FastSimMuonIsoScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				//~ leptonFastSimScaleFactor1 = leptonFastSimScaleFactor1 * FastSimMuonIP2DScaleFactorHisto->GetBinContent(FastSimMuonIP2DScaleFactorHisto->GetXaxis()->FindBin(tempPt),FastSimMuonIP2DScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				//~ leptonFastSimScaleFactor1 = leptonFastSimScaleFactor1 * FastSimMuonSIP3DScaleFactorHisto->GetBinContent(FastSimMuonSIP3DScaleFactorHisto->GetXaxis()->FindBin(tempPt),FastSimMuonSIP3DScaleFactorHisto->GetYaxis()->FindBin(tempEta));
	
				if (diMuPt2 > 120.)tempPt = 119.;
				else tempPt = diMuPt2;
				tempEta = fabs(diMuEta2);
					
				leptonFullSimScaleFactor2 = leptonFullSimScaleFactor2 * muonIDScaleFactorHisto->GetBinContent(muonIDScaleFactorHisto->GetXaxis()->FindBin(tempPt),muonIDScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				leptonFullSimScaleFactor2 = leptonFullSimScaleFactor2 * muonIsoScaleFactorHisto->GetBinContent(muonIsoScaleFactorHisto->GetXaxis()->FindBin(tempPt),muonIsoScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				leptonFullSimScaleFactor2 = leptonFullSimScaleFactor2 * muonIP2DScaleFactorHisto->GetBinContent(muonIP2DScaleFactorHisto->GetXaxis()->FindBin(tempPt),muonIP2DScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				leptonFullSimScaleFactor2 = leptonFullSimScaleFactor2 * muonSIP3DcaleFactorHisto->GetBinContent(muonSIP3DcaleFactorHisto->GetXaxis()->FindBin(tempPt),muonSIP3DcaleFactorHisto->GetYaxis()->FindBin(tempEta));
				
				leptonFullSimScaleFactor2 = leptonFullSimScaleFactor2 * muonTrackScaleFactorEtaHisto->GetBinContent(muonTrackScaleFactorEtaHisto->GetXaxis()->FindBin(tempEta));
				leptonFullSimScaleFactor2 = leptonFullSimScaleFactor2 * muonTrackScaleFactorVtxHisto->GetBinContent(muonTrackScaleFactorVtxHisto->GetXaxis()->FindBin(tempNVertices));
				
				//~ leptonFastSimScaleFactor2 = leptonFastSimScaleFactor2 * FastSimMuonIDScaleFactorHisto->GetBinContent(FastSimMuonIDScaleFactorHisto->GetXaxis()->FindBin(tempPt),FastSimMuonIDScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				//~ leptonFastSimScaleFactor2 = leptonFastSimScaleFactor2 * FastSimMuonIsoScaleFactorHisto->GetBinContent(FastSimMuonIsoScaleFactorHisto->GetXaxis()->FindBin(tempPt),FastSimMuonIsoScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				//~ leptonFastSimScaleFactor2 = leptonFastSimScaleFactor2 * FastSimMuonIP2DScaleFactorHisto->GetBinContent(FastSimMuonIP2DScaleFactorHisto->GetXaxis()->FindBin(tempPt),FastSimMuonIP2DScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				//~ leptonFastSimScaleFactor2 = leptonFastSimScaleFactor2 * FastSimMuonSIP3DScaleFactorHisto->GetBinContent(FastSimMuonSIP3DScaleFactorHisto->GetXaxis()->FindBin(tempPt),FastSimMuonSIP3DScaleFactorHisto->GetYaxis()->FindBin(tempEta));
				
				
					
				if (diMuHBHENoiseFilter > 0 && diMuHBHENoiseIsoFilter > 0 && diMuCSCTightHaloFilter > 0 && diMuGoodVertices > 0 && diMuEcalDeadCellTriggerPrimitiveFilter > 0) metFilterSummary = 1;
				//~ if (diMuHBHENoiseFilter > 0 && diMuHBHENoiseIsoFilter > 0 && diMuGoodVertices > 0 && diMuEcalDeadCellTriggerPrimitiveFilter > 0) metFilterSummary = 1;
				else metFilterSummary = 0;
				
				if (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v > 0 || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v > 0 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v > 0 || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v > 0 || HLT_Mu27_TkMu8_v > 0 || HLT_Mu30_TkMu11_v > 0 ) triggerSummary = 1;
				else triggerSummary = 0;
				
	
				
				diMuonFullSimScaleFactor1Branch->Fill();
				diMuonFullSimScaleFactor2Branch->Fill();
				diMuonFullSimScaleFactorErr1Branch->Fill();
				diMuonFullSimScaleFactorErr2Branch->Fill();
				diMuonMetFilterBranch->Fill();
				//~ diMuonFastSimScaleFactor1Branch->Fill();
				//~ diMuonFastSimScaleFactor2Branch->Fill();
				//~ diMuonFastSimScaleFactorErr1Branch->Fill();
				//~ diMuonFastSimScaleFactorErr2Branch->Fill();			
				diMuonMetFilterBranch->Fill();
				diMuonTriggerBranch->Fill();
			}
		}
		else{
			leptonFullSimScaleFactor1 = 1.;
			leptonFullSimScaleFactor2 = 1.;		
			leptonFullSimScaleFactorErr1 = 0.;
			leptonFullSimScaleFactorErr2 = 0.;
			
			//~ leptonFastSimScaleFactor1 = 1.;
			//~ leptonFastSimScaleFactor2 = 1.;		
			//~ leptonFastSimScaleFactorErr1 = 0.;
			//~ leptonFastSimScaleFactorErr2 = 0.;
		
			for (Long64_t i = 0; i < nentries; i++){
				DiMuonTree->GetEntry(i);
				
				if (diMuHBHENoiseFilter > 0 && diMuHBHENoiseIsoFilter > 0 && diMuCSCTightHaloFilter > 0 && diMuGoodVertices > 0 && diMuEcalDeadCellTriggerPrimitiveFilter > 0 && diMuEEBadScFilter > 0) metFilterSummary = 1;
				else metFilterSummary = 0;
				
				if (HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v > 0 || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v > 0 || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v > 0 || HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v > 0 || HLT_Mu27_TkMu8_v > 0 || HLT_Mu30_TkMu11_v > 0 ) triggerSummary = 1;
				else triggerSummary = 0;
				
				
				diMuonFullSimScaleFactor1Branch->Fill();
				diMuonFullSimScaleFactor2Branch->Fill();
				diMuonFullSimScaleFactorErr1Branch->Fill();
				diMuonFullSimScaleFactorErr2Branch->Fill();
				//~ diMuonFastSimScaleFactor1Branch->Fill();
				//~ diMuonFastSimScaleFactor2Branch->Fill();
				//~ diMuonFastSimScaleFactorErr1Branch->Fill();
				//~ diMuonFastSimScaleFactorErr2Branch->Fill();	
				diMuonMetFilterBranch->Fill();
				diMuonTriggerBranch->Fill();
				
			}
			
		}
		
		DiMuonTree->Write(0,TObject::kOverwrite);
		delete f;
	}
}
		
		
	

	
