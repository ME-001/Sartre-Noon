void test()
{
    TFile *file = new TFile("example_Pb.root", "READ");
    TTree *tree = dynamic_cast<TTree*>(file->Get("tree"));

    Long64_t nE = tree->GetEntries();
    printf("Number of entries %lld\n", nE);

    std::vector<Double_t> nESum1;
    std::vector<Double_t> nESum2;
    std:;vector<Double_t> nEta;

    TClonesArray *neutrons1 = new TClonesArray("TLorentzVector");
    TClonesArray *neutrons2 = new TClonesArray("TLorentzVector");
    tree->SetBranchAddress("neutrons1", &neutrons1);
    tree->SetBranchAddress("neutrons2", &neutrons2);

    TClonesArray onXn("TLorentzVector");
    TClonesArray XnXn("TLorentzVector");

    TFile *file2 = new TFile("Test.root", "RECREATE");
    TTree *tre = new TTree("myTree", "Test Tree");

    tre->Branch("onXn", &onXn);
    tre->Branch("XnXn", &XnXn);

    for (Long64_t i = 0; i < nE; i++) {
        tree->GetEntry(i);

        Double_t E1 = 0;
        Double_t E2 = 0;

        onXn.Clear();
        XnXn.Clear();

        if (neutrons1) {
            Int_t size1 = neutrons1->GetEntriesFast();
            printf(" Event %lld: number of neutrons1 %d\n", i, size1);

            for (Int_t j = 0; j < size1; j++) {
                TLorentzVector *neutron1 = dynamic_cast<TLorentzVector*>(neutrons1->At(j));
                if (neutron1) {
                    E1 += neutron1->E();
                    nEta.push_back(neutron1->Eta());
                    new (neutrons2->GetEntriesFast() ? XnXn[XnXn.GetEntriesFast()] : onXn[onXn.GetEntriesFast()]) TLorentzVector(neutron1->Px(),neutron1->Py(),neutron1->Pz(),neutron1->E());
                }
            }
        }

        if (neutrons2) {
            Int_t size2 = neutrons2->GetEntriesFast();
            printf(" Event %lld: number of neutrons2 %d\n", i, size2);

            for (Int_t j = 0; j < size2; j++) {
                TLorentzVector *neutron2 = dynamic_cast<TLorentzVector*>(neutrons2->At(j));
                if (neutron2) {
                    E2 += neutron2->E();
                    nEta.push_back(neutron2->Eta());
                    new (neutrons1->GetEntriesFast() ? XnXn[XnXn.GetEntriesFast()] : onXn[onXn.GetEntriesFast()]) TLorentzVector(neutron2->Px(),neutron2->Py(),neutron2->Pz(),neutron2->E());
                }
            }
        }

        nESum1.push_back(E1 * 0.001);
        nESum2.push_back(E2 * 0.001);
        tre->Fill();
    }

    file->Close();
    

    TH1D *TotE1 = new TH1D("TotE1", "Total Neutron Energy in each Event side 1", 1000, 0, 10000);
    TH1D *TotE2 = new TH1D("TotE2", "Total Neutron Energy in each Event side 2", 1000, 0, 10000);

    
    	for(Int_t i=0; i<nESum1.size();i++)
    	{
        TotE1->Fill(nESum1[i]);
    	}
    	for(Int_t i=0; i<nESum2.size();i++)
    	{
        TotE2->Fill(nESum2[i]);
    	}
    	TCanvas *c1 = new TCanvas("c1","TOtal Neutron Energy per Event side 1",1200,800);
    	TotE1->GetXaxis()->SetTitle("Energy(GeV)");
    	TotE1->SetAxisRange(0,5000,"X");
    	TotE1->GetYaxis()->SetTitle("Counts");
    	TotE1->Draw("L");
    	c1->SetLogy();
    	c1->Write();
    	//c1->SaveAs("Total_Neutron_E_2.76_side1_ss.png");
    	
    	
    	TCanvas *c2 = new TCanvas("c2","TOtal Neutron Energy per Event side 2",1200,800);
    	TotE2->GetXaxis()->SetTitle("Energy(GeV)");
    	TotE2->SetAxisRange(0,5000,"X");
    	TotE2->GetYaxis()->SetTitle("Counts");
    	TotE2->Draw("L");
    	c2->SetLogy();
    	c2->Write();
    	//c2->SaveAs("Total_Neutron_E_2.76_side2_ss.png");
    	
    	
    	
    	TH1D *nEta_hist = new TH1D("nEta", "pseudo rapidity", 1000, -30, 30);
    	for(Int_t i=0; i<nEta.size();i++)
    	{
        nEta_hist->Fill(nEta[i]);
    	}
    	TCanvas *c3 = new TCanvas("c3","Pseudo rapidity ",1200,800);
    	nEta_hist->GetXaxis()->SetTitle("Energy(GeV)");
    	//TotE2->SetAxisRange(0,5000,"X");
    	nEta_hist->GetYaxis()->SetTitle("Counts");
    	nEta_hist->Draw("L");
    	c3->SetLogy();
    	c3->Write();
    	//c3->SaveAs("Pseudo_rapidity.png");
    	
    	file2->Close();
}



