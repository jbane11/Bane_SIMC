



void mc_weighter (string in_name = "", int debug=1){

	if(in_name=="")
	{
		cout <<"Please enter the input root file name(ex: H3_kin7) "<<endl;
		cin >> in_name;
	}


	//Find input file and make a Tree (in_T);
	string in_loc = "/home/jbane/halla_xem/phase_space_fall16/worksim/T2/";
	TFile *in_root = new TFile(Form("%s%s.root",in_loc.c_str(),in_name.c_str()));
	if(in_root==nullptr){ cout<<"Issue with input file ::::" <<endl; return;}
	TTree *in_T = (TTree*)in_root->Get("h1");
	if(debug) cout <<"Input file * "<<in_root<< "  : Tree * "<<in_T <<endl;
	
	//Only use the needed branches to save RAM
	in_T->SetBranchStatus("*",0);
	in_T->SetBranchStatus("hsyptar",1);
	in_T->SetBranchStatus("hsxptar",1);
	in_T->SetBranchStatus("hsdelta",1);
	in_T->SetBranchStatus("hsytar",1);
	in_T->SetBranchStatus("e_recon",1);
	in_T->SetBranchStatus("th_recon",1);
	in_T->SetBranchStatus("z_recon",1);


	//how many events in input tree
	unsigned int in_entries = in_T->GetEntries();
	if(debug) cout << "Number of entries:: " <<in_entries<<endl;




	






	in_root->Close();

	
}//end of code

