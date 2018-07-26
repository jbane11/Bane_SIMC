

void merge(int run){
	TString rootfile="";
	TChain ch("h9040");
	for(int i=0 ; i<=7 ;i++){
		rootfile = Form("rootfiles/mc%d_%d.root",run,i);
		if( access( rootfile, F_OK ) != -1 ){
			cout << i <<endl;
			ch.Add(rootfile);
		}
	}
	ch.Merge(Form("rootfiles/mc%d.root",run));
}



