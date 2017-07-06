//-----c++------
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <map>

//------root--------
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TRint.h"
#include "TGraph.h"
#include "TMultiGraph.h"

using namespace std;
#define NEUTRONMASS 939565.4133
#define PROTONMASS 938272.0813

double LDM(int A, int Z){
	int coeff;
	if((A-Z)%2==0 && Z%2==0)//even even
		coeff = -1;
	else if((A-Z)%2==1 && Z%2==1)//odd odd
		coeff = 1;
	else//rest
		coeff = 0;

	double t_volume = 15.409*A;
	double t_surface = -16.873*pow(A,2./3.);
	double t_asym = -0.695*Z*(Z-1)/pow(A,1./3.);
	double t_coulomb = -22.435*pow(2*Z-A,2)/A;
	double t_pairing = -11.155/pow(A, 0.5)*coeff;

	return(t_volume+t_surface+t_asym+t_coulomb+t_pairing);
}

int main(int argc, char *argv[]){
	//visualization
	TRint *theApp = new TRint("Rint", &argc, argv);
//*******************************************************************************************
//********************************1st exercise***********************************************
//*******************************************************************************************
	//MASS
	//--read data
/*	ifstream file("aud16.dat");
	map<int,double> BE;
	while ( file.good() ) {
		string oneline;
		getline (file,oneline);
		istringstream is(oneline);
		double gabage,e;
		int z,a;
		if (oneline[0] =='#') continue;
		else {
			is >> z >> a >> e >> gabage >> gabage >> gabage;	
			BE[1000*a+z] = e; 
			//cout << e << endl;
		}
	}
	//--plot data
	cout << BE.size()<< endl;

	TH2F* histo = new TH2F("histo","histo",200,0,200,200,0,200);
	for(map<int,double>::iterator it=BE.begin();it!=BE.end();it++){
		//find the decay products
		cout << it->first << " ";
		bool condition = false;
		int iz = it->first%1000;
		int ia = (it->first - iz)/1000;
		int search_id = (ia -1)*1000+iz-1;
		cout << search_id << " ";
		//check if 1 proton decay is bound or not
		if(BE.find(search_id)!=BE.end()){
			if ((BE.find(search_id)->second - it->second)>0){
				//histo->Fill(ia-iz, iz);
				condition = true;
			}
		}
		search_id = (ia -2)*1000+iz-2;
		cout << search_id << endl;
		//check if 2 proton decay is bound or not
		if(BE.find(search_id)!=BE.end()){
			if ((BE.find(search_id)->second - it->second)<0 && condition == true){
				histo->Fill(ia-iz, iz);
			}
		}
		
	}
	histo->Draw("colz");
*/



//*******************************************************************************************
//********************************2nd exercise***********************************************
//*******************************************************************************************
	//LDM
/*	TH2F* histo = new TH2F("histo","histo",200,0,200,200,0,200);
	double sn[200][200]={0};
	for(int i=0;i<100;i++){//Z
		for(int j=i;j<300;j++){//A
			sn[i][j] = LDM(j,i);
		}
	}
	for(int i=0;i<200;i++){//Z
		bool condition = false;
		for(int j=i;j<200;j++){//A
			//if(sn[i][j]==0 || sn[i][j+1]==0) continue;
			if(sn[i][j]-sn[i][j+1]>0 && condition==false){
				histo->Fill(j-i,i);
				condition = true;
			}
		}
	}
	histo->Draw("colz");
*/




//*******************************************************************************************
//********************************3rd exercise***********************************************
//*******************************************************************************************

	//read radius file
/*	ifstream file("rms13.dat");
	map<int, double> radius;
	while ( file.good() ) {
		string oneline;
		getline (file,oneline);
		istringstream is(oneline);
		double rms,gabage;
		int z,n,id;
		if (oneline[0] =='#') continue;
		else {
			is >> z>>n >> gabage >> rms >> gabage >> gabage;
			id = z*1000+n;
			radius[id] = rms;	
		}
	}
	file.close();
	TGraph* gr[100];
	TMultiGraph* mg = new TMultiGraph();
	int start_isotope = 40;
	for(int ii=start_isotope;ii<80;ii++){
		int counter = 0;
		double x[100], y[100];	
		for(map<int, double>::iterator it=radius.begin();it!=radius.end();it++){
			if(it->first/1000==ii && radius.find(it->first-1)!=radius.end()){
				x[counter] = it->first%1000;
				y[counter] = it->second - radius.find(it->first-1)->second;
				counter++;
			}
		}
		if (counter>0) {
			gr[ii-start_isotope] = new TGraph(counter, x, y);
			mg->Add(gr[ii-start_isotope]);
		}
	}
	mg->Draw("ac");
*/

//*******************************************************************************************
//********************************4th exercise***********************************************
//*******************************************************************************************


	ifstream file("toiee.dat");
	TH1F* histo = new TH1F("histo","histo",200,0,10);
 	//all the levels
	map<int,double> Eexcitation;
	while ( file.good() ) {
		string oneline;
		getline (file,oneline);
		istringstream is(oneline);
		//getline (is,oneline,' ');
		string column = "";
		double gabage,Ex;
		int a,z,J;
		if (oneline[0] =='#') continue;
		else {
			is >> gabage >> a >> z >> J >> gabage >>gabage >> Ex>>column;
			if(J==6||J==8){
				string::iterator it=column.end();
				it--;
				if(*it != '+') continue; // only positive needed
				Eexcitation[a*10000+z*10+J]=Ex;
			}
		}
	}
	file.close();
	for(map<int,double>::iterator it=Eexcitation.begin();it!=Eexcitation.end();it++){
		cout << "here" << endl;
		if(Eexcitation.find(it->first+2)!=Eexcitation.end()){
			int ID = it->first;
			histo->Fill(Eexcitation.find(it->first+2)->second - it->second);
		
		}
	}
	histo->Draw();



	//end of visualization
	theApp->Run(kTRUE);
	delete theApp;
	return 0;




}
