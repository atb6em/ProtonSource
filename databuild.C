#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include <iostream>
#include <cstring>
#include <string>
#include <cmath>
#include <cstdio>
#include <cstdlib>
using namespace std;
// run via root -l "databuild.C(\"file.txt\",n)"
long long int line_count(TString file);
TFile *databuild(TString datafile = "1_25.txt", Int_t time_mult=1, Int_t get=0, Int_t print=0) {

 typedef struct{
  double A;
  double V;
 }Detector;

  double time; //All double for optimization
  double Fil_A;
  double Fil_V;
  double Fil_cntrl;
  double Bias_V;
  double Bias_cntrl;
  double Temp;
  double Emiss_A;
  double Count_p;
  double Count_h; //H2+
  double mbar;
  double mbar2;
  double magnet_T;
  double magnet_V;
  double HVacc;
  double kPressure;
  double kTemp;
  Detector detect_p;
  Detector detect_h; //H2+

  const Int_t nvars = 21;

  double* arry[nvars]={&time, &Fil_A, &Fil_V, &Fil_cntrl, &Bias_V, &Bias_cntrl, &Temp, &Emiss_A, &Count_p, &Count_h, &mbar, &mbar2, &magnet_T, &magnet_V, &HVacc, &detect_p.V, &detect_h.V, &detect_p.A, &detect_h.A, &kPressure, &kTemp};

  TString filename = datafile.Copy();
  filename.ReplaceAll(".txt","");
  filename.Append(".root");
  printf("%s being created...\n",filename.Data());
	cout<<"time multiplier set as second argument:"<<time_mult<<endl;
  
  TFile hfile(filename,"RECREATE");
  
   TH1F *thist[nvars];

   TTree *tree = new TTree("T","30kV Proton Source Data Tree");
   TTree *ttree = new TTree("T_time","30kV Proton Source Data Time Evolution");   

   TBranch *branch[nvars];
   TBranch *tbranch[nvars];

   char varnames[nvars][25] = {"Time", "Filament A", "Filament V", "Filament Control A",  "Bias V", "Bias Control V", "Temp", "Emission #muA", "Proton Count", "H_{2} Count", "Pressure mbar", "Pressure2 mbar", "Magnetic Field T", "Magnet V", "HV_{acc} V", "Proton Detector V", "H_{2} Detector V", "Proton Detector A", "H_{2} Detector A", "Keller Pressure bar", "Keller Temperature C"};
   
   char vars[nvars][25] = {"Time", "FilamentA", "FilamentV", "Filament_Cntrl",  "BiasV", "Bias_Cntrl", "Temp", "Emission#muA", "ProtonCount", "H2Count", "Pressure_mbar", "Pressure2_mbar", "Mag_FieldT", "MagnetV", "HVacc_V", "Proton_DetectorV", "H2DetectorV", "P_DetectorA", "H2DetectorA", "KellerP_bar", "KellerT_C"};
   long long int lines = line_count(datafile)-1;
   long long int nbins = lines;
   float txmin=0,txmax=lines;
   char *thistname=new char[nvars];
   char *timeevolv=new char[nvars];
   char *branchtyp=new char[nvars];

	 char *timemultS=new char[10];
	 sprintf(timemultS,"%d",time_mult);

   for(int i=0;i<nvars;i++){
    thistname=vars[i];

		if(time_mult==1)sprintf(timeevolv,"%s Time Evolution;Time s;%s",varnames[i],varnames[i]);
		else if(time_mult>0)sprintf(timeevolv,"%s Time Evolution;Time %ss;%s",varnames[i],timemultS,varnames[i]);
		else if(time_mult<0){
			cout<<"time multipler cannot be negative"<<endl;
			exit(EXIT_FAILURE);
		}
     sprintf(branchtyp,"%s/D",vars[i]);

     thist[i] = new TH1F(thistname,timeevolv,nbins,txmin,txmax);
     branch[i]= tree->Branch(vars[i],arry[i],branchtyp);
     tbranch[i]=ttree->Branch(vars[i],"TH1F",&thist[i],32000,0);
   }

  FILE *fp2 = fopen(datafile, "r");
	
	if(!fp2){
		cout<<"file not found"<<endl;
		exit(EXIT_FAILURE);
	}

  char header[10000];
  fgets(header,10000,fp2); //Read Header
  

 if(fp2!=NULL){
  for(long long int i=0; i<lines;i++){
    char buffer[10000];
    fgets(buffer,10000,fp2);
    sscanf(buffer,"%lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg", &time, &Fil_A, &Fil_V, &Fil_cntrl, &Bias_V, &Bias_cntrl, &Temp, &Emiss_A, &Count_p, &Count_h, &mbar, &mbar2, &magnet_T, &magnet_V, &HVacc, &detect_p.V, &detect_h.V, &detect_p.A, &detect_h.A, &kPressure, &kTemp);
		Count_h=Count_h/time_mult;
		Count_p=Count_p/time_mult;
  for(int j=0;j<nvars;j++){
    thist[j]->SetBinContent(i,*arry[j]);
  }
  tree->Fill();
  ttree->Fill();
  }
 }
   if (print){
     tree->Print();
     ttree->Print();
   }

   tree->Write();
   for(int i=0;i<nvars;i++){
	ttree->SetBranchAddress(vars[i],&thist[i]);
	thist[i]->Write();
   }
   //ttree->Write();
   fclose(fp2);
   //hfile.Write();
   hfile.Close();
   return 0;
}

long long int line_count(TString file){	//Counts the number of lines in a file
 FILE *fp;
 fp = fopen(file, "r");
 long long int a, lines=0;
 if(fp!=NULL){
  while(!feof(fp)){
   a = fgetc(fp);
   //printf("%c",a);
   if(a=='\n'){
    lines++;
   }
  }
 }
 //printf("\n");
 fclose(fp);
 return lines;
}
