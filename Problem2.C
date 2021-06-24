//Amperiadou Dimitra, AEM:4386

//Problem 2

//We use ROOT framework for drawing and displaying histograms and plots
#include<iostream>
#include<cstdlib>
#include<vector>
#include <math.h>
#include"TROOT.h"
#include"TApplication.h"
#include"TH1F.h"
#include"TCanvas.h"
#include"TGraph.h"
#include"TFile.h"
#include"TLegend.h"
#include"TF1.h"
#define pi M_PI

using namespace std;


int main(){
  int aem=4386,N=20,max=10,min=-10,sumM;
  int sizex=451,sizey=451; //grid size
  int pxc=(sizex-1)/2,pyc=(sizey-1)/2; //position of the central particle
  int a[sizex][sizey];//initializing the grid
  int radius=200; //the radius of the circle
  //pointers to the current position of the particle that performs random walk
  int px,py;
  //the angle for the position on the circle of radius 200
  double xc,yc,L[N/2],M[N][N/2],avgM[N/2],r,r2,theta;
  //variable which tells us weather or not to break the random walk
  bool br=false;
  vector<double>x;//vectors for the plot of dla
  vector<double>y;
  srand(aem);
  //L and M matrix for the calculation of fractal dimension
  for(int j=0; j<N/2;j++){
    L[j]=(j+1)*10;
    for(int i=0;i<N;i++){
      M[i][j]=0;
    }

  }
  

  
  //we need 20 dlas in order to calculate the average of M for each L
  for(int k=0;k<N;k++){
  
    for(int i=0;i<sizex;i++){
      for(int j=0;j<sizey;j++){
	a[i][j]=0;
      }
    }
    a[pxc][pyc]=1; //we place the initial particle in the center of the grid
 
  
  
    //beginning of the process
    while(true){
      //initial position of the particle on the periphery of the circle
      //with center the central particle and radius 200
      theta=(double) rand()*2*pi/(RAND_MAX);
      px=(int) (pxc+radius*sin(theta));
      py=(int) (pyc+radius*cos(theta));
    
      //the particle starts its random walk
      //until it finds a particle to collide
      while(true){
     
	//if particle hits the grid,
	//it goes back to the periphery of the circle
	if(abs(px+1-pxc)>=pxc || abs(py+1-pyc)>=pyc ||
	   abs(px-1-pxc)>=pxc || abs(py-1-pyc)>=pyc){
	  theta=(double) rand()*2*pi/(RAND_MAX);
	  px=(int) (pxc+radius*sin(theta));
	  py=(int) (pyc+radius*cos(theta));
	}
      
	//random steps
	r=(double) rand()/(RAND_MAX);
	if(r<=0.25) px=px+1;

	if(r>0.25 && r<=0.5) px=px-1;
      
    
	if(r>0.5 && r<=0.75) py=py+1;
      
	if(r>0.75) py=py-1;

	//check if there is any neighbor particle
	//we scan all the 8 positions starting by the point (px-1,py-1)
	for(int i=px-1;i<px+2;i++){
	  for(int j=py-1;j<py+2;j++){	 
	    if(a[i][j]==1){//collision happens here
	      br=true;
	      a[px][py]=1;
	  
	    }//if
	  }//for columns  
	}//for rows
     
	if(br==true) {
	  br=false;
	  break; //stop random walk,the particle had the collision
	}
      }//while for random walk

      //the whole process stops when the growing
      //aggregate from the center touches the circumference of the circle

      r2=(px-pxc)*(px-pxc)+(py-pyc)*(py-pyc);
    
      if(r2==radius*radius) break;
    
    }//while for the number of the particles
    //that will perform random walk until the collision

  
    //calculating M matrix for the fractal dimension
    //set a center near the center of the circle
    xc=(double) (max-min)*rand()/(RAND_MAX)-max;
    yc=(double) (max-min)*rand()/(RAND_MAX)-max;
    xc=(int) xc+pxc;
    yc=(int) yc+pyc;
    for(int l=0;l<N/2;l++){
      for(int n=xc-L[l]/2;n<xc+L[l]/2;n++){
	for(int m=yc-L[l]/2;m<yc+L[l]/2;m++){
	  if(a[n][m]==1) M[k][l]+=1;
	}

      }
    }
  


  
    if (k==0){//draw only 1st DLA in order to decrease computational time
      for(int i=0;i<sizex;i++){
	for(int j=0;j<sizey;j++){
	  if(a[i][j]==1){
	    x.push_back(i);
	    y.push_back(j);

	  }

	}

      }
  
   
      TCanvas *c1=new TCanvas("Diffusion Limited Aggregation");
      TGraph *gr=new TGraph(x.size(),&y[0],&x[0]);
      gr->GetXaxis()->SetTitle("");
      gr->GetYaxis()->SetTitle("");
      gr->Draw("A*");
      TFile *file = new TFile("Problem2.root","RECREATE");
      file->WriteTObject(c1);
      file->Close();
  

      delete gr;
      delete c1;
    }

    cout<<"N="<<k;

  }//for N (runs)
  
  

  //find the average of M for each dimension L of the squares
  for(int j=0;j<N/2;j++){
    sumM=0;
    for(int i=0;i<N;i++){
      sumM=sumM+M[i][j];
    }

    avgM[j]=sumM/N;

  }

  for(int j=0;j<N/2;j++){
    L[j]=log(L[j]);
    avgM[j]=log(avgM[j]);
   

 
  }

  TCanvas *c2=new TCanvas("LogM-logL");
  c2->SetGrid();
  TGraph *gr2=new TGraph(N/2,L,avgM);
  gr2->SetMarkerColor(4);
  gr2->SetMarkerSize(1.5);
  gr2->SetMarkerStyle(21);
  gr2->GetXaxis()->SetTitle("logL");
  gr2->GetYaxis()->SetTitle("logM");
  TF1 *fg2 = new TF1("fg2","pol1",L[0],L[N/2-1]); //we use polynomial 1st degree for linear fitting, ROOT uses least squares to calculate parameters
  gr2->Draw("AP");
  gr2->Fit(fg2);
  double po=fg2->GetParameter(0);
  double p1=fg2->GetParameter(1);
  TLegend *leg1 = new TLegend(0.6, 0.8, 0.89, 0.89);
  leg1->AddEntry(fg2, TString::Format("y=%f *x+%f",p1,po),"l");
  leg1->Draw();
  TFile *file = new TFile("Problem2.root","UPDATE");
  file->WriteTObject(c2);
  file->Close();
  
  
  

  
  return 0;
}
  
  

  
  
