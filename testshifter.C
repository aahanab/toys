//gently modified from Sara Kurdi's Shifter code, for use in stand-alone handling of distortion maps.
class Shifter {
public:
  Shifter(TString truthfilename, TString correctionfilename);
  TVector3 Shift(TVector3 position);
  TVector3 ShiftForward(TVector3 position); //only shift with forward histogram
  TVector3 ShiftBack(TVector3 position); //
  TFile *forward, *back, *average;
  bool hasTruth, hasCorrection;
  TH3F *hX, *hY, *hZ, *hR, *hPhi, *hXave, *hYave, *hZave, *hRave, *hPhiave, *hXBack, *hYBack, *hZBack;  
};

Shifter::Shifter(TString truthfilename, TString correctionfilename=""){
  //load a 'truth' distortion map and, optionally, a map of a measured correction to those distortions
  //this code is currently set up to load a particular correction map that doesn't have distortions
  // in X,Y, and Z components, but rather only in R, R*Phi, and Z components.
  
  //single event distortion file
  hasTruth=false;//assume the file doesn't load correctly, until we prove otherwise.
  forward=NULL;
  hX=NULL;hY=NULL;hZ=NULL;
  if (truthfilename!=""){
    forward=TFile::Open(truthfilename,"READ"); 
    if (forward!=NULL){
      hX=(TH3F*)forward->Get("hIntDistortionX");
      hY=(TH3F*)forward->Get("hIntDistortionY");
      hZ=(TH3F*)forward->Get("hIntDistortionZ");

      //not strictly needed, but handy:
      hR=(TH3F*)forward->Get("hIntDistortionR");
      hPhi=(TH3F*)forward->Get("hIntDistortionP");

    }
  }
  if (hX!=NULL && hY!=NULL && hZ!=NULL){
    hasTruth=true;
  }

   //single event distortion file
  hasCorrection=false;//assume the file doesn't load correctly, until we prove otherwise.
  average=NULL;
  hX=NULL;hY=NULL;hZ=NULL;
  if (correctionfilename!=""){
    //average=TFile::Open(correctionfilename,"READ");
     average=TFile::Open("/phenix/u/hpereira/sphenix/work/g4simulations/distortion_maps_rec/Distortions_full_realistic_micromegas_all-coarse.root","READ");
    if (forward!=NULL){
      hXave=(TH3F*)forward->Get("hIntDistortionX");
      hYave=(TH3F*)forward->Get("hIntDistortionY");
      hZave=(TH3F*)forward->Get("hIntDistortionZ");
       hZave=(TH3F*)average->Get("hIntDistortionZ");
  
  hRave=(TH3F*)average->Get("hIntDistortionR");
  hPhiave=(TH3F*)average->Get("hIntDistortionP");
    }
  }
  if (hX!=NULL && hY!=NULL && hZ!=NULL){
    hasCorrection=true;
  }

}

TVector3 Shifter::ShiftForward(TVector3 position){
double x, y, z, xshift, yshift, zshift;
  const double mm = 1.0;
  const double cm = 10.0;
  TVector3 shiftposition;

  x= position.X();
  y= position.Y();
  z= position.Z();

  double r=position.Perp();
  double phi=position.Phi();
  if(position.Phi() < 0.0){
    phi = position.Phi() + 2.0*TMath::Pi(); 
  }

  //distort coordinate of stripe
  xshift=0;
  yshift=0;
  zshift=0;
  if (hasTruth){
    xshift=hX->Interpolate(phi,r,z);
    yshift=hY->Interpolate(phi,r,z);
    zshift=hZ->Interpolate(phi,r,z);
  }
  
  //remove average distortion
  if (hasCorrection){
    double raveshift=hRave->Interpolate(phi,r,z);
    double paveshift=hPhiave->Interpolate(phi,r,z); //hugo confirms the units are cm
    double cosphi=cos(phi);
    double sinphi=sin(phi);
    xshift-=raveshift*cosphi-paveshift*sinphi;
    yshift-=raveshift*sinphi+paveshift*cosphi;
    
    zshift-=hZave->Interpolate(phi,r,z);
  }

  TVector3 forwardshift(x+xshift,y+yshift,z+zshift);

  return forwardshift;
}

TVector3 Shifter::ShiftBack(TVector3 forwardshift){
double x, y, z, xshift, yshift, zshift;
  const double mm = 1.0;
  const double cm = 10.0;
  TVector3 shiftposition;

  x= forwardshift.X();
  y= forwardshift.Y();
  z= forwardshift.Z();

  double rforward=forwardshift.Perp();
  double phiforward=forwardshift.Phi();
  if(forwardshift.Phi() < 0.0){
    phiforward += 2.0*TMath::Pi();
  }
  
  double xshiftback=-1*hXBack->Interpolate(phiforward,rforward,z);
  double yshiftback=-1*hYBack->Interpolate(phiforward,rforward,z);
  double zshiftback=-1*hZBack->Interpolate(phiforward,rforward,z);
    
  shiftposition.SetXYZ(x+xshiftback,y+yshiftback,z+zshiftback);

  return shiftposition;
}

TVector3 Shifter::Shift(TVector3 position){
  
  return ShiftBack(ShiftForward(position));
}



 void testshifter() {
   float phi, r , z;
   r=20.0;
  
  Shifter *test= new Shifter ("/sphenix/user/rcorliss/distortion_maps/2021.04/apr07.average.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root
");
  TH1F * hXposition= new TH1F("hXposition","Shift in x Position", 100,0,80.0);// 100 bins, from 0.0 to 80.0
 // TH1F * hRposition= new TH1F("hRposition","Shift in Radial Posiition" , 100, 20.0, 30.0); 

   for( r=20.0; r<=30; r++){ 
   TVector3 pos(r, 10, 50);
   TVector3 newpos=test->ShiftForward(pos);
   float deltax= newpos.X()-pos.X();
  //float deltar=newpos.X()-pos.X();
   hXposition->Fill(deltax);
    // hRposition->Fill(deltar);
   }
   
  hXposition->Draw();
 //hRposition->Draw(); 
   
   
   
   
  } //func
