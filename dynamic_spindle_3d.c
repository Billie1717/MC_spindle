/******************************************************/
/******* FLUID MEMBRANE    DT                  ********/
/******* July 1 2003                           ********/
/******************************************************/
 
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//#define L 80  // it must be an even number

#define MCSTEPS 1000 //2000
#define WRITE 1 //5000
#define WRITE_CONF 1 //5000

int seed;
double alpha;
double beta;
double r_nuc;
int counter;
int N_max;
int N_init;
double Pbr0, D;
double r_cut, antangle;
double xCr,yCr,zCr;
double r_circ;
double SpindDist;
double step;
double L;
int counter;
typedef struct{

  double sidex,sidey;
  int N;
  int Nfil;
  int Ntips;
} SYSTEM;

typedef struct{
  double x,y,z;  
  double theta;
  double phi;
  int fil;
  int tip;
 // int start;
  int Nverl;
} PARTICLE;

typedef struct{
  int number;
} FILAMENT;


SYSTEM S;
PARTICLE *Elem;
FILAMENT *Fil;

char readname[100];

char Roname[100];
char Prob_name[100];

//FILE *read,*wr1,*wr2,*prob;

void set_all(void);
void dissolve_fil(void);
int Interaction_1p(int);
int Interaction_alpha(int);
int Interaction_antpar(int);
int Interaction_antparINFO(int);
int Annihilate(int);
void add_mon(void);
int find_tip(int);
void MC_fil(int, int);
void MC_branch(int);
void painter(void);
double magnitude(double,double,double);
double normalRandom(double, double);
void crossProduct(double,double,double,double,double,double);
int liesInCylinder(double,double,double,double,double,double,double,double,double);
double lattice;

main(int argc, char* argv[]){
  int i,j,kk;
  double rnad;
  FILE *in0;
  counter = 0;
  step = 1;
  L = 500;
  N_max = 50000;
  N_init = 6;
  //r_circ = 5;
//SpindDist = 25.;
  //Pbr0 = 0.1;//0.1; //0.001;
  D = 10;
  in0=fopen("in.ves3d","r");
  fscanf(in0,"%lf %lf %lf %lf %lf %lf %lf %lf %d\n",&alpha,&beta,&r_nuc,&Pbr0,&SpindDist,&r_circ,&r_cut,&antangle,&seed); 
  printf("Parameters alpha:%lf beta:%lf rnuc:%lf Pbr0:%lf SpindDist:%lf rcirc:%lf rcut:%lf angle:%lf seed:%d\n",alpha,beta,r_nuc,Pbr0,SpindDist,r_circ,r_cut,antangle,seed); 
  set_all();
  
   for (i=1;i<=MCSTEPS;i++){painter(); //This is where the loop is {{if (j==S.Nfil)printf("Nfils: %d, Ns: %d \n",S.Nfil, S.N);}
   for (j=1;j<=S.Nfil;j++) {rnad = drand48();MC_fil(1+(int)(rnad*S.Nfil), i);} //pick a random filament printf("Filament %d",1+(int)(drand48()*S.Nfil));
   for (kk=1;kk<=S.N;kk++) {MC_branch(1+(int)(drand48()*S.N));} //pick a random bead
       //N is keeping track of how many particles are in our system
    //if (i == 1 || i%WRITE_CONF==1 ) {painter();} //&& i>10
  }
//if (i>46)printf("Nfil chosen: %d, Ns: %d \n",1+(int)(rnad*S.Nfil), S.N);
  //write();
  //painter3();
}

/*---------------------------*/
void set_all(void){
  int i, t;
  double dN_init;
  srand48(seed);
    
  S.sidex=(double)(L+2);
  S.sidey=(double)(L+2);
  dN_init = (double)N_init;
  S.N = N_init*2;
  S.Ntips =  N_init*2;
  S.Nfil = N_init*2;
  Elem=(PARTICLE *) malloc((N_max)*sizeof(PARTICLE));
  Fil=(FILAMENT *) malloc((N_max)*sizeof(FILAMENT));
  for (t=1;t<=N_init;t++){ 
      
      Elem[t].fil = t; //we start with N_init separate filaments
      if (t < 5){
          Elem[t].theta = (2*M_PI/4)*t;
          Elem[t].phi = M_PI/2;}//with polarity between 0 and pi
      if (t > 4){
          Elem[t].theta = 0;
          Elem[t].phi = (M_PI)*(t-5);}
      //printf("t: %d, phi: %.1f, theta: %.1f",t, Elem[t].theta, Elem[t].phi);
      Elem[t].x = r_circ*cos(Elem[t].theta)*sin(Elem[t].phi)-SpindDist/2;
      Elem[t].y = r_circ*sin(Elem[t].theta)*sin(Elem[t].phi);
      Elem[t].z = r_circ*cos(Elem[t].phi);
      Elem[t].tip = 1; //start with all being tips
      //Elem[t].start = 1; //start with them all being starts
      //S.Ntips+=1;
      Fil[Elem[t].fil].number = 1; //the number of particles in filament 'Elem[t].fil'
  }
  for (t=N_init+1;t<=N_init*2;t++){ 
      if (t < 11){
          Elem[t].theta = (2*M_PI/4)*(t-6);
          Elem[t].phi = M_PI/2;}//with polarity between 0 and 2pi
      if (t > 10){
          Elem[t].theta = 0;
          Elem[t].phi = (M_PI)*(t-11);}
      Elem[t].fil = t; //we start with N_init separate filaments
      //Elem[t].theta = (2*M_PI/dN_init)*t; //with polarity between 0 and 2pi
      Elem[t].x = r_circ*cos(Elem[t].theta)*sin(Elem[t].phi)+SpindDist/2;
      Elem[t].y = r_circ*sin(Elem[t].theta)*sin(Elem[t].phi);
      Elem[t].z = r_circ*cos(Elem[t].phi);
      Elem[t].tip = 1; //start with all being tips
      //Elem[t].start = 1; //start with them all being starts
      //S.Ntips+=1;
      Fil[Elem[t].fil].number = 1; //the number of particles in filament 'Elem[t].fil'
  }
  //S.Nfil = t;
	
}

/*---------------------------*/  
void MC_fil(int ff, int MCstep){
  int tt,nn,nnn,nnnn, touch;
  double newtheta, newphi;
  double mu, sigma, rrr, rr;
  double s;
  int touchAntPar;
  int h;
  counter = 0;
  rrr = drand48();
  rr = drand48();
  nn = S.N+1; //the new counter atom 
  tt = find_tip(ff);
  
    

  if (abs(Fil[ff].number)>3000){
    printf("Filament blown up in MC_fil ff %d, is size %d",ff,Fil[ff].number);exit(-1);}
  if (tt==0){
    printf("Not found tip  =  %d,  %d",Elem[tt].fil,ff);return;}//exit(-1);}
 
  touchAntPar=Interaction_antpar(tt);
      if (touchAntPar > 0){
          //printf("got here, antipar destroy. ");
          //if (Elem[tt].x>SpindDist){
             // h = Fil[Elem[tt].fil].number-1; // height of cylinder in units of steps
  //Fil[Elem[nn].fil].number =2
              //x1 = Elem[tt].x-h*cos(Elem[tt].theta)*sin(Elem[tt].phi); //getting the other end of the filament coords
              //y1 = Elem[tt].y-h*sin(Elem[tt].theta)*sin(Elem[tt].phi);
              //z1 = Elem[tt].z-h*cos(Elem[tt].phi);
              //theta1 = fmod(Elem[tt].theta*(180./M_PI),360.);
              //phi1 = fmod(Elem[tt].phi*(180./M_PI),360.);
              //Interaction_antparINFO(tt);
            //  printf("tip %d filament %d with x,y,z: %.1f, %.1f, %.1f, theta,phi: %.1f, %.1f, h: %d \n",tt, Elem[tt].fil, Elem[tt].x, Elem[tt].y, Elem[tt].z, Elem[tt].theta, Elem[tt].phi, h);
             // printf("ALT tip %d filament %d with x,y,z: %.1f, %.1f, %.1f, theta,phi: %.1f, %.1f, h: %d \n",tt, Elem[tt].fil, Elem[tt].x, Elem[tt].y, Elem[tt].z, Elem[tt].theta, Elem[tt].phi, h);
            //  printf("end with x,y,z: %.1f, %.1f, %.1f \n",Elem[tt].x-h*cos(Elem[tt].theta)*sin(Elem[tt].phi), Elem[tt].y-h*sin(Elem[tt].theta)*sin(Elem[tt].phi), Elem[tt].z-h*cos(Elem[tt].phi));
             // exit(-1);}
          //printf("got to antipar annhiliation");
          //if (ff == 897){printf("fil to blow anti parallel anhiliated and is size = %d, num fils: %d\n",Fil[ff].number,S.Nfil);}
          Annihilate(ff);Annihilate(Elem[touchAntPar].fil);return;}
  if (rrr<(r_nuc/(r_nuc+alpha+beta))){ //nucleation
      newtheta = drand48()*(2*M_PI); //at the moment, just picking a random angle between 0 and 2pi
      newphi  = drand48()*(M_PI);
      
      
      if (rr>0.5){
          s = -1;}
      if (rr<0.5){
          s = 1;}

      Elem[nn].x = r_circ*cos(newtheta)*sin(newphi)+s*SpindDist/2;
      Elem[nn].y = r_circ*sin(newtheta)*sin(newphi);
      Elem[nn].z = r_circ*cos(newphi);

      touch=Interaction_1p(nn);
      if (touch!=1){  // accept
      //Elem[nn].x = r_circ*cos(newtheta);
      //Elem[nn].y = r_circ*sin(newtheta);
      //Elem[nn].z = 0;
      Elem[nn].fil =  S.Nfil+1; //making new filament
      //if (ff == 897){printf("fil 897 nucleated and is size = %d, num fils: %d\n",Fil[ff].number,S.Nfil);}
      //if (ff == 55){printf("fil 55 nucleated and is size = %d, num fils: %d\n",Fil[ff].number,S.Nfil);}
      Elem[nn].theta = newtheta;
      Elem[nn].phi = newphi;
      Elem[nn].tip = 1;
  
      
      Fil[Elem[nn].fil].number =1; //the number of particles in this new filament is 1
      S.Ntips+=1;
      S.N =S.N+1;
      S.Nfil +=1;
      }
      }     
  if (rrr>(r_nuc/(r_nuc+alpha+beta)) && rrr<((alpha+r_nuc)/(alpha+beta+r_nuc))){ //adding one monomer to the filament
      Elem[nn].x = Elem[tt].x+step*cos(Elem[tt].theta)*sin(Elem[tt].phi);
      Elem[nn].y = Elem[tt].y+step*sin(Elem[tt].theta)*sin(Elem[tt].phi);
      Elem[nn].z = Elem[tt].z+step*cos(Elem[tt].phi);
      touch=Interaction_1p(nn); //=Interaction_alpha(nn); 
      if (touch != 1){
      //if (ff == 897){printf("fil to blow elongated and is size = %d, num fils: %d\n",Fil[ff].number,S.Nfil);}
      Elem[nn].fil =  Elem[tt].fil;
      Elem[nn].theta = Elem[tt].theta;
      Elem[nn].phi = Elem[tt].phi;
      Elem[nn].tip = 1;
      Elem[tt].tip = 0;
      S.N += 1;
      
      Fil[Elem[tt].fil].number+=1;
          }
      } 
      //return ;
  if (rrr>(alpha+r_nuc)/(alpha+beta+r_nuc)) {
      //if (ff == 897){printf("fil to blow anhiliated and is size = %d, num fils: %d\n",Fil[ff].number,S.Nfil);}
      //if (MCstep >=86){
         // counter+=1;}
      Annihilate(ff);
  }
  //if (rrr<(alpha+r_nuc)/(alpha+beta+r_nuc)) {
  
      
  //}
 
    
  if (fabs(Elem[nn].x)>L){
    printf("Reached edge of box in x-direction Elem[k].x  =  %.2f",Elem[nn].x);exit(-1);}
  if (fabs(Elem[nn].y)>L){
    printf("Reached edge of box in y-direction");exit(-1);}
  if (S.Ntips==0){
    printf("No tips left!");exit(-1);}
return ;
}

void MC_branch(int oo){
  int nn,nnn,touch,rr;
  double distCentre,newtheta1, newphi1;
  double Pbranch;
  double normrand, normrand2;
  double mu, sigma;
  double minAng;
  nn = S.N+1; //the new counter atom
  nnn = S.N+2;
  if (abs(Fil[Elem[oo].fil].number)>3000){
    printf("Filament blown up in MC_branch oo %d, is size %d",oo,Fil[Elem[oo].fil].number);exit(-1);}
  if (Elem[oo].x<0){
      distCentre  = sqrt((Elem[oo].x+SpindDist/2)*(Elem[oo].x+SpindDist/2)+Elem[oo].y*Elem[oo].y+Elem[oo].z*Elem[oo].z);}
  if (Elem[oo].x>0){
      distCentre  = sqrt((Elem[oo].x-SpindDist/2)*(Elem[oo].x-SpindDist/2)+Elem[oo].y*Elem[oo].y+Elem[oo].z*Elem[oo].z);}
  Pbranch = Pbr0*0.3; //*exp(-distCentre/D);
  //printf("branching at oo, %d Elem[oo].tip %d",oo,Elem[oo].tip);
  if (Elem[oo].tip != 1){ //don't branch on a tip
  if (drand48()<(Pbranch/(1+Pbranch))){
      mu = 0;
      sigma = 0.3;
      rr = drand48();
      minAng = M_PI/10;
      if (rr < 0.5){minAng = -M_PI/10;}

      normrand = normalRandom(mu, sigma);
      normrand2 = normalRandom(mu, sigma);
      newtheta1 = Elem[oo].theta+normrand; //here, need to ensure between 0 and 2pi
      newphi1 = Elem[oo].phi+normrand2;
      if (newtheta1<0){newtheta1=newtheta1+2*M_PI;}
      if (newphi1<0){newphi1=newphi1+M_PI;}
      if (newtheta1>2*M_PI){newtheta1=newtheta1-2*M_PI;}
      if (newphi1>M_PI){newphi1=newphi1-M_PI;}
      
      //newtheta = fmod(newtheta,2*M_PI);
      //newphi = fmod(newphi,2*M_PI);
      
      Elem[nn].x = Elem[oo].x+step*cos(newtheta1)*sin(newphi1);
      Elem[nn].y = Elem[oo].y+step*sin(newtheta1)*sin(newphi1);
      Elem[nn].z = Elem[oo].z+step*cos(newphi1);
      //Elem[nn].x = Elem[oo].x+step*cos(newtheta);
      //Elem[nn].y = Elem[oo].y+step*sin(newtheta);
      Elem[nnn].x = Elem[nn].x+step*cos(newtheta1)*sin(newphi1);
      Elem[nnn].y = Elem[nn].y+step*sin(newtheta1)*sin(newphi1);
      Elem[nnn].z = Elem[nn].z+step*cos(newphi1);
      //Elem[nnn].x = Elem[nn].x+step*cos(newtheta);
      //Elem[nnn].y = Elem[nn].y+step*sin(newtheta);
      //Elem[nnn].z = 0;
      
      touch=Interaction_1p(nnn); 
      //if (fabs(newtheta-Elem[oo].theta)>M_PI/10){
      if (touch!=1){  // reject, once we have added this tip to this filament, if it is too close to something else, inactivate it
        //Elem[nn].tip = 0;
        //S.Ntips-=1;
         //return ;
      
      Elem[nn].fil =  S.Nfil+1; //making new filament
      Elem[nn].theta = newtheta1;
      Elem[nn].phi = newphi1;
      Elem[nn].tip = 0;
      
      Elem[nnn].fil =  Elem[nn].fil; //making new filament
      Elem[nnn].theta = Elem[nn].theta;
      Elem[nnn].phi = Elem[nn].phi;
      Elem[nnn].tip = 1;
     // Elem[nn].start = 1;
      Fil[Elem[nn].fil].number =2; //the number of particles in this new filament is 1 
      //if (Elem[oo].fil == 897){printf("fil to blow as a mother branched and is size = %d, num fils: %d\n",Fil[Elem[oo].fil].number,S.Nfil);}
      //if (Elem[nn].fil == 897){printf("fil to blow is daughter branch and is size = %d, num fils: %d\n",Fil[Elem[oo].fil].number,S.Nfil);}
      S.Ntips+=1;
      S.N += 2;
      S.Nfil +=1;
      }
      return;
  }
  }
}

int Interaction_antpar(int tt){ //input is the tip
  int j;
  double x1,y1,z1,cyltouch,x2,y2,z2;
  double xv,yv,zv;
  double h;
  double theta1, theta2, phi1, phi2;
  double dthet, dphi;
  h = Fil[Elem[tt].fil].number-1; // height of cylinder in units of steps
  x1 = Elem[tt].x-h*cos(Elem[tt].theta)*sin(Elem[tt].phi); //getting the other end of the filament coords
  y1 = Elem[tt].y-h*sin(Elem[tt].theta)*sin(Elem[tt].phi);
  z1 = Elem[tt].z-h*cos(Elem[tt].phi);
  
  theta1 = Elem[tt].theta;
  phi1 = Elem[tt].phi;

      
  for (j=1;j<=S.N;j++){ //Elem[k].Nverl
      if (Elem[j].fil != Elem[tt].fil){ 
          
         
          theta2 = Elem[j].theta;
          phi2 = Elem[j].phi;
          dthet = fabs(theta1-theta2)*(180/M_PI);
          dphi = fabs(phi1-phi2)*(180/M_PI);
          if (dthet >180){
              dthet  = 360 - dthet;}

          if (dthet > antangle){
              if (dphi > antangle/4){
              
                cyltouch = liesInCylinder(Elem[tt].x,Elem[tt].y,Elem[tt].z,x1,y1,z1,Elem[j].x,Elem[j].y,Elem[j].z);
                if (cyltouch == 1){
                    

                 return j;  // contact-reject
              }
              }
          }
      }
       
  }  

  return 0;
}

/*int Interaction_antparINFO(int tt){ //input is the tip
  int j;
  double x1,y1,z1,cyltouch,x2,y2,z2;
  double xv,yv,zv;
  double h;
  double theta1, theta2, phi1, phi2;
  double dthet, dphi;
  h = Fil[Elem[tt].fil].number-1; // height of cylinder in units of steps
  x1 = Elem[tt].x-h*cos(Elem[tt].theta)*sin(Elem[tt].phi); //getting the other end of the filament coords
  y1 = Elem[tt].y-h*sin(Elem[tt].theta)*sin(Elem[tt].phi);
  z1 = Elem[tt].z-h*cos(Elem[tt].phi);
  xv = Elem[tt].x - x1;
  yv = Elem[tt].y - y1;
  zv = Elem[tt].z - z1;
  
  theta1 = Elem[tt].theta;
  phi1 = Elem[tt].phi;
  theta1 = atan(yv/xv); //fmod(Elem[tt].theta*(180./M_PI),360.);
  if (xv<0){
      theta1 = theta1+M_PI;
      if (yv<0){
          theta1 = theta1+M_PI;
      }
  }
  if (theta1<0){
      theta1 = theta1+M_PI*2;}
  if (theta1>M_PI*2){
      theta1 = theta1-M_PI*2;}
  if (phi1<0){
      phi1 = phi1+M_PI;}
  if (phi1>M_PI){
      phi1 = phi1-M_PI;}
    
    
  theta1 = theta1 * (180/M_PI);
  phi1 =  acos(zv/h)*180/M_PI;  //fmod(Elem[tt].phi*(180./M_PI),360.);
      
  for (j=1;j<=S.N;j++){ //Elem[k].Nverl
      if (Elem[j].fil != Elem[tt].fil){ 
          
          h = Fil[Elem[j].fil].number-1; // height of cylinder in units of steps
  //Fil[Elem[nn].fil].number =2
          x2 = Elem[j].x-h*cos(Elem[j].theta)*sin(Elem[j].phi); //getting the other end of the filament coords
          y2 = Elem[j].y-h*sin(Elem[j].theta)*sin(Elem[j].phi);
          z2 = Elem[j].z-h*cos(Elem[j].phi);
          xv = Elem[j].x - x1;
          yv = Elem[j].y - y1;
          zv = Elem[j].z - z1;
          
          theta2 = atan(yv/xv); //fmod(Elem[tt].theta*(180./M_PI),360.);
          if (xv<0){
              theta2 = theta2+M_PI;
              if (yv<0){
                  theta2 = theta2+M_PI;
              }
          }
          if (theta2<0){
              theta2 = theta2+M_PI*2;}
          if (theta2>M_PI*2){
              theta2 = theta2-M_PI*2;}
          if (phi2<0){
              phi2 = phi2+M_PI;}
          if (phi2>M_PI){
              phi2 = phi2-M_PI;}
          
          theta2 = theta2 * (180/M_PI);
          phi2 =  acos(zv/h)*180/M_PI;  //fmod(Elem[tt].phi*(180./M_PI),360.);
          
          
          theta2 = Elem[j].theta;
          phi2 = Elem[j].phi;
          dthet = fabs(theta1-theta2)*(180/M_PI);
          dphi = fabs(phi1-phi2)*(180/M_PI);
          if (dthet >180){
              dthet  = 360 - dthet;}
          
          //dthet = fabs(fmod(theta1-theta2,360));
          //dphi = fabs(fmod(theta1-theta2,360));
          if (dthet > antangle){
              if (dphi > antangle/2){
              
                cyltouch = liesInCylinder(Elem[tt].x,Elem[tt].y,Elem[tt].z,x1,y1,z1,Elem[j].x,Elem[j].y,Elem[j].z);
                if (cyltouch == 1){
                    printf("INSIDE tip %d filament %d theta,phi: %.1f, %.1f, ALT theta,phi: %.1f, %.1f  \n",tt, Elem[tt].fil, Elem[tt].theta*(180/M_PI), Elem[tt].phi*(180/M_PI),theta1, phi1);
                    printf("INSIDE crossed %d filament %d theta,phi: %.1f, %.1f,ALT theta,phi: %.1f, %.1f  \n",j, Elem[j].fil, Elem[j].theta*(180/M_PI), Elem[j].phi*(180/M_PI), theta2, phi2);

                 return 1;  // contact-reject
              }
              }
          }
      }
       
  }  

  return 0;
}*/


int Interaction_1p(int k){
  int j;
  double dx,dy,dz,r2,rs,dd;
  
  for (j=1;j<=S.N;j++){ //Elem[k].Nverl
      if (j != k){
        dx=Elem[k].x-Elem[j].x;
        dy=Elem[k].y-Elem[j].y;
        dz=Elem[k].z-Elem[j].z;
   
        r2=(dx*dx+ dy*dy+ dz*dz);
   
    if (r2<0.99) return 1;  // contact-reject
      } 
  }  

  return 0;
}

/*int Interaction_alpha(int k){
  int j;
  double dx,dy,dz,r2,rs,dd;
  
  for (j=1;j<=S.N;j++){ //regecting the next one along instead of that alpha
      if (j != k){
        dx=Elem[k].x+step*cos(Elem[k].theta)-Elem[j].x;
        dy=Elem[k].y+step*cos(Elem[k].theta)-Elem[j].y;
        dz=Elem[k].z+step*cos(Elem[k].theta)-Elem[j].z;
   
        r2=(dx*dx+ dy*dy+ dz*dz);
   
    if (r2<0.99) return 1;  // contact-reject
      } 
  }  

  return 0;
}*/

int Annihilate(int k){
    
  int j, jj, b, bb,bbb;
  int fixedFilno;
  int lastindex;
  int index, anfil;
  //printf("Annihilate")
  //k is now the label of the filament to be deleted
  fixedFilno = Fil[k].number; //number of particles in the filament which oo is in
  //if (k == 897){printf("inside annhiliate size = %d, num fils: %d\n",Fil[k].number,S.Nfil);}
  anfil = k; //Elem[k].fil; //label of the filament we are annhilating
  //printf("Annihilate fil  =  %d, nmon =  %d\n",anfil,fixedFilno);
  if (Fil[k].number==0){printf("Number in filament is 0--> S/T WRONG");exit(-1);}
// we want to eliminate a total of Fil[Elem[k].fil].number particles, all of which have the label Elem[k].fil
// we are going to replace the coordinates of these with the other atoms in our system
// (may not work if there are more atoms in that filament than there are in all the other filaments combined)
  bb = 0;
  bbb=0;
  lastindex = 1;
  for (j=1;j<=S.N;j++){
      if (Elem[j].fil == anfil){ 
          bbb+=1;
      }
  }
  for (b=0;b<S.N;b++){ //b is the index at the end of our list (S.N long) that we are moving to fill the gap left by disappearing filament element
      if (Elem[S.N-b].fil != anfil){ //make sure the element we are moving isn't part of this filament
          for (jj=lastindex;jj<=S.N;jj++){
              //if (bb == fixedFilno){break;} //we have been through enough bbs and 
              if (Elem[jj].fil == anfil){ //getting the indices of the ones we are taking out
                  Elem[jj].x=Elem[S.N-b].x; //indexes[j] & S.N+1-j --------->S.N-b will NOT be in the filament
                  Elem[jj].y=Elem[S.N-b].y;
                  Elem[jj].z=Elem[S.N-b].z;
                  Elem[jj].fil=Elem[S.N-b].fil;
                  Elem[jj].theta=Elem[S.N-b].theta;
                  Elem[jj].phi=Elem[S.N-b].phi;
                  Elem[jj].tip=Elem[S.N-b].tip;
                //  Elem[jj].start=Elem[S.N-b].start;
                  bb+=1; //how many beads have we replaced
                  lastindex = jj; //next time we will start from one along jj, otherwise we will 
                  break;} //every time this IF is satisfied, we need to move bb along
              
            } 
          }
  }

  
  if (bb!=fixedFilno){ printf("bb= %d bbb= %d fixedFilno= %d\n filament index= %d\n", bb,bbb, fixedFilno, anfil);} //exit(-1);//}
  for (j=1;j<=S.N;j++){
      if (Elem[j].fil == S.Nfil){ //replacing the beads in the last index filament with the index of the one that's just been destroyed
          Elem[j].fil = anfil;
          Fil[anfil].number = Fil[S.Nfil].number;
      }
  }
  //if (anfil == 119){printf("inside annhiliate size = %d, num fils: %d\n",Fil[anfil].number,S.Nfil);}
  S.N = S.N-bb;
  S.Nfil = S.Nfil-1;
  S.Ntips=S.Ntips-1;
    
  return 0;
}

/*--------------------------------------------------*/ 
void painter(void){
  int p;
  char iiiii[100];
  char iiii[100];
  FILE *o, *h;
  //printf("S.N: %d, Elem[S.N].x, %.2f, Elem[S.N].y, %.2f\n",S.N, Elem[S.N].x ,Elem[S.N].y) ;
  //printf("S.N-1: %d, Elem[S.N-1].x, %.2f, Elem[S.N-1].y, %.2f\n",S.N-1, Elem[S.N-1].x ,Elem[S.N-1].y) ;
  sprintf(iiiii,"out3D_CN_Beta_%.2f_Alpha_%.2f_rnuc_%.2f_dSpind_%.2f_rcut_%.2f_antang_%.1f_seed_%d.xyz",beta,alpha,r_nuc, SpindDist, r_cut, antangle,seed);

  o=fopen(iiiii,"a+");

  fprintf(o,"%d\n",(S.N));
  fprintf(o,"Atoms\n");
  
 for (p=1;p<=S.N;p++){ //printf("p=%d N=%d\n",p,Elem[p].Ng);
 
   //if (p<=S.N) fprintf(o,"%d %f %f %d %d %d %d\n",
	//		     Elem[p].fil,Elem[p].x,Elem[p].y,Elem[p].z,Fil[Elem[p].fil].number,S.Nfil,S.N,Elem[p].tip);
    //else fprintf(o,"%d %f %f %d %d %d %d\n", Elem[p].fil,Elem[p].x,Elem[p].y,Elem[p].z,Fil[Elem[p].fil].number,S.Nfil,S.N,Elem[p].tip);
  // if (p<=S.N) fprintf(o,"%d %f %f %f %f %f %d\n",
//			     Elem[p].fil,Elem[p].x,Elem[p].y,Elem[p].z,Elem[p].theta,Elem[p].phi,Fil[Elem[p].fil].number);
   // else fprintf(o,"%d %f %f %f %f %f %d\n", Elem[p].fil,Elem[p].x,Elem[p].y,Elem[p].z,Elem[p].theta,Elem[p].phi,Fil[Elem[p].fil].number);
    if (p<=S.N) fprintf(o," %d %f %f %f %f %f %d %d\n",
			     Elem[p].fil,Elem[p].x,Elem[p].y,Elem[p].z,Elem[p].theta*180/M_PI,Elem[p].phi*180/M_PI,Fil[Elem[p].fil].number, counter);
    else fprintf(o," %d %f %f %f %f %f %d %d\n", Elem[p].fil,Elem[p].x,Elem[p].y,Elem[p].z,Elem[p].theta*180/M_PI,Elem[p].phi*180/M_PI,Fil[Elem[p].fil].number, counter);
       
    
  }

  fclose(o);
    
    
  /*sprintf(iiii,"FilNumbers_Beta_%.2f_Alpha_%.2f.dat",beta,alpha);

  h=fopen(iiii,"a+");

  //fprintf(h,"%d\n",(S.N));
  //fprintf(h,"Atoms\n");
  
 for (p=1;p<=S.Nfil;p++){ //printf("p=%d N=%d\n",p,Elem[p].Ng);
 
   if (p!=S.Nfil)fprintf(h,"%d ",Fil[p].number);
   if (p == S.Nfil)fprintf(h,"%d \n",Fil[p].number);
    //else fprintf(o,"%d %f %f %d %d %d %d\n", Elem[p].fil,Elem[p].x,Elem[p].y,Fil[Elem[p].fil].number,S.Nfil,S.N,Elem[p].tip);
       
    
  }

  fclose(h);*/

}


double normalRandom(double mu, double sigma){
  double U1, U2, W, mult;
  static double X1, X2;
  static int call = 0;
  if (call == 1){
      call = !call;
      return (mu + sigma * (double) X2);
    }
 
  do{
      U1 = -1 + ((double) rand () / RAND_MAX) * 2;
      U2 = -1 + ((double) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
    }
  while (W >= 1 || W == 0);
 
  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;
 
  call = !call;
 
  return (mu + sigma * (double) X1);
}


int find_tip(int ff){
    int b;
    
    for (b=1;b<=S.N;b++){
        if (Elem[b].fil == ff && Elem[b].tip == 1){
            return b;}
    }
    return 0;
}

void crossProduct(double xA,double yA,double zA,double xB,double yB,double zB)
{
    xCr = yA * zB - zA * yB;
    yCr = zA * xB - xA * zB;
    zCr = xA * yB - yA * xB;
    
}

int liesInCylinder(double xA,double yA,double zA,double xB,double yB,double zB,double xp,double yp,double zp)
{
    double mx,my,mz,ey,ex,ez,EE,MM,v3x,v3y,v3z;
    double w1, w2, d;
    double VV;
    
    
    crossProduct(xA,yA,zA,xB,yB,zB);
    mx = xCr;
    my = yCr;
    mz = zCr;
    ex = xB - xA;
    ey = yB - yA;
    ez = zB - zA;
    
    EE = magnitude(ex,ey,ez);
    MM = magnitude(mx,my,mz);
       
    v3x = xp-xA;
    v3y = yp-yA;
    v3z = zp-zA;
    VV = magnitude(v3x,v3y,v3z);
    w1 = v3x*ex/(EE*VV) + v3y*ey/(EE*VV) + v3z*ez/(EE*VV);
    crossProduct(ex,ey,ez,v3x,v3y,v3z);
    d = magnitude(xCr,yCr,zCr)/EE;
    
    v3x = xp-xB;
    v3y = yp-yB;
    v3z = zp-zB;
    VV = magnitude(v3x,v3y,v3z);
    w2 = v3x*ex/(EE*VV) + v3y*ey/(EE*VV) + v3z*ez/(EE*VV);

      
    if (w1>=0 && w2 <=0 && d < r_cut) {return 1;} //we are inside this cylinder!
    
    return 0;
    
    
}

double magnitude(double xA,double yA,double zA){
    
    return xA*xA+yA*yA+zA*zA;
    }