
//**********************************************************************Bibliotecas****************************************************************

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <omp.h> 

//*********************************************************************definições***************************************************************

#define mu (1.0) //Mobilidade
#define n (15) //Número de partículas no anel
#define fator ((1.0)/(1.0*n))
#define d0 (1.0) 
#define A0 ((0.25*n*d0*d0)/tan((M_PI)/(1.0*n))) // Área de Equilíbrio
#define plot_gnuplot 200
#define t0 (0.0) //Tempo Inicial
#define tf (10.0) //Tempo Final 1000000
#define COLOR 7 //Cor dos anéis


#if !(_SVID_SOURCE || _XOPEN_SOURCE)
double drand48(void) {
    return rand() / (RAND_MAX + 1.0);
}

long int lrand48(void) {
    return rand();
}

long int mrand48(void) {
    return rand() > RAND_MAX / 2 ? rand() : -rand();
}

void srand48(long int seedval) {
    srand(seedval);
}
#endif

//********************************************************************Estrutura***********************************************************************************************

//Definindo o anel como uma estrutura

typedef struct {

//*************************************************
//                                                * 
//                                                *
double x;      //posição em x                     *
double y;      //posição em y                     *
double vx;     //componente x da velocidade       *  
double vy;     //componente y da velocidade       *
double fx;    //componente x da força             *
double fy;    //componente y da força             *
double nx;     //componente x da autopropulsão    *
double ny;     //componente y da autopropulsão    *        
int color;     //cor da partícula                 * 
double theta;  //orientação da autopropulsão      * 
double x_virtual;
double y_virtual;
double nx_right;
double nx_left;
double ny_up;
double ny_down;
double A; //Área do anel
double x_aux;
double y_aux;
//*************************************************
} particle;

//define RING's struct
typedef struct {

//*********************************************
//                                            * 
//                                            *
particle PART[n+1];       //Anel da partícula *

//*********************************************
} ring;

//********************************************************************Funções**********************************************************************
double box_muller(double *z);
void VIEW_GNUPLOT(ring *RING,double t,double L,double PE,double G,double PHI,double Kadh,int N,double p0,double rc);
void initialization(ring *RING,double L,int N,double rc);
void Motion_evol(ring *RING,int ai,double L,int N,double PE,double G,double dt);
void forces(ring *RING,int **list,int *nviz,int ai,double L,double Kadh,double radh,int N,double Ka,double rc,double Kl,double Kc);
void PBC(ring *RING,int i,int j,double L,double *z1,double *z2,double *z3,double *z4,int N);
void fadh(double dr,double rc,int ai, int aj, int bi, int bj, double *force2,ring *RING,double Kadh);
void frep(double dr,double rc,int ai, int aj, int bi, int bj, double *force1,ring *RING,double Kc);
void updating_list(ring *RING,int ai,double P0,int *nviz,int **list,double L, int N);
void RBC(ring *RING,int i,int j,double L);

//**********************************************************************Programa******************************************************************

int main(int argc, char *argv[]){
  
   int i,j,k,count,count2,samp,ai,aj,bi,bj,N,delay,color;
   double radh,t,L,G,Kl,Ka,Kc,Kadh;
   double PE,p0,PHI,P0,rc,dt;

   //Parâmetros de entrada
   PE = atof(argv[1]); //Intensidade da atividade (Peclet)
   dt = 0.005; //passo de tempo
   G = atof(argv[2]);  //Intensidade do acoplamento
   PHI = atof(argv[3]); //Fração de empacotamento (densidade)
   N = atoi(argv[4]); //Número de anéis
   Ka = atof(argv[5]); //Intensidade da força de área
   Kl = atof(argv[6]); //Intensidade da Constante de mola do anel
   Kc = atof(argv[7]); //Intensidade da força de repulsão
   Kadh = atof(argv[8]); //Intensidade da força de adesão
   p0 = atof(argv[9]); //Tensão cortical
   radh = atof(argv[10]); //Alcance da força
   rc = (p0*sqrt(A0)/(1.0*n));
   samp = 2; //Semente de número aleatório
   srand48(samp);
   P0 = (rc*n); //Perímetro de equilibrio do anel
   L = sqrt((N*(A0+0.125*n*M_PI*rc*rc))/(PHI)); //Tamanho do sistema
   PHI = (N*(A0+0.125*n*M_PI*rc*rc))/(L*L);

   int **list = (int**)malloc(N * sizeof(int*));//indice da j-esima PARTicula na i-esima celula
   int *nviz = (int*)malloc(N * sizeof(int)); //number of PARTicles into box i
   ring *RING = (ring*)malloc((N+1) * sizeof(ring));
   
   for(i=0;i<N;i++){ 
      list[i] = (int*) malloc(N * sizeof(int));
      nviz[i] = 0;
      for(j=0;j<N;j++) list[i][j] = 0;
   }
  
   //Inicialização
   initialization(RING,L,N,rc);
   count = 0;
   
   //Começando a dinâmica
   
   printf("COMECOU\n");
   clock_t begin = clock();
   for(t=t0; t<tf; t+=dt) {
   
      //Gnuplot Movie (caso queira assistir video em tempo real)
      // if((count%plot_gnuplot)==0) VIEW_GNUPLOT(RING,t,L,PE,G,PHI,Kadh,N,p0,rc);
                  

      for(ai=0;ai<(N-1);ai++) {
         
         //Formando a lista de vizinhos
         updating_list(RING,ai,P0,nviz,list,L,N);
    
         //cálculo das forças 
         forces(RING,list,nviz,ai,L,Kadh,radh,N,Ka,rc,Kl,Kc);
      }
       
      ai = N-1;     
      //cálculo das forças 
      forces(RING,list,nviz,ai,L,Kadh,radh,N,Ka,rc,Kl,Kc); 
      RING[N].PART[n].x = 0.0;
      RING[N].PART[n].y = 0.0;  
      RING[N].PART[n].x_virtual = 0.0;
      RING[N].PART[n].y_virtual = 0.0; 
      for(ai=0;ai<N;ai++){
	 	 
	 //Integração das equações de movimento
         Motion_evol(RING,ai,L,N,PE,G,dt);
         //Calculando centro de massa do sistema
	 RING[N].PART[n].x += RING[ai].PART[n].x;
	 RING[N].PART[n].y += RING[ai].PART[n].y;
	 RING[N].PART[n].x_virtual += RING[ai].PART[n].x_virtual;
	 RING[N].PART[n].y_virtual += RING[ai].PART[n].y_virtual;
	 nviz[ai]=0;

	 
	} 
	        
      //Calculando centro de massa do sistema	        	 
      RING[N].PART[n].x = (1.0/(1.0*N))*RING[N].PART[n].x;
      RING[N].PART[n].y = (1.0/(1.0*N))*RING[N].PART[n].y;
      RING[N].PART[n].x_virtual = (1.0/(1.0*N))*RING[N].PART[n].x_virtual;
      RING[N].PART[n].y_virtual = (1.0/(1.0*N))*RING[N].PART[n].y_virtual;
      count++;
      
   }
   printf("TERMINOU: ");
   clock_t end = clock();
   double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
   printf("%f", time_spent);
   return(0);
}   



void initialization(ring *RING,double L,int N,double rc){
  
   int ai,aj,k,count,c,Nc;
   double theta,theta2,dtheta,phi,grau,r,dx,dy,prod,fat,h1,h2,x,y,vx,vy,color,rr,XCM,YCM,x0,y0;
   double Fn_red,Fn_green,nx,ny,norm,alpha,radius,R,ds,li,dx_cm,dy_cm,xj,yj,xcm,ycm;

   Nc = (int)(sqrt(N));
   ds = (1.0*d0);
   R = ((0.5*n*ds)/M_PI);
   dtheta = (ds/R);
   x0 = (R + (0.5*d0) + (0.3*L));
   y0 = (R + (0.5*d0) + (0.3*L));
   theta =  0.0;
   k = 0;
   count = 1;   
   //srand48(time(NULL)); 
   for(ai=0;ai<N;ai++){
      
      

	 
      RING[ai].PART[n].color = COLOR;
      alpha = 2*M_PI*drand48();
      RING[ai].PART[n].theta= alpha;
      RING[ai].PART[n].nx = cos(alpha);
      RING[ai].PART[n].ny = sin(alpha);
      RING[ai].PART[n].vx =  RING[ai].PART[n].nx;
      RING[ai].PART[n].vy =  RING[ai].PART[n].ny;
      RING[ai].PART[n].x = 0.0;
      RING[ai].PART[n].y = 0.0;

      for(aj=0;aj<n;aj++) {
	 
	 RING[ai].PART[aj].x = R*cos(theta)+ x0 + (ai-(ai/Nc)*Nc)*((2*R)+(0.5*rc));
	 RING[ai].PART[aj].y = R*sin(theta)+ y0 + (ai/Nc)*((2*R)+(0.5*rc));
	 RING[ai].PART[aj].vx = RING[ai].PART[n].vx;
	 RING[ai].PART[aj].vy = RING[ai].PART[n].vy;
	 RING[ai].PART[aj].x_virtual = RING[ai].PART[aj].x;
	 RING[ai].PART[aj].y_virtual = RING[ai].PART[aj].y;
	 theta2 = 2*M_PI*drand48();     
	 RING[ai].PART[aj].theta = RING[ai].PART[n].theta;
	 RING[ai].PART[aj].color = RING[ai].PART[n].color;
	 RING[ai].PART[aj].fx = 0.0;
	 RING[ai].PART[aj].fy = 0.0;
	 RING[ai].PART[aj].nx_right = 0.0;
	 RING[ai].PART[aj].nx_left = 0.0;
	 RING[ai].PART[aj].ny_up = 0.0;
	 RING[ai].PART[aj].ny_down = 0.0;
	 RING[ai].PART[n].x += RING[ai].PART[aj].x;
	 RING[ai].PART[n].y += RING[ai].PART[aj].y;
         
	 k++;
	 count++;
	 theta+=dtheta;
      }
      
      RING[ai].PART[n].x = RING[ai].PART[n].x*fator;
      RING[ai].PART[n].y = RING[ai].PART[n].y*fator;
      RING[ai].PART[n].x_virtual = RING[ai].PART[n].x;
      RING[ai].PART[n].y_virtual = RING[ai].PART[n].y;    
   }  
   
  
   RING[N].PART[n].x = 0.0;
   RING[N].PART[n].y = 0.0;
   RING[N].PART[n].color = COLOR;
   for(ai=0;ai<N;ai++){
      RING[N].PART[n].x += RING[ai].PART[n].x;
      RING[N].PART[n].y += RING[ai].PART[n].y;

   }
   
   RING[N].PART[n].x = (1.0/(1.0*N))*RING[ai].PART[n].x;
   RING[N].PART[n].y = (1.0/(1.0*N))*RING[ai].PART[n].y;
   RING[N].PART[n].x_virtual = RING[N].PART[n].x;
   RING[N].PART[n].y_virtual = RING[N].PART[n].y;
}


void VIEW_GNUPLOT(ring *RING,double t,double L,double PE,double G,double PHI,double Kadh,int N,double p0,double rc){
   
   int ai,aj,k,color;
   double x,y,vx,vy,rr,v;

   //Gnuplot
   // printf("unset key\n");
   // printf("unset xtics\n");
   // printf("unset ytics\n");
   // printf("set border lt 5\n");
   // printf("set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb\"black\" behind\n");
   // printf("set xrange [0.0:%lf]\n",L);
   // printf("set yrange [0.0:%lf]\n",L); 
   // printf("set size square \n");
   // printf("set title \"t: %.0lf, n:%d, N:%d, PHI:%.1lf, PE:%.1lf, G:%.0lf, p0:%.1lf, Kadh:%.1lf\" tc lt 5 \n",t,n,N,PHI,PE,G,p0,Kadh);
   // printf("plot \"-\" u 1:2:3:4 w circles lc variable fs solid 0.8\n");

   //print N RINGS
  
   for(ai=0;ai<(N);ai++){
      for(aj=0;aj<(n);aj++){ 
         x = RING[ai].PART[aj].x;
         y = RING[ai].PART[aj].y;
         vx = RING[ai].PART[aj].vx;
         vy =  RING[ai].PART[aj].vy;
         rr = 0.5*rc;
	      color = RING[ai].PART[aj].color;
	      //printf("%lf %lf %lf %lf %lf %d\n",x,y,vx/v,vy/v,rr,color);
	      //printf("%lf %lf %lf %d %lf %lf\n",x,y,rr,color,vx,vy);
      }
      
   }
   //printf("e\n");
}

void PBC(ring *RING,int ai,int aj,double L,double *z1,double *z2,double *z3,double *z4,int N){
   
   
   double x,y;
   int k;
   x = RING[ai].PART[aj].x;
   y = RING[ai].PART[aj].y;
   k= ai*N + aj;
   if(x > L){
      
      RING[ai].PART[aj].x = RING[ai].PART[aj].x - L;
      *z1 = 1.0;
   }
   
   else *z1 = 0.0;
   if(x < 0.0){
      
      RING[ai].PART[aj].x = L + RING[ai].PART[aj].x;
      *z2 = 1.0;
   }
   
   else *z2 = 0.0;
   if(y > L){
      
      RING[ai].PART[aj].y = RING[ai].PART[aj].y- L;
      *z3 = 1.0;
   }
   
   else *z3 = 0.0;
   
   if(y < 0.0){
      
      RING[ai].PART[aj].y = L + RING[ai].PART[aj].y;
      *z4 = 1.0;
   }
   
   else *z4 = 0.0;
}


void frep(double dr,double rc,int ai, int aj, int bi, int bj, double *force1,ring *RING,double Kc){
    
   *force1 = -(Kc)*(dr-rc);
   
}


void fadh(double dr,double rc,int ai, int aj, int bi, int bj, double *force2,ring *RING,double Kadh){
   
   *force2 = -(Kadh)*(dr-rc);
   
}

void Motion_evol(ring *RING,int ai,double L,int N,double PE,double G,double dt){
   
   int k,count,lx,ly,ind,aj;
   double r,prod,dtheta,x,y,x_left,y_up,y_down,x_right,dx,dy;
   double random,prod_CM,VCM;
   box_muller(&random);
   prod_CM = ((RING[ai].PART[n].nx*RING[ai].PART[n].vy)-(RING[ai].PART[n].ny*RING[ai].PART[n].vx));
   VCM = sqrt(RING[ai].PART[n].vx*RING[ai].PART[n].vx + RING[ai].PART[n].vy*RING[ai].PART[n].vy);
   if(VCM >0.0) dtheta = (dt*G*(prod_CM/VCM)) + sqrt(2.0*dt)*random;
   else dtheta =  sqrt(2.0*dt)*random;
   
   //updating the theta to theta + dtheta
   RING[ai].PART[n].theta += dtheta;

   RING[ai].PART[n].x_virtual = 0.0;
   RING[ai].PART[n].y_virtual = 0.0;
   RING[ai].PART[n].x = 0.0;
   RING[ai].PART[n].y = 0.0;
   RING[ai].PART[n].vx = 0.0;
   RING[ai].PART[n].vy = 0.0;
  
   for(aj=0;aj<n;aj++){
      
      //Repulsive Boundary Condition
      //RBC(RING,ai,aj,L);
                  
      //Active Overdamped Langevin Equation
      RING[ai].PART[aj].vx = (PE*(RING[ai].PART[n].nx) + (RING[ai].PART[aj].fx));
      RING[ai].PART[aj].vy = (PE*(RING[ai].PART[n].ny) + (RING[ai].PART[aj].fy));
      
      //updating the r to r + dr
      RING[ai].PART[aj].x += (dt*RING[ai].PART[aj].vx);
      RING[ai].PART[aj].y += (dt*RING[ai].PART[aj].vy);
      
      //Periodic Boundary Conditions
      PBC(RING,ai,aj,L,&x_right,&x_left,&y_up,&y_down,N);
      
      RING[ai].PART[aj].x_aux = RING[ai].PART[aj].x;
      RING[ai].PART[aj].y_aux = RING[ai].PART[aj].y;
      
      //Number of times the ring left the box
      RING[ai].PART[aj].nx_right += x_right;
      RING[ai].PART[aj].nx_left += x_left;
      RING[ai].PART[aj].ny_up +=y_up;
      RING[ai].PART[aj].ny_down +=y_down;
      
      //Virtual positions
      RING[ai].PART[aj].x_virtual = (RING[ai].PART[aj].x +(RING[ai].PART[aj].nx_right*L) - (RING[ai].PART[aj].nx_left*L));
      RING[ai].PART[aj].y_virtual = (RING[ai].PART[aj].y +(RING[ai].PART[aj].ny_up*L) - (RING[ai].PART[aj].ny_down*L));
      
      //Calculating virtual CM position of the ring
      RING[ai].PART[n].x_virtual+=RING[ai].PART[aj].x_virtual;
      RING[ai].PART[n].y_virtual+=RING[ai].PART[aj].y_virtual;
           
      //Calculating  CM position of the ring
      RING[ai].PART[n].x += RING[ai].PART[aj].x;
      RING[ai].PART[n].y += RING[ai].PART[aj].y;
       
       
      //Calculating virtual CM velocity of the ring
      RING[ai].PART[n].vx += RING[ai].PART[aj].vx;
      RING[ai].PART[n].vy += RING[ai].PART[aj].vy;
      
      //zerando as forças
      RING[ai].PART[aj].fx = 0.0;
      RING[ai].PART[aj].fy = 0.0;
 
      
   }
   
   
   //Updating the (nx,ny) to t+dt
   RING[ai].PART[n].nx = cos(RING[ai].PART[n].theta);
   RING[ai].PART[n].ny = sin(RING[ai].PART[n].theta);
   
   //CM position's ring
   RING[ai].PART[n].x = RING[ai].PART[n].x*fator;
   RING[ai].PART[n].y = RING[ai].PART[n].y*fator;
   PBC(RING,ai,n,L,&x_right,&x_left,&y_up,&y_down,N);
   RING[ai].PART[n].x_virtual = RING[ai].PART[n].x_virtual*fator;
   RING[ai].PART[n].y_virtual = RING[ai].PART[n].y_virtual*fator;
   RING[ai].PART[n].vx = RING[ai].PART[n].vx*fator;
   RING[ai].PART[n].vy = RING[ai].PART[n].vy*fator;
   RING[ai].PART[n].x_aux = 0.0;
   RING[ai].PART[n].y_aux = 0.0;
   
   for(aj=0;aj<n;aj++){	
      
      dx = RING[ai].PART[aj].x_aux-RING[ai].PART[(aj+1)%n].x_aux;
      dy =  RING[ai].PART[aj].y_aux-RING[ai].PART[(aj+1)%n].y_aux;
      if((fabs(dx)) >(0.5*L)) RING[ai].PART[(aj+1)%n].x_aux += ((dx)/fabs(dx))*L;
      if((fabs(dy))>(0.5*L))  RING[ai].PART[(aj+1)%n].y_aux += ((dy)/fabs(dy))*L;
      RING[ai].PART[n].x_aux+= RING[ai].PART[(aj+1)%n].x_aux;
      RING[ai].PART[n].y_aux+= RING[ai].PART[(aj+1)%n].y_aux;
   }
     
   //CM position's ring after correction
   RING[ai].PART[n].x = (RING[ai].PART[n].x_aux*fator);
   RING[ai].PART[n].y = (RING[ai].PART[n].y_aux*fator);
   //PBC(RING,ai,n,L,&x_right,&x_left,&y_up,&y_down,N);
   
}


void updating_list(ring *RING,int ai,double P0,int *nviz,int **list,double L, int N){
   
   int bi;
   double dX,dY,dR;
  
   list[ai][nviz[ai]] = ai;
   nviz[ai]++;
   for(bi=(ai+1);bi<N;bi++){
      
      dX = RING[ai].PART[n].x - RING[bi].PART[n].x;
      dY = RING[ai].PART[n].y - RING[bi].PART[n].y;
      dX = dX-rint(dX/L)*L;
      dY = dY-rint(dY/L)*L;
      dR = sqrt((dX*dX) + (dY*dY)); 
      if(dR<(0.45*P0)){ 
	 
	 list[ai][nviz[ai]] = bi;
	 list[bi][nviz[bi]] = ai;
	 nviz[ai]++;
	 nviz[bi]++;
	 
      }     
   }
}


void forces(ring *RING,int **list,int *nviz,int ai,double L,double Kadh,double radh,int N,double Ka,double rc,double Kl,double Kc){
   
   int bi,bj,Viz,aj,ai_aux,bi_aux,bj_aux,color[n*N];
   double dx,dy,dr,f,dx_cm,dy_cm,xj,yj,xcm,ycm,theta,dr_aux,f_ext; 
   double dx_right,dx_left,dy_right,dy_left,dr_right,dr_left;
   RING[ai].PART[n].A = 0.0;
   xj = 0.0;
   yj = 0.0;
   #pragma omp parallel for num_threads(4)
   for(aj=0;aj<n;aj++){

      //dx
      dx = RING[ai].PART[aj].x - RING[ai].PART[(aj+n-1)%n].x;
      //dy
      dy = RING[ai].PART[aj].y - RING[ai].PART[(aj+n-1)%n].y;
      
      //Minimum Image Convention
      dx = dx -rint(dx/L)*L;
      dy = dy -rint(dy/L)*L;
   
      //particles' positions
      xj = RING[ai].PART[aj].x;
      yj = RING[ai].PART[aj].y;
      //center of mass
      xcm = RING[ai].PART[n].x;
      ycm = RING[ai].PART[n].y;
      dx_cm = (xj-xcm);
      dy_cm= (yj-ycm);
  
      if((fabs(dx_cm)>(0.5*L)) && (dx_cm>0.0)) dx_cm = dx_cm-L;
      if((fabs(dx_cm)>(0.5*L)) && (dx_cm<0.0)) dx_cm = dx_cm+L;
      if((fabs(dy_cm)>(0.5*L)) && (dy_cm>0.0)) dy_cm = dy_cm-L;
      if((fabs(dy_cm)>(0.5*L)) && (dy_cm<0.0)) dy_cm = dy_cm+L;
      if((fabs(dx)>(0.5*L)) && (dx>0.0)) dx = dx-L;
      if((fabs(dx)>(0.5*L)) && (dx<0.0)) dx = dx+L;
      if((fabs(dy)>(0.5*L)) && (dy>0.0)) dy = dy-L;
      if((fabs(dy)>(0.5*L)) && (dy<0.0)) dy = dy+L;
      //Calculating area
      RING[ai].PART[n].A += 0.5*((dx_cm*dy)-(dx*dy_cm));

   }
      //Loop sobre as partículas do anel ai
      #pragma omp parallel for num_threads(4)
      for(aj=0;aj<n;aj++){
       //Loop sobre os vizinhos do anel ai	 
        for(Viz=0;Viz<nviz[ai];Viz++){
	 //indice do anel vizinho  
         bi = list[ai][Viz];
	 //loop sobre as partículas do anel vizinho bi  
         for(bj=aj;bj<n;bj++){
           
	    //External interations
	    if(ai!=bi){
	       
	       dx = RING[ai].PART[aj].x - RING[bi].PART[bj].x;
	       dy = RING[ai].PART[aj].y - RING[bi].PART[bj].y;
	       dx = dx-rint(dx/L)*L;
	       dy = dy-rint(dy/L)*L;
	       dr = sqrt((dx*dx) + (dy*dy));
	       
	       //Repulsive force
	       if(dr<rc){
		  
                  frep(dr,rc,ai,aj,bi,bj,&f,RING,Kc);
		  RING[ai].PART[aj].fx +=(f*dx)/(dr);
		  RING[bi].PART[bj].fx -=(f*dx)/(dr);
		  RING[ai].PART[aj].fy +=(f*dy)/(dr);
		  RING[bi].PART[bj].fy -=(f*dy)/(dr);
		 
	          
	  }	       	       
	       //Adhesion force
	       if((dr>(rc) && (dr<radh))){
		  
                  fadh(dr,rc,ai,aj,bi,bj,&f,RING,Kadh);
		  RING[ai].PART[aj].fx +=(f*dx)/(dr);
		  RING[bi].PART[bj].fx -=(f*dx)/(dr);
		  RING[ai].PART[aj].fy +=(f*dy)/(dr);
		  RING[bi].PART[bj].fy -=(f*dy)/(dr);
				  
	       }  
	    }
	    	    
	    //Internal interactions
	    if(ai==bi){
	       if(aj!=bj){
	       if(fabs(aj-bj)!=1){
		  
		  dx = RING[ai].PART[aj].x - RING[bi].PART[bj].x;
		  dy = RING[ai].PART[aj].y - RING[bi].PART[bj].y;
		  dx = dx-rint(dx/L)*L;
		  dy = dy-rint(dy/L)*L;
		  dr = sqrt((dx*dx) + (dy*dy));
		  //Repulsive force
		  if(dr<rc){
		     
                     frep(dr,rc,ai,aj,bi,bj,&f,RING,Kc);
		     RING[ai].PART[aj].fx +=(f*dx)/(dr);
		     RING[bi].PART[bj].fx -=(f*dx)/(dr);
		     RING[ai].PART[aj].fy +=(f*dy)/(dr);
		     RING[bi].PART[bj].fy -=(f*dy)/(dr);
		     
		  }
	        }  
	     }
	    }
	  }
        }  
      
      //bond force: forming the ring
           
      dx_right =  RING[ai].PART[(aj+1)%n].x - RING[ai].PART[aj].x;
      dy_right =  RING[ai].PART[(aj+1)%n].y - RING[ai].PART[aj].y;
      dx_right = dx_right-rint(dx_right/L)*L;
      dy_right = dy_right-rint(dy_right/L)*L;
      dr_right = sqrt((dx_right*dx_right)  +  (dy_right*dy_right));
  
      dx_left = RING[ai].PART[aj].x - RING[ai].PART[(aj+n-1)%n].x;
      dy_left = RING[ai].PART[aj].y - RING[ai].PART[(aj+n-1)%n].y;
      dx_left = dx_left-rint(dx_left/L)*L;
      dy_left = dy_left-rint(dy_left/L)*L;
      dr_left = sqrt((dx_left*dx_left)  +  (dy_left*dy_left));
     
      f = -Kl*(dr_right-rc);
      RING[ai].PART[aj].fx += -((f*dx_right)/(dr_right));
      RING[ai].PART[aj].fy += -((f*dy_right)/(dr_right));
      
      f = -Kl*(dr_left-rc);
      RING[ai].PART[aj].fx += ((f*dx_left)/(dr_left));
      RING[ai].PART[aj].fy += ((f*dy_left)/(dr_left));
          
      //Area force: cortical tension
      RING[ai].PART[aj].fx += (-0.5*Ka)*(RING[ai].PART[n].A - A0)*(dy_right + dy_left);
      RING[ai].PART[aj].fy += (0.5*Ka)*(RING[ai].PART[n].A - A0)*(dx_right + dx_left);

}
}
double box_muller(double *z){
   
   double rand1,ranrc_ext;
   rand1 = drand48();
   ranrc_ext = drand48();
   *z = sqrt(-2*log(rand1))*sin(2*M_PI*ranrc_ext);
   
}

void RBC(ring *RING,int ai,int aj,double L){ 
   
   double r,ld,ewall,rwall,fwall;
   ewall = 100.0;
   rwall = 1.0;
   //WALL 1 : (0,y)
   r = RING[ai].PART[aj].x;
   fwall = -(ewall/rwall)*((r/rwall)-1.0);
   if( RING[ai].PART[aj].x < 1.0)  RING[ai].PART[aj].fx +=  fwall;
   //WALL 2 : (x,L)
   r = L-RING[ai].PART[aj].y;
   fwall = -(ewall/rwall)*((r/rwall)-1.0);
   if((L-RING[ai].PART[aj].y) < 1.0) RING[ai].PART[aj].fy += -fwall;	
   //WALL 3 : (L,y)
   r = (L-RING[ai].PART[aj].x);
   fwall = -(ewall/rwall)*((r/rwall)-1.0);
   if((L-RING[ai].PART[aj].x) < 1.0) RING[ai].PART[aj].fx += -fwall;
   //WALL 4 : (x,0)
   r = RING[ai].PART[aj].y;
   fwall = -(ewall/rwall)*((r/rwall)-1.0);
   if(RING[ai].PART[aj].y < 1.0) RING[ai].PART[aj].fy += fwall;
   
}

