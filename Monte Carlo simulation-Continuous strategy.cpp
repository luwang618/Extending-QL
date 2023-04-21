// Extending Q-learning to continuous and mixed strategy games based on spatial reciprocity
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
//#include <windows.h>
using namespace std;
// define priority classes
#define NORMAL_PRIORITY_CLASS       0x00000020
#define IDLE_PRIORITY_CLASS         0x00000040
#define HIGH_PRIORITY_CLASS         0x00000080
#define REALTIME_PRIORITY_CLASS     0x00000100

// define parameters
#define L           200      /* lattice size                   */
#define SIZE        (L*L)    /* number of sites                */
#define MC_STEPS    50000   /* run-time in MCS     */
#define OUT_STEPS   1000   /* Last 5000 steps  */
#define NAMEOUT     "K4b075r5Q2"
#define RANDOMIZE   3145215
double eta,gamma,epsilon,T,S;// 
double Coo;
//
void prodgraph(void);      /* creates host graph                    */
void initial(void);        /* initial state                         */
double payoff(int y);
void update_stra(int x);
int is_max(double x, double y); 
double tongji(void);

ofstream outfile1;
ofstream outfile2;

/*************************** RNG procedures ****************************************/
#define NN 624
#define MM 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)
static unsigned long mt[NN]; /* the array for the state vector  */
static int mti=NN+1; /* mti==NN+1 means mt[NN] is not initialized */
void sgenrand(unsigned long seed)
{int i;
 for (i=0;i<NN;i++) {mt[i] = seed & 0xffff0000; seed = 69069 * seed + 1;
                     mt[i] |= (seed & 0xffff0000) >> 16; seed = 69069 * seed + 1;
  }
  mti = NN;
}
void lsgenrand(unsigned long seed_array[])
{ int i; for (i=0;i<NN;i++) mt[i] = seed_array[i]; mti=NN; }
double genrand() 
{
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    if (mti >= NN) 
    {
        int kk;
        if (mti == NN+1) sgenrand(4357); 
        for (kk=0;kk<NN-MM;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+MM] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<NN-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(MM-NN)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[NN-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[NN-1] = mt[MM-1] ^ (y >> 1) ^ mag01[y & 0x1];
        mti = 0;
    }  
    y = mt[mti++]; y ^= TEMPERING_SHIFT_U(y); y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C; y ^= TEMPERING_SHIFT_L(y);
    return y;  
}
double randf(){ return ( (double)genrand() * 2.3283064370807974e-10 ); }
long randi(unsigned long LIM){ return((unsigned long)genrand() % LIM); }

/********************** END of RNG ************************************/
								  
struct Agent
{
    double Strat;
	double Q_matrix[5][5];
	int Choose;
	int Neighbours[4];  
};

struct Agent Player[SIZE];		
 
void initial(void)
{
	for (int i=0; i<SIZE; i++)
	{
		Player[i].Strat=randf();
		Player[i].Choose=0;		
		memset(Player[i].Q_matrix,0,sizeof(Player[i].Q_matrix));
	}                          
}

void prodgraph(void)
{
	int  iu, ju;
	long int  player1,player2;
	int i,j;
	// set up an initial square lattice
	for(i=0; i<L; i++)                     
	{
		for(j=0; j<L; j++)
		{ 
			// the first player
			player1 = L * j + i;             
			// and its four nearest neighbors
			iu = i + 1;         
			ju = j;     
			if (iu==L) iu = 0;
			player2 = L * ju + iu;  
			Player[player1].Neighbours[0] = player2;

			iu = i;             
			ju = j + 1; 
			if (ju==L) ju = 0;
			player2 = L * ju + iu;  
			Player[player1].Neighbours[1] = player2;

			iu = i - 1;         
			ju = j;     
			if (i==0) iu = L - 1;
			player2 = L * ju + iu;  
			Player[player1].Neighbours[2] = player2;

			iu = i;             
			ju = j - 1; 
			if (j==0) ju = L - 1;
			player2 = L * ju + iu;  
			Player[player1].Neighbours[3] = player2;			
		}
	}

}


double payoff(int y)
{	double pay = 0,pi=0,strat,strat3;
	int player3;
	strat = Player[y].Strat;      	    	
	for(int j=0;j<4;j++)
    {
	 player3=Player[y].Neighbours[j];
     strat3=Player[player3].Strat;
     pi=S*strat+T*strat3+(1-S-T)*strat*strat3;
     pay= pay+pi;
    }	
	return pay;
}


int is_max(int x ,int y)	 
{ double max;     
  int a[5],k,jj,suiji;	
  max=0;
  k=0; 
  jj=0;
  memset(a,0,sizeof(int)*5);
  max=Player[x].Q_matrix[y][0];
  a[k]=0;
  for(int i=1; i<5; i++) 
	{  if(Player[x].Q_matrix[y][i]>max)
		  {		   
		   max=Player[x].Q_matrix[y][i];
		   memset(a,0,sizeof(int)*5); 
		   k=0; 
		   a[k]=i;	 
		  }
		  else if (Player[x].Q_matrix[y][i]==max)
		  {
		  	k++; 
			a[k]=i;		  	
		  }		  
	}
      if(k>0)
     {
      k++;
	  jj=(int)randi(k);
	  suiji=a[jj];	
     }
     else
     suiji=a[k]; 
  
     return suiji;
}


void update_stra(int x)
{
	int action, state1, state2;
	double P1 = 0;
    int i,j,Neig;

	state1=Player[x].Choose;
	if(randf()<epsilon) 
	{
	 action=(int)randi(5); 
     if(action!=0)
	 {
	  Neig=Player[x].Neighbours[action-1];
	  Player[x].Strat = Player[Neig].Strat;
	 }		
	} 
	else
	{  	action = is_max(x,state1); 
	
		  if (action!=0)	
		   {
		     Neig=Player[x].Neighbours[action-1];
		     Player[x].Strat = Player[Neig].Strat;
		   }
        
	}
    Player[x].Choose=action;

    state2 = action;
	
	P1 = payoff(x);		
	int temp = is_max(x,state2);
	Player[x].Q_matrix[state1][action] = (1-eta) * Player[x].Q_matrix[state1][action] + eta * (P1 + gamma * Player[x].Q_matrix[state2][temp]);

}


double tongji()
{	Coo=0;
	for(int i=0; i<SIZE; i++)
	{
	 Coo=Coo+Player[i].Strat;		
	}
}

int main()
{
	int steps;
    double fc,fc_sum,Afc;

	outfile1.open("frequency.txt");
	outfile2.open("average.txt");
    if(!outfile1){
	 cout<<"can not open";
	 abort();
	}  
	if(!outfile2){
	 cout<<"can not open";
	 abort();
	}
	
	sgenrand(RANDOMIZE);
	prodgraph();
	
	epsilon=0.02; 
	gamma = 0.8;
    eta = 0.8;
    S=0;
	for(T=1.2;T<1.25;T=T+0.05)
    {
	  fc_sum=0;
	  initial();  

	for (steps=0; steps<MC_STEPS; steps++)
	{		    	    
	    for (int i = 0; i < SIZE; i++)
		{
		update_stra((int)randi(SIZE));
		}
	    tongji();	    		    						
		fc =(double) Coo/SIZE;	
		cout<<steps<<'\t'<<fc<<endl;	    
	    
		if(steps > MC_STEPS - OUT_STEPS - 1)
			{fc_sum+=fc;
		     Afc=(double)fc_sum/OUT_STEPS;
		    }     		        	   	 
	}
   outfile2<<T<<'\t'<<Afc<<endl;
  }	     		
	outfile1.close();
	outfile2.close();

	return 0;
}


