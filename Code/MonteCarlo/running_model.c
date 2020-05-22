/* This code runs a Monte Carlo simulation of NBALLS for a maximum of NSTEPS on piecewise linear uneven terrain.
The output is a text (.dat) file that contains the number of steps before failure for each ball and the failure mode.*/

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

/* Simulation parameters */
#define NSTATE 6
#define NSTEPS 500
#define NBALLS 1000000 /* I'm plotting with 10^6 to ensure that all distributions have definitely converged */

/* Ball and physical parameters */
#define vscale 3.1320919526732 /*sqrt(rballd*gd). To non-dimensionalize the simulations.*/
#define rballd 1.0;
const double rball = 1.0;
#define md 63.0 /* From data */
#define MId 10.92 /* From data */
#define gr 1.0
#define m 1.0
const double MI = MId/md;

// Control variables
const double rt = 0.99;
const double rn = 0.63;

/* Failure parameters */
# define xtol 1.71 
# define tol M_PI/6  /* Chosen based on static stance at 'step length'*/

/* Intial conditions */
double vx0 = 3/vscale;
double vy0 = 0.8/(vscale);

// Terrain parameters for pregenerated terrain
const double space = 0.1;
#define hmean 0.0
#define scaleh 0.060/2

// Defining failure modes and initializing to SUCCESS
enum success_t {SUCCESS, HEADHIT, SPEED, TRIP};
enum success_t bouncesuccess = SUCCESS;

int flightdynamics (double *state, double* terrain, double counter, int txcount);
int intersect (double y, double* t, int i, double xstart, double xlift, double vxlift, double vylift, double* out);
void impulsecalc (double *state, double* out);
void collisiondynamics (double *state, double thland);
int jumpdynamics (double *state, double *impulsevecnorm, double thland);
void initialize (double vxin, double vyin, double *vec);
double absol (double num);
void generateterrain(double* terrain, int size);
double getval(double a0, double a1, double a2, double a3, double a4, double x);
double roots(double y, double vx, double vy, double xt, double yt, double m1, double m2, double xlift);
void rpoly_ak1(double op[], int *degree, double zeror[], double zeroi[]);

int main()
{
    double thland, bcounter;
    int i, j, k, l, s,txcount, size = (int)(NSTEPS*10*(double)rball/space);
    unsigned int result;
    time_t seed;
    FILE *output;
  
    output = fopen ("failures.dat","w");
   
     /* Big for loop */
    
    //Initializing seed for random number generator
    seed = time(NULL);
    srand(seed);

    int* bounce = (int*)malloc(NBALLS*sizeof(int)); // Number of bounces per ball
    int* fail = (int*)malloc(NBALLS*sizeof(int)); // Failure mode of each ball*/
    
    for (i=0; i<NBALLS; i++)
    {
        bouncesuccess = SUCCESS;

        //Initializing terrain, ball physical state and push-off vector
        double *terrain = (double *) malloc (size*sizeof(double));
        double *nextstate = (double*) malloc((NSTATE+4)*sizeof(double));
        double *impulsevecnorm = (double *)malloc(3*sizeof(double));
        
        // Generate piecewise-flat terrain
        generateterrain(terrain,size);
    
        bcounter = 0.0; //x-position of ball contact point
        txcount = 0; //index of terrain pt before ball
        
        initialize (vx0, vy0, nextstate);
        
        //Calculate impulse required for steady motion of flat ground
        impulsecalc (nextstate, impulsevecnorm);
		
        for (j=0; j<NSTEPS; j++)
        {
            
            // Propagate ball through flght phase of running
            result = flightdynamics (nextstate, terrain, bcounter, txcount);
            if (result==0) // if failure
            {
                break;
            }
            // Extract variables required for collision
            thland = nextstate[6];
            bcounter = nextstate[7];
            txcount = (int) nextstate[8];
			
            collisiondynamics (nextstate, thland);
			
            result = jumpdynamics (nextstate, impulsevecnorm, thland);
			
            if (result==0) //if failure
            {
                break;
            }
        }
		bounce[i] = j;
		fail[i] = bouncesuccess;
		
        free(terrain);
        free(nextstate);
        free(impulsevecnorm);
    }
    // Print output
    for (i = 0; i<NBALLS; i++)
    {
        fprintf(output,"%d, %d, \n", bounce[i], fail[i]);
    }
    free(bounce);
    free(fail);

    fclose(output);
    return 0;
}

int flightdynamics (double *state, double* terrain, double counter, int txcount)
//counter is ball contact point, tcount is grid pt before ball contact point
{
    double *hit = malloc(6*sizeof(double));
	double *t = malloc(3*sizeof(double));
    double xlift = state[0], ylift = state[1];
    double vxlift = state[2], vylift = state[3];
    double philift = state[4], wlift = state[5];
    double tflight, thland, landx,ballh;
    double b0, b1, b2, x, yt;
    int i, result=1, flag = 0;
    
    i = txcount;

    //Defining coefficients of parabolic flight path
    b0 = ylift - xlift*vylift/vxlift-gr*xlift*xlift/(2*vxlift*vxlift);
    b1 = vylift/vxlift + gr*xlift/(vxlift*vxlift);
    b2 = -gr/(2*vxlift*vxlift);
   
    
    while (flag == 0)
    {
        x = (double)i*space;
        ballh = b0 + b1*(x+rball) + b2*(x+rball)*(x+rball);
        yt = terrain[i];
        // Searching for terrain patch where ball lands
        if (ballh < yt+rball)
        {
            flag = intersect(ylift, terrain, i, x, xlift, vxlift, vylift, hit);
        }
        i++;
    }
    
   if (flag == 1) //found terrain patch
   {
       landx = hit[0];
       thland = hit[1];
       
       state[0] = hit[2]; /*xland */
       state[1] = hit[3]; /*yland  */
       tflight = (state[0] - xlift)/vxlift;
       state[3] = vylift - gr*tflight; /* vyland */
       state[2] = vxlift; /* vxland */
       state[4] = philift + wlift*tflight; /* philand */
       state[5] = wlift; /* wland */
       state[6] = thland;
       state[7] = landx;
       state[8] = hit[5];
	   state[9] = hit[4];
       
       if ( (state[4] > tol) || (state[4] < -tol) )
       {
           bouncesuccess = HEADHIT;
           result = 0;
       }
	
   }
    else
    {
        result = 0;
    }
    
    free(t);
	free(hit);
    return result;
}

void generateterrain(double* terrain, int size)
{
    // Generates piecewise flat terrain using random numbers to pick heights of terrain grid points
    int i, j;
	
    for (j = 0; j < (int) round(rball/space); j++)
    {
        terrain[j] = 0;
    }
	for (i = (int) round(rball/space); i< size; i++)
	{
		terrain[i] = (hmean + scaleh*((2.0*rand()/RAND_MAX)-1.0));
	}
}

int intersect (double ylift, double* t, int i, double xstart, double xlift, double vxlift, double vylift, double* out)
{
    // Check for intersection of ball's flight path with terrain patch

    double dx=space,xk,xk1,yk,dyk,yk1,denom,disc,sl2,A,B,C,sn,cs,sl,rootmax,m0,m01,xb1,thtrip=0,sl1,xb=0,lastxb = xstart + 10*rball,roottrue=0,sltrue=0,ytrue=0,b0,b1,b2;
    int flag = 0, j=0,k,counttrue=i;
  
    b0 = ylift - xlift*vylift/vxlift-gr*xlift*xlift/(2*vxlift*vxlift);
    b1 = vylift/vxlift + gr*xlift/(vxlift*vxlift);
    b2 = -gr/(2*vxlift*vxlift);

    do
    {
        xk = xstart+j*space;
        xk1 = xk+space;
        yk = t[i+j];
        yk1 = t[i+j+1];
        dyk = yk1-yk;
        denom = sqrt(dx*dx+dyk*dyk);
        
        sn = dyk/denom;
        cs = dx/denom;
        sl = dyk/dx;
        if (i+j > 1)
        {
            sl1 = (yk - t[i+j-1])/dx;
        }
        else
        {
            sl1 = 0;
        }
        
        
        A = b2;
        B = b1-2*rball*sn*b2-sl;
        C = b0-b1*rball*sn+b2*rball*rball*sn*sn+sl*xk-yk-rball*cs;
        disc = B*B-4*A*C;
        
        
        if (disc>=0)
        {
            rootmax = (-B-sqrt(disc))/(2*A);

            if (rootmax < xk1 && rootmax > xk)
            {
                xb = rootmax - rball*sn;
                if (xb > xlift)
                {
                    lastxb = xb;
                    roottrue = rootmax;
                    sltrue = sl;
                    thtrip = sltrue;
                    ytrue = yk+(roottrue-xk)*sltrue+rball*cs;
                    counttrue = i+j;
                    flag = 1;
                }
            }
            else if (rootmax <= xk)
            {
                if (sl == 0 && sl1 == 0)
                {
                    m0 = 0;
                }
                else
                {
                    m0 = roots(ylift,vxlift,vylift,xk,yk,sl1,sl,xlift);
                }   
                xb = xk - rball*m0/pow(m0*m0+1,0.5);
                if (xb > xlift)// && m0 != -100)
                {
                    lastxb = xb;
                    roottrue = xk;
                    sltrue = m0;
                    thtrip = sl;
                    ytrue = yk+rball/pow(1+m0*m0,0.5);
                    counttrue = i+j;
                    flag = 1;
                }
               
            }
        }
        else
        {
            sl2 = (t[i+j+2]-yk1)/space;
                if (sl == 0 && sl1 == 0)
                {
				    m0 = 0;
                }
                else
                {
                    m0 = roots(ylift,vxlift,vylift,xk,yk,sl1,sl,xlift);
                }
	            xb = xk - rball*m0/pow(m0*m0+1,0.5);
				
				if (sl == 0 && sl2 == 0)
                {
                    m01 = 0;
                }
                else
                {
                    m01 = roots(ylift,vxlift,vylift,xk1,yk1,sl,sl2,xlift);
                }
	            xb1 = xk1 - rball*m01/pow(m01*m01+1,0.5);
				
				if ((xb<xb1 && m0 != -100) && sl1 > sl)
                {
                    lastxb = xb;
                    roottrue = xk;
                    sltrue = m0;
                    thtrip = sl;
                    ytrue = yk+rball/pow(1+m0*m0,0.5);
                    flag = 1;
                    counttrue = i+j;
                }
				else if (m01 != -100 && sl > sl2)
				{
		            lastxb = xb1;
		            roottrue = xk1;
		            sltrue = m01;
					thtrip = sl2;
		            ytrue = yk1+rball/pow(1+m01*m01,0.5);
		            flag = 1;
                    counttrue = i+j+1;
				}
				
			
        }
        j++;
    } while (flag == 0);
    
    xstart = xk1;
    
    if (flag == 1)
    {
        for (k = 0; k < (int) round(0.75*rball/space); k++)
        {
            xk = xstart+k*space;
            xk1 = xk+space;
            yk = t[i+j+k];
            yk1 = t[i+j+k+1];
            dyk = yk1-yk;
            denom = sqrt(dx*dx+dyk*dyk);
            
            sn = dyk/denom;
            cs = dx/denom;
            sl = dyk/dx;
            sl1 = (yk - t[i+j+k-1])/space;
            
            A = b2;
            B = b1-2*rball*sn*b2-sl;
            C = b0-b1*rball*sn+b2*rball*rball*sn*sn+sl*xk-yk-rball*cs;
            disc = B*B-4*A*C;
            
            if (disc>=0)
            {
                rootmax = (-B-sqrt(disc))/(2*A);
                if (rootmax < xk1 && rootmax > xk)
                {
                    xb = rootmax - rball*sn;
                    if (xb < lastxb)
                    {
                        lastxb = xb;
                        roottrue = rootmax;
                        sltrue = sl;
                        ytrue = yk+(roottrue-xk)*sltrue+rball*cs;
                        thtrip = sltrue;
                        counttrue = i+j+k;
                    }
                }
                else if (rootmax <= xk && sl1 > sl)
                {
                    if (sl == 0 && sl1 == 0)
                    {
                        m0 = 0;
                    }
                    else
                    {
                        m0 = roots(ylift,vxlift,vylift,xk,yk,sl1,sl,xlift);
                    }
                    if (m0 != -100)
                    {
                        xb = xk - rball*m0/pow(m0*m0+1,0.5);
                        if (xb < lastxb)
                        {
                            lastxb = xb;
                            roottrue = xk;
                            sltrue = m0;
                            ytrue = yk+rball/pow(1+m0*m0,0.5);
                            thtrip = sl;
                            counttrue = i+j+k;
                           
                        }
                    }
                }
            }
            else
                {                    
                    sl2 = (t[i+j+2]-yk1)/space;
                    if (sl == 0 && sl1 == 0)
                    {
                        m0 = 0;
                    }
                    else
                    {
                        m0 = roots(ylift,vxlift,vylift,xk,yk,sl1,sl,xlift);
                    }
                        
                        xb = xk - rball*m0/pow(m0*m0+1,0.5);

                        if (sl == 0 && sl2 == 0)
                        {
                            m01 = 0;
                        }
                        else
                        {
                            m01 = roots(ylift,vxlift,vylift,xk1,yk1,sl,sl2,xlift);
                        }
                        
                        xb1 = xk1 - rball*m01/pow(m01*m01+1,0.5);
                    
                        if (m0 != -100 && sl1 > sl)
                        {
                            
                            if (xb < lastxb)
                            {
                                lastxb = xb;
                                roottrue = xk;
                                sltrue = m0;
                                thtrip = sl;
                                ytrue = yk+rball/pow(1+m0*m0,0.5);
                                counttrue = i+j+k;
                            }
                        }
                        else
                        {
                            if (xb1 < lastxb && m01 != -100 && sl > sl2)
                            {
                                lastxb = xb1;
                                roottrue = xk1;
                                sltrue = m01;
                                thtrip = sl2;
                                ytrue = yk1+rball/pow(1+m01*m01,0.5);
                                counttrue = i+j+k+1;
                            }
                        }
                }
        }
    }
    if (flag == 1)
    {
        out[0]= roottrue;
        out[1] = atan(sltrue);
        out[2] = lastxb;
        out[3] = ytrue;
		out[4] = thtrip;
        out[5] = (double)counttrue;
    }

    return flag;
}

double roots(double y, double vx, double vy, double xt, double yt, double m1, double m2, double x)
{
    // Find roots of 4th order polynomial corresponding to intersection of parabolic trajectory and terrain patch

    double mnew,x0,m0,y0,xc,mind;
    int l=0,k, degree = 4,i=0,f=0;
    double op[101] = {0},zeror[100] = {0},zeroi[100]={0},d[4]={0};
	double *root=malloc(4*sizeof(double)),*xpos=malloc(4*sizeof(double)),*ypos=malloc(4*sizeof(double));
    
    mind = 0.01; //tolerance

    // Coefficients of 4th order polynomial

    op[4] = -1 + (vy*pow(x - xt,3))/pow(vx,3) + pow(x - xt,4)/(4.*pow(vx,4)) +
    pow(y,2) - (2*vy*(x - xt)*(y - yt))/vx - 2*y*yt + pow(yt,2) +
    (pow(x - xt,2)*(pow(vy,2) - y + yt))/pow(vx,2);
	
    op[3] = (-2*vy)/vx - (2*(x - xt))/pow(vx,2);
	
    op[2] = (pow(x - xt,2)*(-1 + pow(x,2) - 2*x*xt + pow(xt,2)))/(2.*pow(vx,4)) +
    (vy*(2*pow(x,3) + xt - 6*pow(x,2)*xt - 2*pow(xt,3) +
         x*(-1 + 6*pow(xt,2))))/pow(vx,3) +
    (pow(vy,2)*(-1 + 2*pow(x,2) - 4*x*xt + 2*pow(xt,2)) -
     (1 + 2*pow(x,2) - 4*x*xt + 2*pow(xt,2))*(y - yt))/pow(vx,2) -
    (4*vy*(x - xt)*(y - yt))/vx + (-2 + 4*pow(y,2) - 8*y*yt + 4*pow(yt,2))/2.;
	
    op[1] = (-2*vy)/vx - (2*(x - xt))/pow(vx,2);
    
    op[0] = pow(-1 + pow(x,2) - 2*x*xt + pow(xt,2),2)/(4.*pow(vx,4)) +
    (vy*(pow(x,3) + xt - 3*pow(x,2)*xt - pow(xt,3) +
         x*(-1 + 3*pow(xt,2))))/pow(vx,3) +
    (pow(vy,2)*(-1 + pow(x,2) - 2*x*xt + pow(xt,2)) -
     (1 + pow(x,2) - 2*x*xt + pow(xt,2))*(y - yt))/pow(vx,2) -
    (2*vy*(x - xt)*(y - yt))/vx + pow(y - yt,2);
	
    // Using Jenkins-Traub algorithm to evaluate roots
	rpoly_ak1(op, &degree, zeror, zeroi);
	
	for (k = 0; k<degree;k++)
	{
        if (fabs(zeroi[k]) < 0.00034526698300124393200064010223) //round-off error
        {
            zeroi[k]=0;
        }
		if (zeroi[k] == 0 && (zeror[k] > m2 && zeror[k] < m1))
		{
			m0=zeror[k];
			x0 = xt - rball*(m0/pow(m0*m0+1,0.5));
            y0 = (-gr*((x0-x)*(x0-x))/(2*vx*vx))+((x0-x)*vy/vx)+y;//Re-write this
			if (x0 > x+0.002 && y0 > 0.5)
			{
				root[l] = m0;
				xpos[l] = x0;
				ypos[l]= yt + rball/pow(m0*m0+1,0.5);
                d[l] = fabs(1-((x0-xt)*(x0-xt)+(y0-yt)*(y0-yt)));
				l++;
			}
		}	
	}
	mnew = -100;
    xc = xpos[0]+100;
	
	if (l==0)
    {
        do
        {
            if (zeroi[i]!=0 && (zeror[i] > m2 && zeror[i] < m1) && fabs(zeroi[i])<0.1)
            {
                mnew = zeror[i];
                f=1;
            }
            i++;
        }while(f==0 && i < 4);
    }
    else if (l ==1)
	{
		mnew = root[0];
	}
	else
	{
		for (k = 0; k<l;k++)
		{
			if(d[k] <= mind)
			{
				if (xpos[k]<xc)
                {
                    xc = xpos[k];
                    mnew = root[k];
                }
            }
		}
	}

	free(root);
	free(xpos);
    free(ypos);
	return mnew;
}

double getval(double a0, double a1, double a2, double a3, double a4, double x)
{
	double out = a0+a1*x+a2*x*x+a3*x*x*x+a4*x*x*x*x;
	return out;
}


void collisiondynamics (double *state, double thland)
{
    /* Look at Big_ball_impulse_gen_MI_map.nb for derivation of transition equations */
    double vxland = state[2], vyland = state[3],  wland = state[5], vxp, vyp;
	
    double et = rt; // for open-loop runners, this would differ
	
    // Collision equations
    state[2] = ((-m*rball*rball*(-1 + rn) + MI*(-rn + et))*vxland +
                2*MI*rball*(-1 + et)*wland*cos(thland) +
                (m*rball*rball*(1 + rn) + MI*(rn + et))*(vxland*cos(2*thland) + vyland*sin(2*thland)))/
    (2*(MI + m*rball*rball));
    
    state[3] = 1/(2*(MI + m*rball*rball))*
    ((-m*rball*rball*(-1 + rn) + MI*(-rn + et))*vyland +
     2*MI*rball*(-1 + et)*wland*sin(thland) +
     (m*rball*rball*(1 + rn) + MI*(rn + et))*(-vyland*cos(2*thland) +  vxland*sin(2*thland)));
    
    state[5] = ((MI + m*rball*rball*et)*wland +
                m*rball*(-1 + et)*(vxland*cos(thland) + vyland*sin(thland)))
    /(MI + m*rball*rball);

    
}


int jumpdynamics (double *state, double* impulsevecnorm, double thland)
{
    /* Look at Big_ball_impulse_gen_MI_map.nb for derivation of transition equations */
    double vtimp = impulsevecnorm[0], vnimp = impulsevecnorm[1], vximp, vyimp, wimp=impulsevecnorm[2], thtrip = state[9];
    vximp = vtimp*cos(thland)-vnimp*sin(thland);
    vyimp = vtimp*sin(thland)+vnimp*cos(thland);
	
    int result = 1;
    
    state[2] += vximp;

    state[3] += vyimp;
	
    state[5] += wimp;
	
   if (state[2] <= 0.01)
   {
       bouncesuccess = SPEED;
       result = 0;
   }
    else if ( state[3]/state[2] < thtrip )
    {
        bouncesuccess = SPEED;
        result = 0;
    }
    return result;
}

void impulsecalc (double *state, double* out1)

{
    if (rball == 0)
    {
        out1[0]= (1-rt)*state[2];
        out1[1] = (1 - rn)*state[3]; /* vyimp */
        out1[2] = 0; /* wimp */
        
    }
    else
    {
        out1[0]= (1 - rt)*state[2] - rball*rt*state[5] +
        rball*(MI*state[5] + m*rball*((-1 + rt)*state[2] + rball*rt*state[5]))/
        (MI + m*rball*rball);
        out1[1]= (1 - rn)*state[3];
        out1[2]= state[5] - (MI*state[5] + m*rball*((-1 + rt)*state[2] + rball*rt*state[5]))/
        (MI + m*rball*rball);
    }
}

void initialize(double vxin, double vyin, double *out)
{
    /* format is [x,y,vx,vy,phi,w] */
    out[0] = 0;
    out[1] = rball;
    out[2] = vxin;
    out[3] = vyin;
    out[4] = 0;
    out[5] = 0;
    out[6] = 0;
    out[7] = 0;
    out[8] = 0;
	out[9] = 0;
}

double absol (double num)
{
    double out = num;
    if (num < 0)
    {
        out = -out;
    }
    return out;
}