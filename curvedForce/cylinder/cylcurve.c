#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
//----------------------------------------------------------------------
#define nx  600
#define ny 	401
#define q		  9
//----------------------------------------------------------------------
const int time			= 200000;
const int noOfSnaps	= 		 6;
const int dispFreq	= 	 20;
const double D			= 	15.0;
const double tau		=		0.56;
const double rho0		= 	 1.0;
const double uInlet	=	 0.085;
const double xc			=	  75.0;
const double yc			=	 201.0;
//----------------------------------------------------------------------
int main(void)
{
	int i,j,a,a1,ts,ia,ja,cnt;
	static int ex[q],ey[q],kb[q];
	static int isn[nx+2][ny+2],wb[nx+2][ny+2];

	double tmp1, tmp2, tmp3, invC, rhoAvg, feq, f_neq, fx_t, fy_t, Cd, Cl;
	double dudx, dudy, dvdx, dvdy,ii,jj,chi,i_,j_,ia_,ja_;//,Delta[1000]
	double fx[2], fy[2],Q_xx[q],Q_yy[q],Q_xy[q],ubfx,ubfy,uwx,uwy;
	static double f[q][nx+2][ny+2], ft[q][nx+2][ny+2];
	static double wt[q], ux[nx+2][ny+2], uy[nx+2][ny+2], rho[nx+2][ny+2];

	char prefix[]="snap_",type[]=".dat",filename[15],solstr[5];
	int solnumber=0;
	double Delta[1000] = {0.0};
	invC=3.0;
	FILE *soln, *fid;
	//----------------------------------------------------------------------
	ex[0] =	0; ey[0] = 0;
	ex[1] = 1; ey[1] = 0;
	ex[2] = 0; ey[2] = 1;
	ex[3] =-1; ey[3] = 0;
	ex[4] = 0; ey[4] =-1;
	ex[5] = 1; ey[5] = 1;
	ex[6] =-1; ey[6] = 1;
	ex[7] =-1; ey[7] =-1;
	ex[8] = 1; ey[8] =-1;
	//----------------------------------------------------------------------
	for (a=0; a<9; a++)
	{
		if (a==0) 				{wt[a] = 4.0/9.0 ;}
		if (a>=1 && a<=4)	{wt[a] = 1.0/9.0 ;}
		if (a>=5 && a<=8)	{wt[a] = 1.0/36.0;}
	}
	//----------------------------------------------------------------------
	for (a=0; a<q; a++)
	{
		for (a1=a; a1<q; a1++)
		{
			if ( ex[a]+ex[a1]==0 && ey[a]+ey[a1]==0)
			{
				kb[a]	= a1;
				kb[a1]=  a;
			}
		}
	}
	//----------------------------------------------------------------------
	for (i=0; i<=nx+1; i++)
	{
		for (j=0; j<=ny+1; j++)
		{
			for (a=0; a<9; a++)
			{
				tmp1		= uInlet*ex[a];
				tmp2		= uInlet*uInlet;
				f[a][i][j] = wt[a]*rho0*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2);
			}
		}
	}
	//----------------------------------------------------------------------
	fx[0]  = 0.0;
	fy[0]  = 0.0;

	fid = fopen("tRhofxfy.dat","w");
	fprintf(fid,"Variables=time,rho,fx,fy\n");
	//----------------------------------------------------------------------
	for (i=1; i<=nx; i++)
	{
		for (j=0; j<=ny+1; j++)
		{
			isn[i][j] = 0;
			wb[i][j] = 0;
			tmp1 = pow((pow(i-xc,2.0) + pow(j-yc,2.0)),0.5);

			if (tmp1 <= 0.5*D)	isn[i][j] = 1;
			else								isn[i][j] = 0;
		}
	}

	cnt=-1;
	for (i=1; i<=nx; i++) //BC
	{
		for (j=1; j<=ny; j++)
		{
			if (isn[i][j] == 0)
			{
				for (a=0; a<q; a++)
				{
					ia = i+ex[a];
					ja = j+ey[a];

					if (isn[ia][ja]==1) //Cylinder
					{
						wb[ia][ja] = 1;
						i_=i;
						j_=j;
						ia_=ia;
						ja_=ja;
						ii=(i_+ia_)/2.0;
						jj=(j_+ja_)/2.0;
						//printf("%s\n","A" );
						while(fabs( (ii-xc)*(ii-xc) + (jj-yc)*(jj-yc) - D*D/4.0) > 0.0001)
						{
							//printf("%s\n","in" );
							if( (ii-xc)*(ii-xc) + (jj-yc)*(jj-yc) - D*D/4.0 <= 0.0)
							{
								ia_=ii;
								ja_=jj;
							}
							else
							{
								i_=ii;
								j_=jj;
							}

							ii=(i_+ia_)/2.0;
							jj=(j_+ja_)/2.0;
						}

						cnt=cnt+1;
						Delta[cnt] = sqrt((i-ii)*(i-ii) + (j-jj)*(j-jj))/sqrt((i-ia)*(i-ia) + (j-ja)*(j-ja));
						// chi[cnt][1] = i;
						// chi[cnt][2] = j;
						// chi[cnt][3] = ia;
						// chi[cnt][4] = ja;
						// printf("%10.6f\t%d\t%d\t%d\t%d\t%d\t%d\t%10.6f\t%10.6f\n",Delta[cnt],cnt,a,i,j,ia,ja,ii,jj);
					}
				}
			}
		}
	}

	//return 0;
	//----------------------------------------------------------------------
	for (ts=0; ts<=time; ts++)
	{

		rhoAvg = 0.0;

		for (i=1; i<=nx; i++)
		{
			for (j=1; j<=ny; j++)
			{
				tmp1 = 0.0;
				tmp2 = 0.0;
				tmp3 = 0.0;

				for (a=0; a<q; a++)
				{
					tmp1 += f[a][i][j];
					tmp2 += f[a][i][j]*ex[a];
					tmp3 += f[a][i][j]*ey[a];
				}

				rho[i][j] = tmp1;
				ux[i][j]	= tmp2/tmp1;
				uy[i][j]	= tmp3/tmp1;

				rhoAvg	 += tmp1;
			}
		}

		rhoAvg /=(nx*ny);
		//----------------------------------------------------------------------
		for (i=1; i<=nx; i++)
		{
			for (j=1; j<=ny; j++)
			{
				for (a=0; a<q; a++)
				{
					tmp1				= ux[i][j]*ex[a]  + uy[i][j]*ey[a] ;
					tmp2				= ux[i][j]*ux[i][j] + uy[i][j]*uy[i][j];
					feq					= wt[a]*rho[i][j]*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2);
					ft[a][i][j]	= f[a][i][j] - (f[a][i][j]-feq)/tau; //collision
				}
			}
		}
		//----------------------------------------------------------------------
		for (i=1; i<=nx; i++) //Streaming post-collision
		{
			for (j=1; j<=ny; j++)
			{
				for (a=0; a<q; a++)
				{
					ia = i+ex[a];
					ja = j+ey[a];

					//if (ia<1 )	{ ia = nx;  }
					//if (ia>nx)	{ ia = 1 ;  }

					f[a][ia][ja] = ft[a][i][j];
				}
			}
		}

		fx[1]		=0.0;
		fy[1]		=0.0;
		uwx=0.0;
		uwy=0.0;
		cnt=-1;
		for (i=1; i<=nx; i++) //BC
		{
			for (j=1; j<=ny; j++)
			{
				if (isn[i][j] == 0)
				{
					for (a=0; a<q; a++)
					{
						ia = i+ex[a];
						ja = j+ey[a];

						if (isn[ia][ja]==1) //Cylinder
						{
							cnt=cnt+1;
							if(Delta[cnt]>=0.0 && Delta[cnt]<0.5)
							{
								chi=(2.0*Delta[cnt]-1.0)/(tau-2.0);
								ubfx = ux[i-ex[a]][j-ey[a]];
								ubfy = uy[i-ex[a]][j-ey[a]];
							}
							if(Delta[cnt]>=0.5 && Delta[cnt]<1.0)
							{
								chi=(2.0*Delta[cnt]-1.0)/(tau+0.5);
								ubfx = 0.5*(2.0*Delta[cnt]-3.0)/Delta[cnt]*ux[i][j] + 1.5*uwx/Delta[cnt];
								ubfy = 0.5*(2.0*Delta[cnt]-3.0)/Delta[cnt]*uy[i][j] + 1.5*uwy/Delta[cnt];
							}

							tmp1				= ux[i][j]*ex[a]  + uy[i][j]*ey[a] ;
							tmp2				= ux[i][j]*ux[i][j] + uy[i][j]*uy[i][j];
							feq					= wt[a]*rho[i][j]*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2);
							ft[kb[a]][ia][ja] = ft[a][i][j] - chi*(ft[a][i][j] - feq) + wt[a]*rho[i][j]*3.0*(ex[a]*(ubfx-ux[i][j]-2.0*uwx) + ey[a]*(ubfy-uy[i][j]-2.0*uwy));

							tmp1				= -ux[i][j]*ex[a]  - uy[i][j]*ey[a] ;
							tmp2				= ux[i][j]*ux[i][j] + uy[i][j]*uy[i][j];
							feq					= wt[a]*rho[i][j]*(1.0 + 3.0*tmp1 + 4.5*tmp1*tmp1 - 1.5*tmp2);
							ft[a][ia][ja] = ft[kb[a]][i][j] - chi*(ft[kb[a]][i][j] - feq) + wt[a]*rho[i][j]*3.0*(-ex[a]*(ubfx-ux[i][j]-2.0*uwx) - ey[a]*(ubfy-uy[i][j]-2.0*uwy));

							fx[1] += ex[a]*(ft[a][ia][ja] + ft[kb[a]][ia-ex[a]][ja-ey[a]])*(1.0-wb[ia-ex[a]][ja-ey[a]]);
							fy[1] += ey[a]*(ft[a][ia][ja] + ft[kb[a]][ia-ex[a]][ja-ey[a]])*(1.0-wb[ia-ex[a]][ja-ey[a]]);

							f[kb[a]][i ][j ] = ft[a		 ][i ][j ];
							f[a    ][ia][ja] = ft[kb[a]][ia][ja];
							//
							// fx[1] += ex[a]*2.0*(-ft[kb[a]][ia][ja] + ft[a][i][j]);
							// fy[1] += ey[a]*2.0*(-ft[kb[a]][ia][ja] + ft[a][i][j]);
						}
					}
				}
			}
		}


		//----------------------------------------------------------------------
		for(j=1,i=1;i<=nx;i++)
		{
			f[5][i][j]=f[8][i][j-1];
			f[2][i][j]=f[4][i][j-1];
			f[6][i][j]=f[7][i][j-1];
		}

		for(j=ny,i=1;i<=nx;i++)
		{
			f[8][i][j]=f[5][i][j+1];
			f[4][i][j]=f[2][i][j+1];
			f[7][i][j]=f[6][i][j+1];
		}

		for(j=1;j<=ny;j++)
		{
			//Inlet
			i=1;
			ux[i][j] = uInlet;
			uy[i][j] = 0.0;
			rho[i][j] = (f[0][i][j] + f[2][i][j] + f[4][i][j] + 2.0*(f[3][i][j] + f[6][i][j] + f[7][i][j]))/(1-ux[i][j]);
			tmp2 = ux[i][j]*ux[i][j] + uy[i][j]*uy[i][j];
			dudx = (-3.0*ux[i][j] + 4.0*ux[i+1][j] - 1.0*ux[i+2][j])/2.0; //dudx
			if(j == 0)
			{
				dudy = (-3.0*ux[i][j] + 4.0*ux[i][j+1] - 1.0*ux[i][j+2])/2.0;
			}
			else if(j == ny)
			{
				dudy = (3.0*ux[i][j] - 4.0*ux[i][j-1] + 1.0*ux[i][j-2])/2.0;
			}
			else
			{
				dudy = (ux[i][j+1] - ux[i][j-1])/2.0;
			}

			dvdx = (-3.0*uy[i][j] + 4.0*uy[i+1][j] - 1.0*uy[i+2][j])/2.0;
			if(j == 0)
			{
				dvdy = (-3.0*uy[i][j] + 4.0*uy[i][j+1] - 1.0*uy[i][j+2])/2.0;
			}
			else if(j == ny)
			{
				dvdy = (3.0*uy[i][j] - 4.0*uy[i][j-1] + 1.0*uy[i][j-2])/2.0;
			}
			else
			{
				dvdy = (uy[i][j+1] - uy[i][j-1])/2.0;
			}
			for(a=0;a<9;a++)
			{
				Q_xx[a] = ex[a]*ex[a] - 1.0/invC;
				Q_yy[a] = ey[a]*ey[a] - 1.0/invC;
				Q_xy[a] = ex[a]*ey[a];
				f_neq = -wt[a]*invC*tau*(Q_xx[a]*dudx + Q_xy[a]*dudy + Q_xy[a]*dvdx + Q_yy[a]*dvdy);
				tmp1 = ex[a]*ux[i][j] + ey[a]*uy[i][j];
				feq = wt[a]*rho[i][j]*(1+tmp1*invC + 0.5*tmp1*tmp1*invC*invC - 0.5*invC*tmp2);
				f[a][i][j] = feq + rho[i][j]*f_neq;
			}

			//Outlet
			i = nx;
			rho[i][j] = 1.0;
			ux[i][j]=(f[0][i][j]+f[2][i][j]+f[4][i][j]+2*(f[1][i][j]+f[5][i][j]+f[8][i][j]))/rho[i][j] - 1;
			tmp2 = ux[i][j]*ux[i][j] + uy[i][j]*uy[i][j];
			dudx = (3.0*ux[i][j] - 4.0*ux[i-1][j] + 1.0*ux[i-2][j])/2.0; //dudx
			if(j == 0)
			{
				dudy = (-3.0*ux[i][j] + 4.0*ux[i][j+1] - 1.0*ux[i][j+2])/2.0;
			}
			else if(j == ny)
			{
				dudy = (3.0*ux[i][j] - 4.0*ux[i][j-1] + 1.0*ux[i][j-2])/2.0;
			}
			else
			{
				dudy = (ux[i][j+1] - ux[i][j-1])/2.0;
			}

			dvdx = (3.0*uy[i][j] - 4.0*uy[i-1][j] + 1.0*uy[i-2][j])/2.0; //dvdx
			if(j == 0)
			{
				dvdy = (-3.0*uy[i][j] + 4.0*uy[i][j+1] - 1.0*uy[i][j+2])/2.0;
			}
			else if(j == ny)
			{
				dvdy = (3.0*uy[i][j] - 4.0*uy[i][j-1] + 1.0*uy[i][j-2])/2.0;
			}
			else
			{
				dvdy = (uy[i][j+1] - uy[i][j-1])/2.0;
			}
			for(a=0;a<9;a++)
			{
				Q_xx[a] = ex[a]*ex[a] - 1.0/invC;
				Q_yy[a] = ey[a]*ey[a] - 1.0/invC;
				Q_xy[a] = ex[a]*ey[a];
				f_neq = -wt[a]*invC*tau*(Q_xx[a]*dudx + Q_xy[a]*dudy + Q_xy[a]*dvdx + Q_yy[a]*dvdy);
				tmp1 = ex[a]*ux[i][j] + ey[a]*uy[i][j];
				feq = wt[a]*rho[i][j]*(1+tmp1*invC + 0.5*tmp1*tmp1*invC*invC - 0.5*invC*tmp2);
				f[a][i][j] = feq + rho[i][j]*f_neq;
			}
		}
		// for(i=1,j=1;j<=ny;j++)
		// {
		// 	rho[i][j]=(f[0][i][j]+f[2][i][j]+f[4][i][j]+2*(f[6][i][j]+f[3][i][j]+f[7][i][j]))/(1-uInlet);
		// 	f[1][i][j]=f[3][i][j]+((2.0/3.0)*rho[i][j]*uInlet);
		// 	f[5][i][j]=f[7][i][j]-(0.5*(f[2][i][j]-f[4][i][j]))+((1.0/6.0)*rho[i][j]*uInlet);
		// 	f[8][i][j]=f[6][i][j]+(0.5*(f[2][i][j]-f[4][i][j]))+((1.0/6.0)*rho[i][j]*uInlet);
		// }
		//
		// for(i=nx,j=1;j<=ny;j++)
		// {
		// 	ux[i][j]=(f[0][i][j]+f[2][i][j]+f[4][i][j]+2*(f[1][i][j]+f[5][i][j]+f[8][i][j]))/rho0 - 1;
		// 	f[3][i][j]=f[1][i][j]-((2.0/3.0)*rho0*ux[i][j]);
		// 	f[6][i][j]=f[8][i][j]-0.5*(f[2][i][j]-f[4][i][j])-((1.0/6.0)*rho0*ux[i][j]);
		// 	f[7][i][j]=f[5][i][j]+0.5*(f[2][i][j]-f[4][i][j])-((1.0/6.0)*rho0*ux[i][j]);
		// }
		//----------------------------------------------------------------------
		//----------------------------------------------------------------------
		fx_t 	= 0.5*(fx[0]  +  fx[1]);
		fy_t 	= 0.5*(fy[0]  +  fy[1]);

		Cd = fx_t/(0.5*rho0*uInlet*uInlet*D);
		Cl = fy_t/(0.5*rho0*uInlet*uInlet*D);
		//----------------------------------------------------------------------
		if (ts % dispFreq == 0)
		{
			fprintf(fid,"%d\t%10.6f\t%10.6f\t%10.6f\n", ts, rhoAvg, Cd, Cl);
			printf(    "%d\t%10.6f\t%10.6f\t%10.6f\n",	ts, rhoAvg, Cd, Cl);
		}
		//----------------------------------------------------------------------
		fx[0] = fx[1];
		fy[0] = fy[1];
		//----------------------------------------------------------------------
		if(ts<=time && ts % (time/(noOfSnaps-1)) == 0)
		{
			solnumber++;

			strcpy(filename,prefix);
			sprintf(solstr,"%d",solnumber);
			strcat(filename,solstr);
			strcat(filename,type);
			soln = fopen(filename,"w");

			fprintf(soln,"Variables=x,y,u,v,rho,region,wb\n");
			fprintf(soln,"Zone I= %d,J= %d\n\n",nx,ny);

			for (j=1; j<=ny; j++)
			{
				for (i=1; i<=nx; i++)
				{
					fprintf(soln,"%d %d %12.8f %12.8f %12.8f %d %d\n",i,j,ux[i][j],uy[i][j],rho[i][j],isn[i][j],wb[i][j]);
				}
				fprintf(soln,"\n");
			}
			fclose(soln);
			printf("snap %d recorded at time %d\n",solnumber,ts);
		}
		//----------------------------------------------------------------------
	}//Time loop Ends

	fclose(fid);
	return 0;
}
