#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
//----------------------------------------------------------------------
#define nx 100
#define ny 30
#define q 9
//----------------------------------------------------------------------
const int time = 10000;
const int noOfSnaps = 5;
const int dispFreq = 100;
const double tau = 0.8;
const double invTau = 1.0 / tau;
const double rho0 = 1.0;
const double uWall = 0.08;
//----------------------------------------------------------------------
int kdf(int i, int j);

int main(void)
{
  int i, j, a, a1, ts, ia, ja;
  static double f[q][nx + 2][ny + 2], ft[q][nx + 2][ny + 2], wt[q], fneq[q], ux_a[ny + 2];
  static double ux[nx + 2][ny + 2], uy[nx + 2][ny + 2], rho[nx + 2][ny + 2], sigma[2][2], u[2];
  double tmp1, tmp2, tmp3, feq, rhoAvg;
  static int isn[nx + 2][ny + 2], ex[q], ey[q], kb[q], ci[q][2];

  char prefix[] = "snap_", type[] = ".dat", filename[15], solstr[5];
  int solnumber = 0;
  FILE *soln;
  //----------------------------------------------------------------------
  for (j = 1; j <= ny; j++)
  {
    ux_a[j] = uWall * (j - 0.5) / 30.0;
  }
  //----------------------------------------------------------------------
  ex[0] = 0;
  ey[0] = 0;
  ex[1] = 1;
  ey[1] = 0;
  ex[2] = 0;
  ey[2] = 1;
  ex[3] = -1;
  ey[3] = 0;
  ex[4] = 0;
  ey[4] = -1;
  ex[5] = 1;
  ey[5] = 1;
  ex[6] = -1;
  ey[6] = 1;
  ex[7] = -1;
  ey[7] = -1;
  ex[8] = 1;
  ey[8] = -1;

  for (int a = 0; a < q; a++)
  {
    ci[a][0] = ex[a];
    ci[a][1] = ey[a];
  }

  //----------------------------------------------------------------------
  for (a = 0; a < 9; a++)
  {
    if (a == 0)
    {
      wt[a] = 4.0 / 9.0;
    }
    if (a >= 1 && a <= 4)
    {
      wt[a] = 1.0 / 9.0;
    }
    if (a >= 5 && a <= 8)
    {
      wt[a] = 1.0 / 36.0;
    }
  }
  //----------------------------------------------------------------------
  for (a = 0; a < q; a++)
  {
    for (a1 = a; a1 < q; a1++)
    {
      if (ex[a] + ex[a1] == 0 && ey[a] + ey[a1] == 0)
      {
        kb[a] = a1;
        kb[a1] = a;
      }
    }
  }
  //----------------------------------------------------------------------
  for (i = 1; i <= nx; i++)
  {
    for (j = 0; j <= ny + 1; j++)
    {
      isn[i][j] = 0;
      if (j == 0)
      {
        isn[i][j] = 1;
      }
      if (j == ny + 1)
      {
        isn[i][j] = 2;
      }
    }
  }
  //----------------------------------------------------------------------
  for (i = 1; i <= nx; i++)
  {
    for (j = 0; j <= ny + 1; j++)
    {
      for (a = 0; a < 9; a++)
      {
        f[a][i][j] = wt[a] * rho0;
      }
    }
  }
  //----------------------------------------------------------------------
  for (ts = 0; ts <= time; ts++)
  {
    rhoAvg = 0.0;

    for (i = 1; i <= nx; i++)
    {
      for (j = 1; j <= ny; j++)
      {
        tmp1 = 0.0;
        tmp2 = 0.0;
        tmp3 = 0.0;

        for (a = 0; a < q; a++)
        {
          tmp1 += f[a][i][j];
          tmp2 += f[a][i][j] * ex[a];
          tmp3 += f[a][i][j] * ey[a];
        }

        rho[i][j] = tmp1;
        ux[i][j] = tmp2 / tmp1;
        uy[i][j] = tmp3 / tmp1;
        rhoAvg += rho[i][j];

        for (a = 0; a < q; a++)
        {
          tmp1 = ux[i][j] * ex[a] + uy[i][j] * ey[a];
          tmp2 = ux[i][j] * ux[i][j] + uy[i][j] * uy[i][j];
          feq = wt[a] * rho0 * (1.0 + 3.0 * tmp1 + 4.5 * tmp1 * tmp1 - 1.5 * tmp2);
          ft[a][i][j] = f[a][i][j] - (f[a][i][j] - feq) / tau; // collision step
        }
      }
    }

    rhoAvg /= (nx * ny);
    //----------------------------------------------------------------------
    for (i = 1; i <= nx; i++) // Streaming post Collision
    {
      for (j = 0; j <= ny + 1; j++)
      {
        if (isn[i][j] == 0)
        {
          for (a = 0; a < q; a++)
          {
            ia = i + ex[a];
            ja = j + ey[a];

            if (ia < 1)
            {
              ia = nx;
            }
            if (ia > nx)
            {
              ia = 1;
            }

            f[a][ia][ja] = ft[a][i][j];
          }
        }
      }
    }
    //----------------------------------------------------------------------
    for (i = 1; i <= nx; i++) // BC
    {
      for (j = 1; j <= ny; j++)
      {
        if (isn[i][j] == 0)
        {
          for (a = 0; a < q; a++)
          {
            ia = i + ex[a];
            ja = j + ey[a];

            if (ia < 1)
            {
              ia = nx;
            }
            if (ia > nx)
            {
              ia = 1;
            }

            if (isn[ia][ja] == 1) // Lower Wall
            {
              f[kb[a]][i][j] = f[a][ia][ja];
            }

            if (isn[ia][ja] == 2) // Upper Wall
            {
              f[kb[a]][i][j] = f[a][ia][ja] + 6.0 * wt[a] * uWall * ex[kb[a]];
            }
          }
        }
      }
    }

    for (int ii = 0; ii < 2; ii++)
    {
      for (int jj = 0; jj < 2; jj++)
      {

        i = nx / 2;
        j = ny / 4;
        u[0] = ux[i][j];
        u[1] = uy[i][j];

        for (a = 0; a < q; a++)
        {
          tmp1 = ux[i][j] * ex[a] + uy[i][j] * ey[a];
          tmp2 = ux[i][j] * ux[i][j] + uy[i][j] * uy[i][j];
          feq = wt[a] * rho[i][j] * (1.0 + 3.0 * tmp1 + 4.5 * tmp1 * tmp1 - 1.5 * tmp2);
          fneq[a] = f[a][i][j] - feq;
          // printf("%d %10.6f\n", a, fneq[a]);
        }

        // tmp1 = 0.0;
        // for (a = 0; a < q; a++)
        // {
        //   tmp1 = tmp1 + (ci[a][ii] - u[ii]) * (ci[a][jj] - u[jj]) * f[a][i][j];
        // }

        // sigma[ii][jj] = -(1.0 / 6.0) * invTau * rho[i][j] * kdf(ii, jj) - (1.0 - 0.5 * invTau) * tmp1;

        tmp1 = 0.0;
        for (a = 0; a < q; a++)
        {
          tmp1 = tmp1 + fneq[a] * (ci[a][ii] * ci[a][jj]); // - 0.5 * (ci[a][0] * ci[a][0] + ci[a][1] * ci[a][1]) * kdf(ii, jj));
        }

        sigma[ii][jj] = (1.0 - 0.5 * invTau) * tmp1;

        // printf("%d %d %d %d", kdf(0, 0), kdf(0, 1), kdf(1, 0), kdf(1, 1));
        // return 0;
        // sigma[ii][jj] = -(1.0 - 0.5 * invTau) * tmp1;
      }
    }
    sigma[0][0] = sigma[0][0] - (1.0 / 3.0) * rho[i][j];
    sigma[1][1] = sigma[1][1] - (1.0 / 3.0) * rho[i][j];
    //----------------------------------------------------------------------
    if ((ts % dispFreq) == 0)
    {
      printf("ts = %d\t rhoAvg = %10.6f\t stress=%12.8f%12.8f%12.8f%12.8f\n", ts, rhoAvg, sigma[0][0], sigma[0][1], sigma[1][0], sigma[1][1]);
    }
    //----------------------------------------------------------------------
    if (ts <= time && ts % (time / (noOfSnaps - 1)) == 0)
    {
      solnumber++;

      strcpy(filename, prefix);
      sprintf(solstr, "%d", solnumber);
      strcat(filename, solstr);
      strcat(filename, type);
      soln = fopen(filename, "w");

      fprintf(soln, "Variables=x,y,u,u_a,v,rho,region\n");
      fprintf(soln, "Zone I= %d,J= %d\n\n", nx, ny);

      for (j = 1; j <= ny; j++)
      {
        for (i = 1; i <= nx; i++)
        {
          fprintf(soln, "%d %d %12.8f %12.8f %12.8f %12.8f %d\n", i, j, ux[i][j], ux_a[j], uy[i][j], rho[i][j], isn[i][j]);
        }
        fprintf(soln, "\n");
      }
      fclose(soln);
      printf("snap %d recorded\n", solnumber);
    }
    //----------------------------------------------------------------------
  } // Time loop Ends

  return 0;
}

int kdf(int i, int j)
{
  if (i == j)
    return 1;
  else
    return 0;
}
