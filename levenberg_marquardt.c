/* Levenberg-Marquardt Algorithm for curve fitting
   Can fit any function given the input function and its partial derivatives

   y_dat = data set to fit
   np = number of data points in y_dat
   params = array containing initial parameter guesses
   k = number of parameters in fitting function
   function = pointer to fitting function
   partials = partials of fitting function wrt parameters
   plot = an integer which determines whether or not to plot the fit ( fit is plotted for any non-zero value )
*/

/* gcc -lm -g -o levenberg_marquardt levenberg_marquardt.c */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "levenberg_marquardt.h"

double jac[MAX_POINTS][MAX_PARAMS];
double jacT[MAX_PARAMS][MAX_POINTS];
double H[MAX_PARAMS][MAX_PARAMS];
double d[MAX_PARAMS][MAX_PARAMS];
double r[MAX_POINTS];
double delta[MAX_PARAMS];
double *p;
double *f;
double *goodness_of_fit;
double *p_errors;
double lambda;
const double lambda_0 = pow(10,-2);		// initial value of lambda

int final;

double dotProd(double a[], double b[], int len);

fitData(double y_dat[], int len_y, double params_guess[], int k, double (*function)(int, double[]), double (*partials[])(int, double[]), int plot)
{
  
  int i,j;
  int iteration = 0;
  int deg_freedom = len_y-k-1;

  goodness_of_fit = (double*)calloc(NUM_TESTS, sizeof(double));
  p_errors = (double*)calloc(k, sizeof(double));

  final = 0;
  lambda = lambda_0;

  /* f = fitting function
     p = parameters */
  f = (double*)calloc(len_y, sizeof(double));
  p = (double*)calloc(len_y, sizeof(double));
  for (i=0;i<k;i++)
    {
      p[i] = params_guess[i];
    }
   
  // compute the simple average of the data set
  double y_bar = 0.0;	    
  for (i=0;i<len_y;i++)
    {
      y_bar += y_dat[i];
    }
  y_bar /= len_y;
  
 
  /* Procedure for each iteration */
  while (iteration < MAX_ITERATIONS && final == 0)
    {

      for (i=0;i<len_y;i++)
	{
	  f[i] = (*function)(i, p);
	  r[i] = y_dat[i] - f[i];
	}
      updateJac(len_y, k, partials);
      solveDelta(y_dat, len_y, k);
      checkConvergence(y_dat, len_y, k, function, iteration);
      iteration++;
    }
  // Update arrays containing fit statistics
  goodnessOfFit(goodness_of_fit, y_dat, len_y, k, y_bar);
  paramErrors(p_errors, y_dat, len_y, k, function);
  //printStats(goodness_of_fit, p_errors, deg_freedom, k, iteration);

  if (plot != 0)
    {
      // plot fit
      FILE *pipe = popen("gnuplot -persist", "w");
      fprintf(pipe, "reset\n");
      //fprintf(pipe, "set terminal png size 1280,480\n");
      //fprintf(pipe, "set output '/home/acheeseman/2014/Plots/fit_BLR4.png'\n");
      fprintf(pipe, "set style fill solid 0.5\n");
      fprintf(pipe, "set style line 2 lt 1 lw 2 pt 3 lc rgb 'blue'\n");
      fprintf(pipe, "set xtics out nomirror\n");
      fprintf(pipe, "set mxtics 10\n");
      fprintf(pipe, "plot '-' w boxes notitle, '-' w l ls 2 t 'Gaussian Fit'\n");
      for (i=0;i<len_y;i++)
	{
	  fprintf(pipe, "%d\t%g\n", i, y_dat[i]);
	}
      fprintf(pipe, "e\n");
      for (i=0;i<len_y;i++)
	{
	  fprintf(pipe, "%d\t%g\n", i, f[i]);
	}
      fprintf(pipe, "e\n");
      //fprintf(pipe, "pause mouse\n");
      fclose(pipe);
    }
  
  free(f);
  free(p);
  free(p_errors);
  free(goodness_of_fit);
  return 0;
}

/* Updates the Jacobian, its transpose, their product, and its
   diagonal for a Gaussian function */
updateJac(int len_y, int k, double (*partials[])(int, double[]))
{
 
  int i,j,m;
   for (i=0;i<len_y;i++)
    {
      for (j=0;j<k;j++)
	{
	  jac[i][j] = (*partials[j]) (i, p);
	  jacT[j][i] = (*partials[j]) (i, p);
	}
      //printf("%g\n", jac[i][1]);
    }
  
  // compute Hessian
  memset(H, 0, sizeof(H));
  for (i=0;i<k;i++)
    {
      for (j=0;j<k;j++)
	{
	  for (m=0;m<len_y;m++)
	    {
	      H[i][j] += jacT[i][m]*jac[m][j];
	    }
	}
    }
  // diag(H)
  for (i=0;i<k;i++)
    {
      for (j=0;j<k;j++)
	{
	  if ( i==j )
	    {
	      d[i][j] = H[i][j];
	    }
	  else
	    {
	      d[i][j] = 0.0;
	    }
	}
    }
  
  return 0;
}

/* Solves for delta (amount to increment the parameter by) from:
   [jTj + lambda*d]delta = jacT[y_dat - f(p)]                   */
solveDelta(double y_dat[], int len_y, int k)
{
  int i,j,m;
  double A[k][k];
  double b[k];
  double C[k][k+1];

  memset(A, 0, sizeof(A));
  memset(b, 0, sizeof(b));
  // Calculate A & b
   for (i=0;i<k;i++)
    {
      for (j=0;j<k;j++)
	{
	  A[i][j] = H[i][j] + (lambda*d[i][j]);
	  //printf("%g, %g\n", H[i][j], d[i][j]);
	}
    }
  for (i=0;i<k;i++)
    {
      for (j=0;j<len_y;j++)
	{
	  b[i] += jacT[i][j]*r[j];
	}
    }

  // Gaussian Elimination
  int max;
  double t;
  for (i=0;i<k;i++)
    {
      C[i][k] = b[i];
      for (j=0;j<k;j++)
	{
	  C[i][j] = A[i][j];
	}
    }
  for (i=0;i<k;i++) {
    max = i;
    for (j=i+1;j<k;j++)
      if (C[j][i]>C[max][i])
	max = j;
    for (j=0;j<k+1;j++) {
      t = C[max][j];
      C[max][j] = C[i][j];
      C[i][j] = t;
    }
    for (j=k;j>=i;j--)
      for (m=i+1;m<k;m++)
	C[m][j] -= C[m][i]/C[i][i] * C[i][j];
  }
  for (i=k-1;i>=0;--i)
    {
      C[i][k] = C[i][k]/C[i][i];
      C[i][i] = 1;
      for (j=i-1;j>=0;j--)
	{
	  C[j][k] -= C[j][i]*C[i][k];
	  C[j][i] = 0;
	}
    }
  // set delta
  for (i=0;i<k;i++)
    {
      delta[i] = C[i][k];
     }
  return 0;
}

/*Return X^2 value based on input parameters*/
double chiSquared(double y_dat[], int len_y, double f[])
{
  double x2 = 0;
  x2 = (1./2)*dotProd(y_dat,y_dat,len_y) - dotProd(y_dat,f,len_y) + (1./2)*dotProd(f,f,len_y);
  return x2;
}

/*Check for convergence, then check if X2(p)-X2(p+delta) > epsilon_3*deltaT(lambda*delta + jT(y-f(p)),
and update lambda accordingly                                                                         */
checkConvergence(double y_dat[], int len_y, int k, double (*function)(int, double[]), int iteration)
{
  int i,j;
  double p_new[k];
  double f_new[len_y];
  double x2_old, x2_new;

  const double epsilon_1 = pow(10,-6);	// convergence tolerance for gradient
  const double epsilon_2 = pow(10,-6);	// convergence tolerance for parameters
  const double epsilon_3 = pow(10,-5);	// convergence tolerance for Chi-square

  for (i=0;i<k;i++)
    {
      p_new[i] = p[i] + delta[i];
    }
   for (i=0;i<len_y;i++)
    {
      f_new[i] = (*function)(i,p_new);
    }

   double temp2[k];
   memset(temp2,0,sizeof(temp2));
   for (i=0;i<k;i++)
     {
       for (j=0;j<len_y;j++)
	 {
	   temp2[i] += jacT[i][j]*r[j];
	 }
     }
  
   // calculate Chi-squared values
   x2_old = chiSquared(y_dat, len_y, f);
   x2_new = chiSquared(y_dat, len_y, f_new);
   
   // test for convergence
   double max1 = 0;
   double max2 = 0;
   
   for (i=0;i<k;i++)
     {
       if (fabs(temp2[i]) > max1)
	 {
	   max1 = fabs(temp2[i]);   
	 }
     }
   for (i=0;i<k;i++)
     {
        if ((fabs(delta[i]/p[i])) > max2)
	 {
	   max2 = fabs((delta[i]/p[i]));
     	 }
     }

   if (max1 < epsilon_1 | max2 < epsilon_2 | x2_old/(len_y) < epsilon_3)
     {
         final = 1;
     }
   else
     {
       // test if new parameters are sufficiently better
       double t;
       //printf("epsilon: %g, delta: %g, temp2: %g, k: %d\n", epsilon_3, delta[0], temp2[0], k);
       t = epsilon_3 * dotProd(delta,temp2,k);
      
        if ((x2_old-x2_new) > t)
     {
       for (i=0;i<k;i++)
	 {
	   p[i] = p_new[i];
	 }
       lambda /= 10;
     }
       else
	 {
	   lambda *= 10;
	 }
     }
 return 0;
}

/* Updates an array containing 'goodness of fit' statistics
 arr[0] = SST (Total Sum of Squares)
 arr[1] = SSR (Residual Sum of Squares)
 arr[2] = R^2 Value
 arr[3] = Adjusted R^2 Value */
goodnessOfFit(double arr[], double y_dat[], int len_y, int k, double y_bar)
{
  int i;

  for (i=0;i<len_y;i++)
    {
      arr[SSR] += (y_dat[i] - f[i])*(y_dat[i] - f[i]);
      arr[SST] += (y_dat[i] - y_bar)*(y_dat[i] - y_bar);
    }
  arr[R2] = 1-(arr[SSR]/arr[SST]);
  arr[R2_ADJ] = 1- ((1-arr[R2])*((len_y-1.0)/(len_y-k-1)));
  return 0;
}

/*Updates an array containing 1 sigma errors for parameters of fitting function */
paramErrors(double arr[], double y_dat[], int len_y, int k, double (*function)(int, double[]))
{
  int i,j;
  double f_new[len_y];
  double p_new[k];
  double delta = 1.0/1000000;
  double x2 = chiSquared(y_dat, len_y, f);
  double x2_new;
  double param;


  for (i=0;i<k;i++)
    {
      x2_new = x2;
      memcpy(p_new, p, k*sizeof(double));
      while ((x2_new-x2) < 1)
	{
	  p_new[i] += delta;
	  for (j=0;j<len_y;j++)
	    {
	      f_new[j] = (*function)(j,p_new);
	    }
	  x2_new = chiSquared(y_dat, len_y, f_new);
	  //printf("x2_new = %g\n", x2_new);
	}
      arr[i] = p_new[i] - p[i];
      
    }     
  return 0;
}

/*Prints a 'fit summary' to the screen */
printStats(double goodness_of_fit[], double p_errors[], double deg_freedom, int k, int iteration)
{
  int i;
  printf("\nAfter %d iterations the fit converged.\n\n", iteration);
  printf("Final Parameters are:\n");
  for (i=0;i<k;i++)
    {
    printf("p%d: %g +/- %g\n", i, p[i], p_errors[i]);
  }
  printf("Final sum of squares of residuals: %g\n", goodness_of_fit[SSR]);
 printf("Final explained sum of squares: %g\n", goodness_of_fit[SST] );
  printf("R^2 Value: %g\n", goodness_of_fit[R2]);
  printf("Adjusted R^2 Value: %g\n",goodness_of_fit[R2_ADJ]);
  printf("Degrees of Freedom: %d\n", deg_freedom);
  printf("Reduced Chi-square: %g\n", (goodness_of_fit[SSR])/(double) (deg_freedom));
}

/* Computes the dot product of two vectors, a and b, of length 'len' */
double dotProd(double a[], double b[], int len)
{
  double p = 0;
  if (sizeof(a) == sizeof(b))
    {
      int i;
      for (i=0;i<len;i++)
	{
	  p+= (a[i]*b[i]);
	}
      return p;
    }
  else
    {
      return 0;
    }
}


