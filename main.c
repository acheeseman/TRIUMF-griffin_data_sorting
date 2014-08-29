/* high rate data sorting of long waveforms contained in csv files

Implements Baseline Follower and Moving Window Deconvolution algorithms for pulse height evaluation

To compile: gcc -lm -g -o grifsort2 main.c gnuplot.c levenberg-marquardt.c
To run: ./grifsort2 ../inputfile.csv (can be multiple inputfiles) */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gnuplot.h"
#include "levenberg_marquardt.h"

// file properties - UPDATE THESE to match the properties of the input files you are using
#define RECORD_LENGTH 0.20     // length of each file (in s)
#define NUM_SAMPLES 20000000   // number of samples in each file
#define NUM_LINES_HEADER 19    // number of lines at the top of the file which do not contain waveform data
#define MAX_PULSES 100000      // maximum number of pulses to include in spectrum

// energy gates - if you want to look at pulses in a specific energy range only
#define EN_GATE_LOW 2350
#define EN_GATE_HIGH 2370

// signal processing
#define GAMMA_ENERGY 1332      // energy of the peak to fit - to calibrate the spectrum
#define RISETIME 20            // risetime of the signal (~200ns for HPGe)
#define INTEGRATION_DELAY 50   // how long after threshold crossing (this should maybe be changed to CFD time?) time to begin integration
#define WINDOW_WIDTH 600       // L parameter of MWD
#define INTEGRATION_TIME 530   // K parameter of MWD
#define DECAY_CONSTANT 5000    // M parameter of MWD - also used for fitting exponential baseline
#define GAUSSIAN_PEAK_WIDTH 15 // half the number of samples to include when fitting the peak 

// baseline follower
#define BLR_SPEED 64.0         // speed of the baseline restorer
#define BLR_DELAY 5000         // number of samples to calculate baseline before detecting pulses

// triggering
#define TRIG_THRESHOLD 80      // threshold crossing height to detect a pulse
#define TRIG_DIFF_WIDTH 4      // differentiation width when detecting triggers
#define TRIG_LENGTH 5          // number of samples above threshold height to detect a pulse

// peak finding
#define PEAK_NUM 2             // which peak to find in the spectrum
#define PEAK_THRESHOLD 40      // theshold above which to detect peaks
#define PEAK_RANGE 50          // range in which to search for the peak centroid

// processing states
#define BASELINE 0
#define PULSE 1

// store time of previous, current and next triggers
#define PREV 0
#define CURR 1
#define NEXT 2

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))
#define length(x) sizeof(x)/sizeof(x[0])

// curve fitting
extern double *p;
extern double *goodness_of_fit;
extern double *p_errors;

double gauss(int x, double p[3]);
double (*gauss_dp[3]) (int x, double p[3]);
double dgauss_dp0(int x, double p[3]);
double dgauss_dp1(int x, double p[3]);
double dgauss_dp2(int x, double p[3]);

double exp_1(int x, double p[2]);
double (*exp_1dp[2]) (int x, double p[2]);
double dexp_1dp0(int x, double p[2]);
double dexp_1dp1(int x, double p[2]);

// Gaussian function and partial derivatives for fitting
double gauss(int x, double p[3])
{
  double y;
  y = p[0] * exp(-((x-p[1])*(x-p[1]))/(2.*p[2]*p[2]));
  return y;
}
double dgauss_dp0(int x, double p[3])
{
  double y;
  y = exp(-((x-p[1])*(x-p[1]))/(2.*p[2]*p[2]));
  return y;
}
double dgauss_dp1(int x, double p[3])
{
  double y;
  y = p[0]*((x-p[1])/(p[2]*p[2]))*exp(-((x-p[1])*(x-p[1]))/(2.*p[2]*p[2]));
  return y;
}
double dgauss_dp2(int x, double p[3])
{
  double y;
  y = p[0]*(((x-p[1])*(x-p[1]))/(p[2]*p[2]*p[2]))*exp(-((x-p[1])*(x-p[1]))/(2.*p[2]*p[2]));
  return y;
}

// Exponential function and partial derivatives for fitting
double exp_1(int x, double p[2])
{
  double y;
  y = p[0] * exp(-(1.0/DECAY_CONSTANT)*x) + p[1];
  return y;
}

double dexp_1dp0(int x, double p[2])
{
  double y;
  y = exp(-(1.0/DECAY_CONSTANT)*x);
  return y;
}

double dexp_1dp1(int x, double p[2])
{
  double y;
  y = 1.0;
  return y;
}

// function declarations
float pulseHeightFloat(float input[WINDOW_WIDTH+RISETIME*2], int K, int trig_time);
int peakFind(int input[], int len, int peakNum, int thresh, int width);

// global variables
int pulse_ctr;
int total_pulses;
int pulse_height[MAX_PULSES];
int *pulse_height_histo;
int min, max, bin;

int main(int argc, char *argv[])
{
  int i;
  int peak_centre, peak_amplitude, peak_start;
  double gaussian_peak[2*GAUSSIAN_PEAK_WIDTH+1];
  double params_gauss[3] = {0};  // intitial parameter guesses for gaussian fit

  // zero variables
  min = 10000;
  max = 0;
  pulse_ctr = 0;
  total_pulses = 0;
  memset(pulse_height, 0, sizeof(pulse_height));
	  
  // process all input files
  for (i = 1; i < argc; i++)
    {
      printf("Processing file %s\n", argv[i]);
      process_file(argv[i], INTEGRATION_TIME, DECAY_CONSTANT);
    }
      
  // create pulse height spectrum from all files
  pulse_height_histo = (int*) calloc(max+1, sizeof(int));
  for (i=0;i<total_pulses;i++)
    {
      bin = pulse_height[i];
      if ( bin > 0 )
	{
	  pulse_height_histo[bin]++;
	}
    }  
  // plot spectrum to screen 
  plotHisto(pulse_height_histo, max, 0, max+10, 0);

  //for Co-60 spectrum we want to find the 2nd peak (1332 keV) and fit this
  peak_start = peakFind(pulse_height_histo, max+1, PEAK_NUM, PEAK_THRESHOLD, PEAK_RANGE);
  if ( peak_start < 0) { return -1; }
  else
    {
      peak_centre = 0;
      peak_amplitude = 0;
      for (i=peak_start;i<peak_start+GAUSSIAN_PEAK_WIDTH*2+1;i++)
	{
	  if (pulse_height_histo[i] > peak_amplitude)
	    {
	      peak_amplitude = pulse_height_histo[i];
	      peak_centre = i;
	    }
	}
      // store the peak we want to fin in an array of (GAUSSIAN_PEAK_WIDTH*2+1) samples
      for (i=0;i<length(gaussian_peak);i++)
	{
	  gaussian_peak[i] = (double) pulse_height_histo[peak_centre-GAUSSIAN_PEAK_WIDTH+i];
	}
      free(pulse_height_histo);

      //set array of partial derivatives for Gaussian fitting function
      gauss_dp[0] = dgauss_dp0;
      gauss_dp[1] = dgauss_dp1;
      gauss_dp[2] = dgauss_dp2;
      //set initial parameter guess for Gaussian fit
      params_gauss[0] = peak_amplitude*5.5;
      params_gauss[1] = GAUSSIAN_PEAK_WIDTH;
      params_gauss[2] = 0.5;
      // fit Gaussian function to peak
      fitData(gaussian_peak, 2*GAUSSIAN_PEAK_WIDTH+1, params_gauss, length(params_gauss), gauss, gauss_dp, 1);
      
      //calculate FWHM from fit parameters, write summary to stdout (FWHM = 2.35 * sigma)
      printf("Finished processing %d files\n", argc);
      printf("Total counts in spectrum = %d and total triggers = %d\n", pulse_ctr, total_pulses);
      printf("BLR speed: %f, BLR delay: %d\n", BLR_SPEED, BLR_DELAY);
      fprintf(stdout, "%d, %d, %g\n", INTEGRATION_TIME, DECAY_CONSTANT, (2.35*p[2])/(p[1]-GAUSSIAN_PEAK_WIDTH+peak_centre)*GAMMA_ENERGY);
    }
  return 0;
}

/* Processes each file, detecting pulses and evaluating their pulse height.
Pulse heights for each waveform are stored in the array pulse_height.

filename: the name of the pile to process
K: the integration width for MWD
M: the decay constant for MWD
 */
int process_file(char *filename, int K, int M)
{
  int i, j, k;
  int trigger;  // whether or not there is a trigger at the current location in the signal (1 = trigger, 0 = no trigger)
  int triggers[3] = {0};  // store the time of previous, current and next triggers
  int num_pulses_in_window, baseline;
  int baseline_ctr = 0;  // counts the number of samples of baseline which have been sampled from this file
  float slope, baseline_avg1, baseline_avg2, deficit;
  char buffer[200];
  FILE *inputfile;
  
  //arrays to store waveforms in various stages of processing
  float *waveform_input;
  waveform_input = (float*) calloc(NUM_SAMPLES, sizeof(float));
  float waveform_corrected[WINDOW_WIDTH+RISETIME*2] = {0};
  float waveform_temp;

  if ((inputfile = fopen(filename, "r")) == NULL)
    {
      printf("ERROR: Cannot open input file %s.\n", filename);
      return 1;
    }
  else
    {
      i = 0;
      fseek(inputfile, 0, SEEK_SET); // make sure we are at the beginning of the file
      // skip header of file which does not contain waveform data
      while (i<NUM_LINES_HEADER)
	{
	  fgets(buffer, 200, inputfile);
	  i++;
	}
      // read in waveform from input file
      for (i=0;i<NUM_SAMPLES;i++)
	{
	  // waveforms are in csv format "time, waveform_sample"
	  fscanf(inputfile,"%*g%*c%f\n", &waveform_input[i]);
	  waveform_input[i] *= 10000; // gain
	}
      fclose(inputfile);  
      i = WINDOW_WIDTH; // the window with is lost off the beginning of the signal due to the difference filter
      baseline = 0; // initially we assume the baseline to be zero and then calculate deviations from this.

      int state = BASELINE;
      while(i<NUM_SAMPLES){
	switch(state){
	case BASELINE: 
	  // check if we have exceeded the max pulses to be counted in the spectrum and if we have, exit the loop
	  if (pulse_ctr >= MAX_PULSES) {i = NUM_SAMPLES+1;}
	  
	  // follow baseline of corrected (differentiation and decay correction filters applied) signal
	  deficit = 0;
	  for (j=i-WINDOW_WIDTH;j<i;j++)
	    {
	      deficit += (waveform_input[i] - waveform_input[i-WINDOW_WIDTH]) - (baseline/BLR_SPEED);
	    }
	  waveform_temp = waveform_input[i] - waveform_input[i-WINDOW_WIDTH] + ((1/(float) M) * deficit);

	  // compare corrected signal (waveform_temp) to 0 and adjust the baseline accordingly
	  if (waveform_temp > 0) { baseline++; }
	  else if (waveform_temp < 0) { baseline--; }

	  // fast differentiation of signal for triggering  
	  for (j=0;j<TRIG_LENGTH;j++)
	    {
	      if (waveform_input[i+j] - waveform_input[i-TRIG_DIFF_WIDTH+j] >TRIG_THRESHOLD) { trigger = 1; }
	      else
		{
		  trigger = 0;
		  break;
		}
	    }
	  if ( trigger != 0)
	    {
	      total_pulses++; // count every trigger
	      triggers[PREV] = triggers[CURR];
	      triggers[CURR] = i;
	      // only evaluate pulses after BLR_DELAY so that the baseline error is properly established
	      if (baseline_ctr > BLR_DELAY) { state = PULSE;}
	      else { i += WINDOW_WIDTH+RISETIME*2; }
	    }
	  baseline_ctr++;  // count the number of samples of baseline
	  i++;
	  break;

	case PULSE:
	  //plotLine(waveform_input, i+WINDOW_WIDTH, i-RISETIME);
	  
	  // attempt to sort waveforms with effects of pileup from those without
	  //currently only the pulse height of waveforms with no effects of pilueup are evaluated and included in the spectrum

	  // find location of next trigger to determine if there is pileup
	  for (j=i+TRIG_DIFF_WIDTH+RISETIME*2;j<NUM_SAMPLES;j++)
	    {   
	      for (k=0;k<TRIG_LENGTH;k++)
		{
		  if ( waveform_input[j+k] - waveform_input[j-TRIG_DIFF_WIDTH+k] > TRIG_THRESHOLD) { trigger = 1; }
		  else
		    {
		      trigger = 0;
		      break;
		    }
		}
	      if (trigger != 0)
		{
		  triggers[NEXT] = j;
		  break;
		}
	      // if there is no next trigger we have reached the end of this trace
	      triggers[NEXT] = NUM_SAMPLES;
	    }
	  
	  // Calculate approximate slope of baseline over the differentiation width to determine baseline shift
	  baseline_avg1 = 0;
	  baseline_avg2 = 0;
	  for (j=0;j<50;j++)
	    {
	      baseline_avg1 += waveform_input[i-100+j];
	      baseline_avg2 += waveform_input[i-WINDOW_WIDTH-100+j];
	    }
	  baseline_avg1 /= 50;
	  baseline_avg2 /= 50;
	  slope = (baseline_avg1 - baseline_avg2);
	  //plotLine(waveform_input, i+500, i-WINDOW_WIDTH);

	  if (triggers[NEXT]-triggers[CURR] > WINDOW_WIDTH+RISETIME*2 && slope > -15.0)
	    {
	      // single pulse and no decay in baseline
	      //plotLine(waveform_input, i+4000-600, i-600);
	      for (j=0;j<length(waveform_corrected);j++)
		{
		  deficit = 0;
		  for (k=i-WINDOW_WIDTH;k<i;k++)
		    {
		      deficit += (waveform_input[k] - waveform_input[k-WINDOW_WIDTH]) - (baseline/BLR_SPEED);
		    }
		  waveform_corrected[j] = waveform_input[i] - waveform_input[i-WINDOW_WIDTH] + ((1/(float) M) * deficit);
		  i++;
		}
	      //plotLine(waveform_corrected, length(waveform_corrected), 0);
	      pulse_ctr++;
	      if ((pulse_height[pulse_ctr-1] = (int)(pulseHeightFloat(waveform_corrected, K, 50) +0.5)) < 0)
		{pulse_height[pulse_ctr-1] = 0; }
	      
	      // Find the minimum and maximum pulse heights
	      if (pulse_height[pulse_ctr-1] > max) { max = pulse_height[pulse_ctr-1]; }
	      if (pulse_height[pulse_ctr-1] < min) { min = pulse_height[pulse_ctr -1];}
	    }
	  else if (triggers[NEXT]-triggers[CURR] > WINDOW_WIDTH+RISETIME*2)
	    {
	      // single pulse with decaying baseline
	      //plotLine(waveform_input, i+500, i-600);
	      // fit the decay to determine the true baseline
	      i += WINDOW_WIDTH+RISETIME*2;
	    }
	  else
	    { 
	      // multiple pulses in window
	      // skip to next trigger
	      //plotLine(waveform_input, i+500, i-600);
	      i = triggers[NEXT]; 
	    }
	  state = BASELINE;
	  break;
	}
      }
      free(waveform_input);
      return 0;
    }
}

/*** SIGNAL PROCESSING FUNCTIONS ***/

/* Difference and Decay Correction - floating point */
int differenceFilterFloat(int L, float input[NUM_SAMPLES], float output[NUM_SAMPLES])
{
  int i;
  for (i=L;i<NUM_SAMPLES;i++)
    {
      output[i] = input[i] - input[i-L];
    }
  return 0;
}

int deconvolutionFilterFloat(int L, int M, float input[NUM_SAMPLES], float output[NUM_SAMPLES], int baseline)
{ 
  int i,j;
  float deficit;

  for (i=L;i<NUM_SAMPLES;i++)
    {
      deficit = 0.0;
      for (j=i-L;j<i;j++)
	{
	  deficit += input[j];
	}
      output[i] = input[i] + (1/ (float) M) * (deficit - baseline);
    }
  return 0;
  }

/* Returns pulse height of the filtered waveform 

input: input waveform of which we want to 
*/
float pulseHeightFloat(float input[], int K, int trig_time)
{
  //evaluate the energy for one waveform
  float pulseHeight = 0;
  int j;
  for (j=0;j<K;j++)
    {
      pulseHeight += (1/ ((float) K) ) * (input[trig_time+INTEGRATION_DELAY+j]);
    }
  return pulseHeight;
}

/* Finds the location of the nth  peak of an input array

input: array of data containing the spectrum in which you would like to search for peaks
len: len of the input array
peakNum: which peak you would like to find in the spectrum (ie. peakNum = 2 will return the location of 2nd peak)
thresh: threshold height above which to detect peaks (should be slightly higher than background)
width: number of points in which to search for the centroid of the peak
*/
int peakFind(int input[], int len, int peakNum, int thresh, int width)
{
  int peakList[100] = {0};
  int peakWidthHalf = (int) MAX(1,floor(width*0.5));
  int max;
  int minIndx, maxIndx, peakIndx;
  int i, j, n;
  n = 0;
  
  for (i=500;i<len;i++)
    {
      if (input[i] > thresh)
	{
	  minIndx = i - peakWidthHalf;
	  maxIndx = i + peakWidthHalf;
	  if ( minIndx < 0 ) { minIndx = 0; }
	  if ( maxIndx > len ) { maxIndx = len; }
	  max = 0;
	  for (j=minIndx;j<maxIndx;j++)
	    {
	      if ( input[j] > max ) { 
		max = input[j];
		peakIndx = j; }
	    }
	  if ( peakIndx == i ) { peakList[n++] = i; }
	}
    }
      if (n < peakNum)
	{
	  printf("Could not find peak number %d\n", peakNum);
	  return -1;
	}
  return peakList[peakNum-1];
}
