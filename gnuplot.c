#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gnuplot.h"

/* Plots the input signal as a line plot in gnuplot */
int plotLine(float input[], int len, int start_sample)
{
  int i;
  FILE *pipe = popen("gnuplot -persist", "w");
  // save plot to file
  //fprintf(pipe, "set terminal png\n");
  //fprintf(pipe, "set output '/home/acheeseman/2014/fast_differentiation.png'\n");
  fprintf(pipe, "set style line 1 lw 1 lc rgb 'black'\n");
  fprintf(pipe, "set xlabel 'Time (10 ns samples)'\n");
  fprintf(pipe, "set ylabel 'Charge (ADC Units)'\n");
  fprintf(pipe, "set xtics nomirror\n");
  fprintf(pipe, "set ytics nomirror\n");
  //fprintf(pipe, "set xrange [%d:]\n", start_sample);
  //fprintf(pipe, "set yrange [-0.2:0.2]\n");
  fprintf(pipe, "set zeroaxis\n");
  fprintf(pipe, "plot '-' w l ls 1 notitle\n");
  for(i=start_sample;i<len;i++)
    {
      fprintf(pipe, "%.2f\n", input[i]);
    }
  fprintf(pipe, "e\n");
  fprintf(pipe, "pause mouse\n");
  fclose(pipe);
  return 0;
}

/* Plots two input signals on the same axes in gnuplot */
int plotLine2(float input1[], float input2[], int len1, int len2, int start_1, int start_2)
{
  int i;
  FILE *pipe = popen("gnuplot -persist", "w");
  // save plot to file
  //fprintf(pipe, "set terminal png\n");
  //fprintf(pipe, "set output '/home/acheeseman/high_rate/decay_correction/waveform_%d.png'\n", i);
  fprintf(pipe, "set style line 1 lw 1 lc rgb 'black'\n");
  fprintf(pipe, "set style line 2 lw 1 lc rgb 'red'\n");
  fprintf(pipe, "set xlabel 'Time (10 ns samples)'\n");
  fprintf(pipe, "set xrange [%d:]\n", start_1);
  fprintf(pipe, "plot '-' w l ls 1 notitle, '-' w l ls 2 notitle\n");
  for(i=start_1;i<len1;i++)
    {
      fprintf(pipe, "%d %.2f\n", i, input1[i]);
    }
  fprintf(pipe, "e\n");
  for(i=start_2;i<len2;i++)
    {
      fprintf(pipe, "%d %.2f\n", i, input2[i]);
    }
  fprintf(pipe, "e\n");
  fprintf(pipe, "pause mouse\n");
  fclose(pipe);
  return 0;
}

/* Plots the input data as a box plot */
int plotHisto(int input[], int maximum, int x_start, int x_end, char *filename)
{
  int i;
  FILE *pipe = popen("gnuplot -persist", "w");

  if (x_end < maximum) { maximum = x_end; }

  fprintf(pipe, "reset\n");
  if (filename != NULL)
    {
      // save plot to file
      fprintf(pipe, "set terminal png size 1280,480\n");
      fprintf(pipe, "set output '%s'\n", filename);
    }
  fprintf(pipe, "set style fill solid 0.5\n");
  fprintf(pipe, "set yrange [0:]\n");
  fprintf(pipe, "set xrange [%d:%d]\n", x_start, x_end);
  fprintf(pipe, "set xtics 1000 nomirror out\n");
  fprintf(pipe, "set mxtics\n");
  fprintf(pipe, "set xlabel 'Charge'\n");
  fprintf(pipe, "set ylabel 'Counts'\n");
  fprintf(pipe, "set title 'Pulse Height Spectrum of Co-60 Source Evaluated with L=6us,K=5.3us,M=65us'\n");
  fprintf(pipe, "plot '-' w boxes notitle\n");
  for (i=0;i<=maximum;i++)
    {
      fprintf(pipe, "%d\t%d\n", i, input[i]);
    }
  fprintf(pipe, "e\n"); 
  fclose(pipe);
  return 0;
}
