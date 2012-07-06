//----------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
//----------------------------------------------------------------------
const int AVE = 1<<12; // Number of Samples
const double h_gamma = 0.1; // Damping coefficient
const double T = 1.0; // Aimed Temperature
const double dt = 1.0; // Time step
const int LOOP = 1 << 14; // Total Time steps
const double D = sqrt(2.0*h_gamma*T/dt);
double data[LOOP];
double fft_data[LOOP];
//----------------------------------------------------------------------
double
myrand(void){
  return static_cast<double>(rand())/static_cast<double>(RAND_MAX);
}
//----------------------------------------------------------------------
double
get_gauss(void){
  const double r1 = myrand();
  const double r2 = myrand();
  return sqrt(-2.0*log(r1))*cos(M_PI*2.0*r2);
}
//----------------------------------------------------------------------
void
observe(fftw_complex *in, fftw_complex *out){
 double p = 0.0;
  double ave = 0.0;
  double var = 0.0;
  for(int i=0;i<LOOP;i++){
    double r = get_gauss()*D;
    p += (-h_gamma * p + r)*dt;
    data[i] = p;
    ave += p;
    var += p*p;
  }
  ave /= static_cast<double>(LOOP);
  var = var - ave*ave*static_cast<double>(LOOP);
  var /= static_cast<double>(LOOP-1);
  const double std = sqrt(var);
  for(int i=0;i<LOOP;i++){
    in[i][0] =  (data[i] - ave)/std;
    in[i][1] = 0.0;
  }

  fftw_plan plan = fftw_plan_dft_1d(LOOP, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(plan);
  for(int i=0;i<LOOP;i++){
    double re = out[i][0];
    double im = out[i][1];
    in[i][0] = (re*re + im*im)/static_cast<double>(LOOP);
    in[i][1] = 0.0;
  }
  fftw_plan plan2 = fftw_plan_dft_1d(LOOP, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(plan2);
  for(int i=0;i<LOOP;i++){
    const double re = out[i][0]/static_cast<double>(LOOP);
    fft_data[i] += re;
  }
}
//----------------------------------------------------------------------
int
main(void){
  for(int i=0;i<LOOP;i++){
    fft_data[i] = 0.0;
  }
  fftw_complex *in;
  fftw_complex *out;
  in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*LOOP);
  out  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*LOOP);
  for(int i=0;i<AVE;i++){
    observe(in,out);
  }
  for(int i=0;i<LOOP;i++){
    fft_data[i] /= static_cast<double>(AVE);
    const double t = static_cast<double>(i)*dt;
    printf("%f %f\n",t,fft_data[i]);
  }
  fftw_free(in);
  fftw_free(out);
}
//----------------------------------------------------------------------
