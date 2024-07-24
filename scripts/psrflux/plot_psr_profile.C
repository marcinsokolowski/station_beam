// 
// Script for fitting pulse profile to scattered pulsar profiles like VELA:
// There are two options for the onset of the pulses:
//   - linear (fast up)
//   - gaussian
// 
// this is followed by exponetial decay 
// 

// always verify line :
// gROOT->Reset();
#include <TROOT.h>
#include <TF1.h>
#include <TH1F.h>
#include <TPaveText.h>
#include <TGraph.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h>
#include <TComplex.h>
#include <TFile.h>


int gLog=0;
int gVerb=0;
bool gNormaliseInputData=true;
double gPulsarPeriodMS = 89.328385024;
double gOriginalRMS=1.00;
double gOriginalMean=0.00;
double gSigmaSimulated=0.1120;
double gNoiseStart=0.00;
double gNoiseEnd=1e20;
bool   gUseFitResidualsRMS=false;
char   gInputFileName[256];

#define MAX_ROWS 10000000

Double_t Amplitude = 0.00;

// No scattering tau = 0.050
// Scattering tau = 2 seconds ?
Double_t Pulse_with_linear_onset( Double_t* x, Double_t* y )
{
   Double_t t = x[0];

   Double_t offset = y[0];
   Double_t t0 = y[1];           // start of the pulse
   Double_t t_peak = y[2];
   Double_t peak_flux = y[3];
   Double_t tau = y[4];

   double flux = offset;
   if( t > t0 ){
      if( t < t_peak ){
         double a = peak_flux / ( t_peak - t0 );
         flux = a*(t - t0 );
      }else{
         // exponential decay
         flux = peak_flux*exp( - (t-t_peak) / tau );
      }
   }


   return flux;

   // FRB171020 , 200 Jy ms @ 1.4 GHz 
   //             1500 Jy ms @ 185 MHz and spectral index  = -1
   //             11400 Jy ms @ 185 MHz and spectral index = -2
   //             7600 Jy ms @ 185 MHz and spectral index  = -1.8
/*   Double_t fluence = 7.60;  // Jy x sec assuming spectral index  = -1.8

   Double_t T   = 5;          // seconds - duration
   Double_t tau = y[1];      // 50ms - WARNING, when <0.050 - the integral value is incorrect !!!
   Amplitude = fluence / tau; // assuming pulse : Amplitude x exp( - t / tau )

   if( t < t0 ){
     return 0;
   }

   double pulse = Amplitude* exp( - (t-t0) / tau );
   return pulse;*/
}


Double_t Pulse_with_gauss_onset( Double_t* x, Double_t* y )
{
   Double_t t = x[0];

   Double_t offset = y[0];
//   Double_t t0 = y[1];           // start of the pulse
   Double_t t_peak = y[1];
   Double_t peak_flux = y[2];
   Double_t sigma = y[3];
   Double_t tau = y[4];


   double flux = offset;
//   if( t > t0 ){
//      double gaussian = (1.00/sqrt(sigma*2*TMath::Pi()))*exp(-0.5*(x - t_peak)*(x - t_peak)/(sigma*sigma) )

   double norm = peak_flux*(1.00/sqrt(sigma*2*TMath::Pi()));
   double gaussian = norm*exp(-0.5*(t - t_peak)*(t - t_peak)/(sigma*sigma) );

   if( t < t_peak ){
      flux = gaussian;
   }else{
      // exponential decay
      flux = norm*exp( - (t-t_peak) / tau );
   }
//   }


   return flux;
}

Double_t Pulse_gauss( Double_t* x, Double_t* y )
{
   Double_t t = x[0];

   Double_t offset = y[0];
   Double_t t_peak = y[1];
   Double_t peak_flux = y[2];
   Double_t sigma = y[3];

   double flux = offset;
   double norm = peak_flux*(1.00/sqrt(sigma*2*TMath::Pi()));
   double gaussian = norm*exp(-0.5*(t - t_peak)*(t - t_peak)/(sigma*sigma) );

   flux = gaussian;

   return flux;
}




Double_t HorizontalLine( Double_t* x, Double_t* y )
{
   return y[0];
}

Double_t Line( Double_t* x, Double_t* y )
{
   Double_t a = y[0];
   Double_t b = y[1];

   //a = (483.269-33121.78)/(68.97-841.87);
   //b = 483.269 - a*68.97;

   return (x[0]*a+b);
   // return (x[0]*y[0]);
}

TGraphErrors* DrawGraph( Double_t* x_values, Double_t* y_values, int numVal, 
         long q, TPad* pPad, const char* fit_func_name=NULL, 
         double min_y=-10000, double max_y=-10000,
         const char* szStarName="", const char* fname="default",
         int bLog=0, const char* szDescX=NULL, const char* szDescY=NULL,
         double fit_min_x=-100000, double fit_max_x=-100000, 
         Double_t* y_values_errors=NULL )
{
    int MarkerType = 20;
    int ColorNum = kRed;

    Double_t z,sigma_z,const_part;
    TF1 *line = NULL;
    TF1 *line_draw = NULL;
   TF1 *part1 = NULL;
   TF1 *part2 = NULL;
    Double_t maxX=-100000,maxY=-100000;
    Double_t minX=100000,minY=100000;


    TGraphErrors* pGraph = new TGraphErrors(q);
    for(int i=0;i<numVal;i++){
        if( gVerb ){
           printf("q=%d %f %f\n",(int)q, x_values[i], y_values[i] );
        }

        pGraph->SetPoint( i, x_values[i], y_values[i] );

        if( y_values_errors ){
           pGraph->SetPointError( i, 0.00, y_values_errors[i] );
        }

        if(x_values[i]>maxX)
            maxX = x_values[i];
        if(y_values[i]>maxY)
            maxY = y_values[i];
      
        if(x_values[i]<minX)
            minX = x_values[i];
        if(y_values[i]<minY)
            minY = y_values[i];
    }
    printf("Found min_x=%.2f , max_x=%.2f\n",minX,maxX);
    printf("Found min_y=%.2f , max_y=%.2f\n",minY,maxY);
    Double_t stepX = (maxX-minX)/10.00;
    Double_t stepY = (maxY-minY)/10.00;
    printf("Found  stepX=%.2f , stepY=%.2f\n",stepX,stepY);

    // pPad->DrawFrame(minX-stepX,minY-stepY,maxX+stepX,maxY+stepY);
   if( min_y!=-1 && max_y!=-1 ){
      // pPad->DrawFrame(minX-stepX,min_y,maxX+stepX,max_y);
   }else{
      // pPad->DrawFrame(minX-stepX,0,maxX+stepX,maxY+stepY);
   }
  //    pPad->SetGrid();   


   double r=1.00;
   double norm = 1;
   double mean = 1;
   double sigma = 1;


    pGraph->SetMarkerStyle(MarkerType);
    pGraph->SetMarkerColor(ColorNum);
    if( min_y>-10000 && max_y>-10000 ){
       pGraph->SetMinimum( min_y );
       pGraph->SetMaximum( max_y );
       printf("Set Min/Max = (%.2f,%.2f)\n",min_y,max_y);
    }else{      
       double minVal = minY-stepY;
       if( bLog ){
          if( minVal<0 ){
             minVal = 1;
          }
       }

       pGraph->SetMinimum( minVal );
       pGraph->SetMaximum( maxY+stepY );
       printf("Set Min/Max = (%.2f,%.2f)\n",minVal,(maxY+stepY));
    }
//    pGraph->SetMinimum( -90 );
//    pGraph->SetMaximum( 90 );
    pGraph->SetTitle( szStarName );
// WORKS 
//    pGraph->GetXaxis()->SetLimits(0,1);
    pGraph->Draw("APL");

   if( fit_min_x<=-100000 ){
      fit_min_x = minX;
   }
   if( fit_max_x<=-100000 ){
      fit_max_x = maxX;
   }

   printf("fitting in x range (%.2f,%.2f)\n",fit_min_x,fit_max_x);

   Double_t par[10];
   par[0] = 0;
   par[1] = 11.00;
   par[2] = 0.0;
   par[3] = 0.0;

   int local_func=0;
   if( fit_func_name && strlen(fit_func_name) ){
      if( strcmp( fit_func_name, "pulse" )==0 ){
         printf("Fitting pulsar profile in the range : %.6f - %.6f\n",fit_min_x,fit_max_x);
         line = new TF1("fit_func",Pulse_with_linear_onset,fit_min_x,fit_max_x,5);
         line_draw = new TF1("fit_func2",Pulse_with_linear_onset,minX,maxX,5);
         local_func=1;

         par[0] = 0;
         par[1] = 11.00;
         par[2] = 15.0;
         par[3] = 17.0;
         par[4] = 7.00; // very long decay ...

         if( gNormaliseInputData ){
            par[1] = par[1] / numVal;
            par[2] = par[2] / numVal;
         }else{
            par[0] = 15685190000000.000000;
            par[3] = 15.692*1e12;
         }

         line->SetParName(0,"Offset");
         line->SetParName(1,"#t_{start}");
         line->SetParName(2,"#t_{peak}");
         line->SetParName(3,"Peak_flux");
         line->SetParName(4,"#tau");
      }


      if( strcmp( fit_func_name, "pulse_gauss" )==0 ){
         printf("Fitting pulsar pulse_gauss profile in the range : %.6f - %.6f\n",fit_min_x,fit_max_x);
         line = new TF1("fit_func",Pulse_with_gauss_onset,fit_min_x,fit_max_x,5);
         line_draw = new TF1("fit_func2",Pulse_with_gauss_onset,minX,maxX,5);
         local_func=1;

         par[0] = 0.00; // offset (mean off-pulse)
         par[1] = 0.2; // t_peak  0.53 for the MWA data
         par[2] = 5; // peak flux 
         par[3] = 0.112; // sigma 
         par[4] = 0.3; // very long decay ... - 10 for MWA fit 
         par[5] = 0.5; // tau 

/*         if( gNormaliseInputData ){
            par[1] = par[1] / numVal;
            par[2] = par[2] / numVal;
         }else{
            par[0] = 15685190000000.000000;
            par[3] = 15.692*1e12;
         }*/

         line->SetParName(0,"SNR offset");
         line->SetParName(1,"t_{peak}");
         line->SetParName(2,"f_{peak}");
         line->SetParName(3,"#sigma_{gauss}");
         line->SetParName(4,"#tau");

//         line->SetParName(1,"T_start");
//         line->SetParName(4,"#tau");

      }

      if( strcmp( fit_func_name, "pulse_gauss_only" )==0 ){
         printf("Fitting pulsar pulse_gauss profile in the range : %.6f - %.6f ( %.8f - %.8f )\n",fit_min_x,fit_max_x,fit_min_x,fit_max_x);
         line = new TF1("fit_func",Pulse_gauss,fit_min_x,fit_max_x,4);
         line_draw = new TF1("fit_func2",Pulse_gauss,minX,maxX,4);
         local_func=1;
         
         double phase_max=-1, max_value=-1;
         for(int i=0;i<numVal;i++){
            if( y_values[i] > max_value ){
                max_value = y_values[i];
                phase_max = x_values[i];
            }
         }

         par[0] = 0.00; // offset (mean off-pulse)
         par[1] = phase_max;
         par[2] = int(max_value);
//         par[1] = 0.75; // t_peak  0.53 for the MWA data
//         par[2] = 8; // peak flux 
         par[3] = 0.112; // sigma 
         printf("DEBUG : initial values for pulse_gauss_only fitting : %.4f, %.4f, %.4f, %.4f\n",par[0],par[1],par[2],par[3]);

/*         if( gNormaliseInputData ){
            par[1] = par[1] / numVal;
            par[2] = par[2] / numVal;
         }else{
            par[0] = 15685190000000.000000;
            par[3] = 15.692*1e12;
         }*/

         line->SetParName(0,"SNR offset");
         line->SetParName(1,"t_{peak}");
         line->SetParName(2,"f_{peak}");
         line->SetParName(3,"#sigma_{gauss}");

//         line->SetParName(1,"T_start");
//         line->SetParName(4,"#tau");

      }



      if( strcmp( fit_func_name, "horline" )==0 || fit_func_name[0]=='h' ){
         printf("Fitting horizontal line\n");
         line = new TF1("fit_func",HorizontalLine,fit_min_x,fit_max_x,1);
         line_draw = new TF1("fit_func2",HorizontalLine,minX,maxX,1);
         local_func=1;
      }
      if( strcmp( fit_func_name, "line" )==0 || fit_func_name[0]=='l' ){
         printf("Fitting straight line\n");
         line = new TF1("fit_func",Line,fit_min_x,fit_max_x,2);
         line_draw = new TF1("fit_func2",Line,minX,maxX,2);
         local_func=1;
      }
   }


   if( local_func ){
      line->SetParameters(par);
   }
   if( fit_func_name && strlen(fit_func_name) ){
      printf("fitting function : %s",fit_func_name);

      if( strcmp(fit_func_name,"gaus")==0 || fit_func_name[0]=='g' || fit_func_name[0]=='G' ){
         pGraph->Fit("gaus");
      }

      if( strstr(fit_func_name,"line") || fit_func_name[0]=='l' || fit_func_name[0]=='L'
          || fit_func_name[0]=='h' || fit_func_name[0]=='H' || strstr(fit_func_name,"pulse") || fit_func_name[0]=='p' || strstr(fit_func_name,"pulse_gauss") || strstr(fit_func_name,"pulse_gauss_only")
        ){
         pGraph->Fit("fit_func","R");
         printf("DEBUG : fitted fit_func (%s)\n",fit_func_name);
      }

      if( local_func ){
         line->GetParameters(par);
//         line->Draw("same");
         line_draw->SetParameters(par);
         line_draw->Draw("same");
         line_draw->SetLineColor(kBlack);

         printf("Fitted paramaters:\n");
         printf("\tpar[0] = %.8f\n",par[0]);
         printf("\tpar[1] = %.8f\n",par[1]);
         printf("\tpar[2] = %.8f\n",par[2]);
         printf("\tpar[3] = %.8f\n",par[3]);
         printf("TEST VALUE = %.8f vs. %.8f\n",line->Eval(0.5),line_draw->Eval(0.5));
      }


      double sum_resid=0.00;
      double sum_resid2=0.00;
      int count_resid=0;
      bool bIgnorePeak=true;
      FILE* out_f = fopen("fit.txt","w"); 
      for(int i=0;i<numVal;i++){
         bool bInclude=true;
         double x = x_values[i];
         double y = y_values[i];

         if( bIgnorePeak ){
            double dist_from_peak = fabs( x - par[1] );
            if( dist_from_peak <= 0.2 ){
               printf("WARNING : ignoring phase bin %.8f - too close to the peak at %.8f\n",x,par[1]);
               bInclude=false;
            }
         }

         double fit = line_draw->Eval(x);
         double resid = (y-fit);
         fprintf(out_f,"%.8f %.8f %.8f %.8f\n",x,fit,y,resid);

         if( bInclude ){
            sum_resid += resid;
            sum_resid2 += resid*resid;
            count_resid += 1; 
         }
      }
      fclose(out_f);     
       

      double mean_resid = sum_resid / count_resid;
      double rms_resid = sqrt( sum_resid2/count_resid - mean_resid*mean_resid);
      printf("Checking residuals of the fit:\n");
      printf("\tMean residual = %.6f\n",mean_resid);
      printf("\tRMS  residual = %.6f\n",rms_resid);

      // calculate W10 
      int n_bins = numVal;
      double max_flux = -1000;
      double phase=0.00;
      double delta_phase=(1.00/n_bins); // 0.001;
      double fluence = 0.00;
      double sum_profile = 0.00; // fit 
      int n_bins_fit=0;
      while( phase < 1.00 ){ // was 2.0
         double flux = line_draw->Eval( phase );
         printf("DEBUG : phase = %.3f -> flux = %.4f\n",phase,flux);
         if( flux > max_flux ){
            max_flux = flux;            
         }
         fluence += flux*delta_phase;
         sum_profile += flux;

         phase += delta_phase;
         n_bins_fit += 1;
      }
      printf("DEBUG : max_flux = %.8f\n",max_flux);

      if( n_bins_fit != n_bins ){
         printf("WARNING : n_bins_fit = %d != %d = n_bins -> should stick to the same number of bins, as in the data !!!\n",n_bins_fit,n_bins);
      }

      double sum_profile_data = 0.00;
      for(int i=0;i<numVal;i++){
         if( x_values[i] >= gNoiseEnd ){
            sum_profile_data += y_values[i];
            // printf("DEBUG = %e\n",y_values[i]);
         }
      }
      double sum_profile_data_before_rescale =  sum_profile_data*gOriginalRMS;

      printf("Max flux = %.6f [units]\n",max_flux);
      printf("Mean flux = %.6f [units] (or %.6f) = %.8f / %.8f\n",(fluence/phase),(sum_profile/n_bins),fluence,phase);
      printf("Fluence  = %.6f [units*phase]\n",fluence);
      printf("Sum profile = %.8f\n",sum_profile);
      printf("Sum profile data = %e (before re-scaling = %e)\n",sum_profile_data,sum_profile_data_before_rescale);

      double w10_start=-1;
      double w10_end=-1;

      phase=0.00;
      double prev_phase=-1;
      double prev_ratio=-1;
      while( phase < 2.00 ){
         double flux = line_draw->Eval( phase );
         double ratio = flux/max_flux;
         if( ratio >= 0.1 && prev_ratio < 0.1 ){
            printf("Start phase for W10 = %.6f\n",phase);
            w10_start = phase;
         }
         if( ratio < 0.1 && prev_ratio >= 0.1 ){
            printf("End phase for W10 = %.6f\n",phase);
            w10_end = phase;
         }
         
         prev_phase = phase;
         prev_ratio = ratio;
         phase += delta_phase;
      }

      double w10 = (w10_end - w10_start);
      printf("W10 phase range = %.8f - %.8f\n",w10_start,w10_end);
      printf("W10 = %.6f [ms]\n",(w10*gPulsarPeriodMS));

//      int n_bins = numVal;
      double mean_flux = (sum_profile*gSigmaSimulated)/(gOriginalRMS*n_bins);
      double peak_flux_jy = max_flux*(gSigmaSimulated/gOriginalRMS);

      if( gNormaliseInputData ){
         // if normalised by RMS -> sum_profile is already = sum_profile_ORIGNAL/gOriginalRMS -> no need to divide by gOriginalRMS again
         mean_flux = (sum_profile*gSigmaSimulated)/(n_bins);
         printf("DEBUG : mean_flux = %.8f * %.8f / %d\n",sum_profile,gSigmaSimulated,n_bins);

         // same applies to peak flux :
         peak_flux_jy = max_flux*gSigmaSimulated; // max_flux = max_flux_original/gOriginalRMS

         if( gUseFitResidualsRMS ){
            printf("INFO : using residuals of the FIT to calculate zeta\n");
            mean_flux = (sum_profile*gSigmaSimulated)/(n_bins*rms_resid);
            peak_flux_jy = max_flux*gSigmaSimulated/rms_resid;
         }
      }

      // double norm = peak_flux*(1.00/sqrt(sigma*2*TMath::Pi()));
      double fitted_peak = par[2]/sqrt(2*TMath::Pi()*par[3]);
      printf("Fitted peak flux: SNR = %.8f [no units]\n",fitted_peak);
 
      printf("Mean flux density = %.6f Jy\n",mean_flux);
      printf("Peak flux density = %.6f Jy\n",peak_flux_jy);

      FILE* out2_f = fopen("final_flux.txt","a+");
      fprintf(out2_f,"%s %.8f %.8f\n",gInputFileName,mean_flux,peak_flux_jy);
      fclose(out2_f);

      printf("DEBUG - before re-scaling:\n");
      double mean_flux_no_rescale = (sum_profile_data_before_rescale*gSigmaSimulated)/(n_bins*gOriginalRMS);
      printf("\tDEBUG : Mean flux density without rescaling = %e\n",mean_flux_no_rescale);
      double zeta_debug = gSigmaSimulated/gOriginalRMS;
      printf("\tDEBUG : sum_profile_data*zeta / n_bins = %.8f * %.8f / %d = %.8f\n",sum_profile_data,zeta_debug,n_bins,(sum_profile_data*zeta_debug / n_bins));
      double mean_flux_data = (sum_profile_data*gSigmaSimulated)/(n_bins*gOriginalRMS);
      printf("\tDEBUG : Mean flux density data = %e\n",mean_flux_data);
      printf("\tDEBUG : zeta = %.8f\n",(gSigmaSimulated/gOriginalRMS));
   }   

   pGraph->GetXaxis()->SetTitleOffset(1.00);
   pGraph->GetYaxis()->SetTitleOffset(0.60);
   pGraph->GetXaxis()->SetTitleSize(0.05);
   pGraph->GetXaxis()->SetLabelSize(0.05);
   pGraph->GetYaxis()->SetTitleSize(0.05);
   pGraph->GetYaxis()->SetLabelSize(0.05);


   if( szDescX && szDescX[0] ){
      pGraph->GetHistogram()->SetXTitle( szDescX );
   }else{
      pGraph->GetHistogram()->SetXTitle("Channel");
   }
   if( szDescY && szDescY[0] ){
       pGraph->GetHistogram()->SetYTitle( szDescY );
   }else{
      pGraph->GetHistogram()->SetYTitle("Power");
   }
   // pGraph->GetHistogram()->SetYTitle("sigmaG1_homeo/SigmaG1");

   TLatex lat;
   lat.SetTextAlign(23);
   lat.SetTextSize(0.05);
                                                                                
   char szDesc[255];
   sprintf(szDesc,"sigma/mean = %.8f\n",r);
   
   return pGraph;
}

int ReadResultsFile( const char* fname, Double_t* x_values, Double_t* y_values,
                     int CondCol=-1, int CondValue=-1, int x_col=0, int y_col=0 )
{
   const int lSize=1000;
   char buff[1000];
   
   FILE *fcd = fopen(fname,"r");
   if(!fcd){
       printf("Could not open file : %s\n",fname);
       return -1;
   }

   printf("x_col=%d y_col=%d\n",x_col,y_col);

   Int_t all = 0;
   Double_t fval1,fval2,fval3,fval4,fval5,fval6,mag,x,y;
   long lval1,lval2,lval3,lval4;
   double max_val=-100000;
   double max_x = -100000;

   all=0;
   int ncols=-1;
   while (1) {
      if(fgets(buff,lSize,fcd)==0)
         break;
      if(buff[0]=='#')
         continue;      
// what if nan is in the column that we are not interested in 
//      if(strstr(buff,"nan"))
//         continue;

      //  9.39  377.000000  8.000000 328.757587  77.088256   298.312992 65.223146   44.506926
      // ncols = sscanf(buff,"%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",&mag,&x,&y,&fval1,&fval2,&fval3,&fval4,&fval5);
      Double_t x_val;
      Double_t y_val;
//      ncols = sscanf(buff,"%f\t%f\n",&x_val,&y_val);

      int ncols=0;
      char* ptr=NULL;
      char* search_ptr=buff;
      int col=0;
      while( (ptr = strtok(search_ptr," \t")) ){
         search_ptr = NULL;
         if( gVerb ){
              printf("ptr = %s\n",ptr);
                        }

         if( col == x_col ){
            if( strstr(ptr,"nan") ){
               break;
            }
            x_val = atof( ptr );
            ncols++;
         }
         if( col == y_col ){
            if( strstr(ptr,"nan") ){
               break;
            }
            y_val = atof( ptr );
            ncols++;
         }
         col++;
      }


      if( ncols < 2 ){
         // changed due to stupid plots , single column rows are skiped
         printf("ERROR - less then 2 columns, line skiped !!!\n");
         printf("%s\n",buff);
         continue;
         // y_val = x_val;
         // x_val = all;
      }

      if (ncols < 0) break;
     
     if(CondCol>=0){
          if(lval1!=CondValue)
            continue;                                        
     }
//     printf("%s\n",buff);
     
     x_values[all] = x_val;
     y_values[all] = y_val;
     if( y_val > max_val ){
        max_val = y_val;
     }
     if( x_val > max_x ){
        max_x = x_val;
     }
     if( gVerb || 0 ){
        printf("values : %f %f\n",x_val,y_val);
     }

     all++;
   }
   fclose(fcd);

  

   printf("max_val = %.4f\n",max_val);
   if( max_x > 1 ){
      printf("max_x = %.8f -> normalising to 0 - 1 phase range\n");
      for(int i=0;i<all;i++){
         x_values[i] = x_values[i]/max_x;
      }
    }

//   exit(0);   
   return all;
}  

double calc_rms( Double_t* x_values, Double_t* y_values, int cnt, double start=-1, double end=-1 )
{
   double sum=0,sum2=0;
   double maxX=-1000000000;
   int count=0;

   for(int i=0;i<cnt;i++){
      double x_val = x_values[i];

      if( (start<0 || x_val >= start) && (end<0 || x_val <= end) ){
         sum += y_values[i];
         sum2 += (y_values[i])*(y_values[i]);
         count++;
      }

      if( x_val > maxX ){
         maxX = x_val;
      }
   }

   double rms = sqrt( sum2/count - (sum/count)*(sum/count) );
   double mean = (sum/count);
   printf("MEAN = %e\n",mean);
   printf("RMS  = %e\n",rms);
  
   return rms;
}


double normalise( Double_t* x_values, Double_t* y_values, int cnt, double start=-1, double end=-1, bool bSubtractBaseline=false )
{
   double sum=0,sum2=0;
   double maxX=-1000000000;
   double maxY=-1000000000;
   int count=0;

   for(int i=0;i<cnt;i++){
      double x_val = x_values[i];

      if( (start<0 || x_val >= start) && (end<0 || x_val <= end) ){
         sum += y_values[i];
         sum2 += (y_values[i])*(y_values[i]);
         count++;
      }    

      if( x_val > maxX ){
         maxX = x_val;
      }
      if( y_values[i] > maxY ){
         maxY = y_values[i];
      }
   }

   double rms = sqrt( sum2/count - (sum/count)*(sum/count) );
   double mean = (sum/count);
   printf("MEAN = %e (count = %d)\n",mean,count);
   printf("RMS  = %e (count = %d)\n",rms,count);
   printf("MAX  = %.8f\n",maxY);

   if( bSubtractBaseline ){
      for(int i=0;i<cnt;i++){
         y_values[i] = (y_values[i] - mean)/rms;
//         x_values[i] = x_values[i]/maxX;
      }
   }

   gOriginalMean = mean;
   gOriginalRMS = rms;

   return rms;
}

void plot_psr_profile( const char* basename="sigmaG1_vs_lapSigmaG1_for_root", bool bNormaliseInputData=true, 
                       const char* fit_func_name="pulse", // pulse, pulse_gauss, pulse_gauss_only
                       double noise_start=0, double noise_end=10, 
                       double sigma_simulated=0.1120, // simulated sigma of noise in Jy , sigma_Stokes_I - for the entire duration of the observation !!!
                       bool bIsSigmaSimulPerPhaseBin=false, // if this is true no need to multiply sigma_simulated by sqrt(n_bins)
                       bool bUseFitResidualsRMS=false, // use residuals of the FIT to calculate Sigma_obs which may be slightly different than 1.00 (after normalisation)
                       double min_y=-1000000,  double max_y=-1000000, int bLog=0,
      const char* szDescX="Phase",
      const char* szDescY="Signal to Noise Ratio (SNR)", 
      const char* szTitle=NULL,
      double fit_min_x=-100000, double fit_max_x=-100000,
      int x_col=0, int y_col=1, const char* outpngfile=NULL )
{
   if( !szTitle){
      szTitle = basename;
   }
   gLog = bLog;
   gNormaliseInputData = bNormaliseInputData;
   gNoiseStart = noise_start;
   gNoiseEnd   = noise_end;
   gUseFitResidualsRMS = bUseFitResidualsRMS;
   strcpy(gInputFileName,basename);
   
   // gROOT->Reset();
   // const char* basename = "s_vs_sigma_g_sqr";

   TCanvas* c1 = new TCanvas("c1","xxxx",200,10,700,500);
   c1->SetGridx();
   c1->SetGridy();
   c1->SetFillColor(0);
   c1->SetFillStyle(0);
   if( bLog ){
      c1->SetLogy(1);
   }
   
   gStyle->SetPadTopMargin(0.03);
   gStyle->SetPadBottomMargin(1.2);
   gStyle->SetPadLeftMargin(0.15);
   gStyle->SetPadRightMargin(0.04);
   gStyle->SetTitleXOffset(1.12);
   gStyle->SetTitleYOffset(1.2);
   gStyle->SetFillColor(0);
   gStyle->SetFillStyle(0);
   gStyle->SetOptFit(111);

   // p1 = new TPad("p1","p1",0,0,1,1,-1);
   // p1->Draw();

   


   Double_t* x_value1 = new Double_t[MAX_ROWS];
   Double_t* y_value1 = new Double_t[MAX_ROWS];
   Double_t* y_value1_err = new Double_t[MAX_ROWS];
   Double_t maxX=-100000,maxY=-100000;   
   Double_t minX=100000,minY=100000;

   Int_t ncols;
   Int_t lq1=0,lq2=0,lq3=0,lq5=0,lq9=0,lq25=0;



   lq1 = ReadResultsFile( basename, x_value1, y_value1, -1, -1, x_col, y_col ); 
   
   int n_bins = lq1;
   printf("DEBUG : number of bins = %d\n",n_bins);   
   gSigmaSimulated = sigma_simulated*sqrt(n_bins); // multiply by sqrt(n_bins) to calculated simulated sigma for single time bin in phase !!!
   if( bIsSigmaSimulPerPhaseBin ){ 
      // however if sigma provided is already per phase bin total_obs_time/n_bins seconds -> no need to multiply by number of phase bins:
      gSigmaSimulated = sigma_simulated;
      printf("DEBUG : using sigma_simul as provided (per phase bin) = %.6f [Jy]\n",gSigmaSimulated);
   }else{
      gSigmaSimulated = sigma_simulated*sqrt(n_bins); // multiply by sqrt(n_bins) to calculated simulated sigma for single time bin in phase !!!
      printf("DEBUG : using provided sigma_simul x sqrt(n_bins=%d) = %.6f * %.6f = %.6f [Jy]\n",n_bins,sigma_simulated,sqrt(double(n_bins)),gSigmaSimulated);
   }
   

   double rms = normalise( x_value1, y_value1, lq1, noise_start, noise_end, bNormaliseInputData );
   printf("DEBUG : zeta0 = %e\n",(gSigmaSimulated/gOriginalRMS));
   if( bNormaliseInputData ){
      for(int i=0;i<lq1;i++){
         y_value1_err[i] = 1.00;
      }
   }else{
     y_value1_err = NULL;
   }
   printf("Control RMS after re-scaling:\n");
   double rms_new = calc_rms( x_value1, y_value1, lq1, noise_start, noise_end );

   
   // drawing background graphs here :
   TGraphErrors* pGraph1 = DrawGraph( x_value1, y_value1, lq1, 1, NULL, fit_func_name, min_y, max_y, szTitle,
                                      basename, bLog, szDescX, szDescY, fit_min_x, fit_max_x, y_value1_err );
   
   c1->Update();

   TString szPngName1="images/";
   if( outpngfile ){
      szPngName1 += outpngfile;
   }else{
      szPngName1 += basename;
      szPngName1 += ".png";
   }
   c1->Print(szPngName1.Data());
}
