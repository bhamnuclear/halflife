#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <sys/stat.h>
#include <termios.h>
#include <unistd.h>
#include <limits.h> /*/usr/include/limits.h for INT_MAX*/

#define CHMAX     32768      /*max number of channels in spectra*/
#define AUTO      1          /*set to one to enable auto mode to appear in menu*/
#define MAXCOLS   3          /*max. number data columns in input spectrum*/
#define MAXPTS    16384      /*max. channels in input spectrum*/
#define NUMOPT    13         /*number of options in get_mode() excluding quit*/
#define NUMPARS   4          /*number of parameters to fit*/
#define NUMFUNC   1          /*number of separate functions in fit profile excl. bkgnd*/
#define MITER     70         /*max. number of iterations for fitter routine*/
#define CHLEN     120        /*character length of filename arrays*/
#define CHFLEN    14         /*character length of function name arrays*/
#define CHPLEN    30         /*character length of parameter name arrays*/
#define WN        10         /*window size (in chans) for smoothing for finding spec max./min.*/
#define P_PROG    "xmgrace"  /*plotting program to be called from command line*/
#define P_OPT_BEG ""         /*plotting program options at beginning (before filenames)*/
/*#define P_OPT_END "-world 0 0 2050 1500 -pexec \"xaxis tick major 500\" -pexec \"xaxis tick minor ticks 4\""*/
#define P_OPT_FIN "-world XMIN YMIN XMAX YMAX -pexec \"xaxis tick major XTCKMAJ\" -pexec \"xaxis tick minor ticks XTCKMIN\""
                             /*plotting program options at end (after filenames)
                                XMIN, YMIN, XTICK etc will be substituted with spectrum values*/

/* Carl Wheldon July 2009 */

/* To compile:
gcc halflife.c -Wall -pedantic -o halflife -lm -O2
*/

/*%%%% A program to fit a Gaussian convoluted with an exponential decay %%%%*/

/*structure of the Ortec Maestro header as written to/read from a spectrum*/
struct  maest_header {
    short int q1;           /*must be one*/
    short int q2;           /*MCA/det number*/
    short int q3;           /*segment number*/
    short int q4;           /*ascii second of start time*/
    unsigned int real;      /*real time (increments of 20 ms)*/
    unsigned int lve;       /*live time (increments of 20 ms)*/
    char dt[8];             /*Start date as ASCII DDMMMYY*
                              The * character should be ignored if it is not
                              a "1". If it is a "1", it indicates the data is
                              after the year 2000.*/
    char sttm[4];           /*start time as ASCII HHMM or binary zeros,
                                if not known*/
    short int off;          /*channel offset of data*/
    short int channels;     /*length of data (channels)*/
} maest_header;

/*structure of the Ortec Maestro header as written to/read from a spectrum*/
struct  maest_trailer {
    short int t1;           /*must be 102*/
    short int t2;           /*reserved???*/
    float     g[3];         /*energy calib coeff offset, gain and quadratic term*/  
    char      trailer[496]; /*nothing particularly useful in the rest of the trailer*/
} maest_trailer;

struct ft {
    double pars[NUMPARS], errs[NUMPARS], y[MAXPTS], dy[MAXPTS], x[MAXPTS];
    int    freepars[NUMPARS], nfp, ndp, nad, ad[NUMPARS];
} ft;

/*system function declarations needed for -ansi compile flag*/
double  erf(double x);      /*/usr/include/bits/mathcall.h included via math.h*/
FILE    *popen(const char *command, const char *type);  /*/usr/include/stdio.h*/
int     pclose(FILE *stream);                           /*/usr/include/stdio.h*/

int 	eval(double *pars, double x, double *fit, double *derivs, int ch, int mode);
double  gauss(double x, double fwhm);
double  gauss_exp(double x, double thalf, double fwhm);

int     ask_yn(char *questn, int yn);
int 	ascii_read(char fname[], int *col);
void 	ascii_write(char name[], int numch, int col);
void 	chi2(double *chisq);
void 	col_determ(FILE *file, int *col);
void 	compress2(char inname[], char outname[]);
int 	cswap4(int decim);
int 	cswap2(int decim);
void 	file_status(char *fname, int len);
int 	find_strt_end_zeros(int numch);
int 	fitter(int, double *chisq, int vb);
void 	get_ans(char ans[], int num, int md);
void 	get_line(char ans[], int len);
void    get_line_file(FILE *file, char ans[], int len);
int 	get_mode(int md, int exp);
void 	get_pars(char ans0[], double pars[], int num);
void    get_pars_file(FILE *file, float pars[], int num);
void 	itoa(int n, char s[]);
int     maestro_read(char fname[], int *col);
int 	matinv(double *array, int nip, int npars);
void    ortec_head_srch_skip(char *fname, FILE *file, int *mxchan);
void    plot_data(char ans[], char cmd[], char inname[], char outname[],
            int plt, int mode);
int     plot_pid_get(char cmd[]);
int     plot_prog_find(char cmd[]);
void    plot_scale_set(char ans[], int numch);
void    plot_scale_set_par(char p_par[], double val, char ans[]);
void    prep_write(char inname[], char outname[], char nm[], char ans[],
            int numch, int col);
void 	read_par(double pars[], int offset, int num, char *s);
void 	reverse(char s[]);
void 	scaler(double *fit, double *derivs, int mode);
void 	set_ext(char fname[], char ext[]);
void 	skip_hash(FILE *file);
void    skip_lines(FILE *file, int lns);
void    spec_int(void);
void    spec_max(int mode);
void    spec_prop(double init[]);
void    store_colours();
void 	swapb2(char *buf);
void 	swapb4(char *buf);
void 	write_output(char inname[], char outname[], double init[],
	    double x[], double *fit, double chisq, int mode);

double fit, derivs[NUMPARS], bkgnd, fact[NUMPARS+1];
double x12[2], pi, spc[1], sp_max, sp_sum;
int gdfit, lim[2], pid[40], pcnt;
char par_nm[NUMPARS+1][CHPLEN], func_nm[NUMFUNC+1][CHFLEN], clr[10][12], P_OPT_END[CHLEN];

/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
/* ++++++++++++++++++++++++++++++++++ MAIN ++++++++++++++++++++++++++++++++++ */
/* ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
int main(int argc, char *argv[])
{
    extern double fit, derivs[NUMPARS], bkgnd, fact[NUMPARS+1], spc[1], x12[2];
    extern int gdfit, lim[2], pid[40], pcnt;
    extern char par_nm[NUMPARS+1][CHPLEN], func_nm[NUMFUNC+1][CHFLEN], clr[10][12];
    double  chisq = 0.0, init[NUMPARS+1], bklo = 0.0, bklo2, chisqo = 10000.0;
    int     f, i, md = 0, numch = 0, col = 0, flginp = 0, flginit = 0, r = -1;
    int     aut = 0, bkcnt, flgbkpm, plt = 0, strt = 1;
    char    ans[3*CHLEN] = "", inname[CHLEN] = "", outname[CHLEN] = "";
    char    cmd[3*CHLEN] = "", nofit[CHLEN] = "", nm[CHLEN] = "";
    char    smooth[CHLEN] = "";
    struct  stat statbuf;
    
    pi = (double)3.14159265359;
    
    /*Initially there is no good fit*/
    gdfit = 0;    
    /*store colours*/
    store_colours();

    /*copy parameter names to array*/
    strncpy(par_nm[0], "t_1/2", 5);
    strncpy(par_nm[1], "FWHM", 4);
    strncpy(par_nm[2], "Centroid", 8);
    strncpy(par_nm[3], "Scaling factor", 14);
    strncpy(par_nm[NUMPARS], "Background level (not fitted)", 29);
    
    /*copy function names to function components in eval()*/
    strncpy(func_nm[0], "Gaus-exp conv", 13);
    strncpy(func_nm[NUMFUNC], "Background", 10);
    
    /*initialise background parameter (constant offset) to zero*/
    bkgnd = 0.0;
    /*no. of data points*/
    ft.ndp = 0;
    /*no. of fixed parameters*/
    ft.nfp = 0;
    /*set all parameters to be free*/
    for (i = 0; i < NUMPARS; i++) ft.freepars[i] = 1;

    /*zero arrays*/
    for (i = 0; i < NUMPARS; i++)
    {
	ft.pars[i] = 0.0;
	init[i] = 0.0;
    }
    pcnt = 0;
    for (i = 0; i < 40; i++) pid[i] = 0;

    /*set which parameters have analytic derivatives*/
    for (i = 0; i < NUMPARS; i++) ft.ad[i] = 0;
    ft.ad[3] = 1;
    
    /*initialise converting factors from initial parameters to function parameters*/
    /*convert t_1/2 to -1*tau for use in eval()
      Note this is no longer used as gauss_exp accepts half-life
      fact[0] = -1.0/log(2);*/
    /*convert FWHM to sigma*sqrt(2) for use in eval()
      Note this is no longer used as gauss_exp accepts FWHM
      fact[1] = 0.6005612;*/
    
    /*initialise conversion factors from initial parameters to function parameters*/
    for (i = 0; i < NUMPARS; i++) fact[i] = 1.0;
    fact[NUMPARS] = 1.0;

    /*initialise no. of fixed parameters*/
    ft.nfp = 0;
    for (i = 0; i < NUMPARS; i++)
    {
        /*no. of fixed parameters*/
        if (ft.freepars[i] == 0) ft.nfp++;
        /*set number of analytic derivatives*/
        if (ft.ad[i] == 1) ft.nad++;
    }
    
    /*initialise background parameter (constant offset) to zero*/
    init[NUMPARS] = 0.0;
    bkgnd = init[NUMPARS];
    
    /*zero spectrum limits*/
    x12[0] = 0.0; x12[1] = 0.0;
    lim[0] = 0; lim[1] = 0;
    
    printf("\n\t**** Welcome to program %s ****\nThis fits"
    	   " a Gaussian convoluted with an exponential decay.\n"
	   " The data can be in Maestro (.Spe or .Chn) format\n"
           " or ASCII 1, 2 or 3 column format (i.e. y, x y or x y dy):\n"
           " Lines at any point in the file starting with a # will be skipped.\n\n",
              argv[0]);

    /*argv[i] is the ith argument, i.e. first is the program name*/
    if ( argc > 2 || (argc == 2 && (stat(argv[1], &statbuf))) )
    {
	printf("\nUnrecognised arguments...usage: %s\n"
		" or: %s InputFileName\n",argv[0],argv[0]);
	if (argc == 2) printf(" %s***File %s does not exist%s\n",
            clr[1],argv[1],clr[0]);
        
        printf("    %sExiting...%s\n\n",clr[1],clr[0]);
        return 0;
    }
    else if (argc == 2)
    {
	strcpy(inname,argv[1]);
	printf("Input filename: %s\n",inname);
        /*set input file flag*/
	flginp = 1;
    }
    
    /*check if plotting program exists plt = strlen(cmd) means yes, 0 no*/
    plt = plot_prog_find(cmd);
    
    while (1)
    {
	if ( strt == 1 || (md = get_mode(md,aut)) == 1 )
	{
    	    if (!flginp)
	    {
		while (stat(inname, &statbuf))
                {
                    printf("Type ASCII data filename\n");
		    if (strlen(inname) > 0 ) printf(" [<Enter> for %s, q to exit]\n",inname);
    	    	    get_line(nm, CHLEN);
    	    	    if (strlen(nm) == (int)0 && (nm[0] == ((char)0))) ;
    	    	    else strcpy(inname, nm);
                    if (nm[0] == 'q' && strlen(nm) == (int)1)
                    {
                        printf("    %sExiting...%s\n\n",clr[1],clr[0]);
                        return 0;
                    }
                    else if (stat(inname, &statbuf))
                        printf("%s***File %s does not exist%s\n",clr[1],inname,clr[0]);
                }
	    }
            flginp = 0;

    	    /*check if input file has .fit extension and warn if yes*/
            /*if not equal to NULL, find '.' then copy ext*/
            if ( strrchr(inname,'.') && !strcmp(strrchr(inname,'.'),".fit") )
                printf("%s\tWARNING, INPUT FILE HAS '.fit' EXTENSION!%s\n",clr[1],clr[0]);
            
    	    /*get output filename and set extension*/
    	    strncpy(outname, inname, CHLEN);
    	    set_ext(outname, "_full.fit");
    	    printf("Output filename for fit: %s\n", outname);
    	    file_status(outname, CHLEN);
            strncpy(nofit, outname, CHLEN);
            strncpy(smooth, outname, CHLEN);
            set_ext(smooth,"_func.fit");
            
    	    /*read in ascii data*/
            /*before checking for extension have to check if '.' is present*/
    	    if ( strrchr(inname,'.') && !strcmp( (strrchr(inname,'.')), ".Chn") )
                ft.ndp = maestro_read(inname, &col);
            /*if not assume spectrum is ASCII*/
            else
                ft.ndp = ascii_read(inname, &col);
            /*on file read error exit*/
            if (ft.ndp == -1) return -1;
            
            /*store full spectrum length including zeros*/
            numch = ft.ndp;
            
            /*if plotting prog exists, fill plot options with spectrum parameters
                i.e. X and Y min and max and major tick separation*/
            if (plt && strt) plot_scale_set(ans,numch);
	    strt = 0;
            
            /*check for zeros at start and end of spectrum and set upper and
                lower limits for fit (lim[0] and lim[1] respectively)*/
            ft.ndp = find_strt_end_zeros(ft.ndp);
            
            /*integrate spectrum*/
            spec_int();
                        
            /*if Ortec (.Spe or .Chn) file state (in blue highlighting)
                2 column ASCII data file will be written and used*/
            if (plt && ( (!strcmp( (strrchr(inname,'.')), ".Spe" )) ||
                    ( (!strcmp( (strrchr(inname,'.')), ".Chn" ) ) ) ) )
            {
                printf("%s***ORTEC format file found. Writing 2 col. ASCII data file%s\n"
                    " %sand using this as new input file for full plotting options with %s%s\n",
                        clr[3],clr[0],clr[3],P_PROG,clr[0]);
                /*prepare filenames and write data*/
                prep_write(inname, outname, nm, ans, numch, 2);
            }

            /*ask about opening data in P_PROG to aid choosing initial parameters
                providing P_PROG exists and file is not an ORTEC (.Spe or .Chn) file*/
            plot_data(ans, cmd, inname, outname, plt, 1);       
        }
	/*if auto mode, mode is set first to NUMOPT then straight to md = 2
            to call scaler*/
        else if (md == 2)
	{
            /******************************************************************/
            /*FOR TESTING*/
            /*init[0] = 150.0;
            init[1] = 50.0;
            init[2] = 920.0;
            init[4] = 1.0;
            init[5] = 20.0;*/
            /*FOR TESTING END*/
            /******************************************************************/
	    
            /*read in initial parameter values*/
            for (i = 0; i < NUMPARS+1; i++)
            {
                if (i == 3) continue;
                sprintf(ans, "Enter initial guess for %s\n",par_nm[i]);
	        /*if first pass (flginit = 0) and auto mode enabled don't ask*/
                if (!aut) read_par(init, i, 1, ans);
                if (i < NUMPARS) ft.pars[i] = init[i]*fact[i];
            }    	    

    	    /*the first time, call scaler function to get pars[3]*/
    	    if (!flginit)
            {      
	        ft.pars[3] = 1.0;
                scaler(&fit, derivs, 0);
	        init[3] = ft.pars[3]/fact[3];
            }
	    if (flginit && !aut)
	    {
    	    	read_par(init, 3, 1, "Initial guess for scaling factor");
	    	ft.pars[3] = init[3]*fact[3];
	    }
            /*set bkgnd after calling scaler as otherwise scaler fails*/
	    bkgnd = init[NUMPARS];

            if (aut) printf("%sAuto mode set.%s ",clr[3],clr[0]);
            printf("Initial parameter values:\n");
	    for (i = 0; i < NUMPARS+1; i++) printf("    %s = %.3f\n",par_nm[i],init[i]);
             
            /*set initial parameter flag*/
	    flginit = 1;
            /*for auto mode set mode to 4 for immediate fitting*/
            if (aut) md = 4;
	}
	else if (md == 3)
	{
	    /*fix and free parameters*/
	    while (1)
	    {
                /*reset f each time*/
                f = -1;
	    	printf("Type the number of the parameter to be fixed/freed: * = fixed\n");
                for (i = 0; i < NUMPARS; i++) printf(" %s   %1d) %s\n",
                    (ft.freepars[i] == 0) ? "  * " : "free",i+1,par_nm[i]);
                printf("        0) Finish\n");
		
    	    	get_ans(ans,1,0);
 	    	if (ans[0] >= '0' && ans[0] <= NUMPARS + '0') f = ans[0] - '0';
		
		if (f == 0) break;
		else if (f > 0 && f <= NUMPARS && ft.freepars[f-1] == 0)
		{
		    ft.freepars[f-1] = 1;
		    ft.nfp--;
		}
                /*if free then fix*/
		else if (f > 0 && f <= NUMPARS && ft.freepars[f-1] == 1)
		{
                    /*if there's a good fit, use those parameters*/
                    ft.pars[f-1] /= fact[f-1];
		    sprintf(ans, "Enter value to fix %s\n",par_nm[f-1]);
                    if (gdfit == 1) read_par(ft.pars, f-1, 1, ans);
    	    	    else 
                    {
                        read_par(init, f-1, 1, ans);
		        ft.pars[f-1] = init[f-1];
                    }
		    ft.freepars[f-1] = 0;
		    printf("Parameter %d fixed as = %f\n",f,ft.pars[f-1]);
		    ft.pars[f-1] *= fact[f-1];
		    ft.nfp++;
		}
	    }	    
	}
	else if (md == 4 && flginit)
	{
            if (gdfit == 0 || ask_yn(" Start fit with initial (rather than current) values of parameters? [y/n] (n)", 0))
            {
                /*if yes*/
                printf("***Using initial parameter values to start fit***\n\n");
                for (i = 0; i < NUMPARS; i++)
                {
                    if (ft.freepars[i] != 0) ft.pars[i] = init[i]*fact[i];
                }
            }
            else printf("***Using current values of parameter values to start fit***\n\n");
	    /*perform fit and write output file*/
	    if ( (r = fitter(MITER, &chisq, 1)) == 0 )
	    {
                /*write output file*/
	    	write_output(inname, outname, init, ft.x, &fit, chisq, 0);
                gdfit = 1;

                /*ask about overlaying fit on data in P_PROG
                    providing P_PROG exists and file is not an ORTEC (.Spe) file*/
                if (!aut) plot_data(ans, cmd, inname, outname, plt, 0);
	    }
	    else if (r > 0 && r < 3)
            {
		printf("%sBad fit, try changing initial parameters%s\n",
                        clr[1],clr[0]);
                gdfit = -1;
                write_output(inname, outname, init, ft.x, &fit, chisq, 0);
            }
	    else
            {
                printf("%sTry freeing some parameters%s\n",clr[1],clr[0]);
                gdfit = -1;
            }
            /*if auto mode fit background next*/
            if (aut) md = 6;
	}
	else if ( md == 4 && ! flginit)
	    printf("%s ***Enter initial parameter values before fitting!***%s\n",clr[3],clr[0]);
	else if (md == 5)
	{
	    printf("Data will now be compressed by factor of 2\n"
	    	" Initial guesses for coefficient will also be compressed\n");
	    compress2(inname, outname);
            for (i = 0; i < (NUMPARS - 1); i++) printf("###1init[%d]:%f \n",i,init[i]);
            
            /*change initial guesses for parameters and apply fact[] for use in eval()*/
	    for (i = 0; i < NUMPARS; i++)
            {
                if (ft.ad[i] != 1) init[i] /= 2.0;
                else init[i] *= 2.0;
                ft.pars[i] = init[i]*fact[i];
            }
	    init[NUMPARS] *= 2;
	    bkgnd *= 2;
	}
	else if (md == 6 && r == 0)
	{
            /*store current value of bkgnd for step-by-step comparison
                bklo, and for testing the whether to exit while, bklo2*/
	    bkcnt = 0;
            bklo = bkgnd;
	    bklo2 = bkgnd;
            flgbkpm = 1;
            chisqo = 10000.0;
 	    printf("  Bkgnd     t1/2    Chisq/D.O.F \n");
            while (1)
            {
	        if ( (r = fitter(MITER, &chisq, 0)) == 0 )
	    	{
		    printf("  %.2f    %.2f    %.5f\n",
		        bkgnd,(ft.pars[0]/fact[0]),chisq);
	    	    /*if lower chisq update bkgnd*/
                    if (chisq < chisqo)
                    {
                        bklo = bkgnd;
                        bkcnt = 0;
                    }
                    /*else increment 'no improvement in chisq' counter*/
                    else bkcnt++;
                    
                    /*if no change from original three times set minus flag*/
                    if ( (fabs(bklo2 - bklo) < 0.001) && bkcnt > 2 && flgbkpm == 1)
                    {
		        /*flgbkpm = -1 means subtract to search for background*/
                        flgbkpm = -1;
                        bkgnd = bklo;
                        bkcnt = 0;
                        /*printf("Sign flgbkpm=%d (blko2-bklo:%f - %f)\n",
                            flgbkpm,bklo2,bklo);*/
                    }
		    /*update chisq*/
                    chisqo = chisq;
		}
		else printf(" Bad fit!\n");
/*		chi2(&chisq);
		printf("For bkgnd = %f --> chisq/D.O.F = %.3f\n",bkgnd,chisq);*/
                /*exit loop after four times no change*/
                if (bkcnt > 3)
                {
                    bkgnd = bklo;
                    break;
                }
                /*update bkgnd value for testing*/
                if (flgbkpm == 1 || (flgbkpm == -1 && bkgnd - 0.1 >= 0.0))
                        bkgnd += (double)flgbkpm*0.1;
                /*if can't up-date bkgnd, stop*/
                else break;
                /*printf("bkcnt = %d, bkgnd = %f (%f), %d\n",
                    bkcnt,bkgnd,bklo,flgbkpm);*/
            }
 /*    	    read_par(init, 4, 1, "Enter value for background constant background level (not fitted)\n");*/
	    init[NUMPARS] = bkgnd;
	    /*write output file*/    
	    /*perform fit and write output file*/
	    if ( (r = fitter(MITER, &chisq, 1)) == 0 )
	    {
		/*write output file*/
	    	printf(" ==> Writing output for bkgnd = %f\n",bkgnd);
	    	write_output(inname, outname, init, ft.x, &fit, chisq, 0);
                /*ask about overlaying fit on data in P_PROG
                    providing P_PROG exists and file is not an ORTEC (.Spe) file*/
                plot_data(ans, cmd, inname, outname, plt, 0);
	    }
	    else printf("%s Bad fit! Could not write output file.%s\n",clr[1],clr[0]);
            /*turn off auto mode*/
            if (aut) aut = 0;
    	}
	else if ( md == 6 && r != 0)
        {
	    printf("%s ***Must get a good fit before exploring background level!***%s\n",clr[3],clr[0]);	    
            if (aut) aut = 0;
        }
        else if (md == 7)
        {
            chi2(&chisq);
            write_output(inname, outname, init, ft.x, &fit, chisq, 1);
        }
	else if (md == 8)
	{            
	    /*write output file with more data points if required*/
	    spc[0] = ft.x[1] - ft.x[0];
            sprintf(ans, "Enter desired spacing between adjacent points\n");
            read_par(spc, 0, 1, ans);

            /*in case of no fit calculate chisq*/
            chi2(&chisq);
	    /*write output file*/
            if (spc[0] == ft.x[1] - ft.x[0])
            {
                file_status(nofit, 80);
	        write_output(inname, nofit, init, ft.x, &fit, chisq, 2);
            }
            else
            {
                file_status(smooth, 80);            
	        write_output(inname, smooth, init, ft.x, &fit, chisq, 3);
            }
            /*ask about overlaying fit on data in P_PROG
                providing P_PROG exists and file is not an ORTEC (.Spe) file*/
            plot_data(ans, cmd, inname, outname, plt, 0);
	}
        /*write input to file as 2 or 3 column data*/
	else if (md == 9)
	{
	    while (1)
	    {
	    	printf("Enter number of columns required 2 or 3\n");
    	    	get_ans(ans,1,0);
 	    	if (ans[0] >= '2' && ans[0] <= '3')
		{
		    f = ans[0] - '0';
		    break;
		}
	    }
            /*prepare filenames and write data*/
            prep_write(inname, outname, nm, ans, numch, f);
	}
        /*plot spectrum and overlay fit if gdfit == 1*/
        else if (md == (NUMOPT - 3))
        {
            printf("md = %d (%d)\n",md,gdfit);
            if (gdfit) plot_data(ans, cmd, inname, outname, plt, 2);
            else plot_data(ans, cmd, inname, outname, plt, 3);
        }
        /*kill all or some plot windows*/
        else if (md == (NUMOPT - 2))
        {
            /*ask about closing all plots*/
            sprintf(ans,"%sClose all (%d) %s plot windows [y/n] (y)?%s",
                clr[3],pcnt,P_PROG,clr[0]);
            if (ask_yn(ans,1))
            {
                printf("%sKilling %s plots...%s",clr[2],P_PROG,clr[0]);
                while (pcnt > 0)
                {
                    sprintf(ans,"kill -9 %d",pid[pcnt-1]);
                    if(system(ans)) printf("%sCommand %s failed%s\n",clr[1],ans,clr[0]);
                    pcnt--;
                }
                printf("%sComplete.%s\n\n",clr[2],clr[0]);
            }
            /*above 15 too many plots, so only offer kill all windows option*/
            else if (pcnt < 15)
            {
                /*ask which plot windows to close*/
                printf("Type the %s window number (in order of opening) to close\n",
                    P_PROG);
	        while (pcnt > 0)
	        {
                    /*reset answer each time*/
	    	    f = -1;
                    for (i = 1; i <= pcnt; i++)
                    {
                        if (pcnt < 10) printf(" %1d) pid:  %5d\n",i,pid[i-1]);
                        else printf(" %c%1d) pid:  %5d\n",(i < 10) ? '0' : ((int)(i/10)+'0'),
                                (i < 10) ? i : i-10,pid[i-1]);
                    }
                    if (pcnt < 10) printf(" 0) Finish\n");
                    else printf(" 00) Finish\n");
		
                    if (pcnt < 10) get_ans(ans,1,0);
                    else get_ans(ans,2,0);
                    f = atoi(ans);
		
		    if (f == 0) break;
		    else if (f > 0 && f <= pcnt)
		    {
                        sprintf(ans,"kill -9 %d",pid[f-1]);
                        if(system(ans)) printf("%sCommand %s failed%s\n",clr[1],ans,clr[0]);
                        
                        /*now shuffle the pids down to fill empty array slot*/
                        for (i = f; i <= pcnt; i++) pid[i-1] = pid[i];
                        
                        pcnt--;
 		    }
                }
	    }	    
        }
        /*set or adjust spectrum limits*/
        else if (md == (NUMOPT - 1))
        {
            ft.ndp = find_strt_end_zeros(ft.ndp);
        }
        /*get spectrum properties and fill initial values if auto mode enabled
            then do fit directly*/
        else if (md == NUMOPT)
        {
            /*auto-fill initial values*/
            spec_prop(init);
            
            /*set md = 2 to fill ft.pars using init and to call scaler()*/
            md = 2;
            /*set auto mode flag so get_mode() won't ask for option*/
            aut = 1;
        }
        /*exit*/
	else if (md == 0)
        {
            /*ask about closing all plots*/
            sprintf(ans,"%sClose all (%d) %s plot windows [y/n] (y)?%s",
                clr[3],pcnt,P_PROG,clr[0]);
            if (pcnt > 0 && ask_yn(ans,1))
            {
                printf("%sKilling %s plots...%s",clr[2],P_PROG,clr[0]);
                /*start with most recent first*/
                while (pcnt > 0)
                {
                    sprintf(ans,"kill -9 %d",pid[pcnt-1]);
                    if(system(ans)) printf("%sCommand %s failed%s\n",
                        clr[1],ans,clr[0]);
                    pcnt--;
                }
                printf("%sComplete. Exiting%s\n\n",clr[2],clr[0]);
            }
            break;
        }
    }
    return 0;
} /*END main()*/

/*==========================================================================*/
/*eval: calculate fit using present values of the pars 	    	    	    */
/****************************************************************************/
int eval(double *pars, double x, double *fit, double *derivs, int ch, int mode)
{
/*   back-to-back Gaussians convoluted with an exponential decay
        plus a prompt Gaussian component.
        
    ch can be used for debuggging by only printing out for a particular
    channel range or for a particular value of cnt (=no. of function calls).

     Calculate the fit using present values of the pars
       x is the channel number variable
       pars[0] = the half-life
       pars[1] = the FWHM of the Gaussian component
       pars[2] = the prompt peak centroid position
       pars[3] = the scaling factor for the calculated curve
       pars[5] = the left-to-right-scaling factor*/
    
    static int cnt = 0;
    double  d = 0.001, npars[NUMPARS], ofit = 0.0;
    /*double  ft1 = 0.0, ft2 = 0.0;*/
    int     i, lp = 0, DNUM = 0, pd[NUMPARS];
    
    /*store original parameter values and initialise derivs to zero*/
    for (i = 0; i < NUMPARS; i++)
    {
	npars[i] = pars[i];
	derivs[i] = 0.0;
        pd[i] = 0;
        /*store the derivs that need to be calculated*/
        if (ft.freepars[i] && !ft.ad[i]) pd[DNUM++] = i;
    }
    
    while (lp <= DNUM)
    {
/*      full function*/
        /*gauss-exp convolution right-hand exp. if pars[0] > 0*/
        *fit = pars[3]*( gauss_exp(x - pars[2], pars[0], pars[1]));
        if (mode == 2) return 0;
	
    	/*derivs[3] can be calculated analytically*/
    	if (lp == 0 && mode >= 1) derivs[3] = *fit/pars[3];
        
        /*add the functions*/
        *fit += bkgnd;
	
        if (mode == (NUMFUNC+2))
        {
            *fit = bkgnd;
            return 0;
        }
        
        /* calculate derivs only for mode.ge.1 */
    	if (mode >= 1)
	{
            /*store the original value of fit for comparison later and
                increment the first free parameter with a non-analytical deriv*/
	    if (lp == 0)
            {
                ofit = *fit;               
                pars[pd[lp]] += d;
            }
	                
            if (lp > 0)
            {
                /*workout derivative from previous function evaluation*/
                derivs[pd[lp-1]] = (*fit-ofit)/d;
                /*reset previous parameter*/
                pars[pd[lp-1]] = npars[pd[lp-1]];
                 
                /*set next parameter if not at last parameter*/
                if (lp < DNUM) pars[pd[lp]] += d;
            }
            /*increment counter*/
    	    lp++;
	}
	else return 0;
    }
    *fit = ofit;
    cnt++;
    return 0;
} /*END eval()*/

/*==========================================================================*/
/* gauss: generate a Gaussian profile                                       */
/****************************************************************************/
double gauss(double x, double fwhm)
{
    double  sigma, y;
    
    sigma = fwhm/(2.0*sqrt(2*log(2)));
    y = exp( x*x/(-2.0*sigma*sigma) );
    y /= sigma*sqrt(2.0*pi);
    return y;
} /*END gauss()*/

/*==========================================================================*/
/* gauss_exp: generate a Gaussian convoluted with an exponential profile    */
/****************************************************************************/
double gauss_exp(double x, double thalf, double fwhm)
{
    double  sigma, sign, tau, y;
    
    sigma = fwhm/(2.0*sqrt(2*log(2)));
    tau = thalf/log(2);
    sign = (tau/fabs(tau));    
    /*to avoid infinities before being multiplied by 0 later
        [0 = 1+erf() where erf = -1] check exponent, as exp(40) > 2e17*/
    y = (( x - ((0.5*sigma*sigma)/(tau)) )/(-1.0*tau));
    if ( y > 40) y = 0.0;
    else y = 0.5*exp(y);
    
    y *= (double)1.0 + erf(( x - (sigma*sigma/(tau)) )/(sign*sqrt(2)*sigma));
    /*for unit area of exponential divide by tau and multiply by point spacing
        (note x[1]-x[0] is an approximation for a quadratic calibration
        once a good fit is found this factor can be found for each function
        by taking the scaling factors and dividing by the area using the actual
        point spacing over the full extent of the function. Note this will
        lead to an if statement based on tau being required. Also, for long
        lifetime functions the data points will need to be extrapolated beyond
        channel 2047 using the calibration factors.*/
    y *= (ft.x[1] - ft.x[0])/fabs(tau);
    return y;
} /*END gauss_exp()*/

/*==========================================================================*/
/*ascii_read: read ASCII format data	    	    	    	    	    */
/****************************************************************************/
int ascii_read(char * fname, int *col)
{
    float   rd = 0.0;
    int     chan = 0, i, cnt = 0, hash = 0, res = 0, lchan, mxchan = MAXPTS;
    FILE    *file1;
    
    /*open ascii file*/
    if ((file1= fopen(fname, "r" )) == NULL)
    {
    	printf("Cannot open file: %s \n", fname);
    	return -1;	    	    	    
    }
    
    /*zero spectrum array*/
    for (i = 0; i < MAXPTS; i++)
    {        
    	ft.x[i] = 0.0;
    	ft.y[i] = 0.0;
    	ft.dy[i] = 0.0;
    }
   
    /*check for '.Spe' extension and if ORTEC ascii header is present.
        if so, skip header*/
    ortec_head_srch_skip(fname, file1, &mxchan);
    /*get number of columns in ascii file*/
    col_determ(file1, col);
    if (*col == 0)
    {
	printf("%s***No suitable data in file %s; Exiting...%s\n\n",
                clr[1],fname,clr[0]);
	return -1;
    }
    printf("Ascii %d Column format....\n", *col);
    
    lchan = mxchan;
    for (chan = 0; chan < mxchan; chan++)
    {
        /*skip hashed lines.Needed here if e.g. hashed trailer is present*/
        skip_hash(file1);
	/*only for 1 column data set the x value equal to chan*/
	if (*col == 1)
	{
	    ft.x[chan] = (double)(chan);
	    if (chan == 0)
	    	printf("\t**** Spectrum will start from x = 0 ****\n");
	    res = fscanf(file1, "%f", &rd);	
            ft.y[chan] = (double)rd;
	    ft.dy[chan] = sqrt(fabs(ft.y[chan]));
	}
	else if (*col == 2)
	{
	    res = fscanf(file1, "%f", &rd);
	    ft.x[chan] = (double)rd;
	    res = fscanf(file1, "%f", &rd);
	    ft.y[chan] = (double)rd;
	    ft.dy[chan] = sqrt(fabs(ft.y[chan]));
	}
	else if (*col == 3)
	{
	    res = fscanf(file1, "%f", &rd);
	    ft.x[chan] = (double)rd;
	    res = fscanf(file1, "%f", &rd);
	    ft.y[chan] = (double)rd;
	    res = fscanf(file1, "%f", &rd);
	    ft.dy[chan] = (double)rd;
	}
        if (cnt == 0)
        {
            x12[0] = ft.x[0];
            lim[0] = chan;
            printf("    Data start at data ch#%d (x=%.2f,y=%.2f)\n",
                        chan,ft.x[0],ft.y[0]);
        }
        cnt++;
        switch (res)
	{
	    case 0:
	    {
		if ( (hash = fgetc(file1)) == '#')
		{
    	    	    /*read rest of line*/
    	    	    while ( ( (hash = fgetc(file1)) != '\n' ) && (hash != EOF) )
    	    	    	;
		    
		    chan--;
		    break;
		}
		else if (chan > 1 && hash != '#')
		{		    
		    printf("\nline= %d, ft.y[line]= %f\n", chan,ft.y[chan]);
    	    	    printf("Read error occurred for file: %s,"
			    " but %d data points stored\n", fname,chan);
		    lchan = chan;
	    	    chan = MAXPTS;
		    break;
		}
		else
		{
	    	    printf("\nchannel= %d, ft.y[channel]= %f\n", chan,
		    	  ft.y[chan]);
    	    	    printf("%sRead error occurred for file: %s%s\n",
                        clr[1],fname,clr[0]);
		    /*close ascii file*/
	    	    fclose(file1);
	    	    return -1;
                }
	    }
	    case EOF:
	    {
	    	if (chan < MAXPTS)
		{
/*            	    printf("\nMax. no. of channels expected = %d\n", MAXPTS);*/
	    	    printf("    ...reached EOF...\n"); 
		    if (chan == 0)
		    {
			printf("\n*******Incorrect file format*******\n");
			printf("....Exiting....\n\n");
			return -1;
		    }
/*            	    printf("-->.spe file will still be written......\n"); */
		}
		lchan = chan;
	    	chan = MAXPTS;
    	    	break;
	    }
	    default:
	    {
/*	    	if (chan == 0) printf("Reading ascii spectrum.......\n"); */
	    	break;
    	    }
    	}
    }
    chan = lchan;
    x12[1] = ft.x[chan-1];
    lim[1] = chan-1;
    printf("    Data end at data ch#%d (x=%.2f,y=%.2f)\n",
                        chan-1,ft.x[chan-1],ft.y[chan-1]);

    /*close ascii file*/
    fclose(file1);
    /*return chan as this is number of values read, rather than the
    difference between lim[1] and lim[0] which might be bigger if, 
    e.g. only every 2nd channel was written to file*/
    return chan;
} /*END ascii_read()*/

/*==========================================================================*/
/* ascii_write: write an ASCII format spectrum	    	    	    	    */
/****************************************************************************/
void ascii_write(char name[], int numch, int col)
{
    int i;
    FILE *fasc;
    /*open .txt file*/
    if ( (fasc = fopen(name, "w" )) == NULL)
    {
        printf("Cannot open file: %s \n", name);
        return ; 
    }
   
    /*write .txt file*/
    for (i = 0; i < numch; i++)
    {
    	fprintf(fasc, "%8.1f\t%8.3f", ft.x[i], ft.y[i]);
	if (col == 3) fprintf(fasc, "\t%8.3f\n",ft.dy[i]);
	else fprintf(fasc, "\n");
    }
   
    printf(" written ==> %s, %d chs %d columns.\n", name, numch, col);
    fclose(fasc);
} /*END ascii_write()*/

/*==========================================================================*/
/* ask_yn: 0 means no, 1 means yes. yn == -1 for no default, else yn default*/
/****************************************************************************/
int ask_yn(char *questn, int yn)
{
    char ans[3] = "";
    
    while (1)
    {
    	printf("%s\n", questn);
        /*0 option to get_ans doesn't allow carriage returns*/
        if (yn == -1) get_ans(ans,1,0);
        else get_ans(ans,1,1);
        if (yn != -1 && (ans[0] == '\n')) return yn;
        else if (ans[0] == 'y' || ans[0] == 'Y') return 1;
	else if (ans[0] == 'n' || ans[0] == 'N') return 0;
    }
} /*END ask_yn()*/

/*==========================================================================*/
/* chi2: calculate chisquared	    	    	    	    	    	    */
/****************************************************************************/
void chi2(double *chisq)
{
    static int	npars = NUMPARS;
    extern double fit, derivs[NUMPARS];
    double  diff = 0.0, dat;
    int     i, nip, ndf;
    
    *chisq = 0.0;
   /*no. independent parameters*/
    nip = npars - ft.nfp;
    /*no. degrees of freedom*/
    ndf = ft.ndp - nip;
    
    for (i = (0+lim[0]); i < (ft.ndp+lim[0]); i++)
    {
    	eval(ft.pars, ft.x[i], &fit, derivs, i, 0);
    	diff = ft.y[i] - fit;
    	dat = ft.dy[i] * ft.dy[i];
    	if (dat == 0.f) dat = 1.f;
    	
    	*chisq += diff * diff / dat;
    }
    *chisq /= (double)ndf;
    printf(" %d indept. pars    %d degrees of freedom\n", nip, ndf);
    printf("Chisq/D.O.F. = %.4f (Chisq = %.4f)\n",*chisq,*chisq*(double)ndf);
    printf("\t(c.f. full spectrum with no fixed pars %.3f)\n",
            (*chisq*(double)ndf/(double)(ndf-ft.nfp)));  
} /*END chi2()*/

/*==========================================================================*/
/* col_determ: find if a file has 1, 2, 3 or 4 column format      	    */
/****************************************************************************/
void col_determ(FILE *file1, int *col)
{
    int     i = 0, hash = 0;
    fpos_t  pos;
    
    /*skip comments lines starting with #*/
    skip_hash(file1);

    /*Try to decide if spectrum is 1, 2 or 3 column ascii format*/
    *col = 0;
    i = 0;
    
    /*store current file position in pos*/
    fgetpos(file1, &pos);
    
    while (1)
    {
    	/*read until first digit of first number*/
	while ( isdigit(hash = fgetc(file1)) == 0 ) ;
	
	*col = 1;
    	/*found first digit of first number. read rest of digits*/
    	/*ignore decimal points*/
	while ( (isdigit(hash = fgetc(file1)) != 0) || (hash == '.') ) ;

    	/*push the last character back on the stream*/
    	ungetc(hash, file1);
    	/*now test characters until end of line*/
    	while ( (hash = fgetc(file1)) != '\n' )
    	{
    	    /*if blank space,increment blank counter and continue*/
    	    if ( (hash == ' ') || (hash == '\t') || (hash == '\r')) i++;
    	    /*check how many columns of numbers are present*/
            else if ( (isdigit(hash)) == 0 )
            {
                printf("%s***Input is not a valid data file. Illegal characters found.%s\n",
                        clr[1],clr[0]);
                *col = 0;
                return ;
            }
    	    /*if character is a digit*/
    	    /*check how many columns of numbers are present*/
    	    else if ( (isdigit(hash)) != 0 )
    	    {
    		/*if number increment col flag*/
    		if ( i > 0)
		{
		    /*brackets necessary round *col to increment what it
		    	points to not the address!*/
		    (*col)++;
		    /*printf("col. = %d\n", *col); */
		}
		
    		/*otherwise if > MAXCOLS columns set 1 col. flag and hope for best*/
    		if ( (i != 0) && (*col > MAXCOLS) )
    		{
    		    *col = 1;
    		    break;
    		}
    		/*reset blank counter to zero*/
    		i = 0;
    		/*read rest of digits*/
		while ( (isdigit(hash = fgetc(file1)) != 0) || (hash == '.') ) ;
		
    		/*push the last character back on the stream*/
    		ungetc(hash, file1);
    	    }
    	}
    	break;
    }	    
    /*reset file position to start*/
    fsetpos(file1, &pos);   
} /*END col_determ()*/

/*==========================================================================*/
/*compress2: compress data by factor of 2	     	    	    	    */
/****************************************************************************/
void compress2(char inname[], char outname[])
{
    static int	fact = 2;
    int    i;
    char   ext[10] = "", next[10] = "-", nname[CHLEN] = "";
     
    /*do compression*/
    for (i = 0; i < (int)(ft.ndp/2); i++)
    {
	ft.x[i] = (int)(ft.x[(2*i)]/2);
	ft.y[i] = ft.y[(2*i)] + ft.y[(2*i)+1];
	ft.dy[i] = sqrt( (ft.dy[(2*i)]*ft.dy[(2*i)]) +
		(ft.dy[(2*i)+1]*ft.dy[(2*i)+1]) );
    }
    printf("\t ....done. %d channels\n",i);

    ft.ndp = i;
    /*get output filename extension*/
    strcpy(nname, inname);
    if ( (strrchr(nname,'.')) ) strcpy( ext, (strrchr(nname,'.')) );
    /*if ext is Chn or Spe change to .txt*/
    if ( (!strcmp( (strrchr(inname,'.')), ".Spe" )) ||
                    ( (!strcmp( (strrchr(inname,'.')), ".Chn" ) ) ) )
        strcpy( ext, ".txt");
    
    /*now add the compression factor to the input and output filenames*/
    itoa(fact, next+1);
    strcat(next, ext);
    set_ext(nname, next);
    strcpy(outname, nname);
    set_ext(outname, ".fit");
    if (ask_yn(" Write compressed data to file? [y/n] (y)",1))
    {
   	printf("Output filename for compressed data: %s\n", nname);
    	file_status(nname, CHLEN);
        ascii_write(nname, i, 2);
    	strcpy(inname, nname);	
    }
    fact *= 2;
}/*END compress2()*/

/*==========================================================================*/
/* cswap2: swap bits of a 2 byte number     	    	    	    	    */
/****************************************************************************/
int cswap2(int decim)
{
    int i, j, max = 17;
    char bin[17];
    int swapped = 0;
    
    for (i = 1; i >= 0; --i)
    {
	for (j = 0; j < (max-1)/2; ++j)
	{
    	    if (decim & 0x8000) bin[i*8+j] = '1';
	    else bin[i*8+j] = '0';
	    
	    decim <<= 1;
/*	    printf(" i*8+j = %d \n", i*8+j); */
	}
    }
    bin[max-1] = '\0';
/*    printf("bin = %s \n", bin); */
    for (i = 0; i < max-1; ++i)
    {
	if (bin[max-2-i] == '1')
	{
	    swapped = swapped + (int)pow( (double)2, (double)(i));
	}
    }
    
    return swapped;
} /*END cswap2()*/

/*==========================================================================*/
/* cswap4: swap bits of a 4 byte number	    	    	    	    	    */
/****************************************************************************/
int cswap4(int decim)
{
    int i, j, max = 33;
    char bin[33];
    int swapped = 0;
    
    for (i = 3; i >= 0; --i)
    {
	for (j = 0; j < (max-1)/4; ++j)
	{
    	    if (decim & 0x80000000) bin[i*8+j] = '1';
	    else bin[i*8+j] = '0';
	    
	    decim <<= 1;
/*	    printf(" i*8+j = %d \n", i*8+j); */
	}
    }
    bin[max-1] = '\0';
/*    printf("bin = %s \n", bin); */
    for (i = 0; i < max-1; ++i)
    {
	if (bin[max-2-i] == '1')
	    swapped = swapped + (int)pow( (double)2, (double)(i));
    }
    
    return swapped;
} /*END cswap()*/

/*==========================================================================*/
/* file_status: check file status   	    	    	    	    	    */
/****************************************************************************/
void file_status(char *fname, int len)
{
    char ans[3] = "";
    struct stat statbuf;      /*need include files sys/stat.h and sys/types.h*/
    
    while ( (stat(fname, &statbuf) == 0) )
    {
    	printf("\n*****Output file %s exists. Overwrite [y/n]? (y)\n", fname);
	get_ans(ans,1,1);
	if ((ans[0] == '\n') || (ans[0] == 'y' || ans[0] == 'Y')) break;
	else if (ans[0] == 'n' || ans[0] == 'N')
	{
	    printf("Enter new file name inc. extension (eg %s):\n",
		    strrchr(fname,'.'));
	    get_line(fname, len);
    	}
    }
} /*END file_status()*/

/*==========================================================================*/
/*find_strt_end_zeros: check y values for zeros at spectrum start and end   */
/****************************************************************************/
int find_strt_end_zeros(int numch)
{
    int     i = 0, flgstrt = 0, flgend = 0;
    char    ans[3*CHLEN] = "";
    
    /*check start y value of data for zeros*/
    while (i < numch && ft.y[i] == 0.0f) i++;
    
    if (i < numch && i > 0)
    {
        printf("\n\tData start at data ch#%d (x=%.2f,y=%.2f)\n",
                        i,ft.x[i],ft.y[i]);
        /*here, if set, preserve current limit for display*/
        if (lim[0] == 0) lim[0] = i;
        flgstrt = 1;
    }
    /*else leave limit as is - set in spectrum read function*/
    /*else lim[0] = 0;*/
    
    /*check end y value of data for zeros*/
    i = numch - 1;
    while (i > 0 && ft.y[i] == 0.0f) i--;

    if (i > lim[0] && i < (numch -1) )
    {
        if(!flgstrt) printf("\n");
        printf("\tData end at data ch#%d (x=%.2f,y=%.2f)\n",
                        i,ft.x[i],ft.y[i]);
        /*here, if set, preserve current limit for display*/
        if (lim[1] == (numch - 1)) lim[1] = i;
        flgend = 1;
    }
    /*else leave limit as is - set in spectrum read function*/
    /*else lim[1] = numch - 1;*/
    
    while (1)
    {
        if (flgstrt || flgend)
        {
            printf("\n%sZeros found at the %s%s%s of the spectrum.%s\n",
                clr[3],flgstrt ? "start":"",
                (flgstrt && flgend) ? " and ":"",flgend ? "end":"",clr[0]);
        }
        else printf("\n%sNo zeros found at the start or end of the spectrum.%s\n",
                clr[3],clr[0]);
            
        sprintf(ans,"%sOmit channels from fit or set own limits? [y/n] (y)%s",
                clr[3],clr[0]);
        
        if (ask_yn(ans,1))
        {
            sprintf(ans,"%sAdjust these limits [y/n] (n)?%s",clr[0],clr[0]);
            printf("%sCurrent number of data channels for fitting: %d (%d --> %d)%s\n"
                "%sch#%d (x=%.2f,y=%.2f) --> ch#%d (x=%.2f,y=%.2f)%s\n",
                clr[2],(lim[1]-lim[0]+1),lim[0],lim[1],clr[0],
                clr[2],lim[0],ft.x[lim[0]],ft.y[lim[0]],
                lim[1],ft.x[lim[1]],ft.y[lim[1]],clr[0]);
            
            while (ask_yn(ans,0))
            {
                /*set x12[] to lim[] for printing*/
                x12[0] =lim[0];
                x12[1] =lim[1];
                read_par(x12, 0, 2, "Enter lower and upper limits for data chan#");
                lim[0] = (int)(x12[0]+0.5);
                lim[1] = (int)(x12[1]+0.5);
                x12[0] = ft.x[lim[0]];
                x12[1] = ft.x[lim[1]];
                printf("%sNumber of data channels for fitting: %d (%d --> %d)%s\n"
                    "%sch#%d (x=%.2f,y=%.2f) --> ch#%d (x=%.2f,y=%.2f)%s\n",
                    clr[2],(lim[1]-lim[0]+1),lim[0],lim[1],clr[0],
                    clr[2],lim[0],ft.x[lim[0]],ft.y[lim[0]],
                    lim[1],ft.x[lim[1]],ft.y[lim[1]],clr[0]);
            }
            return lim[1] - lim[0] + 1;
    	    
            /*now remove zeros*/
/*     	    for (i = 0; (i < numch) && ((i + j) < numch); i++)
    	    {
	        while ( (ft.y[i+j] == 0.0f) && ((i + j) < numch) ) j++;
	
	            ft.x[i] = ft.x[i+j];
	            ft.y[i] = ft.y[i+j];
	            ft.dy[i] = ft.dy[i+j];
    	    }   
    	    if ( ((i + j) >= numch) && (ft.y[i+j] == 0.0f) ) i--;
	
	    return i;*/
        }
        else
        {
            lim[0] = 0.0;
            lim[1] = numch;
            printf("%sNumber of data channels for fitting: %d%s\n",
                clr[2],numch,clr[0]);
            return numch;
        }
    }
} /*END find_strt_end_zeros()*/

/*==========================================================================*/
/*fitter: fit parameters. vb == 0 suppresses output info except warnings    */
/****************************************************************************/
int fitter(int maxits,  double *chisq, int vb)
{
    static int npars = NUMPARS;

    extern double fit, derivs[NUMPARS];
    double ddat, alpha[NUMPARS][NUMPARS], array[NUMPARS][NUMPARS];
    double  r1, diff, beta[NUMPARS], b[NUMPARS], delta[NUMPARS];
    double  chisq1, flamda, dat, ers[NUMPARS];
    int    conv, nits, test, i, j, k, l, m, nextp[NUMPARS], ndf, nip;

    /* this subroutine is a modified version of 'CURFIT', in Bevington */
    /* see page 237. From RadWare by D.C. Radford at ORNL*/
  
    *chisq = 0.0;
    
    for (i = 0; i < npars; i++)
    {
    	derivs[i] = 0.0;
    	for (j = 0; j < npars; j++) array[i][j] = (double)0.0;
    }
    fit = 0.0;
  
    /*no. independent parameters*/
    nip = npars - ft.nfp;
    /*no. degrees of freedom*/
    ndf = ft.ndp - nip;
    if (ndf < 1)
    {
    	printf("No Degrees-of-Freedom.\n");
    	return 3;
    }
    if (nip < 2)
    {
    	printf("Too many fixed parameters.\n");
    	return 4;
    }
    k = 0;
    for (j = 0; j < npars; ++j)
    {
    	if (ft.freepars[j]) nextp[k++] = j;
    }
    if (k != nip)
    {
    	printf("No. of independent pars != sum(freepars)!!!\n");
    	return 5;
    }
    flamda = 0.001f;
    nits = 0;
    test = 0;
    derivs[0] = 1.0f;
    for (i = 0; i < npars; ++i)
    {
    	ft.errs[i] = 0.0f;
    	b[i] = ft.pars[i];
    }
    /* evaluate fit, alpha & beta matrices, & chisq */
    NEXT_ITERATION:
    for (j = 0; j < nip; ++j)
    {
    	beta[j] = (double)0.0;
    	for (k = 0; k <= j; ++k) alpha[k][j] = 0.0;
    }
    chisq1 = 0.0f;
    for (i = (0+lim[0]); i < (ft.ndp+lim[0]); ++i)
    {
    	eval(ft.pars, ft.x[i], &fit, derivs, i, 1);
    	diff = ft.y[i] - fit;
    	dat = ft.dy[i] * ft.dy[i];
    	if (dat < 1.0f) dat = 1.0f;
	
    	ddat = (double)dat;
    	chisq1 += diff * diff / dat;
    	for (l = 0; l < nip; ++l)
	{
      	    j = nextp[l];
      	    beta[l] += diff * derivs[j] / dat;
      	    for (m = 0; m <= l; ++m)
	    {
	    	alpha[m][l] += (double) derivs[j] * (double) derivs[nextp[m]] / ddat;
      	    }
    	}
    }
    chisq1 /= (double)ndf;
    
    /* invert modified curvature matrix to find new parameters */
    INVERT_MATRIX:
    array[0][0] = flamda + 1.0f;
    for (j = 1; j < nip; ++j)
    {
    	for (k = 0; k < j; ++k)
	{
      	    if (alpha[j][j] * alpha[k][k] == 0.0)
	    {
	    	printf("Cannot - diag. element %i or %i eq. to 0.0\n", j, k);
	    	return 1;
      	    }
      	    array[k][j] = alpha[k][j] / sqrt(alpha[j][j] * alpha[k][k]);
      	    array[j][k] = array[k][j];
    	}
    	array[j][j] = flamda + 1.0f;
    }
    matinv(array[0], nip, npars);
    if (!test)
    {
    	for (j = 0; j < nip; ++j)
	{
      	    if (alpha[j][j] * alpha[j][j] == 0.0)
	    {
	    	printf("Cannot - diag. element %i eq. to 0.0\n", j);
	    	return 1;
      	    }
      	    delta[j] = 0.f;
      	    for (k = 0; k < nip; ++k)
	    {
	    	delta[j] += beta[k] * array[k][j] / sqrt(alpha[j][j] * alpha[k][k]);
      	    }
    	}
    	/* if chisq increased, increase flamda & try again */
    	*chisq = 0.0f;
    	for (l = 0; l < nip; ++l)
	{
      	    j = nextp[l];
      	    b[j] = ft.pars[j] + delta[l];
    	}
    	for (i = (0+lim[0]); i < (ft.ndp+lim[0]); ++i)
	{
      	    eval(b, ft.x[i], &fit, derivs, i, 0);
      	    diff = ft.y[i] - fit;
      	    dat = ft.dy[i] * ft.dy[i];
      	    if (dat == 0.f) dat = 1.f;
      	    
	    *chisq += diff * diff / dat;
    	}
    	*chisq /= (double)ndf;
    	if (*chisq > chisq1 && flamda < 2.f)
	{
      	    flamda *= 10.f;
      	    goto INVERT_MATRIX;
    	}
    }
    
    /* evaluate parameters and errors */
    /* test for convergence */
    conv = 1;
    for (j = 0; j < nip; ++j)
    {
    	if (array[j][j] < 0.0) array[j][j] = 0.0;
	
    	ers[j] = sqrt(array[j][j] /  alpha[j][j]) * sqrt(flamda + 1.0f);
    	if ((r1 = delta[j], fabs(r1)) >= ers[j] / 1e3f) conv = 0;	
    }
    if (!test)
    {
    	for (j = 0; j < npars; ++j)
	{
      	    ft.pars[j] = b[j];
    	}
    	flamda /= 10.f;
    	++nits;
    	if (! conv && nits < maxits)
	{
	    goto NEXT_ITERATION;
	}

    	/* re-do matrix inversion with flamda=0 to calculate errors */
    	flamda = 0.f;
    	test = 1;
    	goto INVERT_MATRIX;
    }

    /* list data and exit */
    for (l = 0; l < nip; ++l) ft.errs[nextp[l]] = ers[l];
    
    if (vb) printf(" %i indept. pars    %i degrees of freedom\n", nip, ndf);
    if (conv)
    {
    	if (vb) printf("%s%i iteration(s),  Chisq/D.O.F. = %.4f (Chisq = %.4f)%s\n",
            clr[2],nits,*chisq,*chisq*ndf,clr[0]);
     	return 0;
    }
    printf(" Failed to converge after %i iteration(s),  %sChisq/D.O.F. = %.3f%s\n"
    	"  WARNING - do not believe quoted errors.\n",nits,clr[3],*chisq,clr[0]);
    return 2;
} /*END fitter()*/	

/*==========================================================================*/
/* get_ans: get answer which can be carriage return if md = 1. 	    	    */
/****************************************************************************/
void get_ans(char ans[], int num, int md)
{
    int     i;
    struct  termios newt, oldt;
    
    /*This code is a modified version of a function in the RadWare
        software suite, by D.C. Radford*/
    while (1)
    {
    	tcgetattr(0, &oldt);
	newt = oldt;
    	newt.c_lflag &= ~ICANON;
    	newt.c_cc[VMIN] = 1;
    	newt.c_cc[VTIME] = 0;
	/*handle sigs*/
/*    	newt.c_lflag |= ISIG;*/
    	tcsetattr(0, TCSANOW, &newt);
    	i = 0;
    	while( (ans[i++] = (char)getchar()) != '\n' && i < num) ;
    	
	tcsetattr(0, TCSANOW, &oldt);
    	if (ans[i-1] != '\n') printf("\n");
	else if (md == 0 && ans[0] == '\n') continue;
	
    	ans[i] = '\0';
    	return ;
    }
} /*END get_ans()*/

/*==========================================================================*/
/* get_line: read a line into ans, of max length len   	    	    	    */
/****************************************************************************/
void get_line(char ans[], int len)
{
    int     c, i = 0;

    memset(ans,'\0',sizeof(char)*len);
    /*note that scanf() leaves a carriage return in the keyboard buffer*/
    while ( (c = getchar()) != '\n' && c != EOF && i < len) ans[i++] = c;
    
    ans[i] = '\0';
} /*END get_line()*/

/*==========================================================================*/
/* get_line_file: read a line from a file into s, return a length           */
/****************************************************************************/
void get_line_file(FILE *file, char ans[], int len)
{
    int     c, i = 0;

    memset(ans,'\0',sizeof(char)*len);
    /*note that scanf() leaves a carriage return in the keyboard buffer*/
    while ( (c = fgetc(file)) != '\n' && c != EOF && i < len) ans[i++] = c;
    
    /*remove any trailing '\r' carriage returns that are often prior to
      '\n' line feeds*/
    if (ans[i-1] == '\r') i--;
    
    ans[i] = '\0';
} /*END get_line_file()*/

/*==========================================================================*/
/* get_mode: get mode from user input	    	    	    	    	    */
/****************************************************************************/
int get_mode(int md, int aut)
{
    char    ans[10] = "";
    
    /*if auto flag is not set, print options. If it is set,
        get_mode() just returns the value of md that was passed*/
    while(1 && !aut)
    {
    	printf("\n 1) Enter filename and read data\n");
    	printf(" 2) Enter initial parameters\n");
    	printf(" 3) Fix or free parameters\n");
    	printf(" 4) Perform fit and write output to file\n");
    	printf(" 5) Compress data by factor of 2\n");
    	printf(" 6) Explore background value as function of chi^2\n");
    	printf(" 7) Print current values of parameters to screen\n");
    	printf(" 8) Write output with current values of parameters\n"
                "\t\t\t    or different point spacing\n");
    	printf(" 9) Write input to file as 2 (x y) or 3 (x y dy) column data\n");
    	if (AUTO) printf(" a) Auto-fill initial parameters and perform fit\n");
    	printf(" p) Display spectrum %s\n",gdfit ? "and fit" : "");
    	if (pcnt) printf(" k) Kill all or selected %s plot windows\n",P_PROG);
    	printf(" s) Set or adjust spectrum limits for fitting\n");
    	printf(" 0) Quit\n");
	
	get_ans(ans,1,0);
	
 	if (ans[0] >= '0' && ans[0] <= '9')
	{
	    md = ans[0] - '0';
	    break;
	}
	else if (ans[0] == 'p' || ans[0] == 'P')
	{
	    md = NUMOPT-3;
	    break;
	}
	else if (pcnt && (ans[0] == 'k' || ans[0] == 'K'))
	{
	    md = NUMOPT-2;
	    break;
	}
	else if (ans[0] == 's' || ans[0] == 'S')
	{
	    md = NUMOPT-1;
	    break;
	}
	else if (ans[0] == 'a' || ans[0] == 'A')
	{
	    md = NUMOPT;
	    break;
	}
    }
    return md;
} /*END get_mode()*/

/*==========================================================================*/
/* get_pars: extract comma or space separated numbers from string ans0	    */
/****************************************************************************/
void get_pars(char ans0[], double pars[], int num)
{
    int     i, j = 0, k, minus, p = 0;
    char    ans1[CHLEN] = "";

/*    get_line(ans0);*/
    
    for (i = 0; i < num; i++)
    {
	k = 0;
	minus = 0;
	memset(&ans1,'\0',sizeof(ans1[0])*CHLEN);
	while ( ! isdigit(ans0[j]) )
	{
	    if (ans0[j] == '-') minus = 1;
	    j++;
	}
	if (minus == 1)
	{
	    j -= 1;
	    ans0[j] = '-';
	}	
	while( ans0[j] != '\n' && ans0[j] != ' ' && ans0[j] != ','
		&& j < strlen(ans0) ) ans1[k++] = ans0[j++];
	
	if (k > 0) p++;
	j++;
	pars[i] = (double)atof(ans1);
    }
    printf("\n");
} /*END get_pars()*/

/*==========================================================================*/
/* get_pars_file: extract separated numbers from string ans0 read from file */
/****************************************************************************/
void get_pars_file(FILE *file, float pars[], int num)
{
    int     i, j = 0, k, minus;
    char    ans0[CHLEN] = "", ans1[40] = "";

    get_line_file(file, ans0, CHLEN);
    
    for (i = 0; i < num; i++)
    {
	k = 0;
	minus = 0;
    	memset(ans1,'\0',sizeof(ans1));
	while ( ! isdigit(ans0[j]) )
	{
	    if (ans0[j] == '-') minus = 1;
	    j++;
	}
	if (minus == 1)
	{
	    j -= 1;
	    ans0[j] = '-';
	}	
	while( ans0[j] != '\n' && ans0[j] != ' ' && ans0[j] != ','
		&& j < strlen(ans0) ) ans1[k++] = ans0[j++];
	
	j++;
	pars[i] = atof(ans1);
    }
} /*END get_pars_file()*/

/*===========================================================================*/
/* convert integer n to string */
/*****************************************************************************/
void itoa(int n, char s[])
{
    int i = 0, sign;
    
    /*record sign and make n positive if needed*/
    if ((sign = n) < 0) n = -n;
	
    /* generates digits in reverse order */
    do {
        /*get next digit then delete it by dividing by 10*/
    	s[i++] = n % 10 + '0';
    } while ((n /= 10) > 0);
    if (sign < 0) s[i++] = '-';
    
    s[i] = '\0';
    reverse(s);
} /*END itoa*/

/*===========================================================================*/
/* maestro_read: read the maestro format spectrum                            */
/*****************************************************************************/
int maestro_read(char fname[], int *col)
{
    int *counts, i = 0;
    char dt[10] = "";
    FILE *file1;
        
    /*opens read only Maestro file*/
    if ( (file1 = fopen(fname, "r")) == NULL)
    {
    	printf("Cannot open file: %s\n", fname);
	return -1;
    }
    
    /*zero spectrum array*/
    for (i = 0; i < MAXCOLS; i++)
    {        
    	ft.x[i] = 0.0;
    	ft.y[i] = 0.0;
    	ft.dy[i] = 0.0;
    }
    /*default is one column for Maestro Chn binary format*/
    *col = 1;
    
    /*clear the header*/
    memset(&maest_header, 0, sizeof(maest_header));
    /*construct Maestro header*/
    fread(&maest_header, sizeof(maest_header), 1, file1);
        
/*    printf(" maest_header.q1 = %d\n maest_header.q2 = %d \n"
	   " maest_header.q3 = %d\n maest_header.q4 = %d \n"
	   " maest_header.real = %d\n maest_header.lve = %d \n"
	   " maest_header.dt = %s\n maest_header.sttm = %s \n"
	   " maest_header.off = %d\n maest_header.channels = %d\n",
	    maest_header.q1, maest_header.q2, maest_header.q3,
	    maest_header.q4, maest_header.real, maest_header.lve,
            maest_header.dt, maest_header.sttm, maest_header.off,
            maest_header.channels);*/
    
    /*unix byte swapping option*/
    if (maest_header.channels > CHMAX && cswap4(maest_header.channels) > CHMAX )
    {
    	printf("Unrecognised format. Exiting.....\n");
	fclose(file1);
    	return -1;
    }
    
    if (maest_header.channels > CHMAX || 
            maest_header.channels < (i = cswap2(CHMAX)) )
    {
	/*swap the bytes in the headers*/
    	maest_header.channels = cswap2(maest_header.channels);
	maest_header.real = cswap4(maest_header.real);
	maest_header.lve = cswap4(maest_header.lve);
        
	printf(".......SWAPPING BYTES read from file.......\n");
    	/*allocate sufficient memory for spectrum*/
    	counts = (int *) malloc( maest_header.channels*sizeof(int) );
    	/*read the data*/
	fread(counts, maest_header.channels*sizeof(int), 1, file1);
	/*read the trailer*/
	fread(&maest_trailer, sizeof(maest_trailer), 1, file1);
    	maest_trailer.g[0] = cswap4(maest_trailer.g[0]);
    	maest_trailer.g[1] = cswap4(maest_trailer.g[1]);
    	maest_trailer.g[2] = cswap4(maest_trailer.g[2]);

      	/*fill spectrum array*/
    	for (i = 0; i < (int) maest_header.channels; i++)
	    swapb4( (char *) (counts + i) );	
    } /*end of byte swapping loop for unix*/    
    else
    {
    	/*allocate sufficient memory for spectrum*/
    	counts = (int *) malloc( maest_header.channels*sizeof(int) );
	/*read the data*/
	fread(counts, maest_header.channels*sizeof(int), 1, file1);
	/*read the trailer*/
	fread(&maest_trailer, sizeof(maest_trailer), 1, file1);
    }
    /*fill spectrum array*/
    for (i = 0; i < maest_header.channels; i++)
    {
        ft.x[i] = (double)i;
        ft.y[i] = (double)*(counts + i);
        ft.dy[i] = sqrt(ft.y[i]);
        /*if (i < 110 && i > 100) printf("ft.y[%d]=%f\n",i,ft.y[i]);*/
    }
    x12[0] = ft.x[0];
    x12[1] = ft.x[maest_header.channels-1];
    lim[0] = (int)(ft.x[0]+0.5);
    lim[1] = (int)(ft.x[maest_header.channels-1]+0.5);

    /*print real and live times to screen and convert from 20ms units to sec*/
    /*Copy day and month*/
    strncpy(dt, maest_header.dt, 5);
    /*check year, if maest_header.dt[7] == 1 year is 2000+YY else 1900+YY*/
    if (maest_header.dt[7] == '1') strncpy(dt+5, "20", 2);
    else strncpy(dt+5, "19", 2);
    strncpy(dt+7, maest_header.dt+5, 2);
    dt[9] = '\0';
    printf("Spectrum info: %s at %c%c:%s\n"
           "               Real time: %d s\n"
           "               Live time: %d s\n",dt,maest_header.sttm[0],
            maest_header.sttm[1],maest_header.sttm+2,
            (int)(maest_header.real*0.02),(int)(maest_header.lve*0.02));
    
    /*print potentially useful trailer information*/
    printf("Maestro energy calibration coeffs: %f %f %f\n",
          maest_trailer.g[0],maest_trailer.g[1],maest_trailer.g[2]);
    
    free(counts);
    fclose(file1);
    return maest_header.channels;
} /*END maestro_read()*/

/*==========================================================================*/
/*matinv: invert matrix     	    	    	    	    	    	    */
/****************************************************************************/
int matinv(double *array, int norder, int dim)
{
    double amax, save, d1;
    int i, j, k, ik[1000], jk[1000];

    /*From RadWare by D.C. Radford at ORNL*/
    
    for (k = 0; k < norder; ++k)
    {
    	/* find largest element array(i,j) in rest of matrix */
    	amax = 0.f;
    	while (1)
	{
      	    for (i = k; i < norder; ++i)
	    {
	    	for (j = k; j < norder; ++j)
		{
	    	    if (fabs(amax) - (d1 = array[i + j*dim], fabs(d1)) <= 0.f)
		    {
	    	    	amax = array[i + j*dim];
	    	    	ik[k] = i;
	    	    	jk[k] = j;
	    	    }
	    	}
      	    }
      	    if (amax == 0.f) return 0;
      	    /* interchange rows and columns to put amax in array(k,k) */
      	    i = ik[k];
      	    if (i < k) continue;
      	    if (i > k)
	    {
	    	for (j = 0; j < norder; ++j)
		{
	    	    save = array[k + j*dim];
	    	    array[k + j*dim] = array[i + j*dim];
	    	    array[i + j*dim] = -save;
	    	}
      	    }
      	    j = jk[k];
      	    if (j >= k) break;
    	}
    	if (j > k)
	{
      	    for (i = 0; i < norder; ++i)
	    {
	    	save = array[i + k*dim];
	    	array[i + k*dim] = array[i + j*dim];
	    	array[i + j*dim] = -save;
      	    }
    	}
    	/* accumulate elements of inverse matrix */
    	for (i = 0; i < norder; ++i)
	{
      	    if (i != k) array[i + k*dim] = -array[i + k*dim] / amax;
    	}
    	for (i = 0; i < norder; ++i)
	{
      	    for (j = 0; j < norder; ++j)
	    {
	    	if (i != k && j != k)
	    	    array[i + j*dim] += array[i + k*dim] * array[k + j*dim];
      	    }
    	}
    	for (j = 0; j < norder; ++j)
	{
      	    if (j != k) array[k + j*dim] /= amax;
    	}
    	array[k + k*dim] = 1.f / amax;
    }
    /* restore ordering of matrix */
    for (k = norder-1; k >= 0; --k)
    {
    	j = ik[k];
    	if (j > k)
	{
      	    for (i = 0; i < norder; ++i)
	    {
	    	save = array[i + k*dim];
	    	array[i + k*dim] = -array[i + j*dim];
	    	array[i + j*dim] = save;
      	    }
    	}
    	i = jk[k];
    	if (i > k)
	{
      	    for (j = 0; j < norder; ++j)
	    {
	    	save = array[k + j*dim];
	    	array[k + j*dim] = -array[i + j*dim];
	    	array[i + j*dim] = save;
      	    }
    	}
    }
    return 0;
} /*END matinv()*/

/*==========================================================================*/
/*ortec_head_srch_skip: search for ASCII ORTEC header and if found, skip    */
/****************************************************************************/
void ortec_head_srch_skip(char *fname, FILE *file, int *mxchan)
{
    float   rlt[2];
    int     hash = 0;
    char    ans[CHLEN] = "";
    fpos_t  pos;
    
    if ( strrchr(fname,'.') && !strcmp(strrchr(fname,'.'),".Spe") )
    {
        printf("%sFound '.Spe' extension...checking for ORTEC ASCII header.%s\n",
                clr[3],clr[0]);
    	/*store current file position*/
    	fgetpos(file, &pos);
        if ( (hash = fgetc(file)) == '#' || isdigit(hash) != 0 || hash == '+'
                || hash == '-') fsetpos(file, &pos);
        else if (hash == '$')
        {
            printf("%sFound ORTEC ASCII header...skipping to data.%s\n",clr[3],clr[0]);
            /*skip_lines(file, 11);*/
            skip_lines(file, 7);
            
            /*get and print date*/
            get_line_file(file, ans, CHLEN);
            printf("%sSpectrum date: %s%s\n",clr[2],ans,clr[0]);
            
            /*get and print real/live time*/
            skip_lines(file, 1);
            rlt[0] = 0.0; rlt[1] = 0.0;
            get_pars_file(file, rlt, 2);
            printf("    %sLive time: %.2f s, Real time %.2f s%s\n",
                clr[2],rlt[0],rlt[1],clr[0]);
            
            /*read start and end channel numbers from file*/
            skip_lines(file, 1);
            get_pars_file(file, rlt, 2);
            /*add one as read loop for data has < not <=*/
            *mxchan = (int)(rlt[1]) + 1;
            printf("    %sFirst chan: %d, last chan %d%s\n",
                clr[2],(int)rlt[0],(int)rlt[1],clr[0]);
        }
        else
        {
            printf("%sUnknown header characters found...\n"
                    "\t....continuing and hoping for the best.%s\n",clr[1],clr[0]);
    	    fsetpos(file, &pos);
        }
    }
} /*END ortec_head_srch_skip()*/

/*==========================================================================*/
/*plot_data: plot data or data and fit                                      */
/****************************************************************************/
void plot_data(char ans[], char cmd[], char inname[], char outname[],
        int plt, int mode)
{
    if (mode == 1) sprintf(ans, "\n Display data with %s? [y/n] (y)",P_PROG);
    else if (! mode) sprintf(ans, "\n Display fit over data with %s? [y/n] (y)",P_PROG);
    
    /*if mode == 3 don't ask, just plot data only*/
    /*if mode == 2 don't ask, just plot fit over data*/
    
    if (plt && (strcmp( (strrchr(inname,'.')), ".Spe"))
              && (strcmp( (strrchr(inname,'.')), ".Chn")) && (mode > 1 || ask_yn(ans,1)))
    {   
        /*printf("plot_data: mode = %d\n",mode);*/
        if (mode > 1) mode -= 2;
        printf("plot_data: mode = %d\n",mode);
        if (mode) sprintf(ans, " %s &\n",inname);
        else sprintf(ans, " %s %s %s %s &\n",P_OPT_BEG, inname, outname, P_OPT_END);
        
        strcpy(cmd+plt, ans);
        printf("Running command: %s\n\t ....Please wait....\n",cmd);
        if (system(cmd)) printf("\t==>command failed\n");
        else
        {
            pid[pcnt++] = plot_pid_get(ans);
        }
    }
    return ;
} /*END plot_data()*/

/*==========================================================================*/
/*plot_pid_get: return pid of most recent instance of P_PROG as int         */
/****************************************************************************/
int plot_pid_get(char cmd[])
{
    FILE    *fp;
    
    /*write pgrep command to cmd[]; -n option gets newest process only*/
    sprintf(cmd,"pgrep -n %s",P_PROG);
    /*execute command*/
    fp = popen(cmd, "r");
    /*read output (pid) into cmd[]*/
    fgets(cmd, 2*CHLEN, fp);
    pclose(fp);
    /*convert pidf to int and return value*/
    return atoi(cmd);
} /*END plot_pid_get()*/

/*==========================================================================*/
/*plot_prog_find: find path to to the plotting program/check if exists      */
/****************************************************************************/
int plot_prog_find(char cmd[])
{
    FILE    *fp;
    
    /*write which command to cmd[]*/
    sprintf(cmd,"which %s",P_PROG);
    /*look for plotting program P_PROG*/
    fp = popen(cmd, "r");
    /*put output, i.e. path if exists, into cmd[]*/
    fgets(cmd, 2*CHLEN, fp);
    pclose(fp);
    /*remove trailing \n by setting to zero*/
    cmd[strcspn(cmd, "\n")] = 0;
    printf("\n****plotting cmd = %s\n\n",cmd);
    /*check if P_PROG can be executed (X_OK)*/
    if (access(cmd, X_OK) == 0)
    {
        printf("%sUsing plotting program: %s%s\n",clr[2],cmd,clr[0]);
        if ( (int)(strlen(P_OPT_BEG)) || (int)(strlen(P_OPT_END)) )
        {
            printf("    %swith option(s)%s",clr[2],clr[0]);
            if ( (int)(strlen(P_OPT_BEG)) ) printf("%s: %s%s%s\n",
                                              clr[2],P_OPT_BEG,
                                              (int)(strlen(P_OPT_END)) ? " and" : "",clr[0]);
            if ( (int)(strlen(P_OPT_END)) ) printf("%s %s: %s%s\n",
                                              (int)(strlen(P_OPT_BEG)) ? "                 " : "",
                                              clr[2],P_OPT_END,clr[0]);
        }
        return (int)strlen(cmd);
    }
    else 
    {
        printf("%sCheck plotting program definition (P_PROG), currently set to: %s%s\n\n",
                clr[1],P_PROG,clr[0]);
        return 0;
    }
}  /*END plot_prog_find()*/

/*==========================================================================*/
/*plot_scale_set: P_PROG_END string with scale parameters 	    	    */
/****************************************************************************/
void plot_scale_set(char ans[],int numch)
{
    int i;
    
    printf("Setting plotting scale parameters:\n");
    /*get spectrum maximum with mode = 0 for silent*/
    spec_max(0);

    plot_scale_set_par("XMIN",ft.x[0],ans);
    plot_scale_set_par("YMIN",ft.y[0],ans);
    /*XMAX rounded up to nearest 50*/
    plot_scale_set_par("XMAX",(double)(50*(int)(ft.x[numch-1]/50.0 + 1.0)),ans);
    /*YMAX plus twice the uncertainty rounded up to nearest 10*/
    plot_scale_set_par("YMAX",10*(int)((sp_max+2.0*sqrt(sp_max))/10.0 + 1.0),ans);
    /*XTCKMAJ to a fifth of the maximum rounded down to nearest 10*/
    if (ft.x[numch-1] > 100.0) i = 20.0*(int)(ft.x[numch-1]/100.0);
    else if (ft.x[numch-1] > 30.0) i = 10;
    else if (ft.x[numch-1] > 15.0) i = 5;
    else i = 5;
    plot_scale_set_par("XTCKMAJ",i,ans);
    /*note XTCKMIN is the number of ticks not the separation*/
    plot_scale_set_par("XTCKMIN",5,ans);
    /*copy remainder of string*/
    plot_scale_set_par("FFFF",0,ans);
    printf("P_OPT_END: %s\n",P_OPT_END);
    return ;
} /*END plot_scale_set()*/

/*==========================================================================*/
/*plot_scale_set_par: set min, max and tick intervals for plotting 	    */
/****************************************************************************/
void plot_scale_set_par(char p_par[], double val, char ans[])
{
    static int i = 0, j = 0;
    
    while ( i < CHLEN && j < strlen(P_OPT_FIN) && strncmp(P_OPT_FIN+j,p_par,strlen(p_par)) )
    {
        /*copy characters while no match*/
        P_OPT_END[i++] = P_OPT_FIN[j++];
        /*if (i > 0 && j > 0)
            printf("i=%d j= %d (%c:%c)\n",i,j,P_OPT_END[i-1],(char)P_OPT_FIN[j-1]);*/
    }
    if (j >= strlen(P_OPT_FIN)) return ;
    
    /*convert val to string*/
    itoa((int)val,ans);
    strncat(P_OPT_END+i,ans,strlen(ans));
    /*increase offsets ready for next copy*/
    i += strlen(ans);
    j += strlen(p_par);
    return ;
} /*END plot_scale_set_par()*/

/*==========================================================================*/
/*prep_write: prepare to write data to file in ASCII format                 */
/****************************************************************************/
void prep_write(char inname[], char outname[], char nm[], char ans[],
        int numch, int col)
{ 
    strcpy(nm, inname);
    /*put col as char in ans[1]*/
    itoa(col, ans+1);
    ans[0] = '_';
    strncpy(ans+2, "col", 3);
    /*if not equal to NULL, find '.' then copy "_?col"*/
    if ( strrchr(nm,'.') )
    {
        /*check if filename contains _?col already*/
        /*if strncmp is 0 identical*/
        /*if not identical copy extension to end of ans*/
        if ( strncmp( (strrchr(nm,'.')-5), ans, 5) )
            strncpy( ans+5, (strrchr(nm,'.')), 5);	        
        
        /*else identical so copy just extension to ans*/
        else strcpy( ans, (strrchr(nm,'.')));
        
        /*copy ans to nm*/
        strcpy( (strrchr(nm,'.')), ans );
        
        /*check for .Spe/.Chn extension and if found change to .txt*/
        /*strcmp returns zero if identical, i.e. if not equal to NULL*/
        if ( (!strcmp( (strrchr(nm,'.')), ".Spe" )) ||
              (!strcmp( (strrchr(nm,'.')), ".Chn" )) )
            set_ext(nm, ".txt");
    }
    /*if equal to NULL just add _?col to the end*/
    else
    {
        strcat(ans, ".txt");
        /*cat ans to nm*/
        strcat(nm, ans);
    }
    printf("Output filename for %d column data: %s\n",col,nm);
    file_status(nm, CHLEN);
    /*note numch is full spectrum length including zeros. Ignoring limits*/
    ascii_write(nm, numch, col);
    strcpy(inname, nm);
    strcpy(outname, inname);
    set_ext(outname, ".fit");
    printf("New input and output file names are:\n  Input: %s\n  Output: %s\n",
            inname,outname);
    return ;
} /*END prep_write()*/

/*==========================================================================*/
/*read_par: read user input	    	    	    	    	    	    */
/****************************************************************************/
void read_par(double pars[], int offset, int num, char s[])
{
    int     i = 0;
    char    ans[CHLEN] = "";
    
    printf("%s",s);
    printf(" [<Enter> for ");
    for (i = 0; i < num; i++) printf("%f ",pars[offset+i]);
    printf("]\n");
    get_line(ans, CHLEN);
    if (strlen(ans) == (int)0 && (ans[0] == ((char)0))) ;
    else get_pars(ans, pars+offset, num);
}/*END read_par()*/

/*==========================================================================*/
/* reverse: reverse string s in place	    	    	    	    	    */
/****************************************************************************/
void reverse(char s[])
{
    int c, i, j;
	
    for (i = 0, j = strlen(s)-1; i < j; i++, j--)
    {
        c = s[i];
	s[i] = s[j];
	s[j] = c;
    }
} /*END reverse()*/

/*==========================================================================*/
/*scaler: workout approximate scaling factor	    	    	    	    */
/****************************************************************************/
void scaler(double *fit, double *derivs, int mode)
{
    double  step, dsqr, dnsqr, sign, osign;
    double  exp = 0.0, ex = 0.0, sum = 0.0, thr = 0.0, st[16384][2];
    int     i, isl, isu, lp = 1, cnt = 0;

    /*From RadWare by D.C. Radford at ORNL*/
      
    for (i = 0; i < 16384; i++)
    {
	st[i][0] = 0.0;
	st[i][1] = 0.0;
    }
    
    isl = ft.x[0]+lim[0];
    isu = ft.x[ft.ndp - 1+lim[0]];
    for (i = (0+lim[0]); i < (ft.ndp+lim[0]); i++)
    {
	ex = ft.x[i];
	exp = ft.y[i];
	sum += exp;
	st[(int)ex][0] = exp;
    	eval(ft.pars, ft.x[i], fit, derivs, i, 0);
	thr = *fit;
	st[(int)ex][1] = thr;
    }
    /*for mode=1 just calculate curve, don't attempt a fit*/
    if (mode == 1) return ;
    /*sum is total no. of experimental counts*/
    ft.pars[3] = sum;
    step = sum/10.0;
    sign = 0.0;
    while ( (step > 1.0 || lp == 1) && lp < 50)
    {
	cnt++;
	lp++; 
	osign = -1.0*sign;
	dsqr = 0.0;
	dnsqr = 0.0;
	for (i = isl; i <= isu; i++)
	{
	    exp = st[i][0];
	    thr = st[i][1];
	    dsqr += (exp - (thr*ft.pars[3]))*(exp - (thr*ft.pars[3]));
	    dnsqr += (exp-(thr*(ft.pars[3]+1)))*(exp-(thr*(ft.pars[3]+1)));
	}
	if (dnsqr > dsqr) sign = 1.0;
	else sign = -1.0;
	
	if (sign == osign) step *= 0.4;
	
	ft.pars[3] -= sign*step;
    }
    ft.pars[3] += sign*step; 
    printf("  %d iterations for scaling value\n",lp); 
} /*END scaler()*/

/*==========================================================================*/
/* set_ext: set file extension of string fname[] to ext[]   	    	    */
/****************************************************************************/
void set_ext(char fname[], char ext[])
{    
    /*if not equal to NULL, find '.' then copy ext*/
    if ( (strrchr(fname,'.')) ) strcpy( (strrchr(fname,'.')), ext );
    /*if equal to NULL just add ext to the end*/
    else strcat( fname, ext );
} /*END set_ext()*/

/*==========================================================================*/
/* skip_hash: skip lines starting with hash (#) from current file position  */
/****************************************************************************/
void skip_hash(FILE *file)
{
    int     hash = 0;
    fpos_t  pos;
     	   
    /*store position corresponding to start position*/
    fgetpos(file, &pos);
    
    /*test for '\n' and '\r' (Ortec files), necessary as fscanf might not*/
    /*read past it. If not there go back to starting point.*/
    /*in Ortec files there can be more than 1 in a row, hence while loop*/
    while ( (hash = fgetc(file)) == '\n' || hash == '\r' ) fgetpos(file, &pos);
    
    /*set position back to start of line*/
    fsetpos(file, &pos);
    while (1)
    {
    	/*store position corresponding to start of line*/
    	fgetpos(file, &pos);
    	if ( (hash = fgetc(file)) != '#' )
    	{
    	    /*reset position to start of line*/
    	    fsetpos(file, &pos);
    	    break;
    	}
    	/*read rest of line*/
    	while ( ( (hash = fgetc(file)) != '\n' ) && (hash != EOF) )
            ;
    }
} /*END skip_hash()*/

/*==========================================================================*/
/*skip_lines: skip lns lines in file                                        */
/****************************************************************************/
void skip_lines(FILE *file, int lns)
{
    int cnt = 0, hash = 0;
    
    while (cnt < lns)
    {
        while ( ( (hash = fgetc(file)) != '\n' ) && (hash != EOF) )
           ;
        
        if (hash == EOF)
        {
            printf("Found EOF...returning\n");
            return ;
        }
        /*increment counter on each carriage return*/
        cnt++;
     }
} /*END skip_lines()*/

/*==========================================================================*/
/*spec_int: integrate spectrum                                              */
/****************************************************************************/
void spec_int(void)
{
    int i;
    
    sp_sum = 0.0; 
    for (i = (0+lim[0]); i < (ft.ndp+lim[0]); i++) sp_sum += ft.y[i];
    
    printf("%sIntegral: spectrum area =%12.1f cnts%s\n",clr[2],sp_sum,clr[0]);
} /*END spec_int()*/

/*==========================================================================*/
/*spec_max: find location of spectrum maximum                               */
/****************************************************************************/
void spec_max(int mode)
{
    int ch = 0, i;
    
    sp_max = 0.0; 
    for (i = (0+lim[0]); i < (ft.ndp+lim[0]); i++)
    {
        if (sp_max < ft.y[i])
        {
            sp_max = ft.y[i];
            ch = i;
        }
    }
    
    if (mode) printf("%sMaximum: spectrum max:%12.1f ch %12.1f cnts%s\n",
                  clr[2],ft.x[ch],sp_max,clr[0]);
} /*END spec_max()*/

/*==========================================================================*/
/*spec_prop: find spectrum properties: maxima, FHWM and half-life estimates */
/****************************************************************************/
void spec_prop(double init[])
{
    double  mx, sp_mx[3][2], wdth[2];
    double  mn, sp_mn[3][2];
    int     cntz, i, j, k, win, vb = 0;
    
    /*zero sp_mx arrays for storing maxima
      and set sp_mn array for minima to large values*/
    for (i = 0; i < 3; i++)
    {
        sp_mx[i][0] = 0.0;
        sp_mx[i][1] = 0.0;
        sp_mn[i][0] = (double)INT_MAX;
        sp_mn[i][1] = (double)INT_MAX;
    }
    wdth[0] = 0.0;
    wdth[1] = 0.0;
    
    /*printf("sp_mn[i-1][1]=%.1f\n",sp_mn[i-1][1]);*/
    
    /*find three maxima using window of width win channels*/
    win = WN;
    for (i = (0+lim[0]); i < (ft.ndp+lim[0]-win); i++)
    {
        mx = 0.0;
        for (j = 0; j < win; j++) mx += ft.y[i+j];
        
        /*compare to previous lowest maxima and if greater store new values*/
        k = 2;
        if (k > 0 && mx > sp_mx[k][1])
        {
            /*store new integration*/
            sp_mx[k][1] = mx;
            /*store new central channel*/
            sp_mx[k][0] = (double)i + ((double)win/2.0);
            /*now check other maxima and swap order as required*/
            while (k > 0 && mx > sp_mx[k-1][1]) 
            {
                /*shift old values to next element in array*/
                sp_mx[k][1] = sp_mx[k-1][1];
                sp_mx[k][0] = sp_mx[k-1][0];   
                
                sp_mx[k-1][1] = mx;
                sp_mx[k-1][0] = (double)i + ((double)win/2.0);
            
                k--;
            }
        }
    }
    /*print maxima results*/
    if (vb)
    {
        printf("Spectrum maxima (x3) using window size of %d channels for integration.\n",win);
        for (i = 0; i < 3; i++)
            printf("ch = %5f int = %10.1f\n",sp_mx[i][0],sp_mx[i][1]);
    }
    
    /*find three minima using window of width win channels*/
    /*note, a reliable minimum is needed to approximate the background level*/
    win = WN;
    for (i = (0+lim[0]); i < (ft.ndp+lim[0]-win); i++)
    {
        mn = 0.0;
        cntz = 0;
        for (j = 0; j < win; j++)
        {
            mn += ft.y[i+j];
            /*count ft.y[] values if equal to zero*/
            if (ft.y[i+j] == 0.0) cntz++;
        }
        
        /*compare to previous lowest minima and if lower store new values
          with additional check that few than half ft.y[] values are zero*/
        k = 2;
        if (k > 0 && mn < sp_mn[k][1] && mn != 0 && cntz < (int)(win/2.0) )
        {
            /*store new integration*/
            sp_mn[k][1] = mn;
            /*store new central channel*/
            sp_mn[k][0] = (double)i + ((double)win/2.0);
            /*now check other minima and swap order as required*/
            while (k > 0 && mn < sp_mn[k-1][1]) 
            {
                /*shift old values to next element in array*/
                sp_mn[k][1] = sp_mn[k-1][1];
                sp_mn[k][0] = sp_mn[k-1][0];   
                
                sp_mn[k-1][1] = mn;
                sp_mn[k-1][0] = (double)i + ((double)win/2.0);
            
                k--;
            }
        }
    }
    /*print minima results*/
    if (vb)
    {
        printf("Spectrum minima (x3) using window size of %d channels for integration.\n",win);
        for (i = 0; i < 3; i++)
            printf("ch = %5f int = %10.1f\n",sp_mn[i][0],sp_mn[i][1]);
    }
    
    
    /*Now, take maximum and move LEFT to find where integration drops to half*/
    for (i = (int)(sp_mx[0][0] - ((double)win/2.0)); i > lim[0]; i--)
    {
        mx = 0.0;
        for (j = 0; j < win; j++) mx += ft.y[i+j];
        
        if (mx <= (0.5*(sp_mx[0][1]-sp_mn[0][1]) + sp_mn[0][1]) )
        {
            wdth[0] = (double)i + ((double)win/2.0);
            break;
        }
    }
    /*Now, take maximum and move RIGHT to find where integration drops to half*/
    for (i = (int)(sp_mx[0][0] - ((double)win/2.0)); i < (ft.ndp+lim[0]-win); i++)
    {
        mx = 0.0;
        for (j = 0; j < win; j++) mx += ft.y[i+j];
        
        if (mx <= (0.5*(sp_mx[0][1]-sp_mn[0][1]) + sp_mn[0][1]))
        {
            wdth[1] = (double)i + ((double)win/2.0);
            break;
        }
    }
    if (vb) printf("Widths around maximum (ch:%d):\n"
           "    0.5FWHM  left at %.1f ch (FWHM = %.1f ch)\n"
           "    0.5FWHM right at %.1f ch (FWHM = %.1f ch)\n",
          (int)sp_mx[0][0],wdth[0],2.0*(sp_mx[0][0]-wdth[0]),
          wdth[1],2.0*(wdth[1]-sp_mx[0][0]));
    /*set widths to be differences from spectrum max.*/
    wdth[0] -= sp_mx[0][0];
    wdth[1] -= sp_mx[0][0];
    
    if (abs(wdth[1]) > abs(wdth[0]))
    {
        mx = wdth[1];
        wdth[1] = wdth[0];
        wdth[0] = mx;
    }
    
    /*fill initial parameters*/
    /*init[0] is half-life estimate*/
    init[0] = abs(wdth[0]);
    /*init[1] is FWHM estimate*/
    init[1] = 2.0*abs(wdth[1]);
    /*init[2] is centroid estimate*/
    init[2] = sp_mx[0][0];
    /*init[NUMPARS] is background estimate
        (increase win by 50% to make sure bkgnd isn't over estimated*/
    init[NUMPARS] = sp_mn[0][1]/((double)(win*1.5));
} /*END spec_prop()*/

/*==========================================================================*/
/* store_colours: store colours in clr[][] array                            */
/****************************************************************************/
void store_colours()
{
    /*none*/
    strcpy(clr[0], "\033[0;0m");
    /*white on red background*/
    strcpy(clr[1], "\033[0;41m");
    /*green no background*/
    strcpy(clr[2], "\033[0;32m");
    /*white on blue background*/
    strcpy(clr[3], "\033[0;44m");
    /*bold*/
    strcpy(clr[4], "\033[1m");
    /*underline*/
    strcpy(clr[5], "\033[0;4m");
} /*END store_colours()*/  

/*==========================================================================*/
/* swapb2: swap 2 array elements (a 2 byte number) in place 	    	    */
/****************************************************************************/
void swapb2(char *buf)
{
    char c;
    
    c = buf[1]; buf[1] = buf[0]; buf[0] = c;
} /*END swapb2()*/

/*==========================================================================*/
/* swapb4: swap 4 array elements (a 4 byte number) in place 	    	    */
/****************************************************************************/
void swapb4(char *buf)
{
    char c;
    
    c = buf[3]; buf[3] = buf[0]; buf[0] = c;
    c = buf[2]; buf[2] = buf[1]; buf[1] = c;    
} /*END swapb4()*/

/*==========================================================================*/
/* write_output: write output to file and coefficients to screen    	    */
/****************************************************************************/
void write_output(char inname[], char outname[], double init[], double x[],
	double *fit, double chisq, int mode)
{
    double  xx = 0.0, sum;
    int     i, j, nip, ndf;
    FILE    *fout = 0;
    
    /*open output file*/
    if (mode != 1)
    {
    	if ((fout = fopen(outname, "w" )) == NULL)
    	{
    	    printf("Cannot open file: %s \n", outname);
    	    return ;
	}			
    }
    
    /*no. independent parameters*/
    nip = NUMPARS - ft.nfp;
    /*no. degrees of freedom*/
    ndf = ft.ndp - nip;

    /*write each initial guesses to output file*/    
    /*mode == 0 correspond to a good fit*/
    /*mode == 1 correspond to print current values to screen option*/
    /*mode == 2 corresponds to no good fit*/
    printf("\n++++++++++++++++Results of fit++++++++++++++++\n");    
    for (i = 0; i < NUMPARS; i++)
    {
        if ( abs(ft.pars[i]/fact[i]) < 0.01 || abs(ft.pars[i]/fact[i]) >= 1.00e5)
            printf("%c%1d) %s = %.3g",
                (i < 10) ? ' ' : '1',(i < 10) ? i : i-10,par_nm[i],(ft.pars[i]/fact[i]));
        else printf("%c%1d) %s = %.4f",
                (i < 10) ? ' ' : '1',(i < 10) ? i : i-10,par_nm[i],(ft.pars[i]/fact[i]));
        
        /*printf(" %s = %.3f",par_nm[i],(i == 0) ? (f*ft.pars[i]/fact[i]) : (ft.pars[i]/fact[i]));*/
        
        if (mode != 2 || (mode == 2 && ft.freepars[i])) printf(" %s",(ft.freepars[i]) ? "+/- " : "FIXED\n");
        if (ft.freepars[i]) printf("%.3f\n",(ft.errs[i]/fact[i]));
        if (mode == 2 && !ft.freepars[i]) printf(" (no fit)\n");
    }   
    printf("    %s = %.2f\n",par_nm[NUMPARS],bkgnd);
    /*if md = 0, good fit, so print chisq in green*/
    if (mode == 0) printf("\n# %sChisq/D.O.F. = %.5f (Chisq = %.4f, D.O.F. = %d)%s\n",
                      clr[2],chisq,chisq*ndf,ndf,clr[0]);
    if (mode != 0) printf("\n# Chisq/D.O.F. = %.5f (Chisq = %.4f, D.O.F. = %d)\n",chisq,chisq*ndf,ndf);
    printf("++++++++++++++++++++++++++++++++++++++++++++++\n");
    
    /*for mode == 1 just write out result of fit don't write to file*/
    if (mode == 1) return ;
    
    fprintf(fout, "# Input file: %s\n",inname);
    fprintf(fout, "# Integral: spectrum area =%12.3f cnts\n",sp_sum);
    fprintf(fout, "# Output file: %s\n",outname);
    for (i = 0; i < NUMPARS; i++)
        fprintf(fout, "#%c%1d) Initial guess %s = %.2f\n",
            (i < 10) ? ' ' : '1',(i < 10) ? i : i-10,par_nm[i],init[i]);
    
    fprintf(fout, "#++++++++++++++++Results of fit++++++++++++++++\n");

    for (i = 0; i < NUMPARS; i++)
    {
        if ((ft.pars[i]/fact[i]) < 0.01 || (ft.pars[i]/fact[i]) >= 1.00e5)
            fprintf(fout, "#%c%1d) %s = %.3g",
                (i < 10) ? ' ' : '1',(i < 10) ? i : i-10,par_nm[i],(ft.pars[i]/fact[i]));
        else fprintf(fout, "#%c%1d) %s = %.4f",
                (i < 10) ? ' ' : '1',(i < 10) ? i : i-10,par_nm[i],(ft.pars[i]/fact[i]));

        if (mode != 2 || (mode == 2 && ft.freepars[i])) fprintf(fout, " %s",(ft.freepars[i]) ? "+/- " : "FIXED\n");
        if (ft.freepars[i]) fprintf(fout, "%.3f\n",(i == 0) ? abs(ft.errs[i]/fact[i]) : (ft.errs[i]/fact[i]));
        if (mode == 2 && !ft.freepars[i]) fprintf(fout, " (no fit)\n");
    }            

    fprintf(fout, "#    %s = %.2f\n",par_nm[NUMPARS],bkgnd);
    fprintf(fout, "# Chisq/D.O.F. = %.5f (Chisq = %.4f, D.O.F. = %d)\n",chisq,chisq*ndf,ndf);
    fprintf(fout, "#++++++++++++++++++++++++++++++++++++++++++++++\n");   
    fprintf(fout, "#   Chan \t    fit\n");
    
    /*write the fitted points to output file*/ 
    if (mode != 3)
    {
        /*write main fit points*/
        for (i = (0+lim[0]); i < (ft.ndp+lim[0]); i++)
        {
    	    eval(ft.pars, x[i], fit, derivs, i, 0);
    	    fprintf(fout, "%8.3f\t%10.5f\n", x[i],*fit);    
        }
        /*write each function separately in turn*/
    	printf("Integrals:\n");
        printf("    Spectrum      %*s area =%12.1f cnts\n",CHFLEN,"",sp_sum);
        for (j = 0; j < (NUMFUNC+1); j++)
        {
            fprintf(fout, "& #Function %d (%s)\n",j,func_nm[j]);
            sum = 0.0; 
            for (i = (0+lim[0]); i < (ft.ndp+lim[0]); i++)
            {
    	        eval(ft.pars, x[i], fit, derivs, i, 2+j);
    	        fprintf(fout, "%8.3f\t%10.5f\n", x[i],*fit);
                /*integrate the area of the function*/
                sum += *fit;
            }
    	    fprintf(fout, "#End function %d (%s),%*s area =%12.3f cnts\n",
                    j,func_nm[j],CHFLEN-(int)strlen(func_nm[j]),"",sum);
    	    printf("    Function %d (%s),%*s area =%12.1f cnts\n",
                    j,func_nm[j],CHFLEN-(int)strlen(func_nm[j]),"",sum);
        }
    }
    else if (mode == 3)
    {
        i = 0;
        for (xx = x12[0]; xx <= (x12[1] + (spc[0]/2)); xx += spc[0])
        {
            eval(ft.pars, xx, fit, derivs, i++, 0);
    	    fprintf(fout, "%8.3f\t%10.5f\n", xx,*fit);    
        }
    }
    fclose(fout);
    printf(" --> Output written to file: %s\n",outname);
} /*END write_output()*/
