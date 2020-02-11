 #define _CRT_SECURE_NO_WARNINGS

//  This file contains sourcecode distributed as freeware. 
//  The intellectual property of the sourcecode is shown 
//  here to belong to Carlo Ciulla.

// Disclaimer: 

// The website here named www.sourcecodewebsiteCarloCiulla.com 2013 [1] does not intend 
// to convey the meaning of profit making for what pertains to the content
// provided. --->>> Instead, when the content is downloaded, the user(s) are
// kindly invited to donate money to charity organizations involved in 
// helping people in need of food and water. <<<---


// The Novel Re-sampling Locations have been sized to be a fraction of 
// the pixel size. The programs presented here confirm both concepts and 
// implications brought to knowledge through the unifying theory [1].


// Reference:

// [1] Carlo Ciulla "Improved Signal and Image Interpolation in Biomedical Applications: 
// The Case of Magnetic Resonance Imaging (MRI)." Medical Information Science 
// Reference - IGI Global Publisher - March 2009; ISBN: 978 - 160566202 - 2.

//  Project Title: Sinc Function Odd (Sinc Function SRE-based Interpolation)

#include < iostream >
#include < fstream >
#include < ostream >
#include < string >
#include < io.h >
#include < dos.h >
#include < conio.h >
#include < stdlib.h >
#include < sstream >
#include < stdio.h >
#include < iomanip >
#include < istream >
#include < math.h >
#include < ctype.h >

#define SCALE 0.000001
#define NSCALE 0.01
#define PLUS 1          // One is positive solution for the novel recomputed locations, Zero is set for negative solutions 

#define TH 0

using namespace std;


class SINCOSRESIMUL2013 {

	int n1; // matrix size x
	int n2; // matrix size y
    int n8;  // neighborhood size


public:

	int getNofPixelsX(void) { return this->n1; };

	int getNofPixelsY(void) { return this->n2; };

	void setNofPixelsX(int x) { this->n1 = x; };

	void setNofPixelsY(int y) { this->n2 = y; };

	int getN8(void) { return this->n8; };

	void setN8(int x) { this->n8 = x; };


public:

	struct data {

		double **the_fMRI; // pointer to the matrix entry 

		double **the_m_corr_recomp_fMRI; // pointer to the matrix entry

		double **shifted_fMRI; // pointer to the matrix entry
		
		double **fMRI; // pointer to the matrix entry 

		double **Xsre; // pointer to the matrix entry 

		double **SE; // pointer to the matrix entry 

		double **SE_SRE; // pointer to the matrix entry 

		double **X_recomputed; // pointer to the matrix entry 

		double **X_recomputed_scaled; // pointer to the matrix entry

		double **h_x_sre; // pointer to the matrix entry
		
		double **de2_h_xre; // pointer to the matrix entry

		double **Omega; // pointer to the matrix entry

		double **theta2; // pointer to the matrix entry

		double **theta3; // pointer to the matrix entry

		double **theta4; // pointer to the matrix entry 

		double **ICM; // pointer to the matrix entry

	}*pointer; // pointer to the matrices

public:

	SINCOSRESIMUL2013(int x, int y, int n) : n1(x), n2(y), n8(n) { };

	void allocateData();

	~SINCOSRESIMUL2013() { } // destructor

};

void SINCOSRESIMUL2013::allocateData() { // allocate data

	// (1) allocate struct 'data' (begin)
	 pointer = new data;
	
	 pointer->the_fMRI = new double*[this->n1];

	 pointer->the_m_corr_recomp_fMRI = new double*[this->n1];

	 pointer->shifted_fMRI = new double*[this->n1];

	 pointer->fMRI = new double*[this->n1];

	 pointer->Xsre = new double*[this->n1];

	 pointer->SE = new double*[this->n1];

	 pointer->SE_SRE = new double*[this->n1];

	 pointer->X_recomputed = new double*[this->n1];

	 pointer->X_recomputed_scaled = new double*[this->n1];

     pointer->h_x_sre = new double*[this->n1];

	 pointer->de2_h_xre = new double*[this->n1];

	 pointer->Omega = new double*[this->n1];

	 pointer->theta2 = new double*[this->n1];

     pointer->theta3 = new double*[this->n1];

	 pointer->theta4 = new double*[this->n1];

	 pointer->ICM = new double*[this->n1];


	
	 for( int v=0; v <= 1; v++ ) { // (a)

	 pointer->theta3[v] = new double[this->n8 + 1];

	 pointer->theta4[v] = new double[this->n8 + 1];


	 } // (a)

	 for( int v=0; v < this->n1; v++ ) { // (1)

		 pointer->the_fMRI[v] = new double[this->n2];

		 pointer->the_m_corr_recomp_fMRI[v] = new double[this->n2];

		 pointer->shifted_fMRI[v] = new double[this->n2];

		 pointer->fMRI[v] = new double[this->n2];
		 
		 pointer->Xsre[v] = new double[this->n2];

		 pointer->SE[v] = new double[this->n2];

		 pointer->SE_SRE[v] = new double[this->n2];

		 pointer->X_recomputed[v] = new double[this->n2];

		 pointer->X_recomputed_scaled[v] = new double[this->n2];

		 pointer->h_x_sre[v] = new double[this->n2];

		 pointer->de2_h_xre[v] = new double[this->n2];

		 pointer->Omega[v] = new double[this->n2];
 
		 pointer->theta2[v] = new double[this->n2];

		 pointer->ICM[v] = new double[this->n2];

	  } // (1) allocate struct 'data' (end)


	 	// (2) initialize (begin)
		for( int v=0; v <= 1; v++ ) { // (a)

			for( int f=0; f < this->n8 + 1 ; f++ ) { // (b)

			pointer->theta3[v][f] = (double)0.0;

		    pointer->theta4[v][f] = (double)0.0;

			} //(b)

		 } //(a)

		// (2) initialize (begin)
		for( int v=0; v < this->n1; v++ ) { // (a)

			for( int f=0; f < this->n2 ; f++ ) { // (b)

			pointer->the_fMRI[v][f] = (double)0.0;

			pointer->the_m_corr_recomp_fMRI[v][f] = (double)0.0;

			pointer->shifted_fMRI[v][f] = (double)0.0;

			pointer->fMRI[v][f] = (double)0.0;
		 
			pointer->Xsre[v][f] = (double)0.0;

			pointer->SE[v][f] = (double)0.0;

			pointer->SE_SRE[v][f] = (double)0.0;

			pointer->X_recomputed[v][f] = (double)0.0;

			pointer->X_recomputed_scaled[v][f] = (double)0.0;

			pointer->h_x_sre[v][f] = (double)0.0;

			pointer->de2_h_xre[v][f] = (double)0.0;

			pointer->Omega[v][f] = (double)0.0;

		    pointer->theta2[v][f] = (double)0.0;

			pointer->ICM[v][f] = (double)0.0;

		 
			} //(b)

		 } //(a)
		// (2) initialize (end)

} // allocate data


int main ( int argc, char * argv[] ) {

	char outputFile[128]="SINCOSRESIMUL.log";

	FILE * savedata;

	double MAX = 5000000000000000000.0;
	
	if ((savedata = fopen(outputFile,"w"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // processing (begin)


	std::cout << endl;
	std::cout << "Results from SINCOSRESIMUL2013" << endl;
	fprintf(savedata,"%s\n", "Results from SINCOSRESIMUL2013");
	fprintf(savedata,"\n");

	double RESOL = (double)0.87; 
	int n2 = 1;									  // number of pixels y dim
	int n1 = (int)floor( ((double)360.0/RESOL) ); // number of pixels x dim
	double STEPSIZE = ((double)RESOL/360.0);  

	double NEI = 3.0;     // neighborhood
	double FW = 0.87;    // FWHM

	int n7 = (int)floor((double)NEI/2.0);  // number of right/left pixels to include into the filter
	int XNEI = (int)n7;					   // neighborhood
	int n8 = (int)floor((double)(NEI-1.0)/2.0);

	double fw = ((double) FW / (2.0 * sqrt(2.0 * (double)log((double)2.0))) ); 
	// 'fw' converts  FWHM to the standard deviation
	
	std::cout << "fw = " << fw << " n7 = " << n7 << " n8 = " << n8 << endl;
	
	// build 1D Gauss filters for the FILTER function (begin)
	double * g_filt = 0;

	if (( g_filt = (double *) calloc( ((int)NEI+1), sizeof(double)) ) == NULL)
	{ // allocate memory (begin)

		std::cout << "Not enough memory to allocate Image filter data: Exit." << endl;
		fprintf(savedata,"%s\n", "Not enough memory to allocate Image filter data: Exit.");
	
		exit(0);

	}// allocate memory (end)

	double filter;
	for (int i =-n7; i <= n7; i++) {
	
		filter = (double) exp( -((double)i*i/(2.0*fw*fw)) );

		*(g_filt + i + n7) = (double)filter;	// output 1D filter: g_filt
	
		std::cout << "g_filt = " << *(g_filt + i + n7) << endl;

	} // build 1D Gauss filters for the FILTER function (end)

	/// read image file (begin)

	SINCOSRESIMUL2013 SRE(n1,n2,n8);

	SRE.allocateData();

	std::cout << "Object Constructed" << endl;

	std::cout << "X image size: " << SRE.getNofPixelsX() << endl;
	std::cout << "Y image size: " << SRE.getNofPixelsY() << endl;

	std::cout << "Please wait. The results will be saved into: " << outputFile << endl;


	for (int SIM = 1; SIM <= 3; SIM++) { /// loop across the simulations: sin, cos, log (begin)

	double misplacement_X = ((double)RESOL/360.0);

	// initialize data (begin)
	if ( SIM == 1 ) 
	{ // run sin simulation 

	std::cout << endl;
	std::cout << "Results from SINCOSRESIMUL2013 sin(x) are being processed" << endl;
	std::cout << endl;
	fprintf(savedata,"\n");
	fprintf(savedata,"%s\n", "Results from SINCOSRESIMUL2013 sin(x)");
	fprintf(savedata,"\n");

	for (int i1=0; i1 < SRE.getNofPixelsX(); i1++) {// x dim
       	
		for (int i2=0; i2 < SRE.getNofPixelsY(); i2++) { // y dim

		double signal = (double)sin( ( (double)2.0*3.1415926535897*((double)(i1 * misplacement_X)) ) );

		SRE.pointer->fMRI[i1][i2] = (double)signal;
                          
		} // y dim
        
	}  // x dim 


	for (int i1=0; i1 < SRE.getNofPixelsX(); i1++) {// x dim
       	
		for (int i2=0; i2 < SRE.getNofPixelsY(); i2++) { // y dim

			SRE.pointer->the_fMRI[i1][i2] = (double)SRE.pointer->fMRI[i1][i2];

		} // y dim
        
	}  // x dim 
	
	} // run sin simulation
    else if ( SIM == 2) 
	{ // run cos simulation 

	std::cout << endl;
	std::cout << "Results from SINCOSRESIMUL2013 cos(x) are being processed" << endl;
	std::cout << endl;
	fprintf(savedata,"\n");
	fprintf(savedata,"%s\n", "Results from SINCOSRESIMUL2013 cos(x)");
	fprintf(savedata,"\n");

	for (int i1=0; i1 < SRE.getNofPixelsX(); i1++) {// x dim
       	
		for (int i2=0; i2 < SRE.getNofPixelsY(); i2++) { // y dim

		double signal = (double)cos( ( (double)2.0*3.1415926535897*((double)(i1 * misplacement_X)) ) );

		SRE.pointer->fMRI[i1][i2] = (double)signal;
                          
		} // y dim
        
	}  // x dim 


	for (int i1=0; i1 < SRE.getNofPixelsX(); i1++) {// x dim
       	
		for (int i2=0; i2 < SRE.getNofPixelsY(); i2++) { // y dim

			SRE.pointer->the_fMRI[i1][i2] = (double)SRE.pointer->fMRI[i1][i2];

		} // y dim
        
	}  // x dim 
	
     	
	} // run cos simulation  
	else if ( SIM == 3) 
	{ // run log simulation 

	std::cout << endl;
	std::cout << "Results from SINCOSRESIMUL2013 log(x) are being processed" << endl;
	std::cout << endl;
	fprintf(savedata,"\n");
	fprintf(savedata,"%s\n", "Results from SINCOSRESIMUL2013 log(x)");
	fprintf(savedata,"\n");

	for (int i1=0; i1 < SRE.getNofPixelsX(); i1++) {// x dim
       	
		for (int i2=0; i2 < SRE.getNofPixelsY(); i2++) { // y dim

		double signal = ((double)log(i1 + misplacement_X));

		SRE.pointer->fMRI[i1][i2] = (double)signal;
                          
		} // y dim
        
	}  // x dim 


	for (int i1=0; i1 < SRE.getNofPixelsX(); i1++) {// x dim
       	
		for (int i2=0; i2 < SRE.getNofPixelsY(); i2++) { // y dim

			SRE.pointer->the_fMRI[i1][i2] = (double)SRE.pointer->fMRI[i1][i2];

		} // y dim
        
	}  // x dim 
	
	 
	 // run log simulation 
	}// initialize data (end)

	// standardize data (begin)
	double max=-MAX;
	double min=MAX;

	double Std = 0.0;
	double Average = 0.0;
	long int counting = 0;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

			if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) {

			Average += (double) SRE.pointer->fMRI[i1][i2];

			counting++;

			}
			
		} // y dim
        
	}  // x dim 

	Average /= ((double)counting);

	counting = 0;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

			if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) {

			Std = (double)Std + ((double) (Average - SRE.pointer->fMRI[i1][i2]) * 
				                          (Average - SRE.pointer->fMRI[i1][i2]) );

			counting++;

			}

		} // y dim
        
	}  // x dim 

	Std = (double) sqrt( (double) Std / ((double)counting) );


	// standardize (begin)
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

           SRE.pointer->the_fMRI[i1][i2] = (double) NSCALE * ( (double) exp((double) (SRE.pointer->fMRI[i1][i2] - Average) / Std) );
       
		} // y dim
        
	}  // x dim 
	// standardize data (end)

	std::cout << "Image data scaled" << endl;
	
	fprintf(savedata, "%s\t %s\t\t %s\t\t %s\n", "misplacement_X", "SE", "SRE SE", "SE Difference"); 
	fprintf(savedata,"\n");
	std::cout << "Now Computing..." << endl; 

	 //---> major for loop (begin) <---//
	for (int step = 1; step <= 360; step++) {    


	// build theta3 and theta4 (begin)
	double pi = 3.141592;

	for (int k = 1; k <= n8; k++) { // for

	SRE.pointer->theta3[0][n8] = ((double) (2.0 * pi * k) / NEI );  
           
    SRE.pointer->theta4[0][n8] = ((double)  ( ( 2.0 * pi * k) / NEI) * ( ( 2.0 * pi * k) / NEI)  );
	
	} // for
	// build theta3 and theta4 (end)

	// to shift ahead (begin)
	double cos_add = 0.0;

	for (int k = 1; k <= n8; k++) { // for

	cos_add = (double)cos_add + (double) cos( ((double) ( ( 2.0 * pi * k ) / NEI ) * misplacement_X) );

	} // for

	cos_add = ((double) SCALE * cos_add);

	// to shift ahead (end)
	
	 // compute thetas, convolutions, and SRE points (begin)
	 double add1;
	 int k1;
	 double shift_funct;
	 double a, b, c, d, e, f, g, h, j;
	 double a1, b1, c1, d1, e1, f1, g1, h1;
	 double sre_num, sre_den;
	 double ratio, store;

	 double INITGS = (double)misplacement_X; // Initial Guess for the Newton's Iterative Method
	 double convergence;
	 double CONVG = (double)0.001;           // convergence ratio for the Newton's Iterative Method
	 double SMOOTH = (double)0.01;           // convergenge smoothing constraint

	 
     for (int i1=n7; i1 < n1-1-n7; i1++) {// x dim   
	 // n7 here allows changing neighborhood size	
	
		 for (int i2=0; i2 < n2; i2++) { // y dim

			if ( (double)SRE.pointer->the_fMRI[i1][i2] > TH ) { // if major

			add1 = 0;

			k1 = 0;

			for ( int k = -n7; k <= n7; k++ ) {

				 add1 =  (double) add1 + ((double)SRE.pointer->the_fMRI[i1+k][i2] * (*(g_filt + k1))) ;
				

			   k1 = k1 + 1;
            }

			SRE.pointer->theta2[i1][i2] = ((double) ( ( add1 * 2.0) / NEI ));

			shift_funct = (double)0.0;

         
			for ( int s = -n7; s <= n7; s++ ) { 
			// determine the sinc image by convolution also thetas by convolutions (begin)

                      
			shift_funct = (double)shift_funct + ((double)SRE.pointer->the_fMRI[i1+s][i2] * cos_add);
                               
           
			} // determine the sinc image by convolution also thetas by convolutions (end)


			SRE.pointer->shifted_fMRI[i1][i2] = (double)SRE.pointer->fMRI[i1][i2] + 
			
			((double)2.0 * shift_funct / NEI) + ((double)SRE.pointer->theta2[i1][i2] / 2.0);


			// Determine SRE
		    // Now Proceed with Newton's Iterative Method (begin)

			a = (double)0.0;
			b = (double)0.0;
			c = (double)0.0;
			d = (double)0.0;
			e = (double)0.0;
			f = (double)0.0;
			g = (double)0.0;
			h = (double)0.0;
			j = (double)0.0;

			a1 = (double)0.0;
			b1 = (double)0.0;
			c1 = (double)0.0;
			e1 = (double)0.0;
			f1 = (double)0.0;
			g1 = (double)0.0;
			h1 = (double)0.0;

			sre_num = (double)0.0;
			sre_den = (double)0.0;


			/// Initial Guess ///
			SRE.pointer->Xsre[i1][i2] = (double)INITGS * ((double)SRE.pointer->theta2[i1][i2] / 100.0) ;
      
			convergence = 1.0;

			while (convergence >= CONVG) 
			{ // while loop (begin)

				a = ((double)-SRE.pointer->theta2[i1][i2] / 2.0) * ((double)SRE.pointer->theta2[i1][i2]); //

				e = ((double)SRE.pointer->theta2[i1][i2] * (double)SRE.pointer->theta2[i1][i2] / 2.0 );


				for ( int k = 1; k <= n8; k++ ) { // for


				b = (double)b + ((double)SRE.pointer->theta4[0][k] * (double)sin( ((double)( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2])));

				c = (double)c - ((double)SRE.pointer->Xsre[i1][i2] * ((double)SRE.pointer->theta3[0][k] * SRE.pointer->theta4[0][k]) * 
					             (double)cos( ((double)( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2]) ) );

				f = (double)f + ((double)sin( ((double)( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2]) ) * 
					             (double)SRE.pointer->theta3[0][k] * (double)cos( ((double)( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2])));

   		
				g = (double)g + (((double)( (2.0 * pi * k ) * SRE.pointer->Xsre[i1][i2] / NEI )) * (double)SRE.pointer->theta3[0][k]); //


				h = (double)h + ((double)SRE.pointer->Xsre[i1][i2] * SRE.pointer->theta3[0][k] * SRE.pointer->theta3[0][k]);

	
				j = (double)j + ((double)cos( ((double)( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2]) ) * 
					             (double)cos( ((double)( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2]) ) - 
								 (double)sin( ((double)( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2]) ) + 1.0 );


				a1 = (double)a1 - ((double)( SRE.pointer->theta2[i1][i2] / 2.0) * 
					              SRE.pointer->theta2[i1][i2] * SRE.pointer->theta3[0][k] * SRE.pointer->theta3[0][k]);


				b1 = (double)b1 + ((double)SRE.pointer->theta4[0][k] * 
					               (double)cos( ((double)( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2]) ) / (double)SRE.pointer->theta3[0][k]);


				c1 = (double)c1 + ( - (double)SRE.pointer->theta4[0][k] * 
					                  (double)cos( ((double)( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2] ) / (double)SRE.pointer->theta3[0][k] ) + 
									 ((double)( SRE.pointer->Xsre[i1][i2] * 
									  (double)SRE.pointer->theta4[0][k] * 
									  (double)sin( ((double)( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2] ) ) ) ));

	
				e1 = (double)e1 + ((double)SRE.pointer->theta2[i1][i2] * SRE.pointer->theta2[i1][i2] / 2.0 ) * 
					              ((double)SRE.pointer->theta3[0][k] * SRE.pointer->theta3[0][k]);


				f1 = (double)f1 + ((double)cos( ((double)( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2] )) * 
					               (double)cos( ((double)( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2] )) ) - 
								  ((double)sin( ((double)( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2] )) * 
								   (double)sin( ((double)( ( 2.0 * pi * k ) / NEI ) )) + 1.0 ); //

   		
				g1 = (double)g1 + ((double)cos( ((double)( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2] )) * 
					               (double)cos( ((double)( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2] )) ) - 
								  ((double)sin( ((double)( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2] )) + 1.0);

		
				h1 = (double)h1 + SRE.pointer->Xsre[i1][i2] * (double)( - 2.0 * SRE.pointer->theta3[0][k] * 
					                                         ((double)cos( ((double) ( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2] ) ) * 
															  (double)sin( ((double) ( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2] )))- 
															 ((double)SRE.pointer->theta3[0][k] * 
															  (double)cos( ((double) ( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2] ))) );
				} // for

				d = ((double)a * ( b + c ));
		
				sre_num = ((double)( d - ( e * ( f + g + ( h * j ) ) ) ) );

				d1 = ((double)a1 * ( b1 + c1 ));

				sre_den = ((double)( d1 - ( e1 * ( f1 + g1 + h1 ) ) ) );


				if ( (double)sre_num != 0.0 && (double)sre_den != 0.0 ) 
					
					ratio = ((double)sre_num / sre_den);

				else if ( (double)sre_num == 0.0 && (double)sre_den != 0.0 ) 

					ratio = (double)0.0;

				else if ( (double)sre_num != 0.0 && (double)sre_den == 0.0 )

					ratio = (double)0.0;

				else if ( (double)sre_num == 0.0 && (double)sre_den == 0.0 )
					
					ratio = (double)1.0; // de l'Hopital
        

				store = (double)SRE.pointer->Xsre[i1][i2];

				// end of n-th iteration

				SRE.pointer->Xsre[i1][i2] = (double)SRE.pointer->Xsre[i1][i2] - (double)ratio;      
       
				convergence = (double) abs ( (double)store - (double)SRE.pointer->Xsre[i1][i2] );

				if ( (double)convergence >= (double)SMOOTH )  
				{// do not allow abrupt change of the SRE value
        
					SRE.pointer->Xsre[i1][i2] = (double)store; 
            
					break;
				} // do not allow abrupt change of the SRE value


				if ( (double)abs( (double)SRE.pointer->Xsre[i1][i2] ) >= 1.0 ) 
				{ // do not allow the SRE point outside the pixel
           
					break;

        		} // do not allow abrupt change of the SRE value

			} // while loop (end)

			// Determine SRE
		    // Now Proceed with Newton's Iterative Method (end)

			} // if major
			else if ( (double)SRE.pointer->the_fMRI[i1][i2] <= TH ) { // if major

			    SRE.pointer->Xsre[i1][i2] = (double)0.0;
		   
			} // if major

		} // y dim
        
	}  // x dim 
	 // compute thetas, convolutions, and SRE points (end)

	// scale the XSRE to the misplacement (begin)
	max=-MAX;
	min=MAX;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

			if ( (double)SRE.pointer->the_fMRI[i1][i2] > TH ) { // if major
	
		if( SRE.pointer->Xsre[i1][i2] > (double)max ) 
			
			max = (double)SRE.pointer->Xsre[i1][i2];
              
		if( SRE.pointer->Xsre[i1][i2] < (double)min ) 
			
			min = (double)SRE.pointer->Xsre[i1][i2];
		
			} // if major

		} // y dim
        
	}  // x dim 
		
	
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		   if ( (double)SRE.pointer->the_fMRI[i1][i2] > TH ) { // if major

           if ( max == min ) SRE.pointer->Xsre[i1][i2] = (double)0.0;

           else SRE.pointer->Xsre[i1][i2] = ( (double) (misplacement_X) * 
			                                  (double) fabs ( (min - SRE.pointer->Xsre[i1][i2]) / (min - max) ) );
              
		   } else if ( (double)SRE.pointer->the_fMRI[i1][i2] <= TH ) {

			    SRE.pointer->Xsre[i1][i2] = (double)0.0;

		   } // if major

		} // y dim
        
	}  // x dim 


	std::cout << "XSRE Calculated" << endl;
	// scale the XSRE to the misplacement (end)

	 // compute the OMEGA(X_sre-misplacement_X) & SRE.pointer->X_recomputed -> (begin)
	double c_add_x1, c_add_x2, c_add_x3, c_add_x4; 
	double c_add_x5, c_add_x6, c_add_x7, c_add_x8;
	
	double d_add_x1, d_add_x2, d_add_x3, d_add_x4; 
	double d_add_x5, d_add_x6, d_add_x7, d_add_x8;

	double Ein_x_sre, Ein_x_sre_x_mis, SQRT_TERM, condTerm, L;

	 for (int i1=n7; i1 < n1-1-n7; i1++) {// x dim   
	 // n7 here allows changing neighborhood size
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

			if ( (double)SRE.pointer->the_fMRI[i1][i2] > TH ) {  // if major

		// determine Ein coefficients at X_sre  
		c_add_x1 = (double)0.0;
		c_add_x2 = (double)0.0;
		c_add_x3 = (double)0.0;
		c_add_x4 = (double)0.0;
		c_add_x5 = (double)0.0;
		c_add_x6 = (double)0.0;
		c_add_x7 = (double)0.0;
		c_add_x8 = (double)0.0;


		for (int k = 1; k <= n8; k++) { /// for loop

		c_add_x1 =   (double)c_add_x1 + ( SRE.pointer->theta4[0][k] * (double)cos( ((double)( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2] )) );

        c_add_x2 =   (double)c_add_x2 + (double)sin( ((double)( ( 2.0 * pi * k )/ NEI ) * SRE.pointer->Xsre[i1][i2]) ) ;

        c_add_x4 =   (double)c_add_x4 + ( (double)SRE.pointer->theta4[0][k] * (double)sin( ((double)( ( 2.0 * pi *k ) / NEI ) * SRE.pointer->Xsre[i1][i2] )) );

        c_add_x5 =   (double)c_add_x5 + (double)cos( ((double)( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2] ) ) ;

	    c_add_x8 =   (double)c_add_x8 +  (   (double)sin( ((double)( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2]) ) * 
			                                 (double)cos( ((double)( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2]) )   );


		}  /// for loop


		for (int k = 1; k <= n8; k++) { /// for loop

		c_add_x3 = (double)c_add_x3 + (double)c_add_x2 * ( - (double)SRE.pointer->theta3[0][k] * (double)cos( ((double)( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2] )) );

		c_add_x6 = (double)c_add_x6 - (double)c_add_x5 * (   (double)SRE.pointer->theta3[0][k] * (double)sin( ((double)( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->Xsre[i1][i2] )) );

		c_add_x7 = (double)c_add_x7 + ((double) SRE.pointer->theta3[0][k] * c_add_x8);


		}  /// for loop


		c_add_x1 = - ((double)SCALE * c_add_x1);

        c_add_x2 =   ((double)SCALE * c_add_x2);

        c_add_x3 =   ((double)SCALE * c_add_x3); 
          
        c_add_x4 =   ((double)SCALE * c_add_x4); 

        c_add_x5 =   ((double)SCALE * c_add_x5);

	    c_add_x6 =   ((double)SCALE * c_add_x6);

	    c_add_x7 =   ((double)SCALE * c_add_x7);


		// determine Ein coefficients at ( X_sre - misplacement_X )  
        d_add_x1 = (double)0.0;
        d_add_x2 = (double)0.0;
        d_add_x3 = (double)0.0;
        d_add_x4 = (double)0.0;
        d_add_x5 = (double)0.0; 
	    d_add_x6 = (double)0.0;
        d_add_x7 = (double)0.0;
	    d_add_x8 = (double)0.0;


		for (int k = 1; k <= n8; k++) { /// for loop

		d_add_x1 =   (double)d_add_x1 + ( SRE.pointer->theta4[0][k] * (double)cos( ((double)( ( 2.0 * pi * k ) / NEI ) * ( SRE.pointer->Xsre[i1][i2] - misplacement_X ) )) );

        d_add_x2 =   (double)d_add_x2 + (double)sin( ((double)( ( 2.0 * pi * k ) / NEI ) * ( SRE.pointer->Xsre[i1][i2] - misplacement_X ) )) ;

        d_add_x4 =   (double)d_add_x4 + ( SRE.pointer->theta4[0][k] * (double)sin( ((double)( ( 2.0 * pi * k ) / NEI ) * ( SRE.pointer->Xsre[i1][i2] - misplacement_X ) )) );

        d_add_x5 =   (double)d_add_x5 + (double)cos( ((double)( ( 2.0 * pi * k ) / NEI ) * ( SRE.pointer->Xsre[i1][i2] - misplacement_X )) ) ;

		d_add_x8 =   (double)d_add_x8 +  (   (double)sin( ((double)( ( 2.0 * pi * k ) / NEI ) * ( SRE.pointer->Xsre[i1][i2] - misplacement_X )) ) * 
											 (double)cos( ((double)( ( 2.0 * pi * k ) / NEI ) * ( SRE.pointer->Xsre[i1][i2] - misplacement_X )) )  + 
											 ((double)( ( ( 2.0 * pi * k ) / NEI ) * ( SRE.pointer->Xsre[i1][i2] - misplacement_X) ) ) );

		}  /// for loop


		for (int k = 1; k <= n8; k++) { /// for loop

		d_add_x3 = (double)d_add_x3 + (double)d_add_x2 * ( - (double)SRE.pointer->theta3[0][k] * (double)cos( ((double)( ( 2.0 * pi * k ) / NEI ) * ( SRE.pointer->Xsre[i1][i2] - misplacement_X ) )) );

		d_add_x6 = (double)d_add_x6 - (double)d_add_x5 * (   (double)SRE.pointer->theta3[0][k] * (double)sin( ((double)( ( 2.0 * pi * k ) / NEI ) * ( SRE.pointer->Xsre[i1][i2] - misplacement_X ) )) );

		d_add_x7 = (double)d_add_x7 + ((double) SRE.pointer->theta3[0][k] * d_add_x8);

		}  /// for loop


		d_add_x1 = - ((double)SCALE * d_add_x1);

        d_add_x2 =   ((double)SCALE * d_add_x2);

        d_add_x3 =   ((double)SCALE * d_add_x3); 
          
        d_add_x4 =   ((double)SCALE * d_add_x4); 

        d_add_x5 =   ((double)SCALE * d_add_x5);

	    d_add_x6 =   ((double)SCALE * d_add_x6);

	    d_add_x7 =   ((double)SCALE * d_add_x7);


		// determine h and d2_h at xsre       
        double h = 0.0;

		for ( int s = -n7; s <= n7; s++ ) { /// for loop

		// determine the h at sre             
        h = (double)h + ((double)SRE.pointer->the_fMRI[i1+s][i2] * c_add_x5);

		}  /// for loop

		SRE.pointer->h_x_sre[i1][i2] = ((double) ( 2.0 * h / NEI ) + ( SRE.pointer->theta2[i1][i2] / 2.0 )); 
		// determine the h at sre 
        
        SRE.pointer->de2_h_xre[i1][i2] = ((double)SRE.pointer->theta2[i1][i2] * c_add_x1);


		// determine Ein at X_sre
        Ein_x_sre = - ((double) ( SRE.pointer->theta2[i1][i2] / 2.0 ) * SRE.pointer->theta2[i1][i2] * c_add_x4);


        Ein_x_sre = (double)Ein_x_sre + ((double)( SRE.pointer->theta2[i1][i2] * SRE.pointer->theta2[i1][i2] * c_add_x3 ));


        Ein_x_sre = (double)Ein_x_sre - ((double)( SRE.pointer->theta2[i1][i2] * SRE.pointer->theta2[i1][i2] * c_add_x6 ));


		Ein_x_sre = (double)Ein_x_sre - ((double)( SRE.pointer->theta2[i1][i2] * SRE.pointer->theta2[i1][i2] * c_add_x7 ) / 2.0);


		// determine Ein at ( X_sre - misplacement_X )
 	    Ein_x_sre_x_mis = - ((double) ( SRE.pointer->theta2[i1][i2] / 2.0 ) * SRE.pointer->theta2[i1][i2] * d_add_x4);


        Ein_x_sre_x_mis = (double)Ein_x_sre_x_mis + ((double)( SRE.pointer->theta2[i1][i2] * SRE.pointer->theta2[i1][i2] * d_add_x3 ));


        Ein_x_sre_x_mis = (double)Ein_x_sre_x_mis - ((double)( SRE.pointer->theta2[i1][i2] * SRE.pointer->theta2[i1][i2] * d_add_x6 ));

        
	    Ein_x_sre_x_mis = (double)Ein_x_sre_x_mis - ((double)( SRE.pointer->theta2[i1][i2] * SRE.pointer->theta2[i1][i2] * d_add_x7 ) / 2.0);

	
		/// Calculate the Intensity-Curvature Measure (begin) ///
		
		SRE.pointer->ICM[i1][i2] = (double) Ein_x_sre_x_mis / (double) Ein_x_sre;

	    /// Calculate the Intensity-Curvature Measure (end) ///


		// compute the OMEGA(X_sre-misplacement_X)
        if ( (double)Ein_x_sre_x_mis == 0.0 && (double)Ein_x_sre != 0.0 )

        SRE.pointer->Omega[i1][i2] = (double)0.0;

        else if ( (double)Ein_x_sre_x_mis != 0.0 && (double)Ein_x_sre == 0.0 )
	
		SRE.pointer->Omega[i1][i2] = (double)0.0;

        else if ( (double)Ein_x_sre_x_mis == 0.0 && (double)Ein_x_sre == 0.0 )

		SRE.pointer->Omega[i1][i2] = ((double) SRE.pointer->h_x_sre[i1][i2] * SRE.pointer->de2_h_xre[i1][i2]); // de l'Hopital

		else if ( (double)Ein_x_sre_x_mis != 0.0 && (double)Ein_x_sre != 0.0 )

        SRE.pointer->Omega[i1][i2] = ((double)( Ein_x_sre_x_mis / Ein_x_sre ) * SRE.pointer->h_x_sre[i1][i2] * SRE.pointer->de2_h_xre[i1][i2]);

       
		SQRT_TERM = ((double)SRE.pointer->theta2[i1][i2] * SRE.pointer->theta2[i1][i2] * SRE.pointer->theta2[i1][i2] * SRE.pointer->theta2[i1][i2] / 4.0);

		SQRT_TERM = (double)SQRT_TERM * ((double)( ( 2.0 * pi / NEI ) * ( 2.0 * pi / NEI ) * ( 2.0 * pi / NEI ) * ( 2.0 * pi / NEI ) ));  

		SQRT_TERM = (double)SQRT_TERM - ((double)( 4.0 * SRE.pointer->theta2[i1][i2] * SRE.pointer->theta2[i1][i2] * ( 2.0 * pi / NEI ) * ( 2.0 * pi / NEI ) * SRE.pointer->Omega[i1][i2] ));

		condTerm = ((double) ( 2.0 * SRE.pointer->theta2[i1][i2] * SRE.pointer->theta2[i1][i2] * 2.0 * pi * 2.0 * pi / ( NEI * NEI) ) );

		if ( ((double)SQRT_TERM >= 0.0) && ((double)condTerm != 0.0) ) 
		{ // if condition

		if ( PLUS == 1 )

		L = ((double)( - ( 2.0 * SRE.pointer->theta2[i1][i2] * SRE.pointer->theta2[i1][i2] * pi * pi / ( NEI * NEI ) ) + (double)sqrt( (double)SQRT_TERM ) ) ) /  
		    ((double)( 2.0 * SRE.pointer->theta2[i1][i2] * SRE.pointer->theta2[i1][i2] * 2.0 * pi * 2.0 * pi / ( NEI * NEI) ) );

		else if ( PLUS == 0 )

		L = ((double)( - ( 2.0 * SRE.pointer->theta2[i1][i2] * SRE.pointer->theta2[i1][i2] * pi * pi / ( NEI * NEI ) ) - (double)sqrt((double)SQRT_TERM ) ) ) /  
		    ((double)( 2.0 * SRE.pointer->theta2[i1][i2] * SRE.pointer->theta2[i1][i2] * 2.0 * pi * 2.0 * pi / ( NEI * NEI) ) );

		} else if ( ((double)SQRT_TERM < 0.0) && ((double)condTerm != 0.0) ) {	 // if condition
	
		if ( PLUS == 1 )

		L = ((double)( - ( 2.0 * SRE.pointer->theta2[i1][i2] * SRE.pointer->theta2[i1][i2] * pi * pi / ( NEI * NEI ) )  ) ) /  
		    ((double)( 2.0 * SRE.pointer->theta2[i1][i2] * SRE.pointer->theta2[i1][i2] * 2.0 * pi * 2.0 * pi / ( NEI * NEI) ) );

		else if ( PLUS == 0 )

		L = ((double)( - ( 2.0 * SRE.pointer->theta2[i1][i2] * SRE.pointer->theta2[i1][i2] * pi * pi / ( NEI * NEI ) )  ) ) /  
		    ((double)( 2.0 * SRE.pointer->theta2[i1][i2] * SRE.pointer->theta2[i1][i2] * 2.0 * pi * 2.0 * pi / ( NEI * NEI) ) );

		} else if ( ((double)condTerm == 0.0) ) { // if condition

		L = (double)0.0;
    
		} else if ( ((double)( - ( 2.0 * SRE.pointer->theta2[i1][i2] * SRE.pointer->theta2[i1][i2] * pi * pi / ( NEI * NEI ) ) - (double)sqrt((double)SQRT_TERM ) ) ) == 0.0 
			
					&& ((double)condTerm == 0.0) ) { // if condition

		L = (double)1.0;  // de l'Hopital
    
		} // if condition
		
		if ( (double)abs((double)L) < 1.0 )
    
			SRE.pointer->X_recomputed[i1][i2] = ((double)( SRE.pointer->Xsre[i1][i2] - ( 3.0 * (double)acos ( (double)L ) / ( 2.0 * pi ) ) ));  
    
		else
    
			SRE.pointer->X_recomputed[i1][i2] = (double)SRE.pointer->Xsre[i1][i2];		

			} // if major
			else if ( (double)SRE.pointer->the_fMRI[i1][i2] <= TH ) { // if major

			    SRE.pointer->X_recomputed[i1][i2] = (double)0.0;
			} // if major

		} // y dim
        
	}  // x dim 
	// compute the OMEGA(X_sre-misplacement_X) & SRE.pointer->X_recomputed -> (end)


	std::cout << "Re-computed locations Calculated" << endl;


	/// compute max and min values of recomputed locations (begin)
	/// scale the re-computed locations (begin)
	max=-MAX;
	min=MAX;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		if ( (double)SRE.pointer->the_fMRI[i1][i2] > TH ) { 
	

		if( SRE.pointer->X_recomputed[i1][i2] > (double)max ) max = (double)SRE.pointer->X_recomputed[i1][i2];
              
		if( SRE.pointer->X_recomputed[i1][i2] < (double)min ) min = (double)SRE.pointer->X_recomputed[i1][i2];
	
		
		}

		} // y dim
        
	}  // x dim 

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		   if ( (double)SRE.pointer->the_fMRI[i1][i2] > TH ) {
	
           if ( max == min ) SRE.pointer->X_recomputed_scaled[i1][i2] = (double)0.0;

           else  SRE.pointer->X_recomputed_scaled[i1][i2] = ((double) ( RESOL ) * 
			                                                 (double) fabs ( ( min - SRE.pointer->X_recomputed[i1][i2] ) / ( min - max )) );
                    
		   } else if ( (double)SRE.pointer->the_fMRI[i1][i2] <= TH ) {

				 SRE.pointer->X_recomputed_scaled[i1][i2] = (double)0.0;

		   }

		} // y dim
        
	}  // x dim 

	/// compute max and min values of recomputed locations (end)
	/// scale the re-computed locations (end)


	// motion correct the shifted image of the given misplacement,
	// and motion correct the shifted image of the recomputed location,
	// and compute SE (square error) -> (begin)

	// motion correction with recomputed locations (begin)
	double cos_X_SRE_recomp;

	for (int i1=n7; i1 < n1-1-n7; i1++) {// x dim   
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

			if ( (double)SRE.pointer->the_fMRI[i1][i2] > TH ) {  // if major
           
		cos_X_SRE_recomp = (double)0.0;

           
		for ( int k = 1; k <= n8; k++ ) { // for
                    
		cos_X_SRE_recomp = (double)cos_X_SRE_recomp + 
			               (double)cos( ((double)( ( 2.0 * pi * k ) / NEI ) * SRE.pointer->X_recomputed_scaled[i1][i2] ));
   
        } // for
                   
        cos_X_SRE_recomp = ((double) SCALE * cos_X_SRE_recomp);

        shift_funct = 0.0;
           

		for ( int s = -n7; s <= n7; s++ ) 
		{ // determine the sinc image by convolution 
                      
           shift_funct = (double)shift_funct + ((double)SRE.pointer->the_fMRI[i1+s][i2] * cos_X_SRE_recomp) ;

        } // determine the sinc image by convolution 
           
           
		SRE.pointer->the_m_corr_recomp_fMRI[i1][i2] = (double)SRE.pointer->fMRI[i1][i2] + 
		
		((double)2.0 * shift_funct / NEI) + ((double)SRE.pointer->theta2[i1][i2]/2.0);
   
			} else if ( (double)SRE.pointer->the_fMRI[i1][i2] <= TH ) {  // if major

				   SRE.pointer->the_m_corr_recomp_fMRI[i1][i2] = (double)0.0;

			}  // if major


 		} // y dim
        
	}  // x dim 
	// motion correction with recomputed locations (end)


	double sum_SE = 0.0;
	double sum_SE_SRE = 0.0;

	// calculate Square Errors (begin)
	long int counter = 0;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

			if ( (double)SRE.pointer->the_fMRI[i1][i2] > TH ) {
			
			counter++;

			SRE.pointer->SE[i1][i2] = (double) sqrt( (double)( SRE.pointer->fMRI[i1][i2] - SRE.pointer->shifted_fMRI[i1][i2] ) * 
				                                             ( SRE.pointer->fMRI[i1][i2] - SRE.pointer->shifted_fMRI[i1][i2] ) ) ;  

			
			SRE.pointer->SE_SRE[i1][i2] = (double) sqrt( (double)( SRE.pointer->fMRI[i1][i2] - SRE.pointer->the_m_corr_recomp_fMRI[i1][i2] ) * 
				                                                 ( SRE.pointer->fMRI[i1][i2] - SRE.pointer->the_m_corr_recomp_fMRI[i1][i2] ) ) ;  

			sum_SE = ((double) sum_SE + SRE.pointer->SE[i1][i2]);

			sum_SE_SRE = ((double) sum_SE_SRE + SRE.pointer->SE_SRE[i1][i2]);

			} else if ( (double)SRE.pointer->the_fMRI[i1][i2] <= TH ) { }

		} // y dim
        
	}  // x dim 

	std::cout << "SE without SRE is: " << ((double)sum_SE/( counter ))  << endl;
	std::cout << endl;
	std::cout << "SE with SRE is: " << ((double)sum_SE_SRE/( counter ))  << endl;
	std::cout << endl;
	// calculate Square Errors (end)


	// motion correct the shifted image of the given misplacement,
	// and motion correct the shifted image of the recomputed location,
	// and compute SE (square error) -> (end)

	double DIFF = ((double) ((double)sum_SE_SRE/( counter )) - ((double)sum_SE/( counter )) );

	fprintf(savedata, "%e\t %e\t %e\t %e\n", misplacement_X, ((double)sum_SE/( counter )), ((double)sum_SE_SRE/( counter )), DIFF);

	misplacement_X = ((double)misplacement_X + STEPSIZE);
            
	}  //---> major for loop (end) <---//

	fprintf(savedata,"\n");
	fprintf(savedata,"%s\n", "End of Computation...");

	} /// loop across the simulations: sin, cos, log (end)
	
	std::cout << endl;
	std::cout << "End of Computation..." << endl;
	std::cout << endl;

	fclose(savedata);
	free(g_filt);
	delete SRE.pointer;
	SRE.~SINCOSRESIMUL2013();
	
	} // processing (end)

	system("java processCompleted");
	
	return 0;
} // end of main 