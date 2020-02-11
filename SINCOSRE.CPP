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

#define INTERVAL 2000
#define PLUS 0         // One is positive solution for the novel recomputed locations, Zero is set for negative solutions 
#define INFERIOR 500
#define SUPERIOR 1500

#define TH 0

using namespace std;


class SINCOSRE2013 {

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

	SINCOSRE2013(int x, int y, int n) : n1(x), n2(y), n8(n) { };

	void allocateData();

	void save();

	void OnFourierAnalysis(int m_Xpixels, int m_Ypixels, char imageFilename[]);

	~SINCOSRE2013() { } // destructor

};

void SINCOSRE2013::allocateData() { // allocate data

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

void SINCOSRE2013::save() { // saveImages

	FILE * savedata;
	char outputFile[128];
	
	sprintf(outputFile, "%s","fMRI.img");

	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2 ; f++ ) 
	
		fwrite(&pointer->fMRI[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)


	
	sprintf(outputFile, "%s","Xsre.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2 ; f++ ) 
	
		fwrite(&pointer->Xsre[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)
	
	
	sprintf(outputFile, "%s","shifted_fMRI.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2 ; f++ ) 
			
		fwrite(&pointer->shifted_fMRI[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)



	sprintf(outputFile, "%s","SE.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2 ; f++ ) 
	
		fwrite(&pointer->SE[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)


	
	sprintf(outputFile, "%s","h_x_sre.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2 ; f++ ) 
	
		fwrite(&pointer->h_x_sre[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	

	sprintf(outputFile, "%s","de2_h_xre.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2 ; f++ ) 
	
		fwrite(&pointer->de2_h_xre[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)
	

	sprintf(outputFile, "%s","the_m_corr_recomp_fMRI.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2 ; f++ ) 
	
		fwrite(&pointer->the_m_corr_recomp_fMRI[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)


	
	sprintf(outputFile, "%s","SE_SRE.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2 ; f++ ) 
	
		fwrite(&pointer->SE_SRE[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)


	sprintf(outputFile, "%s","X_recomputed.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2 ; f++ ) 
	
		fwrite(&pointer->X_recomputed[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)



	sprintf(outputFile, "%s","X_recomputed_scaled.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2 ; f++ ) 
	
		fwrite(&pointer->X_recomputed_scaled[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)


	sprintf(outputFile, "%s","Omega.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2 ; f++ ) 
	
		fwrite(&pointer->Omega[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

	
	sprintf(outputFile, "%s","ICM.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->ICM[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

} // saveImages


void SINCOSRE2013::OnFourierAnalysis(int m_Xpixels, int m_Ypixels, char imageFilename[]) 
{
	
	int NofXpixels = m_Xpixels, NofYpixels = m_Ypixels;

	double *kSpaceR = 0, *kSpaceI = 0, *Signal = 0, *magimageFT = 0;
	double *bin_value = 0, *histogram_FT = 0, *spectr_energy = 0;

	int i, j, index;
	int dx, dy;
	int ds, dp; 
	int k2, k3, w, t;
	
	double pi = 3.141592;
	double MAX = 5000000000000000000.0;

	FILE * logfile;
	char logfilename[128];

	sprintf(logfilename, "%s%s%s", "Fourier-T-" , imageFilename, ".log");
	
	if ((logfile = fopen(logfilename,"w+"))==NULL)
	{

		std::cout << "Cannot open Fourier Analysis log file, Now Exit..." << endl;
	
		// clear memory here not necessary here
  	    // because this if statement takes place 
	    // after blocks of memory have not been allocated. 

		exit(0);
	
	} else { // allocate memory (begin) 

	if ((kSpaceR = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");
		exit(0);

	}

	if ((kSpaceI = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");

		// FIFO memory deallocation method
		free(kSpaceR);
		exit(0);

	}

	if ((Signal = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{
	
		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");
	
		// FIFO memory deallocation method
 	   free(kSpaceR);
	   free(kSpaceI);

	   exit(0);

	}
	
	if (( magimageFT = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{

		fprintf(logfile,"%s\n", "Not enough memory to allocate Real Image data: Exit");
	
		// FIFO memory deallocation method
 	   free(kSpaceR);
	   free(kSpaceI);
	   free(Signal);

	   exit(0);

	}

	} // allocate memory (end)

	std::cout << "Direct Fourier Transform Processing Started" << endl;
	std::cout << endl;
    fprintf(logfile,"%s\n", "Direct Fourier Transform Processing Started");

	//// read image data and initialize pointers
	FILE *image;

	char m_ImageFile[128];
	
	sprintf(m_ImageFile, "%s", imageFilename);

	std::cout << "File to process: " << m_ImageFile << endl;
	std::cout << "Pixels X Size: " << NofXpixels << endl;   
	std::cout << "Pixels Y Size: " << NofYpixels << endl;

    if ((image = fopen(m_ImageFile,"rb+"))==NULL)
	{
	
	 fprintf(logfile, "%s%s\n", "Cannot open Image File: ", imageFilename);

	 // FIFO memory deallocation method
 	 free(kSpaceR);
     free(kSpaceI);
	 free(Signal);

	 exit(0);

	} else { // read data and initialize pointers

		double number = 0.0;

		for (i=0; i<NofXpixels; i++)
		{ 
			for (j=0; j<NofYpixels; j++)
			{

				index = ((j*NofXpixels)+i);

				fread(&number,sizeof(double),1,image);

				*(Signal+index)= (double)number;
				
				*(kSpaceR+index) = (double) 0.0;

				*(kSpaceI+index) = (double) 0.0;

				*(magimageFT+index) = (double) 0.0;

			}

		}


		fclose(image);

	}// read data and initialize pointers


	std::cout << "Data Initialized" << endl;
	std::cout << "Fourier Transformation..." << endl;
	std::cout << "Please wait..." << endl;

	double phase, complexR, complexI;
	
	///// Fourier Transform  (begin)//////
	///calculate k-space data (begin)
	for (i=0; i<NofXpixels; i++)
	{ 

		for (j=0; j<NofYpixels; j++)
		{

	
			dx = ((int) i - NofXpixels/2);
		    dy = ((int) j - NofYpixels/2);

	
			k2 = ((int)(dx*NofYpixels)+dy);

			w = ((j*NofXpixels)+i);


			for (int s=0; s<NofXpixels; s++)
			{ ///calculate k-space data 
				for (int p=0; p<NofYpixels; p++)
				{ 
					

		     		ds = ((int) s - NofXpixels/2);
		            dp = ((int) p - NofYpixels/2);
			
					k3 = ((int)(dp*NofXpixels)+ds);  

					t = ((p*NofXpixels)+s);

				
					phase = ((double) 2.0 * pi * k2 * k3 / (NofXpixels*NofYpixels) );

					complexR = (double) cos( (double)phase ) + (double) sin( (double)phase ); 

                    complexI = -(double) sin( (double)phase ) + (double) cos( (double)phase ); 

				
					*(kSpaceR+w) += (double) *(Signal+t) * (double) complexR;

					*(kSpaceI+w) += (double) *(Signal+t) * (double) complexI;

					
			}

		}///calculate k-space data 


		}
	} ///calculate k-space data (end)
	///// Fourier Transform  (end)//////

	std::cout << "End of Fourier Transformation" << endl;
	std::cout << "Now Saving k-space data" << endl;

	double savedata = 0.0;
	FILE * pf;
	char filename[128];

	sprintf(filename, "%s%s", "K-SpaceR-", m_ImageFile);

    fprintf(logfile, "%s\t%s\n", "Now Saving K-Space Signal (Real) in File: ", filename);

    if ((pf = fopen(filename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save K-Space Signal");

	 // FIFO memory deallocation method
 	 free(kSpaceR);
 	 free(kSpaceI);
	 free(Signal);
	 free(magimageFT);

	 exit(0);
	
	} else { // save data


	for (i=0; i<NofXpixels; i++)
	{ ///save k-space data
		for (j=0; j<NofYpixels; j++)
		{

			index = ((j*NofXpixels)+i);

			savedata = (double)*(kSpaceR+index);
          
            fwrite(&savedata,sizeof(double),1,pf);

		}
	} ///save k-space data

	fprintf(logfile,"%s\n", "K-Space Signal (Real) Saved");

	fclose (pf);
	} // save data


	sprintf(filename, "%s%s", "K-SpaceI-", m_ImageFile);

    fprintf(logfile, "%s\t%s\n", "Now Saving K-Space Signal (Imaginary) in File: ", filename);

    if ((pf = fopen(filename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save K-Space Signal");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(Signal);
	 free(magimageFT);

	 exit(0);
	
	} else { // save data


	for (i=0; i<NofXpixels; i++)
	{ ///save k-space data
		for (j=0; j<NofYpixels; j++)
		{

			index = ((j*NofXpixels)+i);

			savedata = (double)*(kSpaceI+index);
          
            fwrite(&savedata,sizeof(double),1,pf);

		}
	} ///save k-space data

	fprintf(logfile,"%s\n", "K-Space Signal (Imaginary) Saved");

	fclose (pf);
	} // save data

	std::cout << "Direct Fourier Transform Processing Completed" << endl;
	std::cout << endl;
    fprintf(logfile,"%s\n", "Direct Fourier Transform Processing Completed");
	

	std::cout << "Start of the Frequency Spectrum Analysis" << endl;
	std::cout << endl;
    fprintf(logfile, "%s\n", "Start of the Frequency Spectrum Analysis");

	// build the frequency magnitudes array 'magimageFT' (begin)
	int signReal = 0;
	int signImag = 0;
	double magnitude;
	
	for (i=0; i<NofXpixels; i++)
	{ 
		for (j=0; j<NofYpixels; j++)
		{

			index = ((j*NofXpixels)+i);
			
			if ( (double)*(kSpaceR+index) >= 0.0 ) signReal = 1;
			else if ( (double)*(kSpaceR+index) < 0.0 ) signReal = -1;

			if ( (double)*(kSpaceI+index) >= 0.0 ) signImag = 1;
			else if ( (double)*(kSpaceI+index) < 0.0 ) signImag = -1;

			magnitude = (double) sqrt( (double)*(kSpaceR+index) * (double)*(kSpaceR+index) +
									   (double)*(kSpaceI+index) * (double)*(kSpaceI+index) );

			*(magimageFT+index) = ((double) (signReal / signImag) * magnitude );

		}
	} 
	// build the frequency magnitudes array 'magimageFT' (end)


	// Frequency Spectral Analysis (begin)
	// save data (begin)
	sprintf(filename, "%s%s", "FT_image-", m_ImageFile);

    fprintf(logfile, "%s\t%s\n", "Now Saving FT Image in File: ", filename);

    if ((pf = fopen(filename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save FT Image");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(Signal);
	 free(magimageFT);

	 exit(0);
	
	} else { // else save data


	for (i=0; i<NofXpixels; i++)
	{ ///save data
		for (j=0; j<NofYpixels; j++)
		{

			index = ((j*NofXpixels)+i);

			savedata = (double)*(magimageFT+index);
          
            fwrite(&savedata,sizeof(double),1,pf);

		}
	} ///save data

	fprintf(logfile,"%s\n", "FT Image Saved");

	fclose (pf);
	} // else save data
	// save data (end)

	double max=-MAX;
	double min=MAX;

	for (i=0; i<NofXpixels; i++)
	{ // (a)
		for (j=0; j<NofYpixels; j++)
		{ // (b)

			index = ((j*NofXpixels)+i);
	
			if( (double)*(magimageFT+index) > (double)max ) max = (double)*(magimageFT+index);
              
			if( (double)*(magimageFT+index) < (double)min ) min = (double)*(magimageFT+index);
		
		} // (b)
        
	}  // (a)

	/// build the Histogram (begin)
	if ((bin_value = (double *) calloc( INTERVAL, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Histogram data: Exit");

		// FIFO memory deallocation method
		free(kSpaceR);
		free(kSpaceI);
		free(Signal);
		free(magimageFT);
		exit(0);

	} else {

		for (int i=0; i<INTERVAL; i++) 

		*(bin_value+i) = ((double) (max - min) * i / INTERVAL);

	} 

	if ((histogram_FT = (double *) calloc( INTERVAL, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Histogram data: Exit");

		// FIFO memory deallocation method
		free(kSpaceR);
		free(kSpaceI);
		free(Signal);
		free(magimageFT);
		free(bin_value);

		exit(0);

	}

	if ((spectr_energy = (double *) calloc( INTERVAL, sizeof(double)) ) == NULL)
	{
   
		fprintf(logfile,"%s\n", "Not enough memory to allocate Histogram data: Exit");

		// FIFO memory deallocation method
		free(kSpaceR);
		free(kSpaceI);
		free(Signal);
		free(magimageFT);
		free(bin_value);
		free(histogram_FT);

		exit(0);

	} else { // else

		for ( int i = 0; i < INTERVAL; i++ )
		{ // for

			*(histogram_FT + i) = (double)0.0;
			*(spectr_energy + i) = (double)0.0;

		} // for


	} // else 

	int k;
	for (i=0; i<NofXpixels; i++)
	{ // (a)
		for (j=0; j<NofYpixels; j++)
		{ // (b)

		   index = ((j*NofXpixels)+i);
	
		   k = 1;
		   while ( k < INTERVAL ) { // while

		   if ( (double)*(magimageFT + index) >= (double)*(bin_value + k-1) && 
			   
			    (double)*(magimageFT + index) < (double)*(bin_value + k) 
			  ){ // if 

			   *(histogram_FT + k) += (double)1.0;
			   *(spectr_energy + k) += ( (double)*(magimageFT + index) * (double)*(magimageFT + index) );

			   break;
			} // if

		   k = k+1;
		   } // while

		} // (b)
        
	}  // (a)

	/// build the Histogram (end)

	// save histogram data in text file (begin)
	FILE *HistogramPointer;

	char HistogramData[200];

	sprintf(filename, "%s%s", "FT_image-", m_ImageFile);
	sprintf(HistogramData, "%s%s", filename, "-Histogram.txt");

	if ((HistogramPointer = fopen(HistogramData,"w"))==NULL)
	{
	 
	fprintf(logfile,"%s%s%s\n", "Cannot Open File: ", HistogramData, " Now Exit. ");
    
	exit(0);
	
	}else{ // save (begin)

	
	for(int i=0; i < INTERVAL; i++)

	fprintf( HistogramPointer, "%20.25lf\n", (double) *(spectr_energy + i) );
	
	fclose(HistogramPointer);

	fprintf(logfile,"%s%s\n", "Data saved into file: ", HistogramData);
    
	}/// save (end)
	// save histogram data in text file (end)

	std::cout << "Frequency Spectrum Analysis Completed" << endl;
	std::cout << endl;
    fprintf(logfile, "%s\n", "Frequency Spectrum Analysis Completed");


	/// Overlay routine (begin)
	double * Overlay = 0;

	if (( Overlay = (double *) calloc( NofXpixels*NofYpixels, sizeof(double)) ) == NULL)
	{ // allocate memory (begin)

		fprintf(logfile,"%s\n", "Not enough memory to allocate Overlay Image data: Exit");
	
		// FIFO memory deallocation method
 		free(kSpaceR);
		free(kSpaceI);
		free(Signal);
		free(magimageFT);
		free(bin_value);
		free(histogram_FT);
		free(spectr_energy);

		exit(0);

	}// allocate memory (end)
	else { // processing (begin)

	///  find max of Signal (begin)
	max=-MAX;
	min=MAX;

	for (i=0; i<NofXpixels; i++)
	{ // (a)
		for (j=0; j<NofYpixels; j++)
		{ // (b)

			index = ((j*NofXpixels)+i);
	
			if( (double)*(Signal+index) > (double)max ) max = (double)*(Signal+index);
              
			if( (double)*(Signal+index) < (double)min ) min = (double)*(Signal+index);
		
		} // (b)
        
	}  // (a)
	

	for (i=0; i<NofXpixels; i++)
	{ ///store data
		for (j=0; j<NofYpixels; j++)
		{

			index = ((j*NofXpixels)+i);

			*(Overlay+index) = (double)*(Signal+index);

			if ( (double)*(magimageFT+index) >= *(bin_value + INFERIOR) &&

				 (double)*(magimageFT+index) < *(bin_value + SUPERIOR)) 
			{ // if true then overlay (begin)
				
				(double)*(Overlay+index) = (double)max;
				
			}// if true then overlay (end)
        }
	} /// store data
	
	// save data (begin)
	sprintf(filename, "%s%s", "Overlay_image-", m_ImageFile);

    fprintf(logfile, "%s\t%s\n", "Now Saving Overlay Image in File: ", filename);

    if ((pf = fopen(filename,"wb+"))==NULL)
	{

	 fprintf(logfile, "%s\n", "Cannot open file to save Overlay Image");

	 // FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(Signal);
	 free(magimageFT);
	 free(bin_value);
	 free(histogram_FT);
	 free(spectr_energy);
	 free(Overlay);

	 exit(0);
	
	} else { // else save data

	for (i=0; i<NofXpixels; i++)
	{ ///save data
		for (j=0; j<NofYpixels; j++)
		{

			index = ((j*NofXpixels)+i);

			savedata = (double)*(Overlay+index);
          
            fwrite(&savedata,sizeof(double),1,pf);

		}
	} ///save data

	fprintf(logfile,"%s\n", "Overlay Image Saved");

	fclose (pf);
	} // else save data
	// save data (end)
	
	} // processing (end)
	//// Overlay routine (end)

	// Frequency Spectral Analysis (end)

	fclose(logfile);
	system( logfilename );

	// FIFO memory deallocation method
	 free(kSpaceR);
	 free(kSpaceI);
	 free(Signal);
	 free(magimageFT);
	 free(bin_value);
	 free(histogram_FT);
	 free(spectr_energy);
	 free(Overlay);
	
} // end of OnFourierAnalysis


int main ( int argc, char * argv[] ) {

	char outputFile[128]="SINCOSRE.log";

	FILE * savedata;

	double MAX = 5000000000000000000.0;

if (argc < 9) { std::cout << endl;
				std::cout << "Please type the image file name" << endl;
				std::cout << "Please make sure that the image format is Analyze 'double': 64 bits real" << endl;
				std::cout << "Please enter the number of pixels along the X direction (integer)" << endl;
				std::cout << "Please enter the number of pixels along the Y direction (integer)" << endl;
				std::cout << "Please enter the pixel size along the X direction (double)" << endl;
				std::cout << "Please enter the pixel size along the Y direction (double)" << endl;
				std::cout << "Please enter the misplacement along the X direction (double)" << endl;
				std::cout << "Please enter the XY rotation angle (double)" << endl;
				std::cout << "Please type y to run the Spectral Power Analysis (y (yes) or n (no))" << endl;
				std::cout << endl;
				exit(0); }

else { // run the program (begin)

	
	if ((savedata = fopen(outputFile,"w"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // processing (begin)

	int n1 = atoi(argv[2]);
	int n2 = atoi(argv[3]);

	double XPixelSize = atof(argv[4]);
	double YPixelSize = atof(argv[5]);

	double x_misplacement_X = atof(argv[6]);

	double atheta = atof(argv[7]);
	
	char imageFileName[128];

	unsigned char Fyn = 'n';

	Fyn = *argv[8];

	sprintf(imageFileName, "%s", argv[1]);

	std::cout << endl;
	std::cout << "The image file name is: " << imageFileName << endl;
	std::cout << "The number of pixels along the X direction is: " << atoi(argv[2]) << endl;
	std::cout << "The number of pixels along the Y direction is: " << atoi(argv[3]) << endl;
	std::cout << "The pixel size along the X direction is: " << atof(argv[4]) << endl;
	std::cout << "The pixel size along the Y direction is: " << atof(argv[5]) << endl;
	std::cout << "The misplacement along the X direction is: " << atof(argv[6]) << endl;
	std::cout << "The XY rotation angle is: " << atof(argv[7]) << endl;
	
	if ( tolower(Fyn) == 'y' )
	std::cout << "The program runs the Spectral Power Analysis: " << Fyn << endl;
	else if ( tolower(Fyn) != 'y' ) std::cout << "The program does not run the Spectral Power Analysis: " << Fyn << endl;
	std::cout << endl;


	fprintf(savedata,"%s%s\n", "The image file name is: " , imageFileName);
	fprintf(savedata,"%s%d\n", "The number of pixels along the X direction is: ", n1);
	fprintf(savedata,"%s%d\n", "The number of pixels along the Y direction is: ", n2);
    fprintf(savedata,"%s%lf\n", "The pixel size along the X direction is: ", XPixelSize);
	fprintf(savedata,"%s%lf\n", "The pixel size along the Y direction is: ", YPixelSize);
	fprintf(savedata,"%s%lf\n", "The misplacement along the X direction is: ", x_misplacement_X);
	fprintf(savedata,"%s%lf\n", "The XY rotation angle is: ", atheta);

	if ( tolower(Fyn) == 'y' )
	fprintf(savedata,"%s%c\n", "The program runs the Spectral Power Analysis: ", Fyn);
	else if ( tolower(Fyn) != 'y' )
	fprintf(savedata,"%s%c\n", "The program does not run the Spectral Power Analysis: ", Fyn);
	fprintf(savedata,"\n");

	double misplacement_X = ((double)1.0 - ( cos( (double)atheta ) + sin( (double)atheta ) ) + x_misplacement_X);
 
	misplacement_X = ((double)misplacement_X/XPixelSize);

	  ///////////////////***********//////////////////////
	  // Above formula scales the misplacement to the   //
	  // pixel size the same way the following formula  //
	  // would do: (min - misplacement)/(min - max)     //  
	  //////////////////***********///////////////////////
	
	double NEI = 3.0;     // neighborhood
	double FW = 0.87;     // FWHM

	int n7 = (int)floor((double)NEI/2.0);  // number of right/left pixels to include into the filter

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

	SINCOSRE2013 SRE(n1,n2,n8);

	SRE.allocateData();

	std::cout << "Object Constructed" << endl;

	std::cout << "X image size: " << SRE.getNofPixelsX() << endl;
	std::cout << "Y image size: " << SRE.getNofPixelsY() << endl;

	FILE * pf;

	if ((pf = fopen(imageFileName,"rb+"))==NULL)
	{

		std::cout << "Cannot open file: " << imageFileName << endl;
		fprintf(savedata,"%s%s\n", "Cannot open file: " , imageFileName );
		exit(0);

	} else { // else

	double number;

	// initialize data (begin)
	for (int i1=0; i1 < SRE.getNofPixelsX(); i1++) {// x dim
       	
		for (int i2=0; i2 < SRE.getNofPixelsY(); i2++) { // y dim
			
		fread(&number,sizeof(double),1,pf);

		SRE.pointer->fMRI[i1][i2] = (double)number;

                          
		} // y dim
        
	}  // x dim 

	fclose (pf);

	} // else 
	/// read image file (end)

	for (int i1=0; i1 < SRE.getNofPixelsX(); i1++) {// x dim
       	
		for (int i2=0; i2 < SRE.getNofPixelsY(); i2++) { // y dim

		SRE.pointer->the_fMRI[i1][i2] = (double)SRE.pointer->fMRI[i1][i2];

		} // y dim
        
	}  // x dim 
	// initialize data (end)
     	
	
	 	std::cout << "Data read in and Initialized" << endl;


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

           else  SRE.pointer->X_recomputed_scaled[i1][i2] = ((double) ( XPixelSize ) * 
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

	fprintf(savedata,"%s%e\n", "SE without SRE is: ", ((double)sum_SE/( counter )) );
	fprintf(savedata,"\n");

	std::cout << "SE with SRE is: " << ((double)sum_SE_SRE/( counter ))  << endl;
	std::cout << endl;

	fprintf(savedata,"%s%e\n", "SE with SRE is: ", ((double)sum_SE_SRE/( counter )) );
	fprintf(savedata,"\n");
	// calculate Square Errors (end)


	// motion correct the shifted image of the given misplacement,
	// and motion correct the shifted image of the recomputed location,
	// and compute SE (square error) -> (end)


	SRE.save(); // save all image data 

	std::cout << "End of Computation..." << endl;
	std::cout << endl;

	fprintf(savedata,"%s\n", "End of Computation...");
	fprintf(savedata,"\n");

	if ( tolower(Fyn) == 'y' ) {  // run Power Spectra Analysis (begin)
	
	/// FT analysis (begin)
	std::cout << "Now Performing FT Analysis and Saving..." << endl;
	std::cout << endl;

	fprintf(savedata,"%s\n", "Now Performing FT Analysis and Saving...");
	fprintf(savedata,"\n");

	SRE.OnFourierAnalysis(n1, n2, "fMRI.img"); 
	// FT Analysis of the original image: "fMRI.img"

	SRE.OnFourierAnalysis(n1, n2, "the_m_corr_recomp_fMRI.img");
	// FT Analysis of the SRE IMAGE: "the_m_corr_recomp_fMRI.img"

	SRE.OnFourierAnalysis(n1, n2, "shifted_fMRI.img");
	// FT Analysis of the NO SRE IMAGE: "shifted_fMRI.img"
	
	} // run Power Spectra Analysis (end)
	/// FT analysis (end)

	fclose(savedata);
	free(g_filt);
	delete SRE.pointer;
	SRE.~SINCOSRE2013();

	
	} // processing (end)

	} // run the program (end)

	system("java processCompleted");

	return 0;
} // end of main 