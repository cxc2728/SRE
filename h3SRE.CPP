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

//  Project Title: parametric Quadratic B-Spline (parametric Quadratic B-Spline SRE-based Interpolation)

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
#define PLUS 1

#define TH 0

using namespace std;

class h3SRE2013 {

	int n1; // matrix size x
	int n2; // matrix size y

public:

	int getNofPixelsX(void) { return this->n1; };

	int getNofPixelsY(void) { return this->n2; };

	void setNofPixelsX(int x) { this->n1 = x; };

	void setNofPixelsY(int y) { this->n2 = y; };

public:

	struct data {

		double **the_fMRI; // pointer to the matrix entry 

		double **the_m_corr_recomp_fMRI; // pointer to the matrix entry

		double **shifted_fMRI; // pointer to the matrix entry
		
		double **fMRI; // pointer to the matrix entry 

		double **Xsre; // pointer to the matrix entry 

		double **SE; // pointer to the matrix entry 

		double **DELTA; // pointer to the matrix entry 

		double **SE_SRE; // pointer to the matrix entry 

		double **X_recomputed; // pointer to the matrix entry 

		double **X_recomputed_scaled; // pointer to the matrix entry

		double **thetaa0; // pointer to the matrix entry

		double **thetaa12; // pointer to the matrix entry 

		double **thetaa32; // pointer to the matrix entry

		double **h_31_recomp; // pointer to the matrix entry

		double **h_32_recomp; // pointer to the matrix entry

		double **ICF; // pointer to the matrix entry

	}*pointer; // pointer to the matrices

public:

	h3SRE2013(int x, int y) : n1(x), n2(y) { };

	void allocateData();

	void save();

	void OnFourierAnalysis(int m_Xpixels, int m_Ypixels, char imageFilename[]);

	~h3SRE2013() { } // destructor

};

void h3SRE2013::allocateData() { // allocate data

	// (1) allocate struct 'data' (begin)
	 pointer = new data;
	
	 pointer->the_fMRI = new double*[this->n1];

	 pointer->the_m_corr_recomp_fMRI = new double*[this->n1];

	 pointer->shifted_fMRI = new double*[this->n1];

	 pointer->fMRI = new double*[this->n1];

	 pointer->Xsre = new double*[this->n1];

	 pointer->SE = new double*[this->n1];

	 pointer->DELTA = new double*[this->n1];

	 pointer->SE_SRE = new double*[this->n1];

	 pointer->X_recomputed = new double*[this->n1];

	 pointer->X_recomputed_scaled = new double*[this->n1];

     pointer->thetaa0 = new double*[this->n1];

	 pointer->thetaa12 = new double*[this->n1];

	 pointer->thetaa32 = new double*[this->n1];

	 pointer->h_31_recomp = new double*[this->n1];

	 pointer->h_32_recomp = new double*[this->n1];

	 pointer->ICF = new double*[this->n1];


	 for( int v=0; v < this->n1; v++ ) { // (1)

		 pointer->the_fMRI[v] = new double[this->n2];

		 pointer->the_m_corr_recomp_fMRI[v] = new double[this->n2];

		 pointer->shifted_fMRI[v] = new double[this->n2];

		 pointer->fMRI[v] = new double[this->n2];
		 
		 pointer->Xsre[v] = new double[this->n2];

		 pointer->SE[v] = new double[this->n2];

		 pointer->DELTA[v] = new double[this->n2];

		 pointer->SE_SRE[v] = new double[this->n2];

		 pointer->X_recomputed[v] = new double[this->n2];

		 pointer->X_recomputed_scaled[v] = new double[this->n2];

		 pointer->thetaa0[v] = new double[this->n2];

		 pointer->thetaa12[v] = new double[this->n2];

		 pointer->thetaa32[v] = new double[this->n2];

		 pointer->h_31_recomp[v] = new double[this->n2];

		 pointer->h_32_recomp[v] = new double[this->n2];

		 pointer->ICF[v] = new double[this->n2];


	  } // (1) allocate struct 'data' (end)

		// (2) initialize (begin)
		for( int v=0; v < this->n1; v++ ) { // (a)

			for( int f=0; f < this->n2 ; f++ ) { // (b)

			pointer->the_fMRI[v][f] = (double)0.0;

			pointer->the_m_corr_recomp_fMRI[v][f] = (double)0.0;

			pointer->shifted_fMRI[v][f] = (double)0.0;

			pointer->fMRI[v][f] = (double)0.0;
		 
			pointer->Xsre[v][f] = (double)0.0;

			pointer->SE[v][f] = (double)0.0;

			pointer->DELTA[v][f] = (double)0.0;

			pointer->SE_SRE[v][f] = (double)0.0;

			pointer->X_recomputed[v][f] = (double)0.0;

			pointer->X_recomputed_scaled[v][f] = (double)0.0;

			pointer->thetaa0[v][f] = (double)0.0;

			pointer->thetaa12[v][f] = (double)0.0;

			pointer->thetaa32[v][f] = (double)0.0;

			pointer->h_31_recomp[v][f] = (double)0.0;

			pointer->h_32_recomp[v][f] = (double)0.0;

			pointer->ICF[v][f] = (double)0.0;

			 } //(b)

		 } //(a)
		// (2) initialize (end)

} // allocate data

void h3SRE2013::save() { // saveImages

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


	
	sprintf(outputFile, "%s","DELTA.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2 ; f++ ) 
	
		fwrite(&pointer->DELTA[v][f],sizeof(double),1,savedata);

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


	sprintf(outputFile, "%s","ICF.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for( int v=0; v < this->n1; v++ ) { // (a)

		for( int f=0; f < this->n2; f++ ) 
	
		fwrite(&pointer->ICF[v][f],sizeof(double),1,savedata);

	} // (a)

	fclose(savedata);

	} // (save)

} // saveImages


void h3SRE2013::OnFourierAnalysis(int m_Xpixels, int m_Ypixels, char imageFilename[]) 
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

	 
} // end of OnFourierAnalysis



int main ( int argc, char * argv[] ) {

	char outputFile[128]="h3SRE.log";
	
	FILE * savedata;

	double MAX = 5000000000000000000.0;

if (argc < 11) { std::cout << endl;
				std::cout << "Please type the image file name" << endl;
				std::cout << "Please make sure that the image format is Analyze 'double': 64 bits real" << endl;
				std::cout << "Before running the program, please make sure that the image is padded of 'n >= 3'";
				std::cout << "number of pixels along X and Y" << endl;
				std::cout << "Please enter the number of pixels along the X direction (integer)" << endl;
				std::cout << "Please enter the number of pixels along the Y direction (integer)" << endl;
				std::cout << "Please enter the pixel size along the X direction (double)" << endl;
				std::cout << "Please enter the pixel size along the Y direction (double)" << endl;
				std::cout << "Please enter the misplacement along the X direction (double)" << endl;
				std::cout << "Please enter the value of the 'a constant' (double)" << endl;
				std::cout << "Please enter the XY rotation angle (double)" << endl;
				std::cout << "Please type y to run the Spectral Power Analysis (y (yes) or n (no))" << endl;
				std::cout << "Please enter that 'n >= 3' number of pixels along X and Y which will pad the image" << endl;
				std::cout << "Some data will be lost if you enter a value of n which is greater than" << endl;
				std::cout << "the actual number of pixels along X and Y which was used to pad the image" << endl;
				std::cout << "before running the program." << endl;
				std::cout << "Please make sure that the value of the 'a constant' is reasonably" << endl;
				std::cout << "not too high neither is too low, otherwise the interpolation error" << endl;
				std::cout << "is not meaningful" << endl;
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
	double the_A_const = atof(argv[7]);

	double theta = atof(argv[8]);

	char imageFileName[128];

	unsigned char Fyn = 'n';

	Fyn = *argv[9];

	int PAD = atoi(argv[10]);

	if ( PAD < 3 ) { // if test PAD
	
		std::cout << "Please make sure to pad the image with 'n >= 3' number of pixels along X and Y" << endl;
	
		fprintf(savedata,"%s\n", "Please make sure to pad the image with 'n >= 3' number of pixels along X and Y");
	
		fclose(savedata);
	
		exit(0);
	
	} // if test PAD

	
	int XNEI = (int)3;

	sprintf(imageFileName, "%s", argv[1]);

	std::cout << endl;
	std::cout << "The image file name is: " << imageFileName << endl;
	std::cout << "The number of pixels along the X direction is: " << atoi(argv[2]) << endl;
	std::cout << "The number of pixels along the Y direction is: " << atoi(argv[3]) << endl;
	std::cout << "The pixel size along the X direction is: " << atof(argv[4]) << endl;
	std::cout << "The pixel size along the Y direction is: " << atof(argv[5]) << endl;
	std::cout << "The misplacement along the X direction is: " << atof(argv[6]) << endl;
	std::cout << "The value of the 'a constant' is: " << atof(argv[7]) << endl;
	std::cout << "The XY rotation angle is: " << atof(argv[8]) << endl;
	
	if ( tolower(Fyn) == 'y' )
	std::cout << "The program runs the Spectral Power Analysis: " << Fyn << endl;
	else if ( tolower(Fyn) != 'y' ) std::cout << "The program does not run the Spectral Power Analysis: " << Fyn << endl;
	std::cout << "The 'n' number of pixels along X and Y which you entered to pad the image is: " << atoi(argv[10]) << endl;
	std::cout << "Some data will be lost if you enter a value of n which is greater than" << endl;
	std::cout << "the actual number of pixels along X and Y which was used to pad the image" << endl;
	std::cout << "before running the program." << endl;
	std::cout << endl;


	fprintf(savedata,"%s%s\n", "The image file name is: " , imageFileName);
	fprintf(savedata,"%s%d\n", "The number of pixels along the X direction is: ", n1);
	fprintf(savedata,"%s%d\n", "The number of pixels along the Y direction is: ", n2);
	fprintf(savedata,"%s%lf\n", "The misplacement along the X direction is: ", x_misplacement_X);
	fprintf(savedata,"%s%lf\n", "The pixel size along the X direction is: ", XPixelSize);
	fprintf(savedata,"%s%lf\n", "The pixel size along the Y direction is: ", YPixelSize);
	fprintf(savedata,"%s%lf\n", "The value of the 'a constant' is: ", the_A_const);
	fprintf(savedata,"%s%lf\n", "The XY rotation angle is: ", theta);
	if ( tolower(Fyn) == 'y' )
	fprintf(savedata,"%s%c\n", "The program runs the Spectral Power Analysis: ", Fyn);
	else if ( tolower(Fyn) != 'y' )
	fprintf(savedata,"%s%c\n", "The program does not run the Spectral Power Analysis: ", Fyn);
	fprintf(savedata,"%s%d\n", "The 'n' number of pixels along X and Y which you entered to pad the image is: ", PAD);
	fprintf(savedata,"%s\n", "Some data will be lost if you enter a value of n which is greater than");
	fprintf(savedata,"%s\n", "the actual number of pixels along X and Y which was used to pad the image");
	fprintf(savedata,"%s\n", "before running the program.");
	fprintf(savedata,"%s\n", "Please make sure that the value of the 'a constant' is reasonably");
	fprintf(savedata,"%s\n", "not too high neither is too low, otherwise the interpolation error is not meaningful.");
	fprintf(savedata,"\n");


	double misplacement_X = ((double)1.0 - ( cos( (double)theta ) + sin( (double)theta ) ) + x_misplacement_X);
 
	misplacement_X = ((double)misplacement_X/XPixelSize);

	   //////////////////***********//////////////////////
	  // Above formula scales the misplacement to the   //
	  // pixel size the same way the following formula  //
	  // would do: (min - misplacement)/(min - max)     //  
	  //////////////////***********//////////////////////

	h3SRE2013 SRE(n1,n2);

	SRE.allocateData();

	std::cout << "Object Constructed" << endl;

	std::cout << "X image size: " << SRE.getNofPixelsX() << endl;
	std::cout << "Y image size: " << SRE.getNofPixelsY() << endl;


	/// read image file (begin)
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

	// finds index of first and last horizonatal lines of the signal (begin)
	long int firstLine = -1;
	long int lastLine = -1;
	
	//first horizonatal line
	  for (int i1=0; i1 < n1; i1++) {// x dim

	   for (int i2=0; i2 < n2; i2++) { // y dim
		 
		   if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) {
		   
			   firstLine = (int)i1;
			   break;
			 			         				        
		   }  
			  
		} // y dim
    
		if ( firstLine != (int)-1) break;
	}  // x dim
	
	std::cout << "first Line of Signal = " << firstLine << endl;

	//last horizonatal line
	   for (int i1=n1-1; i1 >= 0; i1--) {// x dim

		   for (int i2=0; i2 < n2; i2++) { // y dim

		   if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) {
		   
			   lastLine = (int)i1;
			   break;
					   
		   } 	      
			 
		} // y dim
         
		if ( lastLine != (int)-1) break;
	}  // x dim

	std::cout << "last Line of Signal = " << lastLine << endl;
	// finds index of first and last horizonatal lines of the signal (end)
	
	
	// compute thetas (begin)
	double add_up, add_up_1, add_up_2;

	for (int i1=PAD+1; i1 < n1-PAD-1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) { 
	
		if ( i1 == (int)firstLine || i1 == ((int)firstLine + 1) || i1 == ((int)firstLine + 2) )
		{ // calculate the first three horizonatal lines of the signal (begin)

		add_up = ((double) ( SRE.pointer->the_fMRI[i1-3+XNEI][i2] + SRE.pointer->the_fMRI[i1-2+XNEI][i2] +  
			                 SRE.pointer->the_fMRI[i1-1+XNEI][i2] + SRE.pointer->the_fMRI[i1+3+XNEI][i2] +  
							 SRE.pointer->the_fMRI[i1+2+XNEI][i2] + SRE.pointer->the_fMRI[i1+1+XNEI][i2] ));
     
        add_up_1 = ((double) SRE.pointer->the_fMRI[i1-1+XNEI][i2] + SRE.pointer->the_fMRI[i1+1+XNEI][i2] ); 
        
		add_up_2 = ((double) SRE.pointer->the_fMRI[i1-3+XNEI][i2] + SRE.pointer->the_fMRI[i1-2+XNEI][i2] +  
							 SRE.pointer->the_fMRI[i1+3+XNEI][i2] + SRE.pointer->the_fMRI[i1+2+XNEI][i2] ); 

		// calculate the first three horizonatal lines of the signal (end)
		} else 

		if ( i1 == ((int)lastLine) ||  i1 == ((int)lastLine - 1) || i1 == ((int)lastLine - 2) )
		{ // calculate the last three horizonatal lines of the signal (begin)

		add_up = ((double) ( SRE.pointer->the_fMRI[i1-3-XNEI][i2] + SRE.pointer->the_fMRI[i1-2-XNEI][i2] +  
			                 SRE.pointer->the_fMRI[i1-1-XNEI][i2] + SRE.pointer->the_fMRI[i1+3-XNEI][i2] +  
							 SRE.pointer->the_fMRI[i1+2-XNEI][i2] + SRE.pointer->the_fMRI[i1+1-XNEI][i2] )); 
     
        add_up_1 = ((double) SRE.pointer->the_fMRI[i1-1-XNEI][i2] + SRE.pointer->the_fMRI[i1+1-XNEI][i2] ); 
        
		add_up_2 = ((double) SRE.pointer->the_fMRI[i1-3-XNEI][i2] + SRE.pointer->the_fMRI[i1-2-XNEI][i2] +  
							 SRE.pointer->the_fMRI[i1+3-XNEI][i2] + SRE.pointer->the_fMRI[i1+2-XNEI][i2] ); 

		// calculate the last three horizonatal lines of the signal (end)
		} else { // calculate the rest of the lines of the signal (begin)

 		add_up = ((double) ( SRE.pointer->the_fMRI[i1-3][i2] + SRE.pointer->the_fMRI[i1-2][i2] +  
			                 SRE.pointer->the_fMRI[i1-1][i2] + SRE.pointer->the_fMRI[i1+3][i2] +  
							 SRE.pointer->the_fMRI[i1+2][i2] + SRE.pointer->the_fMRI[i1+1][i2] ));
     
        add_up_1 = ((double) SRE.pointer->the_fMRI[i1-1][i2] + SRE.pointer->the_fMRI[i1+1][i2] ); 
        
		add_up_2 = ((double) SRE.pointer->the_fMRI[i1-3][i2] + SRE.pointer->the_fMRI[i1-2][i2] +  
							 SRE.pointer->the_fMRI[i1+3][i2] + SRE.pointer->the_fMRI[i1+2][i2] ); 

		 // calculate the rest of the lines of the signal (end)
		}


		SRE.pointer->thetaa0[i1][i2] = ((double)  ( - 2.0 * the_A_const * SRE.pointer->the_fMRI[i1][i2] * add_up_1 ) + 
			                                      ( 2.0 * the_A_const * SRE.pointer->the_fMRI[i1][i2] * add_up_2 ));
                

		SRE.pointer->thetaa12[i1][i2] = (double)  add_up_1 * ( (double) ( -2.0 * the_A_const * add_up_1 ) +  
			                                                             ( 2.0 * the_A_const * add_up_2 ) );

        
		SRE.pointer->thetaa32[i1][i2] = (double)  add_up *   ( (double) ( -2.0 * the_A_const * add_up_1 ) +  
			                                                             ( 2.0 * the_A_const * add_up_2 ) );

		} else if ( (double)SRE.pointer->fMRI[i1][i2] <= TH ) {

		SRE.pointer->thetaa0[i1][i2] =  (double)0.0;        

		SRE.pointer->thetaa12[i1][i2] = (double)0.0;

		SRE.pointer->thetaa32[i1][i2] = (double)0.0;
		
		}

		} // y dim
        
	}  // x dim 
	// compute thetas (end)

	// compute the h functions (begin)
	double h_31 = (double) SCALE * (double) exp ( - (double)( -2.0 * the_A_const * misplacement_X * misplacement_X + 0.5 * (1.0 + the_A_const ) ) );

	double h_32 = (double) SCALE * (double) exp ( - (double)( the_A_const * misplacement_X * misplacement_X - ( ( 2.0 * the_A_const + 0.5 ) * misplacement_X ) + ( (3.0/4.0) *( 1.0 + the_A_const) ) ) ) ;
	// compute the h functions (end)

	std::cout << "Thetas Calculated" << endl;

	// compute XSRE (begin)
	double a, b;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) { 

         
		a = ((double) SRE.pointer->thetaa32[i1][i2] * ( 2.0 * the_A_const + 0.5 ) / 2.0 );
			
		b = ((double) ( 2.0 * SRE.pointer->thetaa32[i1][i2] * the_A_const / 3.0 ) - ( 2.0 * SRE.pointer->thetaa12[i1][i2] * the_A_const * 2.0 / 3.0 ));
        

		if ( (double)a == 0.0 && (double)b == 0.0 )
		{

		 SRE.pointer->Xsre[i1][i2] = (double)1.0; // de L'Hopital
		
		}
		else if ( (double)a == 0.0 && (double)b != 0.0 )
		{ 
			
         SRE.pointer->Xsre[i1][i2] = (double)0.0;
		 
		} 
		else if ( (double)a != 0.0 && (double)b == 0.0 )
		{ 
			
         SRE.pointer->Xsre[i1][i2] = (double)0.0;
		 
		} else { // else

         SRE.pointer->Xsre[i1][i2] = ((double)a/b); 
		       											  
		} // else
               
  	
		} else if ( (double)SRE.pointer->fMRI[i1][i2] <= TH ) {

		 SRE.pointer->Xsre[i1][i2] = (double)0.0;

		}

		} // y dim
        
	}  // x dim 
	// compute XSRE (end)

	
	// scale the XSRE to the misplacement (begin)
	max=-MAX;
	min=MAX;

 	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) {
	
		if( SRE.pointer->Xsre[i1][i2] > (double)max ) 
			
			max = (double)SRE.pointer->Xsre[i1][i2];
              
		if( SRE.pointer->Xsre[i1][i2] < (double)min ) 
			
			min = (double)SRE.pointer->Xsre[i1][i2];

		}
		
		} // y dim
        
	}  // x dim 

	
  	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		   if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) { 

           if ( max == min ) SRE.pointer->Xsre[i1][i2] = (double)1.0; // de l'Hopital

           else SRE.pointer->Xsre[i1][i2] = ( (double) (misplacement_X) * (max - SRE.pointer->Xsre[i1][i2]) / (max - min) );

		   } else if ( (double)SRE.pointer->fMRI[i1][i2] <= TH ) {

			    SRE.pointer->Xsre[i1][i2] = (double)0.0;
		   }
           
		} // y dim
        
	}  // x dim 
	
	std::cout << "XSRE Calculated" << endl;
	// scale the XSRE to the misplacement (end)


	// compute DELTA (begin)
	double q, k, s, m, f, n;
	double num_DELTA, den_DELTA;
	double q1, k1, s1;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) {  


		q = (double) ( SRE.pointer->Xsre[i1][i2] - misplacement_X ) * ( SRE.pointer->Xsre[i1][i2] - misplacement_X ) * ( SRE.pointer->Xsre[i1][i2] - misplacement_X );
       
        k = (double) ( SRE.pointer->Xsre[i1][i2] - misplacement_X ) * ( SRE.pointer->Xsre[i1][i2] - misplacement_X );

        s = (double) ( SRE.pointer->Xsre[i1][i2] - misplacement_X ) ;


        m = (double) ( SRE.pointer->thetaa32[i1][i2] * the_A_const / 3.0 ) - ( SRE.pointer->thetaa12[i1][i2] * 2.0 * the_A_const / 3.0 ); 
 
        f = (double) ( SRE.pointer->thetaa32[i1][i2] * ( 2.0 * the_A_const + 0.5 ) / 2.0 ) ; 
       
        n = (double) ( SRE.pointer->thetaa0[i1][i2] + ( SRE.pointer->thetaa12[i1][i2] * ( the_A_const + 1.0 ) / 2.0 ) +  ( SRE.pointer->thetaa32[i1][i2] * 3.0 * ( the_A_const + 1.0 ) / 4.0 ) );
           
		 
		num_DELTA = ((double) ( q * m ) - ( k * f ) + ( s * n ));


        q1 = ((double) SRE.pointer->Xsre[i1][i2] *  SRE.pointer->Xsre[i1][i2] * SRE.pointer->Xsre[i1][i2] ) ;
       
        k1 = ((double) SRE.pointer->Xsre[i1][i2] *  SRE.pointer->Xsre[i1][i2] ) ;

        s1 = ((double) SRE.pointer->Xsre[i1][i2] ) ;


        den_DELTA = ((double) ( q1 * m ) - ( k1 * f ) + ( s1 * n ));


		if ( (double)num_DELTA == 0.0 && (double)den_DELTA == 0.0 ) 
			
			SRE.pointer->DELTA[i1][i2] = (double)1.0; // de L'Hopital
		
		else if ( (double)num_DELTA != 0.0 && den_DELTA == 0.0 ) 
			
			SRE.pointer->DELTA[i1][i2] = (double)0.0; 

		else if ( (double)num_DELTA == 0.0 && den_DELTA != 0.0 ) 
			
			SRE.pointer->DELTA[i1][i2] = (double)0.0; 

		else SRE.pointer->DELTA[i1][i2] = ((double) num_DELTA / den_DELTA);


		} else if ( (double)SRE.pointer->fMRI[i1][i2] <= TH ) {

			   SRE.pointer->DELTA[i1][i2] = (double)0.0;

		}

		} // y dim
        
	}  // x dim 
	// compute DELTA (end)

	std::cout << "DELTA Calculated" << endl;	

	/// calculation of the Intensity-Curvature Functional (begin) ///
	double num_ICF;
	double den_ICF;
    for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) {  


        m = (double) ( SRE.pointer->thetaa32[i1][i2] * the_A_const / 3.0 ) - ( SRE.pointer->thetaa12[i1][i2] * 2.0 * the_A_const / 3.0 ); 
 
        f = (double) ( SRE.pointer->thetaa32[i1][i2] * ( 2.0 * the_A_const + 0.5 ) / 2.0 ) ; 
       
        n = (double) ( SRE.pointer->thetaa0[i1][i2] + ( SRE.pointer->thetaa12[i1][i2] * ( the_A_const + 1.0 ) / 2.0 ) +  ( SRE.pointer->thetaa32[i1][i2] * 3.0 * ( the_A_const + 1.0 ) / 4.0 ) );
           
		 
		num_ICF = ((double) SRE.pointer->thetaa0[i1][i2] * misplacement_X);

		den_ICF =  (double) m * misplacement_X * misplacement_X * misplacement_X;
		den_ICF -= (double) f * misplacement_X * misplacement_X;
		den_ICF +=  (double) n * misplacement_X;

		if ( (double)num_ICF == 0.0 && (double)den_ICF == 0.0 ) 
			
			SRE.pointer->ICF[i1][i2] = (double)1.0; // de L'Hopital
		
		else if ( (double)num_ICF != 0.0 && den_ICF == 0.0 ) 
			
			SRE.pointer->ICF[i1][i2] = (double)0.0; 

		else if ( (double)num_ICF == 0.0 && den_ICF != 0.0 ) 
			
			SRE.pointer->ICF[i1][i2] = (double)0.0; 

		else SRE.pointer->ICF[i1][i2] = ((double) num_ICF / den_ICF);


		} else if ( (double)SRE.pointer->fMRI[i1][i2] <= TH ) {

			   SRE.pointer->ICF[i1][i2] = (double)0.0;

		}

		} // y dim
        
	}  // x dim 
    /// calculation of the Intensity-Curvature Functional (end) ///

	// calculate the re-computed locations (begin)
	double ro, f1, f2, r;
	double hsre, p, SQRT_term;	
	
	for (int i1=PAD+1; i1 < n1-PAD-1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) { 

 	    if ( i1 == (int)firstLine || i1 == ((int)firstLine + 1) || i1 == ((int)firstLine + 2) )
		{ // calculate the first three horizonatal lines of the signal (begin)

		f1 = ((double) SRE.pointer->the_fMRI[i1-1+XNEI][i2] + SRE.pointer->the_fMRI[i1+1+XNEI][i2] );  

        f2 = ((double) SRE.pointer->the_fMRI[i1-3+XNEI][i2] + SRE.pointer->the_fMRI[i1-2+XNEI][i2] +  
		               SRE.pointer->the_fMRI[i1-1+XNEI][i2] + SRE.pointer->the_fMRI[i1+3+XNEI][i2] +  
					   SRE.pointer->the_fMRI[i1+2+XNEI][i2] + SRE.pointer->the_fMRI[i1+1+XNEI][i2] );    	

		// calculate the first three horizonatal lines of the signal (end)
		} else 

		if ( i1 == ((int)lastLine) ||  i1 == ((int)lastLine - 1) || i1 == ((int)lastLine - 2) )
		{ // calculate the last three horizonatal lines of the signal (begin)

		f1 = ((double) SRE.pointer->the_fMRI[i1-1-XNEI][i2] + SRE.pointer->the_fMRI[i1+1-XNEI][i2] );  

        f2 = ((double) SRE.pointer->the_fMRI[i1-3-XNEI][i2] + SRE.pointer->the_fMRI[i1-2-XNEI][i2] +  
		               SRE.pointer->the_fMRI[i1-1-XNEI][i2] + SRE.pointer->the_fMRI[i1+3-XNEI][i2] +  
					   SRE.pointer->the_fMRI[i1+2-XNEI][i2] + SRE.pointer->the_fMRI[i1+1-XNEI][i2] );    	

		// calculate the last three horizonatal lines of the signal (end)
		} else { // calculate the rest of the lines of the signal (begin)
 	
	  	f1 = ((double) SRE.pointer->the_fMRI[i1-1][i2] + SRE.pointer->the_fMRI[i1+1][i2] );  

        f2 = ((double) SRE.pointer->the_fMRI[i1-3][i2] + SRE.pointer->the_fMRI[i1-2][i2] +  
		               SRE.pointer->the_fMRI[i1-1][i2] + SRE.pointer->the_fMRI[i1+3][i2] +  
					   SRE.pointer->the_fMRI[i1+2][i2] + SRE.pointer->the_fMRI[i1+1][i2] );    	

		// calculate the rest of the lines of the signal (end)
		}


		a = ((double) ( f2 * ( 2.0 * the_A_const + 0.5) ) * ( f2 * ( 2.0 * the_A_const + 0.5) ));

		r =  (double) SRE.pointer->DELTA[i1][i2];

        s =  ((double) SRE.pointer->the_fMRI[i1][i2] + f1 * (  ( -2.0 * the_A_const * SRE.pointer->Xsre[i1][i2] * SRE.pointer->Xsre[i1][i2] ) + ( 0.5 * ( the_A_const + 1.0 ) ) ));
       
		m =  ((double) f2 * ( ( the_A_const * SRE.pointer->Xsre[i1][i2] * SRE.pointer->Xsre[i1][i2] ) - ( ( 2.0 * the_A_const + 0.5 ) * SRE.pointer->Xsre[i1][i2] ) + ( (3.0/4.0) * ( the_A_const + 1.0) )  ));
    
        hsre = ((double)s + m);

		p = ((double) SRE.pointer->the_fMRI[i1][i2] + ( f1 * 0.5 * ( the_A_const + 1.0 ) ) + ( f2 * (3.0/4.0) * ( the_A_const + 1.0 ) ));

		b = ((double) ( 4.0 * ( - 2.0 * f1 + f2 ) * the_A_const ) * ( ( r * hsre ) - p ));

        SQRT_term = ((double)a + b);

		if ( (SQRT_term) >= 0.0 && ( 2.0 * the_A_const * ( -2.0 * f1 + f2) ) != 0.0 ) { // if (1)

			if ( PLUS == 1 ) 
			ro = ((double) ( f2 * ( 2.0 * the_A_const + 0.5 ) + (double)sqrt( (double)SQRT_term ) ) / ( 2.0 * the_A_const * ( -2.0 * f1 + f2) ));
				
			if ( PLUS == 0 ) 
			ro = ((double) ( f2 * ( 2.0 * the_A_const + 0.5 ) - (double)sqrt( (double)SQRT_term ) ) / ( 2.0 * the_A_const * ( -2.0 * f1 + f2) ));
			
		} //if (1) 
		else if ( (SQRT_term) < 0.0 && ( 2.0 * the_A_const * ( -2.0 * f1 + f2) ) != 0.0 ) { // else (1)

			if ( PLUS == 1 || PLUS == 0 ) 
			ro = ((double) ( f2 * ( 2.0 * the_A_const + 0.5 ) ) / ( 2.0 * the_A_const * ( -2.0 * f1 + f2) ));
				
		} //else (1)
		else if ( ( ((double) ( f2 * ( 2.0 * the_A_const + 0.5 ) + (double)sqrt( (double)SQRT_term ) )) != 0.0 ||
			
				  ((double) ( f2 * ( 2.0 * the_A_const + 0.5 ) - (double)sqrt( (double)SQRT_term ) )) != 0.0 ) &&
			
				  ( 2.0 * the_A_const * ( -2.0 * f1 + f2) ) == 0.0 ) { // else (1)
		
			if ( PLUS == 1 || PLUS == 0 ) 
			ro = (double)0.0;
		
		} 
		else if ( ( ((double) ( f2 * ( 2.0 * the_A_const + 0.5 ) + (double)sqrt( (double)SQRT_term ) )) == 0.0 ||
			
				    ((double) ( f2 * ( 2.0 * the_A_const + 0.5 ) - (double)sqrt( (double)SQRT_term ) )) == 0.0 ) &&
			
				    ( 2.0 * the_A_const * ( -2.0 * f1 + f2) ) == 0.0 ) 
		{ // de L'Hopital
		
			if ( PLUS == 1 || PLUS == 0 ) 
			ro = (double)1.0;
		
		} 		
		//else (1)
		
		SRE.pointer->X_recomputed[i1][i2] = (double)ro;


		} else if ( (double)SRE.pointer->fMRI[i1][i2] <= TH ) {

			   SRE.pointer->X_recomputed[i1][i2] = (double)0.0;

		}


		} // y dim
        
	}  // x dim 
	// calculate the re-computed locations (end)
	
	std::cout << "Re-computed locations Calculated" << endl;
	
	/// compute max and min values of recomputed locations (begin)
	/// scale the re-computed locations (begin)
	max=-MAX;
	min=MAX;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) {

	
		if( SRE.pointer->X_recomputed[i1][i2] > (double)max ) max = (double)SRE.pointer->X_recomputed[i1][i2];
              
		if( SRE.pointer->X_recomputed[i1][i2] < (double)min ) min = (double)SRE.pointer->X_recomputed[i1][i2];
		

		} else if ( (double)SRE.pointer->fMRI[i1][i2] <= TH ) {

			   SRE.pointer->X_recomputed[i1][i2] = (double)0.0;

		}

		} // y dim
        
	}  // x dim 

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) {

           if ( max == min ) SRE.pointer->X_recomputed_scaled[i1][i2] = (double)1.0; // de L'Hopital

           else  SRE.pointer->X_recomputed_scaled[i1][i2] = ((double) ( XPixelSize ) * 
			                                                 (double) fabs( ( min - SRE.pointer->X_recomputed[i1][i2] ) / ( min - max )) );
           
		   
		} else if ( (double)SRE.pointer->fMRI[i1][i2] <= TH ) {

			   SRE.pointer->X_recomputed_scaled[i1][i2] = (double)0.0;

		}

		} // y dim
        
	}  // x dim 


	/// compute max and min values of recomputed locations (end)
	/// scale the re-computed locations (end)

	// compute the h functions (begin)
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

			if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) {

			SRE.pointer->h_31_recomp[i1][i2] = (double) SCALE * (double) exp ( -(double) ( -2.0 * the_A_const * SRE.pointer->X_recomputed_scaled[i1][i2] * SRE.pointer->X_recomputed_scaled[i1][i2] + 0.5 * (1.0 + the_A_const ) )  );


			SRE.pointer->h_32_recomp[i1][i2] = (double) SCALE * (double) exp ( -(double) ( the_A_const * SRE.pointer->X_recomputed_scaled[i1][i2] * SRE.pointer->X_recomputed_scaled[i1][i2] - ( ( 2.0 * the_A_const + 0.5 ) * SRE.pointer->X_recomputed_scaled[i1][i2] ) + ( (3.0/4.0) *( 1.0 + the_A_const) ) ) );

			} else if ( (double)SRE.pointer->fMRI[i1][i2] <= TH ) {


			SRE.pointer->h_31_recomp[i1][i2] = (double) 0.0;


			SRE.pointer->h_32_recomp[i1][i2] = (double) 0.0;


			}

		} // y dim
        
	}  // x dim 
	// compute the h functions (end)

	std::cout << "H functions Calculated" << endl;

	// shift the image of the given misplacement (begin)
	double sum;

	for (int i1=PAD+1; i1 < n1-PAD-1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) {
           
 
	    if ( i1 == (int)firstLine || i1 == ((int)firstLine + 1) || i1 == ((int)firstLine + 2) )
		{ // calculate the first three horizonatal lines of the signal (begin)

        sum = ((double) SRE.pointer->the_fMRI[i1-1+XNEI][i2] + SRE.pointer->the_fMRI[i1+1+XNEI][i2] );     
        

	    add_up = ((double) SRE.pointer->the_fMRI[i1-3+XNEI][i2] + SRE.pointer->the_fMRI[i1-2+XNEI][i2] +  
		                   SRE.pointer->the_fMRI[i1-1+XNEI][i2] + SRE.pointer->the_fMRI[i1+3+XNEI][i2] +  
			    		   SRE.pointer->the_fMRI[i1+2+XNEI][i2] + SRE.pointer->the_fMRI[i1+1+XNEI][i2] );    	

		// calculate the first three horizonatal lines of the signal (end)
		} else 

		if ( i1 == ((int)lastLine) ||  i1 == ((int)lastLine - 1) || i1 == ((int)lastLine - 2) )
		{ // calculate the last three horizonatal lines of the signal (begin)

        sum = ((double) SRE.pointer->the_fMRI[i1-1-XNEI][i2] + SRE.pointer->the_fMRI[i1+1-XNEI][i2] );     
      
        add_up = ((double) SRE.pointer->the_fMRI[i1-3-XNEI][i2] + SRE.pointer->the_fMRI[i1-2-XNEI][i2] +  
		                   SRE.pointer->the_fMRI[i1-1-XNEI][i2] + SRE.pointer->the_fMRI[i1+3-XNEI][i2] +  
					       SRE.pointer->the_fMRI[i1+2-XNEI][i2] + SRE.pointer->the_fMRI[i1+1-XNEI][i2] );    	

		// calculate the last three horizonatal lines of the signal (end)
		} else { // calculate the rest of the lines of the signal (begin)
 
	    sum = ((double) SRE.pointer->the_fMRI[i1-1][i2] + SRE.pointer->the_fMRI[i1+1][i2] );    
       
        add_up = ((double) SRE.pointer->the_fMRI[i1-3][i2] + SRE.pointer->the_fMRI[i1-2][i2] +  
		                   SRE.pointer->the_fMRI[i1-1][i2] + SRE.pointer->the_fMRI[i1+3][i2] +  
				           SRE.pointer->the_fMRI[i1+2][i2] + SRE.pointer->the_fMRI[i1+1][i2] );    	

		// calculate the rest of the lines of the signal (end)
		}

			SRE.pointer->shifted_fMRI[i1][i2] = ((double) SRE.pointer->fMRI[i1][i2] + ( sum * h_31 ) + ( add_up * h_32 ) );

		} else if ( (double)SRE.pointer->fMRI[i1][i2] <= TH ) {

			   SRE.pointer->shifted_fMRI[i1][i2] = (double)0.0;

		}

 		} // y dim
        
	}  // x dim 
	// shift the image of the given misplacement (end)

	std::cout << "The Image is Shifted" << endl;

	// shift the image of the the recomputed location (begin)
	for (int i1=PAD+1; i1 < n1-PAD-1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) {

		if ( i1 == (int)firstLine || i1 == ((int)firstLine + 1) || i1 == ((int)firstLine + 2) )
		{ // calculate the first three horizonatal lines of the signal (begin)

        sum = ((double) SRE.pointer->the_fMRI[i1-1+XNEI][i2] + SRE.pointer->the_fMRI[i1+1+XNEI][i2] );     
        

	    add_up = ((double) SRE.pointer->the_fMRI[i1-3+XNEI][i2] + SRE.pointer->the_fMRI[i1-2+XNEI][i2] +  
		                   SRE.pointer->the_fMRI[i1-1+XNEI][i2] + SRE.pointer->the_fMRI[i1+3+XNEI][i2] +  
			    		   SRE.pointer->the_fMRI[i1+2+XNEI][i2] + SRE.pointer->the_fMRI[i1+1+XNEI][i2] );    	

		// calculate the first three horizonatal lines of the signal (end)
		} else 

		if ( i1 == ((int)lastLine) ||  i1 == ((int)lastLine - 1) || i1 == ((int)lastLine - 2) )
		{ // calculate the last three horizonatal lines of the signal (begin)

        sum = ((double) SRE.pointer->the_fMRI[i1-1-XNEI][i2] + SRE.pointer->the_fMRI[i1+1-XNEI][i2] );     
      
        add_up = ((double) SRE.pointer->the_fMRI[i1-3-XNEI][i2] + SRE.pointer->the_fMRI[i1-2-XNEI][i2] +  
		                   SRE.pointer->the_fMRI[i1-1-XNEI][i2] + SRE.pointer->the_fMRI[i1+3-XNEI][i2] +  
					       SRE.pointer->the_fMRI[i1+2-XNEI][i2] + SRE.pointer->the_fMRI[i1+1-XNEI][i2] );    	

		// calculate the last three horizonatal lines of the signal (end)
		} else { // calculate the rest of the lines of the signal (begin)
 
	    sum = ((double) SRE.pointer->the_fMRI[i1-1][i2] + SRE.pointer->the_fMRI[i1+1][i2] );    
       
        add_up = ((double) SRE.pointer->the_fMRI[i1-3][i2] + SRE.pointer->the_fMRI[i1-2][i2] +  
		                   SRE.pointer->the_fMRI[i1-1][i2] + SRE.pointer->the_fMRI[i1+3][i2] +  
				           SRE.pointer->the_fMRI[i1+2][i2] + SRE.pointer->the_fMRI[i1+1][i2] );    	

		// calculate the rest of the lines of the signal (end)
		}

		SRE.pointer->the_m_corr_recomp_fMRI[i1][i2] = ((double) SRE.pointer->fMRI[i1][i2] + ( sum * SRE.pointer->h_31_recomp[i1][i2] ) + ( add_up * SRE.pointer->h_32_recomp[i1][i2] ) );
		
		} else if ( (double)SRE.pointer->fMRI[i1][i2] <= TH ) {

  			   SRE.pointer->the_m_corr_recomp_fMRI[i1][i2] = (double)0.0;
		
		}

			} // y dim
        
	}  // x dim 

	std::cout << "Image SRE shifted" << endl;
    std::cout << endl;
	// shift the image of the the recomputed location (end)

	// calculate Square Errors (begin)
	double sum_SE = 0.0;
	double sum_SE_SRE = 0.0;
	long int counter = 0;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

			if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) {

			counter++;

		    SRE.pointer->SE[i1][i2] = (double) sqrt( (double)( SRE.pointer->fMRI[i1][i2] - SRE.pointer->shifted_fMRI[i1][i2] ) * 
				                                             ( SRE.pointer->fMRI[i1][i2] - SRE.pointer->shifted_fMRI[i1][i2] ) ) ;  

			
			SRE.pointer->SE_SRE[i1][i2] = (double) sqrt( (double)( SRE.pointer->fMRI[i1][i2] - SRE.pointer->the_m_corr_recomp_fMRI[i1][i2] ) * 
				                                                 ( SRE.pointer->fMRI[i1][i2] - SRE.pointer->the_m_corr_recomp_fMRI[i1][i2] ) ) ;  


			sum_SE = ((double) sum_SE + SRE.pointer->SE[i1][i2]);

			sum_SE_SRE = ((double) sum_SE_SRE + SRE.pointer->SE_SRE[i1][i2]);

			} else if ( (double)SRE.pointer->fMRI[i1][i2] <= TH ) { }

		} // y dim
        
	}  // x dim 

	// ----->>>print and save Square Errors (begin)
	std::cout << "SE without SRE is: " << ((double)sum_SE/( (counter) ))  << endl;
	std::cout << endl;

	fprintf(savedata,"%s%e\n", "SE without SRE is: ", ((double)sum_SE/( (counter) )) );
	fprintf(savedata,"\n");

	std::cout << "SE with SRE is: " << ((double)sum_SE_SRE/( counter ))  << endl;
	std::cout << endl;

	fprintf(savedata,"%s%e\n", "SE with SRE is: ", ((double)sum_SE_SRE/( counter )) );
	fprintf(savedata,"\n");
	//----->>> print and save Square Errors (end)
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
	delete SRE.pointer;
	SRE.~h3SRE2013();
	} // processing (end)

	} // run the program (end)

	system("java processCompleted");

	return 0;
} // end of main 