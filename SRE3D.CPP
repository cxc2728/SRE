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

// Project Title: Trivariate Linear (Trivariate Linear SRE-based Interpolation)


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


#define NSCALE 0.01
#define REVERTEX 1      // use of positive or negative solutions of the Novel Recomputed Locations

#define TH 0

#define INTERVAL 2000	// size of interval of the frequency spectrum analyzer
#define INFERIOR 500
#define SUPERIOR 1500
#define ZPLUS 0 		// 1 indicates positive solutions, else wise 0 indicates negative solutions
#define YPLUS 0         // 1 indicates positive solutions, else wise 0 indicates negative solutions


using namespace std;


class SRE3D2013 {

	int n1; // matrix size x
	int n2; // matrix size y
	int n4; // matrix size z

public:

	int getNofPixelsX(void) { return this->n1; };

	int getNofPixelsY(void) { return this->n2; };

	int getNofPixelsZ(void) { return this->n4; }

	void setNofPixelsX(int x) { this->n1 = x; };

	void setNofPixelsY(int y) { this->n2 = y; };

	void setNofPixelsZ(int z) { this->n4 = z; };

public:

	struct data {

		double *fMRI; // pointer to the matrix entry 

		double *scaledfMRI; // pointer to the matrix entry
		
		double *SRE_Xcoord; // pointer to the matrix entry 
	
		double *SRE_Ycoord; // pointer to the matrix entry 

		double *SRE_Zcoord; // pointer to the matrix entry 

		double *shifted_fMRI; // pointer to the matrix entry 
	
		double *SE; // pointer to the matrix entry 

		double *DELTA; // pointer to the matrix entry 

		double *the_SRE_shifted_fMRI; // pointer to the matrix entry

		double *SE_SRE; // pointer to the matrix entry 

		double *X_recomputed; // pointer to the matrix entry 

		double *Y_recomputed; // pointer to the matrix entry

		double *Z_recomputed; // pointer to the matrix entry

		double *X_recomputed_scaled; // pointer to the matrix entry 

		double *Y_recomputed_scaled; // pointer to the matrix entry

		double *Z_recomputed_scaled; // pointer to the matrix entry


		double *thetaXYZ; // pointer to the matrix entry

		double *F_XYZ; // pointer to the matrix entry

		double *fxyz; // pointer to the matrix entry

		double *F_XY; // pointer to the matrix entry

		double *F_YZ; // pointer to the matrix entry

		double *F_XZ; // pointer to the matrix entry
	
		double *deDELTAdeX; // pointer to the matrix entry

		double *deDELTAdeY; // pointer to the matrix entry

		double *deDELTAdeZ; // pointer to the matrix entry


		double *Recomp_Mx_plus; // pointer to the matrix entry

		double *Recomp_Mx_minus; // pointer to the matrix entry


		double *Recomp_My_plus; // pointer to the matrix entry

		double *Recomp_My_minus; // pointer to the matrix entry


		double *Recomp_Mz_plus; // pointer to the matrix entry

		double *Recomp_Mz_minus; // pointer to the matrix entry

		double *ICF; // pointer to the matrix entry

	}*pointer; // pointer to the matrices


public:

	SRE3D2013(int x, int y, int z) : n1(x), n2(y), n4(z) { };// constructor 
	
	void allocateData();

	void OnFourierAnalysis(int m_Xpixels, int m_Ypixels, char imageFilename[]);

	void save();

	~SRE3D2013() { } // destructor

};

void SRE3D2013::allocateData() { // allocate data


	 // (1) allocate struct 'data' (begin)
	 pointer = new data;

	 int size = ((int)this->n1 * (int)this->n2 * (int)this->n4);
			
	 pointer->fMRI = new double[size];

	 pointer->scaledfMRI = new double[size];

	 pointer->SRE_Xcoord = new double[size];

	 pointer->SRE_Ycoord = new double[size];

	 pointer->SRE_Zcoord = new double[size];

	 pointer->shifted_fMRI = new double[size];

	 pointer->SE = new double[size];

	 pointer->DELTA = new double[size];

	 pointer->the_SRE_shifted_fMRI = new double[size];

	 pointer->SE_SRE = new double[size];

	 pointer->X_recomputed = new double[size];

	 pointer->Y_recomputed = new double[size];

	 pointer->Z_recomputed = new double[size];

	 pointer->X_recomputed_scaled = new double[size];

	 pointer->Y_recomputed_scaled = new double[size];

	 pointer->Z_recomputed_scaled = new double[size];


	 pointer->thetaXYZ = new double[size]; 

	 pointer->F_XYZ = new double[size]; 

	 pointer->fxyz = new double[size];

	 pointer->F_XY = new double[size];

	 pointer->F_YZ = new double[size];

	 pointer->F_XZ = new double[size];
	
	 pointer->deDELTAdeX = new double[size];

	 pointer->deDELTAdeY = new double[size]; 

	 pointer->deDELTAdeZ = new double[size]; 


	 pointer->Recomp_Mx_plus = new double[size]; 

	 pointer->Recomp_Mx_minus = new double[size];


	 pointer->Recomp_My_plus = new double[size];

	 pointer->Recomp_My_minus = new double[size];


	 pointer->Recomp_Mz_plus = new double[size];

	 pointer->Recomp_Mz_minus = new double[size];

	 pointer->ICF = new double[size];


		// (2) initialize (begin)
		for( int v=0; v < size; v++ ) { // (a)
		 
			pointer->fMRI[v] = (double)0.0;

			pointer->scaledfMRI[v] = (double)0.0;

			pointer->SRE_Xcoord[v] = (double)0.0;

			pointer->SRE_Ycoord[v] = (double)0.0;

			pointer->SRE_Zcoord[v] = (double)0.0;

			pointer->shifted_fMRI[v] = (double)0.0;

			pointer->SE[v] = (double)0.0;

			pointer->DELTA[v] = (double)0.0;

			pointer->the_SRE_shifted_fMRI[v] = (double)0.0;

			pointer->SE_SRE[v] = (double)0.0;

			pointer->X_recomputed[v] = (double)0.0;

			pointer->Y_recomputed[v] = (double)0.0;

			pointer->Z_recomputed[v] = (double)0.0;

			pointer->X_recomputed_scaled[v] = (double)0.0;

			pointer->Y_recomputed_scaled[v] = (double)0.0;

			pointer->Z_recomputed_scaled[v] = (double)0.0;


			pointer->thetaXYZ[v] = (double)0.0;

			pointer->F_XYZ[v] = (double)0.0; 

			pointer->fxyz[v] = (double)0.0;

			pointer->F_XY[v] = (double)0.0;

			pointer->F_YZ[v] = (double)0.0;

			pointer->F_XZ[v] = (double)0.0;
	
			pointer->deDELTAdeX[v] = (double)0.0;

			pointer->deDELTAdeY[v] = (double)0.0;

			pointer->deDELTAdeZ[v] = (double)0.0;


			pointer->Recomp_Mx_plus[v] = (double)0.0;

			pointer->Recomp_Mx_minus[v] = (double)0.0;


			pointer->Recomp_My_plus[v] = (double)0.0;

			pointer->Recomp_My_minus[v] = (double)0.0;


			pointer->Recomp_Mz_plus[v] = (double)0.0;

			pointer->Recomp_Mz_minus[v] = (double)0.0;

			pointer->ICF[v] = (double) 0.0;
			

		 } //(a)
		// (2) initialize (end)

} // allocate data


void SRE3D2013::save() { // saveImages

	FILE * savedata;
	char outputFile[128];
	int size = ((int)this->n1 * (int)this->n2 * (int)this->n4);
	
	sprintf(outputFile, "%s","fMRI.img");

	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for (int sl=1; sl <= (int)this->n4; sl++) { // z dim
       	
			for (int i2=0; i2 < (int)this->n2; i2++) { // y dim

				for (int i1=0; i1 < (int)this->n1; i1++) {// x dim

				int v = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)
	
				fwrite(&pointer->fMRI[v],sizeof(double),1,savedata);
                          
			} // y dim
        
		}  // x dim 

	}  // z dim 
		
	fclose(savedata);

	} // (save)

	
	sprintf(outputFile, "%s","SRE_Xcoord.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for (int sl=1; sl <= (int)this->n4; sl++) { // z dim
       	
			for (int i2=0; i2 < (int)this->n2; i2++) { // y dim

				for (int i1=0; i1 < (int)this->n1; i1++) {// x dim

				int v = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)
	
				fwrite(&pointer->SRE_Xcoord[v],sizeof(double),1,savedata);
                          
			} // y dim
        
		}  // x dim 

	}  // z dim 
		
	fclose(savedata);

	} // (save)
	


	sprintf(outputFile, "%s","SRE_Ycoord.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)

	for (int sl=1; sl <= (int)this->n4; sl++) { // z dim
       	
			for (int i2=0; i2 < (int)this->n2; i2++) { // y dim

				for (int i1=0; i1 < (int)this->n1; i1++) {// x dim

				int v = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)
	
				fwrite(&pointer->SRE_Ycoord[v],sizeof(double),1,savedata);
                          
			} // y dim
        
		}  // x dim 

	}  // z dim 

	fclose(savedata);

	} // (save)


	sprintf(outputFile, "%s","SRE_Zcoord.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for (int sl=1; sl <= (int)this->n4; sl++) { // z dim
       	
			for (int i2=0; i2 < (int)this->n2; i2++) { // y dim

				for (int i1=0; i1 < (int)this->n1; i1++) {// x dim

				int v = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)
	
				fwrite(&pointer->SRE_Zcoord[v],sizeof(double),1,savedata);
                          
			} // y dim
        
		}  // x dim 

	}  // z dim 

	fclose(savedata);

	} // (save)
	

	sprintf(outputFile, "%s","shifted_fMRI.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for (int sl=1; sl <= (int)this->n4; sl++) { // z dim
       	
			for (int i2=0; i2 < (int)this->n2; i2++) { // y dim

				for (int i1=0; i1 < (int)this->n1; i1++) {// x dim

				int v = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)
	
				fwrite(&pointer->shifted_fMRI[v],sizeof(double),1,savedata);
                          
			} // y dim
        
		}  // x dim 

	}  // z dim 

	fclose(savedata);

	} // (save)


	sprintf(outputFile, "%s","SE.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for (int sl=1; sl <= (int)this->n4; sl++) { // z dim
       	
			for (int i2=0; i2 < (int)this->n2; i2++) { // y dim

				for (int i1=0; i1 < (int)this->n1; i1++) {// x dim

				int v = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)
	
				fwrite(&pointer->SE[v],sizeof(double),1,savedata);
                          
			} // y dim
        
		}  // x dim 

	}  // z dim 

	fclose(savedata);

	} // (save)

	
	sprintf(outputFile, "%s","DELTA.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for (int sl=1; sl <= (int)this->n4; sl++) { // z dim
       	
			for (int i2=0; i2 < (int)this->n2; i2++) { // y dim

				for (int i1=0; i1 < (int)this->n1; i1++) {// x dim

				int v = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)
	
				fwrite(&pointer->DELTA[v],sizeof(double),1,savedata);
                          
			} // y dim
        
		}  // x dim 

	}  // z dim 

	fclose(savedata);

	} // (save)


	sprintf(outputFile, "%s","the_SRE_shifted_fMRI.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for (int sl=1; sl <= (int)this->n4; sl++) { // z dim
       	
			for (int i2=0; i2 < (int)this->n2; i2++) { // y dim

				for (int i1=0; i1 < (int)this->n1; i1++) {// x dim

				int v = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)
	
				fwrite(&pointer->the_SRE_shifted_fMRI[v],sizeof(double),1,savedata);
                          
			} // y dim
        
		}  // x dim 

	}  // z dim 

	fclose(savedata);

	} // (save)

	
	sprintf(outputFile, "%s","SE_SRE.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for (int sl=1; sl <= (int)this->n4; sl++) { // z dim
       	
			for (int i2=0; i2 < (int)this->n2; i2++) { // y dim

				for (int i1=0; i1 < (int)this->n1; i1++) {// x dim

				int v = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)
	
				fwrite(&pointer->SE_SRE[v],sizeof(double),1,savedata);
                          
			} // y dim
        
		}  // x dim 

	}  // z dim 

	fclose(savedata);

	} // (save)


	sprintf(outputFile, "%s","X_recomputed.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for (int sl=1; sl <= (int)this->n4; sl++) { // z dim
       	
			for (int i2=0; i2 < (int)this->n2; i2++) { // y dim

				for (int i1=0; i1 < (int)this->n1; i1++) {// x dim

				int v = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)
	
				fwrite(&pointer->X_recomputed[v],sizeof(double),1,savedata);
                          
			} // y dim
        
		}  // x dim 

	}  // z dim 

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","Y_recomputed.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)
	
	for (int sl=1; sl <= (int)this->n4; sl++) { // z dim
       	
			for (int i2=0; i2 < (int)this->n2; i2++) { // y dim

				for (int i1=0; i1 < (int)this->n1; i1++) {// x dim

				int v = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)
	
				fwrite(&pointer->Y_recomputed[v],sizeof(double),1,savedata);
                          
			} // y dim
        
		}  // x dim 

	}  // z dim 

	fclose(savedata);

	} // (save)


	sprintf(outputFile, "%s","Z_recomputed.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for (int sl=1; sl <= (int)this->n4; sl++) { // z dim
       	
			for (int i2=0; i2 < (int)this->n2; i2++) { // y dim

				for (int i1=0; i1 < (int)this->n1; i1++) {// x dim

				int v = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)
	
				fwrite(&pointer->Z_recomputed[v],sizeof(double),1,savedata);
                          
			} // y dim
        
		}  // x dim 

	}  // z dim 

	fclose(savedata);

	} // (save)


	sprintf(outputFile, "%s","X_recomputed_scaled.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for (int sl=1; sl <= (int)this->n4; sl++) { // z dim

		for (int i2=0; i2 < (int)this->n2; i2++) { // y dim

				for (int i1=0; i1 < (int)this->n1; i1++) {// x dim

				int v = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)
	
				fwrite(&pointer->X_recomputed_scaled[v],sizeof(double),1,savedata);
                          
			} // y dim
        
		}  // x dim 

	}  // z dim 

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","Y_recomputed_scaled.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for (int sl=1; sl <= (int)this->n4; sl++) { // z dim

		for (int i2=0; i2 < (int)this->n2; i2++) { // y dim

				for (int i1=0; i1 < (int)this->n1; i1++) {// x dim

				int v = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)
	
				fwrite(&pointer->Y_recomputed_scaled[v],sizeof(double),1,savedata);
                          
			} // y dim
        
		}  // x dim 

	}  // z dim 

	fclose(savedata);

	} // (save)

	sprintf(outputFile, "%s","Z_recomputed_scaled.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for (int sl=1; sl <= (int)this->n4; sl++) { // z dim

		for (int i2=0; i2 < (int)this->n2; i2++) { // y dim

				for (int i1=0; i1 < (int)this->n1; i1++) {// x dim

				int v = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)
	
				fwrite(&pointer->Z_recomputed_scaled[v],sizeof(double),1,savedata);
                          
			} // y dim
        
		}  // x dim 

	}  // z dim 


	fclose(savedata);

	} // (save)

		
	sprintf(outputFile, "%s","ICF.img");
	if ((savedata = fopen(outputFile,"wb+"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // (save)


	for (int sl=1; sl <= (int)this->n4; sl++) { // z dim

		for (int i2=0; i2 < (int)this->n2; i2++) { // y dim

				for (int i1=0; i1 < (int)this->n1; i1++) {// x dim

				int v = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)
	
				fwrite(&pointer->ICF[v],sizeof(double),1,savedata);
                          
			} // y dim
        
		}  // x dim 

	}  // z dim 


	fclose(savedata);

	} // (save)

} // saveImages


void SRE3D2013::OnFourierAnalysis(int m_Xpixels, int m_Ypixels, char imageFilename[]) 
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

	char outputFile[128]="SRE3D.log";

	FILE * savedata;

	double MAX = 5000000000000000000.0;

if (argc < 13) { std::cout << endl;
				 std::cout << "Please type the image file name" << endl;
				 std::cout << "Please make sure that the image format is Analyze 'double': 64 bits real" << endl;
				 std::cout << "Before running the program, please make sure that the image is padded of 'k >= 1'";
				 std::cout << "number of pixels along X, Y and Z" << endl;
				 std::cout << "Please enter the number of pixels along the X direction (integer)" << endl;
				 std::cout << "Please enter the number of pixels along the Y direction (integer)" << endl;
				 std::cout << "Please enter the number of pixels along the Z direction (integer)" << endl;
				 std::cout << "Please enter the pixel size along the X direction (double)" << endl;
				 std::cout << "Please enter the pixel size along the Y direction (double)" << endl;
				 std::cout << "Please enter the pixel size along the Z direction (double)" << endl;
				 std::cout << "Please enter the misplacement along the X direction (double)" << endl;
				 std::cout << "Please enter the misplacement along the Y direction (double)" << endl;
				 std::cout << "Please enter the misplacement along the Z direction (double)" << endl;
				 std::cout << "Please enter that 'n = k - 1' number of pixels along X, Y and Z which will pad the image" << endl;
				 std::cout << "Some data will be lost if you enter a value of n which is greater than" << endl;
				 std::cout << "the actual number of pixels along X, Y and Z which was used to pad the image" << endl;
				 std::cout << "before running the program." << endl;
				 std::cout << "Please type y to run the Spectral Power Analysis (y (yes) or n (no))" << endl;
				 std::cout << endl;
				 exit(0); }

else { // run the program (begin)

	
	if ((savedata = fopen(outputFile,"w"))==NULL)
	{

		std::cout << "Cannot open output file, Now Exit..." << endl;

	} else  { // processing (begin)

	char imageFileName[128];

	sprintf(imageFileName, "%s", argv[1]);

	int n1 = atoi(argv[2]);
	int n2 = atoi(argv[3]);
	int n4 = atoi(argv[4]);

	double XPixelSize = atof(argv[5]);
	double YPixelSize = atof(argv[6]);
	double ZPixelSize = atof(argv[7]);

	double x_misplacement_X = atof(argv[8]);
	double y_misplacement_Y = atof(argv[9]);
	double z_misplacement_Z = atof(argv[10]);

	double TOLL = (double)0.0001;			// SRE convergence tolerance
	double SMOOTH = (double)0.1;			// SRE convergence smoothing factor
	double ITERATIONFACTOR = (double)0.01;  // SRE convergence speed factor

	int PAD = atoi(argv[11]);
	

	if ( PAD < 1 ) { // if test PAD
	
		std::cout << "Please make sure to pad the image with 'n >= 1' number of pixels along X, Y and Z" << endl;
	
		fprintf(savedata,"%s\n", "Please make sure to pad the image with 'n >= 1' number of pixels along X, Y and Z");

		fclose(savedata);
	
		exit(0);
	
	} // if test PAD

	unsigned char Fyn = 'n';

	Fyn = *argv[12];	

	std::cout << endl;
	std::cout << "The image file name is: " << imageFileName << endl;
	std::cout << "The number of pixels along the X direction is: " << atoi(argv[2]) << endl;
	std::cout << "The number of pixels along the Y direction is: " << atoi(argv[3]) << endl;
	std::cout << "The number of pixels along the Z direction is: " << atoi(argv[4]) << endl;
	std::cout << "The pixel size along the X direction is: " << atof(argv[5]) << endl;
	std::cout << "The pixel size along the Y direction is: " << atof(argv[6]) << endl;
	std::cout << "The pixel size along the Z direction is: " << atof(argv[7]) << endl;
	std::cout << "The misplacement along the X direction is: " << atof(argv[8]) << endl;
	std::cout << "The misplacement along the Y direction is: " << atof(argv[9]) << endl;
	std::cout << "The misplacement along the Y direction is: " << atof(argv[10]) << endl;
	std::cout << "The 'n' number of pixels along X, Y and Z which you entered to pad the image is: " << atoi(argv[11]) << endl;
	std::cout << "Some data will be lost if you enter a value of n which is greater than" << endl;
	std::cout << "the actual number of pixels along X, Y and Z which was used to pad the image" << endl;
	std::cout << "before running the program." << endl;
	std::cout << endl;
	if ( tolower(Fyn) == 'y' )
	std::cout << "The program runs the Spectral Power Analysis: " << Fyn << endl;
	else if ( tolower(Fyn) != 'y' ) std::cout << "The program does not run the Spectral Power Analysis: " << Fyn << endl;
	std::cout << endl;

	fprintf(savedata,"%s%s\n", "The image file name is: " , imageFileName);
	fprintf(savedata,"%s%d\n", "The number of pixels along the X direction is: ", n1);
	fprintf(savedata,"%s%d\n", "The number of pixels along the Y direction is: ", n2);
	fprintf(savedata,"%s%d\n", "The number of pixels along the Z direction is: ", n4);
	fprintf(savedata,"%s%lf\n", "The pixel size along the X direction is: ", XPixelSize);
	fprintf(savedata,"%s%lf\n", "The pixel size along the Y direction is: ", YPixelSize);
	fprintf(savedata,"%s%lf\n", "The pixel size along the Z direction is: ", ZPixelSize);	
	fprintf(savedata,"%s%e\n", "The misplacement along the X direction is: ", x_misplacement_X);
	fprintf(savedata,"%s%e\n", "The misplacement along the Y direction is: ", y_misplacement_Y);
	fprintf(savedata,"%s%e\n", "The misplacement along the Y direction is: ", z_misplacement_Z);
	fprintf(savedata,"%s%d\n", "The 'n' number of pixels along X and Y which you entered to pad the image is: ", PAD);
	fprintf(savedata,"%s\n", "Some data will be lost if you enter a value of n which is greater than");
	fprintf(savedata,"%s\n", "the actual number of pixels along X, Y and Z which was used to pad the image");
	fprintf(savedata,"%s\n", "before running the program.");
	fprintf(savedata,"\n");
	if ( tolower(Fyn) == 'y' )
	fprintf(savedata,"%s%c\n", "The program runs the Spectral Power Analysis: ", Fyn);
	else if ( tolower(Fyn) != 'y' )
	fprintf(savedata,"%s%c\n", "The program does not run the Spectral Power Analysis: ", Fyn);
	fprintf(savedata,"\n");

    double Mx = ((double)x_misplacement_X);
    double My = ((double)y_misplacement_Y);
	double Mz = ((double)z_misplacement_Z);

    Mx = ((double)Mx/XPixelSize);
    My = ((double)My/YPixelSize);
	Mz = ((double)Mz/ZPixelSize);

	  //////////////////***********//////////////////////
	  // Above formulas scale the misplacement to the  //
	  // pixel size the same way the following formula //
	  // would do: (min - misplacement)/(min - max)    //  
	  //////////////////***********//////////////////////

	double INITGX = (double)Mx;
	double INITGY = (double)My;
	double INITGZ = (double)Mz;

	
	SRE3D2013 SRE(n1,n2,n4);

	SRE.allocateData();

	/// read image file (begin)
	FILE * pf;

	if ((pf = fopen(imageFileName,"rb+"))==NULL)
	{

		std::cout << "Cannot open file: " << imageFileName << endl;
		fprintf(savedata,"%s%s\n", "Cannot open file: " , imageFileName );
		exit(0);

	} else { // else

	double number;

	for (int sl=1; sl <= n4; sl++) { // z dim

		for (int i2=0; i2 < n2; i2++) { // y dim

			for (int i1=0; i1 < n1; i1++) {// x dim

			    int k1 = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)

				fread(&number,sizeof(double),1,pf);
		
				SRE.pointer->fMRI[k1] = (double)number;

				SRE.pointer->scaledfMRI[k1] = (double)number;
                          
			} // y dim
        
		}  // x dim 

	}  // z dim 

      	
    fclose (pf);


	} // else 
	/// read image file (end)

	std::cout << "Image data loaded" << endl;

	/// compute max and min of data (begin)
	double max=-MAX;
	double min=MAX;

	for (int sl=1; sl <= n4; sl++) { // z dim

		for (int i2=0; i2 < n2; i2++) { // y dim

			for (int i1=0; i1 < n1; i1++) {// x dim

			int k1 = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)

			if( SRE.pointer->fMRI[k1] > (double)max ) 
			
			max = (double)SRE.pointer->fMRI[k1];
              
			if( SRE.pointer->fMRI[k1] < (double)min ) 
			
			min = (double)SRE.pointer->fMRI[k1];
		

			} // y dim
        
		}  // x dim
	
	}  // z dim 
	/// compute max and min of data (end)

	// scale (begin)
	for (int sl=1; sl <= n4; sl++) { // z dim

		for (int i2=0; i2 < n2; i2++) { // y dim

			for (int i1=0; i1 < n1; i1++) {// x dim

			    int v = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)

				if ( (double)SRE.pointer->fMRI[v] > TH ) { 

				if ( max == min ) SRE.pointer->scaledfMRI[v] = (double)0.0;

				else SRE.pointer->scaledfMRI[v] = 
				(double) NSCALE * (double) fabs ( (min - SRE.pointer->fMRI[v]) / (min - max) );
             
				} else if ( (double)SRE.pointer->fMRI[v] <= TH ) {
  	
					 SRE.pointer->scaledfMRI[v] = (double)0.0;
			
				}

			} // y dim
        
		}  // x dim 

	}  // z dim 
	// scale (end)


	std::cout << "Image data scaled" << endl;


	int numslice = (int)n4;
      
    int rcxres = (int)n1;
      
    int rcyres = (int)n2;

	int XYZNEI = ((int)PAD + 1);  // to process the image volume consistently 
							      // with the PAD value entered by the user

	int k1, k2, k3, k4, k5, k6, k7, k8;
	double omega, theta1, theta2, theta3;
	double beta1, beta2, beta3;
	double gamma1, gamma2, gamma3;
	double delta1;

	// shift forward (begin)
	for (int sl=PAD; sl <= n4; sl++) { // z dim

		for (int j=0; j < n2; j++) { // y dim

			for (int i=0; i < n1; i++) {// x dim

			int k1 = (j*n1+i) + ((sl-2)*n1*n2); // f(0,0,0)

				if ( (double)SRE.pointer->fMRI[k1] > TH ) {  // major if

                        
					   
                        k2=((j+1)*rcxres+i)+( (sl-2)*rcyres*rcxres); // f(0,1,0) 
			
                        k3=(j*rcxres+(i+1))+( (sl-2)*rcyres*rcxres); // f(1,0,0) 
  			           
                        k4=((j+1)*rcxres+(i+1))+( (sl-2)*rcyres*rcxres); // f(1,1,0) 
            
                        k5=(j*rcxres+i)+( (sl-1)*rcyres*rcxres); // f(0,0,1)
			
                        k6=((j+1)*rcxres+i)+( (sl-1)*rcyres*rcxres); // f(0,1,1) 
					   
                        k7=(j*rcxres+(i+1))+( (sl-1)*rcyres*rcxres); // f(1,0,1) 
  
                        k8=((j+1)*rcxres+(i+1))+((sl-1)*rcyres*rcxres); // f(1,1,1) 
						
              
            omega = ((double) SRE.pointer->scaledfMRI[k8] - SRE.pointer->scaledfMRI[k4] - SRE.pointer->scaledfMRI[k6] - 
			
			                  SRE.pointer->scaledfMRI[k7] + SRE.pointer->scaledfMRI[k5] + SRE.pointer->scaledfMRI[k2] + 
					
					          SRE.pointer->scaledfMRI[k3] - SRE.pointer->scaledfMRI[k1] );    
     
            theta1 = ( -((double) SRE.pointer->scaledfMRI[k1] - SRE.pointer->scaledfMRI[k3] - 
				
				                  SRE.pointer->scaledfMRI[k2] + SRE.pointer->scaledfMRI[k4] ) ); // thetaxy
 
            theta2 = ( -((double) SRE.pointer->scaledfMRI[k1] - SRE.pointer->scaledfMRI[k2] - 
				
				                  SRE.pointer->scaledfMRI[k5] + SRE.pointer->scaledfMRI[k6] ) ); // thetayz 

            theta3 = ( -((double) SRE.pointer->scaledfMRI[k1] - SRE.pointer->scaledfMRI[k3] - 
				
				                  SRE.pointer->scaledfMRI[k5] + SRE.pointer->scaledfMRI[k7] ) ); // thetaxz 

             
            beta1 = (double)Mx * ((double) SRE.pointer->scaledfMRI[k3] - SRE.pointer->scaledfMRI[k1] );  // thetax

            beta2 = (double)My * ((double) SRE.pointer->scaledfMRI[k2] - SRE.pointer->scaledfMRI[k1] );  // thetay

            beta3 = (double)Mz * ((double) SRE.pointer->scaledfMRI[k5] - SRE.pointer->scaledfMRI[k1] );  // thetaz
               
            gamma1 = ((double) Mx * My * theta1);

            gamma2 = ((double) My * Mz * theta2);

            gamma3 = ((double) Mx * Mz * theta3);   

            delta1 = ((double) Mx * My * Mz * omega);
                
                
            SRE.pointer->shifted_fMRI[k1] = (double)SRE.pointer->fMRI[k1] + 
			                               ((double)beta1 + beta2 + beta3 + gamma1 + gamma2 + gamma3 + delta1);
              
				} else if ( (double)SRE.pointer->fMRI[k1] <= TH ) { // major if

					SRE.pointer->shifted_fMRI[k1] = (double)0.0;
			
				} // major if


 			} // x dim
        
		}  // y dim 

	}  // z dim 
	// shift forward (end)

	// calculate SE (begin)
	double SE_sum = (double)0.0;
	long int counter = 0;

	for (int sl=PAD; sl <= n4; sl++) { // z dim

		for (int i2=0; i2 < n2; i2++) { // y dim

			for (int i1=0; i1 < n1; i1++) {// x dim

			 int k1 = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)

				if ( (double)SRE.pointer->fMRI[k1] > TH ) {  // major if

				counter++;

				SE_sum += (double)sqrt( ((double)(SRE.pointer->shifted_fMRI[k1] - SRE.pointer->fMRI[k1] ) * 
					                             (SRE.pointer->shifted_fMRI[k1] - SRE.pointer->fMRI[k1] ) ) ); 

				SRE.pointer->SE[k1] = (double)sqrt( ((double)(SRE.pointer->shifted_fMRI[k1] - SRE.pointer->fMRI[k1] ) * 
					                                         (SRE.pointer->shifted_fMRI[k1] - SRE.pointer->fMRI[k1] ) ) ); 

				}  // major if

			} // y dim
        
		}  // x dim 

	}  // z dim 
	// calculate SE (end)


	std::cout << "SE with no SRE is: " << ((double)SE_sum / (counter))  << endl;
	std::cout << endl;

	fprintf(savedata,"%s%e\n", "SE with no SRE is: ", ((double)SE_sum / (counter)) );
	fprintf(savedata,"\n");

    // calculate SRE (begin)
	std::cout << "Calculating SRE..." << endl;
	fprintf(savedata,"%s\n", "Calculating SRE...");

	double convergence;
	double b1, b2, b3, b4, b5, b6, b7, b8;
	double a1, a2, a3, a4, b, a, c1, c2, c;
	double delta, storeZ, storeY, storeX;
	double d1, d2;

	for (int sl=PAD; sl <= n4; sl++) { // z dim

		for (int j=0; j < n2; j++) { // y dim

			for (int i=0; i < n1; i++) {// x dim

			int k1 = (j*n1+i) + ((sl-2)*n1*n2); // f(0,0,0)

			   if ( (double)SRE.pointer->fMRI[k1] > TH ) {  // major if


					   
               k2=((j+1)*rcxres+i)+( (sl-2)*rcyres*rcxres);     // f(0,1,0) 
			
		       k3=(j*rcxres+(i+1))+( (sl-2)*rcyres*rcxres);     // f(1,0,0) 
  			           
               k4=((j+1)*rcxres+(i+1))+( (sl-2)*rcyres*rcxres); // f(1,1,0) 
            
		       k5=(j*rcxres+i)+( (sl-1)*rcyres*rcxres);         // f(0,0,1)
			
		       k6=((j+1)*rcxres+i)+( (sl-1)*rcyres*rcxres);     // f(0,1,1) 
					   
		       k7=(j*rcxres+(i+1))+( (sl-1)*rcyres*rcxres);     // f(1,0,1) 
  
               k8=((j+1)*rcxres+(i+1))+((sl-1)*rcyres*rcxres);  // f(1,1,1) 



			omega = ((double) SRE.pointer->scaledfMRI[k8] - SRE.pointer->scaledfMRI[k4] - SRE.pointer->scaledfMRI[k6] - 
			
			                  SRE.pointer->scaledfMRI[k7] + SRE.pointer->scaledfMRI[k5] + SRE.pointer->scaledfMRI[k2] + 
					
					          SRE.pointer->scaledfMRI[k3] - SRE.pointer->scaledfMRI[k1] );    
     
            theta1 = ( -((double) SRE.pointer->scaledfMRI[k1] - SRE.pointer->scaledfMRI[k3] - 
				
				                  SRE.pointer->scaledfMRI[k2] + SRE.pointer->scaledfMRI[k4] ) ); // thetaxy
 
            theta2 = ( -((double) SRE.pointer->scaledfMRI[k1] - SRE.pointer->scaledfMRI[k2] - 
				
				                  SRE.pointer->scaledfMRI[k5] + SRE.pointer->scaledfMRI[k6] ) ); // thetayz 

            theta3 = ( -((double) SRE.pointer->scaledfMRI[k1] - SRE.pointer->scaledfMRI[k3] - 
				
				                  SRE.pointer->scaledfMRI[k5] + SRE.pointer->scaledfMRI[k7] ) ); // thetaxz 


			beta1 = ((double) SRE.pointer->scaledfMRI[k3] - SRE.pointer->scaledfMRI[k1] );  // thetax

            beta2 = ((double) SRE.pointer->scaledfMRI[k2] - SRE.pointer->scaledfMRI[k1] );  // thetay

            beta3 = ((double) SRE.pointer->scaledfMRI[k5] - SRE.pointer->scaledfMRI[k1] );  // thetaz


			 // Initialization of the SRE
             SRE.pointer->SRE_Xcoord[k1] = (double)INITGX;
             SRE.pointer->SRE_Ycoord[k1] = (double)INITGY;
             SRE.pointer->SRE_Zcoord[k1] = (double)INITGZ;


			 convergence = (double)1.0;

			 while ( (double)convergence >= (double)TOLL )
			 { // while loop 

				 /// calculation of the SRE Z coordinate (begin)
				 b1 = ((double)( SRE.pointer->scaledfMRI[k1] + ( (double)SRE.pointer->SRE_Ycoord[k1] * beta2 ) ) * (double)theta3 * omega / 2.0);

 				 b2 = ((double)( - SRE.pointer->scaledfMRI[k1] - ( (double)SRE.pointer->SRE_Ycoord[k1] * beta2 ) ) * (double)SRE.pointer->SRE_Ycoord[k1] * omega / 4.0);

				 b3 =  ((double) ( beta1 / 2.0 ) - ( ( (double)SRE.pointer->SRE_Ycoord[k1] * theta1 * omega ) / 2.0 ) ) * 
					   ((double) ( - beta3 + ( ( (double)SRE.pointer->SRE_Ycoord[k1] * theta2 * omega ) / 2.0 ))) ;

 	             b4 = ((double)beta3 * beta1 / 4.0);

				 b5 = ((double)SRE.pointer->SRE_Ycoord[k1] * theta2 * omega * beta1 / 4.0);

				 b6 = ((double)SRE.pointer->SRE_Ycoord[k1] * theta1 * omega * beta3 / 8.0);

				 b7 = ((double)SRE.pointer->SRE_Ycoord[k1] * SRE.pointer->SRE_Ycoord[k1] * theta2 * omega * omega * theta1 / 8.0);

				 b8 =  ( ((double)( theta3 / 4.0 ) * omega ) -  ( (double)SRE.pointer->SRE_Ycoord[k1] * omega / 4.0 ) ) * 
					   ( - (double)SRE.pointer->scaledfMRI[k1] - ( (double)SRE.pointer->SRE_Ycoord[k1] * beta2 / 2.0 ) );


				 b = ((double)b1 + b2 + b3 - b4 + b5 + b6 - b7 + b8);


				 a1 = ((double)( ((double) beta3 / 2.0 ) - ( (double)SRE.pointer->SRE_Ycoord[k1] * omega * theta2 / 2.0 ) ) * theta3 * omega / 2.0);

 				 a2 = ((double) beta3 * SRE.pointer->SRE_Ycoord[k1] * omega / 8.0);

				 a3 = ((double) SRE.pointer->SRE_Ycoord[k1] * SRE.pointer->SRE_Ycoord[k1] * theta2 * omega * omega / 8.0);

 				 a4 =  (  ((double) ( theta3 / 4.0 ) * omega ) -  ( (double)SRE.pointer->SRE_Ycoord[k1] * omega / 4.0 )  )  * 
			           ( - (double)beta3 + ( (double)SRE.pointer->SRE_Ycoord[k1] * theta2 * omega / 2.0 )  );


				 a = ((double)a1 - a2 + a3 + a4);


				 c1 = ((double)( - (double)SRE.pointer->scaledfMRI[k1] - SRE.pointer->SRE_Ycoord[k1] * beta2 )) * 
					  ((double)( ( (double)beta1 / 2.0 ) - ( (double)SRE.pointer->SRE_Ycoord[k1] * theta1 * omega / 4.0 ) ));


 				 c2 = ((double)( ( (double)beta1 / 2.0 ) - ( (double)SRE.pointer->SRE_Ycoord[k1] * theta1 * omega / 2.0 ) )) * 
					  ((double)( - (double)SRE.pointer->scaledfMRI[k1] - SRE.pointer->SRE_Ycoord[k1] * beta2 / 2.0 ));


				 c = ((double)c1 - c2);

				 delta = ((double) b * b - 4.0 * a * c);


				 storeZ = (double)SRE.pointer->SRE_Zcoord[k1];


				 if ( ((double)delta >= 0.0) && ((double)a != 0.0) )
				 { // if

					 if ( ZPLUS == 1 ) 

					 SRE.pointer->SRE_Zcoord[k1] = (double)storeZ + ( (double)ITERATIONFACTOR * ( (double) - b + (double)sqrt( (double)delta ) ) / ( (double)2.0 * a ) );

					 else if ( ZPLUS == 0 ) 

					 SRE.pointer->SRE_Zcoord[k1] = (double)storeZ + ( (double)ITERATIONFACTOR * ( (double) - b - (double)sqrt( (double)delta ) ) / ( (double)2.0 * a ) );
				 				 
				 } // if 
				 else if ( ((double)delta < 0.0) && ((double)a != 0.0) )
				 { // if

					 if ( ZPLUS == 1 ) 

					 SRE.pointer->SRE_Zcoord[k1] = (double)storeZ + ( (double)ITERATIONFACTOR * ( (double) - b ) / ( (double)2.0 * a ) );

					 else if ( ZPLUS == 0 ) 

					 SRE.pointer->SRE_Zcoord[k1] = (double)storeZ + ( (double)ITERATIONFACTOR * ( (double) - b ) / ( (double)2.0 * a ) );
						 
				 } // if 
				 else if ( ((double)a == 0.0) )
				 { // if

					 if ( ZPLUS == 1 ) 

					 SRE.pointer->SRE_Zcoord[k1] = (double)storeZ;

					 else if ( ZPLUS == 0 ) 

					 SRE.pointer->SRE_Zcoord[k1] = (double)storeZ;
				 	
				 } // if
				 /// calculation of the SRE Z coordinate (end)


				 /// calculation of the SRE Y coordinate (begin)
				 b1 = ( - (double)SRE.pointer->scaledfMRI[k1] * theta1 * omega / 4.0 );

				 b2 = (   (double)SRE.pointer->scaledfMRI[k1] * SRE.pointer->SRE_Zcoord[k1] * omega / 4.0 );

 	             b3 = ( - (double)SRE.pointer->scaledfMRI[k1] - ( (double)SRE.pointer->SRE_Zcoord[k1] * beta3 ) );

				 b4 = (((double)  ( - (double)theta1 * omega / 2.0 ) + ((double)SRE.pointer->SRE_Zcoord[k1] * omega / 2.0) + 
					              ((double)beta1) - ((double)SRE.pointer->SRE_Zcoord[k1] * omega * theta3)   ) / 2.0  );

				 b5 = ( ((double)beta2 / 2.0) - ((double)SRE.pointer->SRE_Zcoord[k1] * omega * theta2 / 2.0) ) * 
					  ( (double)beta1  - ((double)SRE.pointer->SRE_Zcoord[k1] * omega * theta3 / 2.0) );

				 b6 = ((double)SRE.pointer->SRE_Zcoord[k1] * beta3 * 
					   ((double)  ( - (double)theta1 * omega / 4.0) + ((double)SRE.pointer->SRE_Zcoord[k1] * omega / 4.0) ) );

		
				 b = ((double)b1 + b2 + ( b3 * b4 ) + ( b5 * 2.0 ) + b6);

				 a1 = ( ((double) beta2 / 2.0) - ((double)SRE.pointer->SRE_Zcoord[k1] * omega * theta2 / 2.0) );

 	             a2 = (  ((double)theta1 * omega / 2.0) - ((double)SRE.pointer->SRE_Zcoord[k1] * omega / 2.0) + 
					  ( - (double)theta1 * omega / 4.0) + ((double)SRE.pointer->SRE_Zcoord[k1] * omega / 4.0) ) ;


				 a = ((double)a1 * a2);

				 c = (double)SRE.pointer->scaledfMRI[k1] * ( ((double)beta1 / 2.0) - ((double)SRE.pointer->SRE_Zcoord[k1] * omega * theta3 / 2.0) );

				 c  = c + ( (double)SRE.pointer->SRE_Zcoord[k1] * beta3 * ( ((double) beta1 / 2.0) - ((double)SRE.pointer->SRE_Zcoord[k1] * theta3 * omega / 2.0 ) ) );


   				 delta = ((double) b * b - ((double)4.0 * a * c));

		
				 storeY = (double)SRE.pointer->SRE_Ycoord[k1];

				 if ( ((double)delta >= 0.0) && ((double)a != 0.0) )
				 { // if

					 if (YPLUS == 1)

					 SRE.pointer->SRE_Ycoord[k1] = (double)storeY + ( (double)ITERATIONFACTOR * ( (double)- b + (double)sqrt((double)delta) ) / ((double)2.0 * a) );

					 else if (YPLUS == 0)

					 SRE.pointer->SRE_Ycoord[k1] = (double)storeY + ( (double)ITERATIONFACTOR * ( (double)- b - (double)sqrt((double)delta) ) / ((double)2.0 * a) );

				 } // if
				 else if ( ((double)delta < 0.0) && ((double)a != 0.0) )
				 { // if

					 if (YPLUS == 1)

					 SRE.pointer->SRE_Ycoord[k1] = (double)storeY + ( (double)ITERATIONFACTOR * ( (double)- b ) / ((double)2.0 * a) );

					 else if (YPLUS == 0)

					 SRE.pointer->SRE_Ycoord[k1] = (double)storeY + ( (double)ITERATIONFACTOR * ( (double)- b ) / ((double)2.0 * a) );

				 } // if
				 else if ( ((double)a == 0.0) )
				 { // if

					 if (YPLUS == 1)

					 SRE.pointer->SRE_Ycoord[k1] = (double)storeY;

					 else if (YPLUS == 0)

					 SRE.pointer->SRE_Ycoord[k1] = (double)storeY;

				 } // if
				 /// calculation of the SRE Y coordinate (end)



				 /// calculation of the SRE X coordinate (begin)
				 storeX = (double)SRE.pointer->SRE_Xcoord[k1];

				 d1 = - ((double)SRE.pointer->scaledfMRI[k1]) - 
					    ((double)SRE.pointer->SRE_Ycoord[k1] * beta2 / 2.0) - 
					    ((double)SRE.pointer->SRE_Zcoord[k1] * beta3) + 
				  	    ((double)SRE.pointer->SRE_Ycoord[k1] * SRE.pointer->SRE_Zcoord[k1] * omega * theta2 / 2.0) ;

	             d2 = ((double) beta1 / 2.0) - ((double)SRE.pointer->SRE_Ycoord[k1] * theta1 * omega / 4.0) - 
					                           ((double)SRE.pointer->SRE_Zcoord[k1] * theta3 * omega / 2.0) + 
											   ((double)SRE.pointer->SRE_Ycoord[k1] * SRE.pointer->SRE_Zcoord[k1] * omega / 4.0);

				 if ( (double)d2 != 0.0 )
				 { // if
				
					SRE.pointer->SRE_Xcoord[k1] = ((double)storeX) + ((double)ITERATIONFACTOR * ( d1 / d2 ));
				 
				 } // if 
				 else if ( (double)d2 == 0.0 ) SRE.pointer->SRE_Xcoord[k1] = ((double)storeX);
				 /// calculation of the SRE X coordinate (end)


				 /// convergence test (begin)
				 convergence = (double)abs ( ((double)storeX - SRE.pointer->SRE_Xcoord[k1]) ) + 
					           (double)abs ( ((double)storeY - SRE.pointer->SRE_Ycoord[k1]) ) + 
							   (double)abs ( ((double)storeZ - SRE.pointer->SRE_Zcoord[k1]) ) ;	    

				 if ( (double)convergence >= (double)SMOOTH ) 
				 // do not allow abrupt change of the SRE value
				 { // if

					  SRE.pointer->SRE_Xcoord[k1] = (double)storeX; 
					  SRE.pointer->SRE_Ycoord[k1] = (double)storeY; 
					  SRE.pointer->SRE_Zcoord[k1] = (double)storeZ; 
    
					  break;

				 } // if
				  /// convergence test (end)

				 // do not allow the SRE point outside the pixel (begin)
				 if ( (double)abs( (double)SRE.pointer->SRE_Xcoord[k1] ) >= (double)XPixelSize )		     
				 break;
     
				 if ( (double)abs( (double)SRE.pointer->SRE_Ycoord[k1] ) >= (double)YPixelSize )
				 break;

				 if ( (double)abs( (double)SRE.pointer->SRE_Zcoord[k1] ) >= (double)ZPixelSize ) 
				 break;
				 // do not allow the SRE point outside the pixel (end)
			

				 } // while loop
			
				 	} else if ( (double)SRE.pointer->fMRI[k1] <= TH ) { // major if

						SRE.pointer->SRE_Xcoord[k1] = (double)0.0;

						SRE.pointer->SRE_Ycoord[k1] = (double)0.0;
							
						SRE.pointer->SRE_Zcoord[k1] = (double)0.0;
			
					} // major if

			} // x dim
        
		}  // y dim 

	}  // z dim 
	// calculate SRE (end)

	// scale SRE (begin)
	//SRE.pointer->SRE_Xcoord (begin)
	max=-MAX;
	min=MAX;

	for (int sl=1; sl <= n4; sl++) { // z dim

		for (int i2=0; i2 < n2; i2++) { // y dim

			for (int i1=0; i1 < n1; i1++) {// x dim

			   int k1 = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)

			   if( SRE.pointer->SRE_Xcoord[k1] > (double)max ) max = (double)SRE.pointer->SRE_Xcoord[k1];
              
		       if( SRE.pointer->SRE_Xcoord[k1] < (double)min ) min = (double)SRE.pointer->SRE_Xcoord[k1];
			
			   
			} // x dim
        
		}  // y dim 

	}  // z dim 

	for (int sl=PAD; sl <= n4; sl++) { // z dim

		for (int i2=0; i2 < n2; i2++) { // y dim

			for (int i1=0; i1 < n1; i1++) {// x dim

			   int k1 = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)

			   if ( (double)SRE.pointer->fMRI[k1] > TH ) {  // major if

			   if ( max == min ) SRE.pointer->SRE_Xcoord[k1] = (double)0.0;

               else SRE.pointer->SRE_Xcoord[k1] = (double) Mx * (double) fabs ( (min - SRE.pointer->SRE_Xcoord[k1]) / (min - max) );
			   
			   } // major if

			} // x dim
        
		}  // y dim 

	}  // z dim 
	//SRE.pointer->SRE_Xcoord (end)


	//SRE.pointer->SRE_Ycoord (begin)
	max=-MAX;
	min=MAX;

	for (int sl=1; sl <= n4; sl++) { // z dim

		for (int i2=0; i2 < n2; i2++) { // y dim

			for (int i1=0; i1 < n1; i1++) {// x dim

			   int k1 = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)

			   if( SRE.pointer->SRE_Ycoord[k1] > (double)max ) max = (double)SRE.pointer->SRE_Ycoord[k1];
              
		       if( SRE.pointer->SRE_Ycoord[k1] < (double)min ) min = (double)SRE.pointer->SRE_Ycoord[k1];
			
			   
			} // x dim
        
		}  // y dim 

	}  // z dim 

	for (int sl=PAD; sl <= n4; sl++) { // z dim

		for (int i2=0; i2 < n2; i2++) { // y dim

			for (int i1=0; i1 < n1; i1++) {// x dim

			   int k1 = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)

			   if ( (double)SRE.pointer->fMRI[k1] > TH ) {  // major if

			   if ( max == min ) SRE.pointer->SRE_Ycoord[k1] = (double)0.0;

               else SRE.pointer->SRE_Ycoord[k1] = (double) My * (double) fabs ( (min - SRE.pointer->SRE_Ycoord[k1]) / (min - max) );
			   
			   } // major if

			} // x dim
        
		}  // y dim 

	}  // z dim 
	//SRE.pointer->SRE_Ycoord (end)

	//SRE.pointer->SRE_Zcoord (begin)
	max=-MAX;
	min=MAX;

	for (int sl=1; sl <= n4; sl++) { // z dim

		for (int i2=0; i2 < n2; i2++) { // y dim

			for (int i1=0; i1 < n1; i1++) {// x dim

			   int k1 = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)

			   if( SRE.pointer->SRE_Zcoord[k1] > (double)max ) max = (double)SRE.pointer->SRE_Zcoord[k1];
              
		       if( SRE.pointer->SRE_Zcoord[k1] < (double)min ) min = (double)SRE.pointer->SRE_Zcoord[k1];
			
			   
			} // x dim
        
		}  // y dim 

	}  // z dim 

	for (int sl=PAD; sl <= n4; sl++) { // z dim

		for (int i2=0; i2 < n2; i2++) { // y dim

			for (int i1=0; i1 < n1; i1++) {// x dim

			   int k1 = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)

			   if ( (double)SRE.pointer->fMRI[k1] > TH ) {  // major if

			   if ( max == min ) SRE.pointer->SRE_Zcoord[k1] = (double)0.0;

               else SRE.pointer->SRE_Zcoord[k1] = (double) Mz * (double) fabs ( (min - SRE.pointer->SRE_Zcoord[k1]) / (min - max) );
			   
			   } // major if

			} // x dim
        
		}  // y dim 

	}  // z dim 
	//SRE.pointer->SRE_Zcoord (end)
	// scale SRE (end)


	std::cout << "Computing Novel Re-sampling locations on the basis of the SRE..." << endl;

	// compute thetaXYZ, fxyz, F_XYZ, F_XY, F_XZ, F_YZ, on the basis of the SRE (begin)
	// & compute the novel re-sampling locations of Mx, My, Mz (begin)
	double theSumofThem, A, B, C;
	double eta1, eta2, eta3, eta4, eta5, eta6;
	double alpha1_star, thetax, thetay, thetaz;
	double beta1_star, beta2_star, beta3_star;
	double gamma1_star, gamma2_star, gamma3_star;
	double delta1_star, F_XYZ_star, epsilon;
	double alpha1, alpha2, alpha3, alpha4;
	

	for (int sl=PAD; sl <= n4; sl++) { // z dim

		for (int j=0; j < n2; j++) { // y dim

			for (int i=0; i < n1; i++) {// x dim

			 int k1 = (j*n1+i) + ((sl-2)*n1*n2); // f(0,0,0)

			   if ( (double)SRE.pointer->fMRI[k1] > TH ) {  // major if

	
				k2=((j+1)*rcxres+i)+( (sl-2)*rcyres*rcxres);     // f(0,1,0) 

                k3=(j*rcxres+(i+1))+( (sl-2)*rcyres*rcxres);     // f(1,0,0) 

                k4=((j+1)*rcxres+(i+1))+( (sl-2)*rcyres*rcxres); // f(1,1,0) 
            
                k5=(j*rcxres+i)+( (sl-1)*rcyres*rcxres);         // f(0,0,1)

                k6=((j+1)*rcxres+i)+( (sl-1)*rcyres*rcxres);     // f(0,1,1) 

                k7=(j*rcxres+(i+1))+( (sl-1)*rcyres*rcxres);     // f(1,0,1) 
  
                k8=((j+1)*rcxres+(i+1))+((sl-1)*rcyres*rcxres);  // f(1,1,1) 


	    	omega = ((double) SRE.pointer->scaledfMRI[k8] - SRE.pointer->scaledfMRI[k4] - SRE.pointer->scaledfMRI[k6] - 
			
			                  SRE.pointer->scaledfMRI[k7] + SRE.pointer->scaledfMRI[k5] + SRE.pointer->scaledfMRI[k2] + 
					
					          SRE.pointer->scaledfMRI[k3] - SRE.pointer->scaledfMRI[k1] );    
     
            theta1 = ( -((double) SRE.pointer->scaledfMRI[k1] - SRE.pointer->scaledfMRI[k3] - 
				
				                  SRE.pointer->scaledfMRI[k2] + SRE.pointer->scaledfMRI[k4] ) ); // thetaxy
 
            theta2 = ( -((double) SRE.pointer->scaledfMRI[k1] - SRE.pointer->scaledfMRI[k2] - 
				
				                  SRE.pointer->scaledfMRI[k5] + SRE.pointer->scaledfMRI[k6] ) ); // thetayz 

            theta3 = ( -((double) SRE.pointer->scaledfMRI[k1] - SRE.pointer->scaledfMRI[k3] - 
				
				                  SRE.pointer->scaledfMRI[k5] + SRE.pointer->scaledfMRI[k7] ) ); // thetaxz 


			if ( (double)omega == 0.0 ) // compute thetaXYZ
                               
            SRE.pointer->thetaXYZ[k1] = (double)0.0;
                
            else if ( (double)omega != 0.0 )    

			SRE.pointer->thetaXYZ[k1] = ((double)(theta1 + theta2 + theta3 ) / omega);


			// compute F_XYZ, F_XY, F_YZ, F_XZ (begin)

			alpha1 = ((double)SRE.pointer->SRE_Xcoord[k1] * SRE.pointer->SRE_Ycoord[k1] * SRE.pointer->SRE_Zcoord[k1]);

            alpha2 = ((double)SRE.pointer->SRE_Xcoord[k1] * SRE.pointer->SRE_Ycoord[k1]);

            alpha3 = ((double)SRE.pointer->SRE_Xcoord[k1] * SRE.pointer->SRE_Zcoord[k1]);

			alpha4 = ((double)SRE.pointer->SRE_Ycoord[k1] * SRE.pointer->SRE_Zcoord[k1]);


			
			beta1 = ((double)SRE.pointer->SRE_Xcoord[k1]/2.0) * 
				    ((double)SRE.pointer->scaledfMRI[k3] - SRE.pointer->scaledfMRI[k1]);  // thetax

            beta2 = ((double)SRE.pointer->SRE_Ycoord[k1]/2.0) * 
				    ((double)SRE.pointer->scaledfMRI[k2] - SRE.pointer->scaledfMRI[k1]);  // thetay

            beta3 = ((double)SRE.pointer->SRE_Zcoord[k1]/2.0) * 
				    ((double) SRE.pointer->scaledfMRI[k5] - SRE.pointer->scaledfMRI[k1]); // thetaz



			gamma1 = ( ((double)(SRE.pointer->SRE_Xcoord[k1] * 
				                 SRE.pointer->SRE_Ycoord[k1]/4.0)) * (double)theta1 );

            gamma2 = ( ((double)(SRE.pointer->SRE_Xcoord[k1] * 
				                 SRE.pointer->SRE_Zcoord[k1]/4.0)) * (double)theta2 );

            gamma3 = ( ((double)(SRE.pointer->SRE_Ycoord[k1] * 
				                 SRE.pointer->SRE_Zcoord[k1]/4.0)) * (double)theta3 );  
			
			
			delta1 = ( ((double)(SRE.pointer->SRE_Xcoord[k1] * 
				                 SRE.pointer->SRE_Ycoord[k1] * 
								 SRE.pointer->SRE_Zcoord[k1]/6.0)) * (double)omega );


			theSumofThem = ((double)SRE.pointer->scaledfMRI[k1] + beta1 + beta2 + beta3 + gamma1 + gamma2 + gamma3 + delta1);

			SRE.pointer->F_XYZ[k1] = ((double)alpha1 * theSumofThem);
			
            SRE.pointer->F_XY[k1]  = ((double)alpha2 * theSumofThem);
 
            SRE.pointer->F_XZ[k1]  = ((double)alpha3 * theSumofThem);

			SRE.pointer->F_YZ[k1] =  ((double)alpha4 * theSumofThem);
			
			// compute F_XYZ, F_XY, F_YZ, F_XZ (end)


			// compute fxyz (begin)

            SRE.pointer->fxyz[k1] = ( ((double)SRE.pointer->scaledfMRI[k1]) + 
				                      ((double)2.0 * beta1) + ((double)2.0 * beta2) + 
									  ((double)2.0 * beta3) + ((double)4.0 * gamma1) + 
									  ((double)4.0 * gamma2) + ((double)4.0 * gamma3) + 
									  ((double)6.0 * delta1) );

			// compute fxyz (end)

			// Calculate the Intensity-Curvature Functional (begin) //

			double E0 = -((double)SRE.pointer->scaledfMRI[k1] * Mx * My * Mz * SRE.pointer->thetaXYZ[k1] * omega / 4.0);

			double thetaX = (double)SRE.pointer->scaledfMRI[k3] - (double)SRE.pointer->scaledfMRI[k1];
			double thetaY = (double)SRE.pointer->scaledfMRI[k2] - (double)SRE.pointer->scaledfMRI[k1];
			double thetaZ = (double)SRE.pointer->scaledfMRI[k5] - (double)SRE.pointer->scaledfMRI[k1];

			double Hxyz = (double)SRE.pointer->scaledfMRI[k1] * Mx * My * Mz + ((double) Mx * Mx * My * Mz * thetaX / 2.0);

			Hxyz += ((double) Mx * My * My * Mz * thetaY / 2.0);

			Hxyz += ((double) Mx * My * Mz * Mz * thetaZ / 2.0);

			Hxyz -= ((double) Mx * Mx * My * My * Mz * theta1 * omega / 4.0);

			Hxyz -= ((double) Mx * My * My * Mz * Mz * theta2 * omega / 4.0);

			Hxyz -= ((double) Mx * Mx * My * Mz * Mz * theta3 * omega / 4.0);

			Hxyz += ((double) Mx * Mx * My * My * Mz * Mz * omega / 8.0);



			double CHIxHxyz = ((double)SRE.pointer->scaledfMRI[k1] * Mx * My * Mz / 2.0) + ((double) Mx * Mx * My * Mz * thetaX / 6.0);

			CHIxHxyz += ((double) Mx * My * My * Mz * thetaY / 4.0);

			CHIxHxyz += ((double) Mx * My * Mz * Mz * thetaZ / 4.0);

			CHIxHxyz -= ((double) Mx * Mx * My * My * Mz * theta1 * omega / 12.0);

			CHIxHxyz -= ((double) Mx * My * My * Mz * Mz * theta2 * omega / 8.0);

			CHIxHxyz -= ((double) Mx * Mx * My * Mz * Mz * theta3 * omega / 12.0);

			CHIxHxyz += ((double) Mx * Mx * My * My * Mz * Mz * omega / 24.0);

			

			double CHIyHxyz = ((double)SRE.pointer->scaledfMRI[k1] * Mx * My * Mz / 2.0) + ((double) Mx * Mx * My * Mz * thetaX / 4.0);

			CHIyHxyz += ((double) Mx * My * My * Mz * thetaY / 6.0);

			CHIyHxyz += ((double) Mx * My * Mz * Mz * thetaZ / 4.0);

			CHIyHxyz -= ((double) Mx * Mx * My * My * Mz * theta1 * omega / 12.0);

			CHIyHxyz -= ((double) Mx * My * My * Mz * Mz * theta2 * omega / 12.0);

			CHIyHxyz -= ((double) Mx * Mx * My * Mz * Mz * theta3 * omega / 8.0);

			CHIyHxyz += ((double) Mx * Mx * My * My * Mz * Mz * omega / 24.0);


			double CHIzHxyz = ((double)SRE.pointer->scaledfMRI[k1] * Mx * My * Mz / 2.0) + ((double) Mx * Mx * My * Mz * thetaX / 4.0);

			CHIzHxyz += ((double) Mx * My * My * Mz * thetaY / 4.0);

			CHIzHxyz += ((double) Mx * My * Mz * Mz * thetaZ / 6.0);

			CHIzHxyz -= ((double) Mx * Mx * My * My * Mz * theta1 * omega / 8.0);

			CHIzHxyz -= ((double) Mx * My * My * Mz * Mz * theta2 * omega / 12.0);

			CHIzHxyz -= ((double) Mx * Mx * My * Mz * Mz * theta3 * omega / 12.0);

			CHIzHxyz += ((double) Mx * Mx * My * My * Mz * Mz * omega / 24.0);


			double EIN = -(double) omega * Hxyz * ( theta1 + theta2 + theta3);

			EIN += (double) omega * ((double)Mx * (Hxyz - CHIxHxyz) + My * (Hxyz - CHIyHxyz) + Mz * (Hxyz - CHIzHxyz));

			SRE.pointer->ICF[k1] = ((double)E0/EIN);

			// Calculate the Intensity-Curvature Functional (end) //

			// compute deDELTAdeX, deDELTAdeY, deDELTAdeZ, DELTA (begin)

			eta1 = ((double)SRE.pointer->SRE_Ycoord[k1] * SRE.pointer->SRE_Zcoord[k1] * 
				            SRE.pointer->thetaXYZ[k1] * SRE.pointer->scaledfMRI[k1] * SRE.pointer->F_XYZ[k1]);

            eta4 = ( (double)3.0 + ((double)SRE.pointer->thetaXYZ[k1] - 
				                            SRE.pointer->SRE_Xcoord[k1] - 
											SRE.pointer->SRE_Ycoord[k1] - 
											SRE.pointer->SRE_Zcoord[k1]) );

            eta2 = ((double)SRE.pointer->SRE_Xcoord[k1] * SRE.pointer->SRE_Ycoord[k1] * 
				            SRE.pointer->SRE_Zcoord[k1] * SRE.pointer->thetaXYZ[k1] * 
							SRE.pointer->scaledfMRI[k1]);
                
            eta3 = ((double)SRE.pointer->F_YZ[k1] * eta4);

			eta5 = ((double)SRE.pointer->F_XZ[k1] * eta4);

            eta6 = ((double)SRE.pointer->F_XY[k1] * eta4);


			if ( ((double)SRE.pointer->F_XYZ[k1] * eta4 ) == 0.0 )
			{ // if

				SRE.pointer->deDELTAdeX[k1] = (double)0.0;

                SRE.pointer->deDELTAdeY[k1] = (double)0.0;

                SRE.pointer->deDELTAdeZ[k1] = (double)0.0;

				SRE.pointer->DELTA[k1] = (double)0.0;

			} // if
			else if ( ((double)SRE.pointer->F_XYZ[k1] * eta4 ) != 0.0 )
			{ // if

			SRE.pointer->deDELTAdeX[k1] = ( ((double)eta1 * eta4) - ((double)eta2 * eta3) - (double)SRE.pointer->F_XYZ[k1] ) / 
					
					                        ((double)SRE.pointer->F_XYZ[k1] * SRE.pointer->F_XYZ[k1] * eta4 * eta4);

			SRE.pointer->deDELTAdeY[k1] = ( ((double)eta1 * eta4) - ((double)eta2 * eta5) - (double)SRE.pointer->F_XYZ[k1] ) / 
				
										    ((double)SRE.pointer->F_XYZ[k1] * SRE.pointer->F_XYZ[k1] * eta4 * eta4);

            SRE.pointer->deDELTAdeZ[k1] = ( ((double)eta1 * eta4) - ((double)eta2 * eta6) - (double)SRE.pointer->F_XYZ[k1] ) / 
				
											((double)SRE.pointer->F_XYZ[k1] * SRE.pointer->F_XYZ[k1] * eta4 * eta4 );
                
            SRE.pointer->DELTA[k1] = ((double)alpha1 * SRE.pointer->thetaXYZ[k1] * SRE.pointer->scaledfMRI[k1] ) / 
				
									 ((double)SRE.pointer->F_XYZ[k1] * eta4 );

			} // if
			// compute deDELTAdeX, deDELTAdeY, deDELTAdeZ, DELTA (end)

			// now computing the novel re-sampling locations (begin)

			thetax = ((double)SRE.pointer->scaledfMRI[k3] - SRE.pointer->scaledfMRI[k1]);

            thetay = ((double)SRE.pointer->scaledfMRI[k2] - SRE.pointer->scaledfMRI[k1]);

            thetaz = ((double)SRE.pointer->scaledfMRI[k5] - SRE.pointer->scaledfMRI[k1]);
               
            alpha1_star = ( ((double)Mx - SRE.pointer->SRE_Xcoord[k1]) * 
				            ((double)My - SRE.pointer->SRE_Ycoord[k1]) * 
						    ((double)Mz - SRE.pointer->SRE_Zcoord[k1]) );

			beta1_star = ( ((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) /2.0) ) * 
				
						   ((double)SRE.pointer->scaledfMRI[k3] - SRE.pointer->scaledfMRI[k1]);
 
            beta2_star = ( ((double)(My - SRE.pointer->SRE_Ycoord[k1])/2.0 ) ) * 
				
				           ((double)SRE.pointer->scaledfMRI[k2] - SRE.pointer->scaledfMRI[k1]); 
 
            beta3_star = ( ((double)(Mz-SRE.pointer->SRE_Zcoord[k1])/2.0 ) ) * 
				
				           ((double)SRE.pointer->scaledfMRI[k5] - SRE.pointer->scaledfMRI[k1]);


			gamma1_star = ( ((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * 
				                     (My - SRE.pointer->SRE_Ycoord[k1])/4.0 )) * (double)theta1;

            gamma2_star = ( ((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * 
				                     (Mz - SRE.pointer->SRE_Zcoord[k1])/4.0 )) * (double)theta2;

            gamma3_star = ( ((double)(My - SRE.pointer->SRE_Ycoord[k1]) * 
				                     (Mz - SRE.pointer->SRE_Zcoord[k1])/4.0 )) * (double)theta3;   

            delta1_star = ( ((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * 
				                     (My - SRE.pointer->SRE_Ycoord[k1]) * 
							         (Mz - SRE.pointer->SRE_Zcoord[k1])/6.0 )) * (double)omega ;

			theSumofThem = ((double)SRE.pointer->scaledfMRI[k1] + 
				                    beta1_star + beta2_star + beta3_star + 
							        gamma1_star + gamma2_star + gamma3_star + delta1_star);

			F_XYZ_star = ((double)alpha1_star * theSumofThem);
                
            epsilon = ((double)1.0 / ( 3.0 + theta1 + theta2 + theta3 - 
					    
				      SRE.pointer->SRE_Xcoord[k1] - SRE.pointer->SRE_Ycoord[k1] - SRE.pointer->SRE_Zcoord[k1] ) );
                
            B = (double)0.0;

            C = (double)0.0;

			A = (double)0.0;

			// compute B term of 2nd order equation to calculate Mx (begin)               

			B = ((double)SRE.pointer->scaledfMRI[k1] * omega) + ((double)thetax * thetaz) + 
				((double)(Mz - SRE.pointer->SRE_Zcoord[k1]) * thetax * omega) + 
				((double)thetax * thetax * thetay) + 
				((double)(My - SRE.pointer->SRE_Ycoord[k1]) * thetax * omega);

			B = B + ((double)(Mz - SRE.pointer->SRE_Zcoord[k1]) * 
				             (My - SRE.pointer->SRE_Ycoord[k1]) * omega * theta2); 

            B = B + ( ((double)(Mz - SRE.pointer->SRE_Zcoord[k1]) * (My - SRE.pointer->SRE_Ycoord[k1]) * omega) - 
				
				      ((double)(Mz - SRE.pointer->SRE_Zcoord[k1]) * theta3) - 
				
				      ((double)theta1 * (My - SRE.pointer->SRE_Ycoord[k1]) ) ) * 
				
				    (  (double)thetaz + ((double)(Mz - SRE.pointer->SRE_Zcoord[k1]) * omega) + 
				 
				        thetay + thetax + ((double)(My-SRE.pointer->SRE_Ycoord[k1]) * omega) );  
			
			B = B + ((double)(My - SRE.pointer->SRE_Ycoord[k1]) * omega * thetay) + 
				    ((double)(Mz - SRE.pointer->SRE_Zcoord[k1]) * omega * thetaz);


			if ( (double)SRE.pointer->F_XYZ[k1] == 0.0 ) B = (double)B;
			else 
			{ // if

				B = B + ( (double) ((double)F_XYZ_star / SRE.pointer->F_XYZ[k1]) * epsilon );
			
			} // if 
			// compute B term of 2nd order equation to calculate Mx (end)  
			
			// compute C term of 2nd order equation to calculate Mx (begin)

			C = ((double)SRE.pointer->scaledfMRI[k1] * (thetax + thetay + thetaz)) + 
				
				((double)SRE.pointer->scaledfMRI[k1] * thetaz) + 
				
			   ( (double)SRE.pointer->scaledfMRI[k1] * omega * 
				
				((double)( (Mz - SRE.pointer->SRE_Zcoord[k1]) + (My - SRE.pointer->SRE_Ycoord[k1]) ) )  );

            C = C + (  ((double) (My - SRE.pointer->SRE_Ycoord[k1])) * (double)thetay * 
				
				       ((double) ( thetax + thetay + thetaz + (Mz - SRE.pointer->SRE_Zcoord[k1]) * omega )) );

			C = C + ( ( ((double) (My - SRE.pointer->SRE_Ycoord[k1]) * (My - SRE.pointer->SRE_Ycoord[k1]) ) * omega * thetay ) + 

				        ((double)(Mz - SRE.pointer->SRE_Zcoord[k1]) * thetaz * 
					  
					             ( thetax + thetay + thetaz + (My - SRE.pointer->SRE_Ycoord[k1]) * omega ) ) );

			C = C + ( ((double)(Mz - SRE.pointer->SRE_Zcoord[k1]) * (Mz - SRE.pointer->SRE_Zcoord[k1]) * omega * thetaz) - 
				      ((double)(Mz - SRE.pointer->SRE_Zcoord[k1]) * (My - SRE.pointer->SRE_Ycoord[k1]) * theta2) *  
					  (thetax + thetay + thetaz + 
					  ((double)(My - SRE.pointer->SRE_Ycoord[k1])) * omega + 
					  ((double)(Mz - SRE.pointer->SRE_Zcoord[k1])) * omega ) );

			if ( (double)SRE.pointer->F_XYZ[k1] == 0.0 ) C = (double)C;
			else 
			{ // if

				 C = C - ( (double)F_XYZ_star / SRE.pointer->F_XYZ[k1] ) * epsilon * 
					      ((double) ( 3.0 + theta1 + theta2 + theta3 - 
					      ((double) My - SRE.pointer->SRE_Ycoord[k1] ) - 
					      ((double) Mz - SRE.pointer->SRE_Zcoord[k1] ) ) );
			} // if 
			// compute C term of 2nd order equation to calculate Mx (end)

			// compute A term of 2nd order equation to calculate Mx (begin) 

			A = ((double)thetax * omega) - ((double)(My - SRE.pointer->SRE_Ycoord[k1]) * theta1 * omega) - 			
				((double)(Mz - SRE.pointer->SRE_Zcoord[k1]) * theta3 * omega) + 			
				((double)omega * omega * 
				((double)(Mz - SRE.pointer->SRE_Zcoord[k1])) * 
				((double)(My - SRE.pointer->SRE_Ycoord[k1])) ); 

			// compute A term of 2nd order equation to calculate Mx (end) 

			// compute Novel re-sampling location of Mx (begin) 
			if ( ((double)A == 0.0) )
			{ // if
	
				SRE.pointer->Recomp_Mx_plus[k1] = (double)0.0;

				SRE.pointer->Recomp_Mx_minus[k1] = (double)0.0;

			} // if 
			else if ( ((double)(B*B - 4.0*A*C) >= 0) && ((double)A != 0.0) )
			{ // if

			SRE.pointer->Recomp_Mx_plus[k1] = ( (double) ( -(double)B + (double)sqrt( ((double)B*B - 4.0*A*C) ) ) / ((double)2.0*A) );

			SRE.pointer->Recomp_Mx_minus[k1] = ( (double) ( -(double)B - (double)sqrt( ((double)B*B - 4.0*A*C) ) ) / ((double)2.0*A) );

			} // if 
			else if ( ((double)(B*B - 4.0*A*C) < 0) && ((double)A != 0.0) )
			{ // if

			SRE.pointer->Recomp_Mx_plus[k1] = ( (double) ( -(double)B ) / ((double)2.0*A) );

			SRE.pointer->Recomp_Mx_minus[k1] = ( (double) ( -(double)B ) / ((double)2.0*A) );


			} // if 
			// compute Novel re-sampling location of Mx (end)

			B = (double)0.0;

            C = (double)0.0;

            A = (double)0.0;

			// compute B term of 2nd order equation to calculate My (begin) 
			B = ((double)SRE.pointer->scaledfMRI[k1] * omega) + 
			    ((double)(Mx - SRE.pointer->SRE_Xcoord[k1])) * thetax * omega + thetay * 
				((double)(thetaz + ((double)(Mz - SRE.pointer->SRE_Zcoord[k1])) * omega + 
				 thetay + thetax + ((double)(Mx - SRE.pointer->SRE_Xcoord[k1])) * omega )) + 
				((double)(Mz - SRE.pointer->SRE_Zcoord[k1]) * thetaz * omega);

			B = B - ((double) omega * theta1 * ((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * (Mx - SRE.pointer->SRE_Xcoord[k1])) ); 

			B = B + ( ((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * theta1)  - 
				      ((double)(Mz - SRE.pointer->SRE_Zcoord[k1]) * theta2) ) * 
				      ((double)thetaz + ((double)(Mz - SRE.pointer->SRE_Zcoord[k1]) * omega + thetay + thetax) );


			B = B + ((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * (Mz - SRE.pointer->SRE_Zcoord[k1]) ) * 				
				    ((double)( 2.0 * omega * theta2 - omega * theta3 + omega * omega * 
					((double)(Mx - SRE.pointer->SRE_Xcoord[k1])) ));

			B = B + (double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * (double)(Mz - SRE.pointer->SRE_Zcoord[k1]) * 				
				   ((double)(thetaz + ((double)(Mz - SRE.pointer->SRE_Zcoord[k1])) * omega + thetay + thetax) * omega);
			
			if ( (double)SRE.pointer->F_XYZ[k1] == 0.0 ) B = (double)B;
			else 
			{ // if
		
				B = B + ( (double)F_XYZ_star / SRE.pointer->F_XYZ[k1] ) * (double)epsilon;

			} // if
			// compute B term of 2nd order equation to calculate My (end) 

			// compute C term of 2nd order equation to calculate My (begin) 

			C = (double)SRE.pointer->scaledfMRI[k1] * ((double)thetax + thetay + thetaz) + 
				( (double)SRE.pointer->scaledfMRI[k1] * thetaz + (double)SRE.pointer->scaledfMRI[k1] * omega * 
			      ( ((double)(Mz - SRE.pointer->SRE_Zcoord[k1])) + ((double)(Mx - SRE.pointer->SRE_Xcoord[k1])) ) );

			C = C + ((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * thetax * thetaz) + 
				
				((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * (Mz - SRE.pointer->SRE_Zcoord[k1]) * thetax * omega) + 
				
				((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * thetax * (thetay + thetax)) + 
				
				((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * (Mx - SRE.pointer->SRE_Xcoord[k1])* omega * omega * thetax);


			C = C + (  ((double)(Mz - SRE.pointer->SRE_Zcoord[k1]) * thetaz) + 
				       ((double)(Mz - SRE.pointer->SRE_Zcoord[k1]) * thetaz) * 
					   ( thetax + thetay + thetaz + ((double)(Mx - SRE.pointer->SRE_Xcoord[k1])) * omega) );


			C = C + ((double) ( (Mz - SRE.pointer->SRE_Zcoord[k1]) * (Mz - SRE.pointer->SRE_Zcoord[k1]) ) * omega * thetaz) - 
				    ((double) (Mx - SRE.pointer->SRE_Xcoord[k1]) * (Mx - SRE.pointer->SRE_Xcoord[k1]) * (Mz - SRE.pointer->SRE_Zcoord[k1]) * theta3 * omega) - 
					((double) (Mx - SRE.pointer->SRE_Xcoord[k1]) * (Mz - SRE.pointer->SRE_Zcoord[k1]) * theta3 *
					  ( thetax + thetay + thetaz + ((double)(Mz - SRE.pointer->SRE_Zcoord[k1]) * omega) ) );

			if ( (double)SRE.pointer->F_XYZ[k1] == 0.0 ) C = (double)C;
			else 
			{ // if

			C = C - ( (double)F_XYZ_star / SRE.pointer->F_XYZ[k1] ) * (double)epsilon * 
				( (double) (3.0 + theta1 + theta2 + theta3 - 
				 ((double)(Mx - SRE.pointer->SRE_Xcoord[k1])) - ((double)(Mz - SRE.pointer->SRE_Zcoord[k1])) ));

			} // if
			// compute C term of 2nd order equation to calculate My (end) 
			
			// compute A term of 2nd order equation to calculate My (begin)

			A = ((double)thetay * omega) - ((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * theta1 * omega) - 	
				((double) (Mz - SRE.pointer->SRE_Zcoord[k1]) * theta2 * omega + omega * omega * 
				          (Mz - SRE.pointer->SRE_Zcoord[k1]) * (Mx - SRE.pointer->SRE_Xcoord[k1]) ); 

			// compute A term of 2nd order equation to calculate My (end)

			// compute Novel re-sampling location of My (begin)
			if ( ((double)A == 0.0) )
			{ // if
	
				SRE.pointer->Recomp_My_plus[k1] = (double)0.0;

				SRE.pointer->Recomp_My_minus[k1] = (double)0.0;

			} // if 
			else if ( ((double)(B*B - 4.0*A*C) >= 0) && ((double)A != 0.0) )
			{ // if

			SRE.pointer->Recomp_My_plus[k1] = ( (double) ( -(double)B + (double)sqrt( ((double)B*B - 4.0*A*C) ) ) / ((double)2.0*A) );

			SRE.pointer->Recomp_My_minus[k1] = ( (double) ( -(double)B - (double)sqrt( ((double)B*B - 4.0*A*C) ) ) / ((double)2.0*A) );

			} // if 
			else if ( ((double)(B*B - 4.0*A*C) < 0) && ((double)A != 0.0) )
			{ // if

			SRE.pointer->Recomp_My_plus[k1] = ( (double) ( -(double)B ) / ((double)2.0*A) );

			SRE.pointer->Recomp_My_minus[k1] = ( (double) ( -(double)B ) / ((double)2.0*A) );


			} // if 
			// compute Novel re-sampling location of My (end)


			B = (double)0.0;

            C = (double)0.0;

            A = (double)0.0;

			// compute B term of 2nd order equation to calculate Mz (begin) 

			B = ((double)SRE.pointer->scaledfMRI[k1] * omega) + 			
				((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * thetax * omega) + 
				((double)(My - SRE.pointer->SRE_Ycoord[k1]) * thetay * omega) + 				
				((double)thetaz * (thetaz + ((double)(My - SRE.pointer->SRE_Ycoord[k1]) * omega) + 
				                   thetay + thetax + ((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * omega) ) );


			B = B - ( ((double)omega * theta1 * 
				      ((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * (My - SRE.pointer->SRE_Ycoord[k1])) ) ) - 
				      ((double)(My - SRE.pointer->SRE_Ycoord[k1])) * 
					   (double)theta2 * ((double)(thetaz + ((double)(My - SRE.pointer->SRE_Ycoord[k1]) * omega) + thetay + thetax)); 



			B = B + ((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * theta2 * (My - SRE.pointer->SRE_Ycoord[k1]) * omega ) - 
				    ((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * (Mx - SRE.pointer->SRE_Xcoord[k1]) * omega * theta3) - 
				    ((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * theta3 * 
					 (thetaz + ((double)(My - SRE.pointer->SRE_Ycoord[k1]) * omega) + thetay + thetax));


			B = B + ((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * (Mx - SRE.pointer->SRE_Xcoord[k1]) * 
				              omega * omega * (My - SRE.pointer->SRE_Ycoord[k1]) ) + 
				    ((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * (My - SRE.pointer->SRE_Ycoord[k1])) *
					((double)(thetaz + ((double)(My - SRE.pointer->SRE_Ycoord[k1]) * omega) + thetay + thetax) * omega);

			if ( (double)SRE.pointer->F_XYZ[k1] == 0.0 ) B = (double)B;
			else 
			{ // if

				 B = B + ( (double)F_XYZ_star / SRE.pointer->F_XYZ[k1] ) * (double)epsilon; 

			} // if
			// compute B term of 2nd order equation to calculate Mz (end) 


			// compute C term of 2nd order equation to calculate Mz (begin)

			C = ((double)SRE.pointer->scaledfMRI[k1] * (thetax + thetay + thetaz)) + 				
				((double)SRE.pointer->scaledfMRI[k1] * thetaz) +
				((double)SRE.pointer->scaledfMRI[k1] * omega * ( (My - SRE.pointer->SRE_Ycoord[k1]) + (Mx - SRE.pointer->SRE_Xcoord[k1]) ) );


			C = C + ((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * thetax * thetaz) + 
				    ((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * (thetax + thetay) * thetax) + 
					((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * (My - SRE.pointer->SRE_Ycoord[k1]) * thetax * omega) + 
					((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * (Mx - SRE.pointer->SRE_Xcoord[k1]) * thetax * omega);


			C = C + ((double)(My - SRE.pointer->SRE_Ycoord[k1])) * thetay * 
				    ((double)( thetax + thetay + thetaz + ((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * omega)) );


			C = C + ((double)( (My - SRE.pointer->SRE_Ycoord[k1]) * (My - SRE.pointer->SRE_Ycoord[k1]) ) * omega * thetay) - 
				    ((double)  (Mx - SRE.pointer->SRE_Xcoord[k1]) * (Mx - SRE.pointer->SRE_Xcoord[k1]) ) * 
					((double)  (My - SRE.pointer->SRE_Ycoord[k1]) * theta1 * omega) - 
					((double)  (Mx - SRE.pointer->SRE_Xcoord[k1]) * (My - SRE.pointer->SRE_Ycoord[k1]) * theta1 * 
					           ( thetax + thetay + thetaz + ((double)(My - SRE.pointer->SRE_Ycoord[k1]) * omega) ) );

			if ( (double)SRE.pointer->F_XYZ[k1] == 0.0 ) C = (double)C;
			else 
			{ // if
			 
				 C = C - ( (double)F_XYZ_star / SRE.pointer->F_XYZ[k1] ) * (double)epsilon * 
					      ((double)(3.0 + theta1 + theta2 + theta3 - 
						  ((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) - (My - SRE.pointer->SRE_Ycoord[k1])) ));
			} // if
			// compute C term of 2nd order equation to calculate Mz (end)


			// compute A term of 2nd order equation to calculate Mz (begin)

			A = ((double)thetaz * omega) - ((double)(My - SRE.pointer->SRE_Ycoord[k1]) * theta2 * omega) - 
				((double)(Mx - SRE.pointer->SRE_Xcoord[k1]) * theta3 * omega) + 
				((double) omega * omega * (Mx - SRE.pointer->SRE_Xcoord[k1]) * (My - SRE.pointer->SRE_Ycoord[k1]) ); 

			// compute A term of 2nd order equation to calculate Mz (end)

			// compute Novel re-sampling location of Mz (begin)
			if ( ((double)A == 0.0) )
			{ // if
	
				SRE.pointer->Recomp_Mz_plus[k1] = (double)0.0;

				SRE.pointer->Recomp_Mz_minus[k1] = (double)0.0;

			} // if 
			else if ( ((double)(B*B - 4.0*A*C) >= 0) && ((double)A != 0.0) )
			{ // if

			SRE.pointer->Recomp_Mz_plus[k1] = ( (double) ( -(double)B + (double)sqrt( ((double)B*B - 4.0*A*C) ) ) / ((double)2.0*A) );

			SRE.pointer->Recomp_Mz_minus[k1] = ( (double) ( -(double)B - (double)sqrt( ((double)B*B - 4.0*A*C) ) ) / ((double)2.0*A) );

			} // if 
			else if ( ((double)(B*B - 4.0*A*C) < 0) && ((double)A != 0.0) )
			{ // if

			SRE.pointer->Recomp_Mz_plus[k1] = ( (double) ( -(double)B ) / ((double)2.0*A) );

			SRE.pointer->Recomp_Mz_minus[k1] = ( (double) ( -(double)B ) / ((double)2.0*A) );


			} // if 
			// compute Novel re-sampling location of Mz (end)
			// now computing the novel re-sampling locations (end)



			   } else if ( (double)SRE.pointer->fMRI[k1] <= TH ) { // major if

						SRE.pointer->Recomp_Mx_plus[k1] = (double)0.0;

						SRE.pointer->Recomp_Mx_minus[k1] = (double)0.0;

						SRE.pointer->Recomp_My_plus[k1] = (double)0.0;

						SRE.pointer->Recomp_My_minus[k1] = (double)0.0;

						SRE.pointer->Recomp_Mz_plus[k1] = (double)0.0;

						SRE.pointer->Recomp_Mz_minus[k1] = (double)0.0;
			
			   } // major if

			} // x dim
        
		}  // y dim 

	}  // z dim 
	// compute thetaXYZ, fxyz, F_XYZ, F_XY, F_XZ, F_YZ, on the basis of the SRE (end)
	// & compute the novel re-sampling locations of Mx, My, Mz (end)



	// scale the novel re-sampling locations of Mx, My, Mz (begin)
	if ( REVERTEX == 1 ) 
	{ // if REVERTEX == 1

	//SRE.pointer->Recomp_Mx_plus (begin)
	max=-MAX;
	min=MAX;

	for (int sl=1; sl <= n4; sl++) { // z dim

		for (int i2=0; i2 < n2; i2++) { // y dim

			for (int i1=0; i1 < n1; i1++) {// x dim

			   int k1 = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)

			   if( SRE.pointer->Recomp_Mx_plus[k1] > (double)max ) max = (double)SRE.pointer->Recomp_Mx_plus[k1];
              
		       if( SRE.pointer->Recomp_Mx_plus[k1] < (double)min ) min = (double)SRE.pointer->Recomp_Mx_plus[k1];
			
			   
			} // x dim
        
		}  // y dim 

	}  // z dim 

	for (int sl=PAD; sl <= n4; sl++) { // z dim

		for (int i2=0; i2 < n2; i2++) { // y dim

			for (int i1=0; i1 < n1; i1++) {// x dim

			   int k1 = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)

			   if ( (double)SRE.pointer->fMRI[k1] > TH ) {  // major if

			   SRE.pointer->X_recomputed[k1] = (double)SRE.pointer->Recomp_Mx_plus[k1]; 

			   if ( max == min ) SRE.pointer->Recomp_Mx_plus[k1] = (double)0.0;

               else SRE.pointer->Recomp_Mx_plus[k1] = (double) XPixelSize * (double) fabs ( (min - SRE.pointer->Recomp_Mx_plus[k1]) / (min - max) );
			   
			   SRE.pointer->X_recomputed_scaled[k1] = (double)SRE.pointer->Recomp_Mx_plus[k1];

			   } // major if

			} // x dim
        
		}  // y dim 

	}  // z dim 
	//SRE.pointer->Recomp_Mx_plus (end)

	//SRE.pointer->Recomp_My_plus (begin)
	max=-MAX;
	min=MAX;

	for (int sl=1; sl <= n4; sl++) { // z dim

		for (int i2=0; i2 < n2; i2++) { // y dim

			for (int i1=0; i1 < n1; i1++) {// x dim

			   int k1 = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)

			   if( SRE.pointer->Recomp_My_plus[k1] > (double)max ) max = (double)SRE.pointer->Recomp_My_plus[k1];
              
		       if( SRE.pointer->Recomp_My_plus[k1] < (double)min ) min = (double)SRE.pointer->Recomp_My_plus[k1];
			
			   
			} // x dim
        
		}  // y dim 

	}  // z dim 

	for (int sl=PAD; sl <= n4; sl++) { // z dim

		for (int i2=0; i2 < n2; i2++) { // y dim

			for (int i1=0; i1 < n1; i1++) {// x dim

			   int k1 = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)

			   if ( (double)SRE.pointer->fMRI[k1] > TH ) {  // major if

			   SRE.pointer->Y_recomputed[k1] = (double)SRE.pointer->Recomp_My_plus[k1];

			   if ( max == min ) SRE.pointer->Recomp_My_plus[k1] = (double)0.0;

               else SRE.pointer->Recomp_My_plus[k1] = (double) YPixelSize * (double) fabs ( (min - SRE.pointer->Recomp_My_plus[k1]) / (min - max) );
			   
			   SRE.pointer->Y_recomputed_scaled[k1] = (double)SRE.pointer->Recomp_My_plus[k1];

			   } // major if

			} // x dim
        
		}  // y dim 

	}  // z dim 
	//SRE.pointer->Recomp_My_plus (end)

	//SRE.pointer->Recomp_Mz_plus (begin)
	max=-MAX;
	min=MAX;

	for (int sl=1; sl <= n4; sl++) { // z dim

		for (int i2=0; i2 < n2; i2++) { // y dim

			for (int i1=0; i1 < n1; i1++) {// x dim

			   int k1 = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)

			   if( SRE.pointer->Recomp_Mz_plus[k1] > (double)max ) max = (double)SRE.pointer->Recomp_Mz_plus[k1];
              
		       if( SRE.pointer->Recomp_Mz_plus[k1] < (double)min ) min = (double)SRE.pointer->Recomp_Mz_plus[k1];
			
			   
			} // x dim
        
		}  // y dim 

	}  // z dim 


	for (int sl=PAD; sl <= n4; sl++) { // z dim

		for (int i2=0; i2 < n2; i2++) { // y dim

			for (int i1=0; i1 < n1; i1++) {// x dim

			   int k1 = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)

			   if ( (double)SRE.pointer->fMRI[k1] > TH ) {  // major if

			   SRE.pointer->Z_recomputed[k1] = (double)SRE.pointer->Recomp_Mz_plus[k1];

			   if ( max == min ) SRE.pointer->Recomp_Mz_plus[k1] = (double)0.0;

               else SRE.pointer->Recomp_Mz_plus[k1] = (double) ZPixelSize * (double) fabs ( (min - SRE.pointer->Recomp_Mz_plus[k1]) / (min - max) );
			   
			   SRE.pointer->Z_recomputed_scaled[k1] = (double)SRE.pointer->Recomp_Mz_plus[k1];

			   } // major if

			} // x dim
        
		}  // y dim 

	}  // z dim 
	//SRE.pointer->Recomp_Mz_plus (end)

	} // if REVERTEX == 1
	else if ( REVERTEX == 0 )
	{ // if REVERTEX == 0

	//SRE.pointer->Recomp_Mx_minus (begin)
	max=-MAX;
	min=MAX;


	for (int sl=1; sl <= n4; sl++) { // z dim

		for (int i2=0; i2 < n2; i2++) { // y dim

			for (int i1=0; i1 < n1; i1++) {// x dim

			   int k1 = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)

			   if( SRE.pointer->Recomp_Mx_minus[k1] > (double)max ) max = (double)SRE.pointer->Recomp_Mx_minus[k1];
              
		       if( SRE.pointer->Recomp_Mx_minus[k1] < (double)min ) min = (double)SRE.pointer->Recomp_Mx_minus[k1];
			
			   
			} // x dim
        
		}  // y dim 

	}  // z dim 


	for (int sl=PAD; sl <= n4; sl++) { // z dim

		for (int i2=0; i2 < n2; i2++) { // y dim

			for (int i1=0; i1 < n1; i1++) {// x dim

			   int k1 = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)

			   if ( (double)SRE.pointer->fMRI[k1] > TH ) {  // major if

			   SRE.pointer->X_recomputed[k1] = (double)SRE.pointer->Recomp_Mx_minus[k1];

			   if ( max == min ) SRE.pointer->Recomp_Mx_minus[k1] = (double)0.0;

               else SRE.pointer->Recomp_Mx_minus[k1] = (double) XPixelSize * (double) fabs ( (min - SRE.pointer->Recomp_Mx_minus[k1]) / (min - max) );
			   
			   SRE.pointer->X_recomputed_scaled[k1] = (double)SRE.pointer->Recomp_Mx_minus[k1];

			   } // major if

			} // x dim
        
		}  // y dim 

	}  // z dim 
	//SRE.pointer->Recomp_Mx_minus (end)

	//SRE.pointer->Recomp_My_minus (begin)
	max=-MAX;
	min=MAX;

	for (int sl=1; sl <= n4; sl++) { // z dim

		for (int i2=0; i2 < n2; i2++) { // y dim

			for (int i1=0; i1 < n1; i1++) {// x dim

			   int k1 = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)

			   if( SRE.pointer->Recomp_My_minus[k1] > (double)max ) max = (double)SRE.pointer->Recomp_My_minus[k1];
              
		       if( SRE.pointer->Recomp_My_minus[k1] < (double)min ) min = (double)SRE.pointer->Recomp_My_minus[k1];
			
			   
			} // x dim
        
		}  // y dim 

	}  // z dim 

	for (int sl=PAD; sl <= n4; sl++) { // z dim

		for (int i2=0; i2 < n2; i2++) { // y dim

			for (int i1=0; i1 < n1; i1++) {// x dim

			   int k1 = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)

			   if ( (double)SRE.pointer->fMRI[k1] > TH ) {  // major if

			   SRE.pointer->Y_recomputed[k1] = (double)SRE.pointer->Recomp_My_minus[k1];

			   if ( max == min ) SRE.pointer->Recomp_My_minus[k1] = (double)0.0;

               else SRE.pointer->Recomp_My_minus[k1] = (double) YPixelSize * (double) fabs ( (min - SRE.pointer->Recomp_My_minus[k1]) / (min - max) );
			   
			   SRE.pointer->Y_recomputed_scaled[k1] = (double)SRE.pointer->Recomp_My_minus[k1];

			   } // major if

			} // x dim
        
		}  // y dim 

	}  // z dim 
	//SRE.pointer->Recomp_My_minus (end)

	//SRE.pointer->Recomp_Mz_minus (begin)
	max=-MAX;
	min=MAX;

	for (int sl=1; sl <= n4; sl++) { // z dim

		for (int i2=0; i2 < n2; i2++) { // y dim

			for (int i1=0; i1 < n1; i1++) {// x dim

			   int k1 = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)

			   if( SRE.pointer->Recomp_Mz_minus[k1] > (double)max ) max = (double)SRE.pointer->Recomp_Mz_minus[k1];
              
		       if( SRE.pointer->Recomp_Mz_minus[k1] < (double)min ) min = (double)SRE.pointer->Recomp_Mz_minus[k1];
			
			   
			} // x dim
        
		}  // y dim 

	}  // z dim 

	for (int sl=PAD; sl <= n4; sl++) { // z dim

		for (int i2=0; i2 < n2; i2++) { // y dim

			for (int i1=0; i1 < n1; i1++) {// x dim

			   int k1 = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)

			   if ( (double)SRE.pointer->fMRI[k1] > TH ) {  // major if

			   SRE.pointer->Z_recomputed[k1] = (double)SRE.pointer->Recomp_Mz_minus[k1];

			   if ( max == min ) SRE.pointer->Recomp_Mz_minus[k1] = (double)0.0;

               else SRE.pointer->Recomp_Mz_minus[k1] = (double) ZPixelSize * (double) fabs ( (min - SRE.pointer->Recomp_Mz_minus[k1]) / (min - max) );
			   
			   SRE.pointer->Z_recomputed_scaled[k1] = (double)SRE.pointer->Recomp_Mz_minus[k1];

			   } // major if

			} // x dim
        
		}  // y dim 

	}  // z dim 
	//SRE.pointer->Recomp_Mz_minus (end)

	} // if REVERTEX == 0
	// scale the novel re-sampling locations of Mx, My, Mz (end)



	// shift forward, SRE -> calculate 'the_SRE_shifted_fMRI' (begin)
	for (int sl=PAD; sl <= n4; sl++) { // z dim

		for (int j=0; j < n2; j++) { // y dim

			for (int i=0; i < n1; i++) {// x dim

			int k1 = (j*n1+i) + ((sl-2)*n1*n2); // f(0,0,0)

				if ( (double)SRE.pointer->fMRI[k1] > TH ) {  // major if

        					   
                        k2=((j+1)*rcxres+i)+( (sl-2)*rcyres*rcxres); // f(0,1,0) 
			
                        k3=(j*rcxres+(i+1))+( (sl-2)*rcyres*rcxres); // f(1,0,0) 
  			           
                        k4=((j+1)*rcxres+(i+1))+( (sl-2)*rcyres*rcxres); // f(1,1,0) 
            
                        k5=(j*rcxres+i)+( (sl-1)*rcyres*rcxres); // f(0,0,1)
			
                        k6=((j+1)*rcxres+i)+( (sl-1)*rcyres*rcxres); // f(0,1,1) 
					   
                        k7=(j*rcxres+(i+1))+( (sl-1)*rcyres*rcxres); // f(1,0,1) 
  
                        k8=((j+1)*rcxres+(i+1))+((sl-1)*rcyres*rcxres); // f(1,1,1) 
						
              
            omega = ((double) SRE.pointer->scaledfMRI[k8] - SRE.pointer->scaledfMRI[k4] - SRE.pointer->scaledfMRI[k6] - 
			
			                  SRE.pointer->scaledfMRI[k7] + SRE.pointer->scaledfMRI[k5] + SRE.pointer->scaledfMRI[k2] + 
					
					          SRE.pointer->scaledfMRI[k3] - SRE.pointer->scaledfMRI[k1] );    
     
            theta1 = ( -((double) SRE.pointer->scaledfMRI[k1] - SRE.pointer->scaledfMRI[k3] - 
				
				                  SRE.pointer->scaledfMRI[k2] + SRE.pointer->scaledfMRI[k4] ) ); // thetaxy
 
            theta2 = ( -((double) SRE.pointer->scaledfMRI[k1] - SRE.pointer->scaledfMRI[k2] - 
				
				                  SRE.pointer->scaledfMRI[k5] + SRE.pointer->scaledfMRI[k6] ) ); // thetayz 

            theta3 = ( -((double) SRE.pointer->scaledfMRI[k1] - SRE.pointer->scaledfMRI[k3] - 
				
				                  SRE.pointer->scaledfMRI[k5] + SRE.pointer->scaledfMRI[k7] ) ); // thetaxz 

          
            beta1 = (double)SRE.pointer->X_recomputed_scaled[k1] * ((double) SRE.pointer->scaledfMRI[k3] - SRE.pointer->scaledfMRI[k1] );  // thetax

            beta2 = (double)SRE.pointer->Y_recomputed_scaled[k1] * ((double) SRE.pointer->scaledfMRI[k2] - SRE.pointer->scaledfMRI[k1] );  // thetay

            beta3 = (double)SRE.pointer->Z_recomputed_scaled[k1] * ((double) SRE.pointer->scaledfMRI[k5] - SRE.pointer->scaledfMRI[k1] );  // thetaz
               
            gamma1 = ((double) SRE.pointer->X_recomputed_scaled[k1] * SRE.pointer->Y_recomputed_scaled[k1] * theta1);

            gamma2 = ((double) SRE.pointer->X_recomputed_scaled[k1] * SRE.pointer->Z_recomputed_scaled[k1] * theta2);

            gamma3 = ((double) SRE.pointer->Y_recomputed_scaled[k1] * SRE.pointer->Z_recomputed_scaled[k1] * theta3);   

            delta1 = ((double) SRE.pointer->X_recomputed_scaled[k1] * SRE.pointer->Y_recomputed_scaled[k1] * SRE.pointer->Z_recomputed_scaled[k1] * omega);
                
                
            SRE.pointer->the_SRE_shifted_fMRI[k1] = (double)SRE.pointer->fMRI[k1] + 
			                                       ((double)beta1 + beta2 + beta3 + gamma1 + gamma2 + gamma3 + delta1);
              
				} else if ( (double)SRE.pointer->fMRI[k1] <= TH ) { // major if

					SRE.pointer->the_SRE_shifted_fMRI[k1] = (double)0.0;
			
				} // major if


 			} // x dim
        
		}  // y dim 

	}  // z dim 
	// shift forward, SRE -> calculate 'the_SRE_shifted_fMRI' (end)


	// calculate SE_SRE (begin)
	double SE_SRE_sum = (double)0.0;
	counter = 0;

	for (int sl=PAD; sl <= n4; sl++) { // z dim

		for (int i2=0; i2 < n2; i2++) { // y dim

			for (int i1=0; i1 < n1; i1++) {// x dim

			int k1 = (i2*n1+i1) + ((sl-1)*n1*n2); // f(0,0,0)

				if ( (double)SRE.pointer->fMRI[k1] > TH ) {  // major if

				counter++;

				SE_SRE_sum += (double)sqrt( ((double)(SRE.pointer->fMRI[k1] - SRE.pointer->the_SRE_shifted_fMRI[k1] ) * 
					                                 (SRE.pointer->fMRI[k1] - SRE.pointer->the_SRE_shifted_fMRI[k1] ) ) ); 

				SRE.pointer->SE_SRE[k1] = (double)sqrt( ((double)(SRE.pointer->fMRI[k1] - SRE.pointer->the_SRE_shifted_fMRI[k1] ) * 
					                                             (SRE.pointer->fMRI[k1] - SRE.pointer->the_SRE_shifted_fMRI[k1] ) ) ); 


				}  // major if                                                         

			} // y dim
        
		}  // x dim 

	}  // z dim 
	// calculate SE_SRE (end)


	std::cout << "SE with SRE is: " << ((double)SE_SRE_sum / (counter))  << endl;
	std::cout << endl;

	fprintf(savedata,"%s%e\n", "SE with SRE is: ", ((double)SE_SRE_sum / (counter)) );
	fprintf(savedata,"\n");

	SRE.save();

	std::cout << "End of Computation..." << endl;
	std::cout << endl;

	fprintf(savedata,"%s\n", "End of Computation...");
	fprintf(savedata,"\n");

	if ( tolower(Fyn) == 'y' ) {  // run Power Spectra Analysis (begin)
	
	/// FFT analysis (begin)
	std::cout << "Now Performing FFT Analysis and Saving..." << endl;
	std::cout << endl;

	fprintf(savedata,"%s\n", "Now Performing FFT Analysis and Saving...");
	fprintf(savedata,"\n");

	FILE *image;
	char filenameToProcess[128];
	double SaveTheData;


	for ( int sl = PAD; sl <= n4 - PAD; sl++ ) 
	{ // for loop (begin)

	int imageNumber = (int)sl;

	sprintf(filenameToProcess, "%s%d%s", "fMRI-", imageNumber , ".img");
	/// save an image slice (begin)
	if ((image= fopen(filenameToProcess,"wb"))==NULL)
	{
	 
	std::cout << "Cannot Open File: " << filenameToProcess << " Now Exit. " << endl;
    fprintf(savedata, "%s%s%s\n", "Cannot Open File: ", filenameToProcess, " Now Exit. ");
    fclose(savedata);
	delete SRE.pointer;
	SRE.~SRE3D2013();
		
	}else{ // save (begin)

	for (int j=0; j < rcyres; j++) {// y dim
       	
			for (int i=0; i < rcxres; i++) { // x dim

				k1=(j*rcxres+i)+( (sl-1)*rcyres*rcxres); // f(0,0,0)

				SaveTheData = (double)SRE.pointer->fMRI[k1];
          
				fwrite(&SaveTheData,sizeof(double),1,image);


			} // x dim

	} // y dim
	fclose(image);
	/// save an image slice (end)
	} // save (end)


	SRE.OnFourierAnalysis(n1, n2, filenameToProcess); 
	// FFT Analysis of the original image: "fMRI.img"


	sprintf(filenameToProcess, "%s%d%s", "the_SRE_shifted_fMRI-", imageNumber, ".img");
	/// save an image slice (begin)
	if ((image= fopen(filenameToProcess,"wb"))==NULL)
	{
	 
	std::cout << "Cannot Open File: " << filenameToProcess << " Now Exit. " << endl;
    fprintf(savedata, "%s%s%s\n", "Cannot Open File: ", filenameToProcess, " Now Exit. ");
	fclose(savedata);
	delete SRE.pointer;
	SRE.~SRE3D2013();
	
	}else{ // save (begin)

	for (int j=0; j < rcyres; j++) {// y dim
       	
			for (int i=0; i < rcxres; i++) { // x dim

				k1=(j*rcxres+i)+( (sl-1)*rcyres*rcxres); // f(0,0,0)

				SaveTheData = (double)SRE.pointer->the_SRE_shifted_fMRI[k1];
          
				fwrite(&SaveTheData,sizeof(double),1,image);


			} // x dim

	} // y dim
	fclose(image);
	/// save an image slice (end)
	} // save (end)


	SRE.OnFourierAnalysis(n1, n2, filenameToProcess);
	// FFT Analysis of the SRE IMAGE: "the_SRE_shifted_fMRI.img"

	sprintf(filenameToProcess, "%s%d%s", "shifted_fMRI-", imageNumber, ".img");
	/// save an image slice (begin)
	if ((image= fopen(filenameToProcess,"wb"))==NULL)
	{
	 
	std::cout << "Cannot Open File: " << filenameToProcess << " Now Exit. " << endl;
    fprintf(savedata, "%s%s%s\n", "Cannot Open File: ", filenameToProcess, " Now Exit. ");
	fclose(savedata);
	delete SRE.pointer;
	SRE.~SRE3D2013();
	
	}else{ // save (begin)

	for (int j=0; j < rcyres; j++) {// y dim
       	
			for (int i=0; i < rcxres; i++) { // x dim

				k1=(j*rcxres+i)+( (sl-1)*rcyres*rcxres); // f(0,0,0)

				SaveTheData = (double)SRE.pointer->shifted_fMRI[k1];
          
				fwrite(&SaveTheData,sizeof(double),1,image);


			} // x dim

	} // y dim
	fclose(image);
	/// save an image slice (end)
	} // save (end)


	SRE.OnFourierAnalysis(n1, n2, filenameToProcess);
	// FFT Analysis of the NO SRE IMAGE: "shifted_fMRI.img"
	
	} // for loop (end)

	} // run Power Spectra Analysis (end)
	/// FFT analysis (end)

	fclose(savedata);
	delete SRE.pointer;
	SRE.~SRE3D2013();
	} // processing (end)

	} // run the program (end)

	system("java processCompleted");
	
	return 0;
} // end of main 