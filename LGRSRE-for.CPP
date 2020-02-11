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

//  Project Title: Cubic Lagrange (Cubic Lagrange polynomial SRE-based Interpolation)

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
#define SREPLUS 0       // positive solution for the SRE, otherwise set 0 for negative solution 
#define NOVELPLUS 0     // positive solution for the NOVEL's third degree equation, otherwise set 0 for negative solution 
#define NOVEL2PLUS 0    // positive solution for the NOVEL's fourth degree equation, otherwise set 0 for negative solution 
#define SQRTNOVELPLUS 0 // positive solution for the NOVEL's final equation, otherwise set 0 for negative solution 


#define TH 0

using namespace std;


class LGRSRE2013 {

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

		double **theta; // pointer to the matrix entry
		
		double **theta1; // pointer to the matrix entry

		double **theta2; // pointer to the matrix entry

		double **theta3; // pointer to the matrix entry

		double **theta4; // pointer to the matrix entry 

		double **theta5; // pointer to the matrix entry

		double **h_LGR1_recomp; // pointer to the matrix entry

		double **h_LGR2_recomp; // pointer to the matrix entry

		double **ICF; // pointer to the matrix entry

	}*pointer; // pointer to the matrices

public:

	LGRSRE2013(int x, int y) : n1(x), n2(y) { };
	
	void allocateData();

	~LGRSRE2013() { } // destructor

};

void LGRSRE2013::allocateData() { // allocate data

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

     pointer->theta = new double*[this->n1];

	 pointer->theta1 = new double*[this->n1];

	 pointer->theta2 = new double*[this->n1];

	 pointer->theta3 = new double*[this->n1];

	 pointer->theta4 = new double*[this->n1];

	 pointer->theta5 = new double*[this->n1];

	 pointer->h_LGR1_recomp = new double*[this->n1];

	 pointer->h_LGR2_recomp = new double*[this->n1];

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

		 pointer->theta[v] = new double[this->n2];

		 pointer->theta1[v] = new double[this->n2];

		 pointer->theta2[v] = new double[this->n2];

		 pointer->theta3[v] = new double[this->n2];

		 pointer->theta4[v] = new double[this->n2];

		 pointer->theta5[v] = new double[this->n2];

		 pointer->h_LGR1_recomp[v] = new double[this->n2];

		 pointer->h_LGR2_recomp[v] = new double[this->n2];

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

			pointer->theta[v][f] = (double)0.0;

		    pointer->theta1[v][f] = (double)0.0;

		    pointer->theta2[v][f] = (double)0.0;

		    pointer->theta3[v][f] = (double)0.0;

		    pointer->theta4[v][f] = (double)0.0;

		    pointer->theta5[v][f] = (double)0.0;

			pointer->h_LGR1_recomp[v][f] = (double)0.0;

			pointer->h_LGR2_recomp[v][f] = (double)0.0;

			pointer->ICF[v][f] = (double)0.0;


			 } //(b)

		 } //(a)
		// (2) initialize (end)

} // allocate data


int main ( int argc, char * argv[] ) {

	char outputFile[128]="LGRSRE-for.log";

	FILE * savedata;

	double MAX = 5000000000000000000.0;

if (argc < 7) { std::cout << endl;
				std::cout << "Please type the image file name" << endl;
				std::cout << "Please make sure that the image format is Analyze 'double': 64 bits real" << endl;
				std::cout << "Before running the program, please make sure that the image is padded of 'n >= 3'";
				std::cout << "number of pixels along X and Y" << endl;
				std::cout << "Please enter the number of pixels along the X direction (integer)" << endl;
				std::cout << "Please enter the number of pixels along the Y direction (integer)" << endl;
				std::cout << "Please enter the pixel size along the X direction (double)" << endl;
				std::cout << "Please enter the pixel size along the Y direction (double)" << endl;
				std::cout << "Please enter that 'n >= 3' number of pixels along X and Y which will pad the image" << endl;
				std::cout << "Some data will be lost if you enter a value of n which is greater than" << endl;
				std::cout << "the actual number of pixels along X and Y which was used to pad the image" << endl;
				std::cout << "before running the program." << endl;
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

	char imageFileName[128];

	int PAD = atoi(argv[6]);

	if ( PAD < 3 ) { // if test PAD
	
		std::cout << "Please make sure to pad the image with 'n >= 3' number of pixels along X and Y" << endl;
	
		fprintf(savedata,"%s\n", "Please make sure to pad the image with 'n >= 3' number of pixels along X and Y");
	
		fclose(savedata);
	
		exit(0);
	
	} // if test PAD

	int XNEI = (int)2;
	sprintf(imageFileName, "%s", argv[1]);

	std::cout << endl;
	std::cout << "The image file name is: " << imageFileName << endl;
	std::cout << "The number of pixels along the X direction is: " << atoi(argv[2]) << endl;
	std::cout << "The number of pixels along the Y direction is: " << atoi(argv[3]) << endl;
	std::cout << "The pixel size along the X direction is: " << atof(argv[4]) << endl;
	std::cout << "The pixel size along the Y direction is: " << atof(argv[5]) << endl;
	std::cout << "The 'n' number of pixels along X and Y which you entered to pad the image is: " << atoi(argv[6]) << endl;
	std::cout << endl;

	fprintf(savedata,"%s%s\n", "The image file name is: " , imageFileName);
	fprintf(savedata,"%s%d\n", "The number of pixels along the X direction is: ", n1);
	fprintf(savedata,"%s%d\n", "The number of pixels along the Y direction is: ", n2);
	fprintf(savedata,"%s%lf\n", "The pixel size along the X direction is: ", XPixelSize);
	fprintf(savedata,"%s%lf\n", "The pixel size along the Y direction is: ", YPixelSize);
	fprintf(savedata,"%s%d\n", "The 'n' number of pixels along X and Y which you entered to pad the image is: ", PAD);
	fprintf(savedata,"\n");

	double STEP_A = (double)0.005;
	double STEP_MX = (double)0.35;

	double x_misplacement_X = (double)0.4;
	double theta = (double)0.01;


	LGRSRE2013 SRE(n1,n2);

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

	std::cout << "Image data scaled" << endl;
	
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
			

	fprintf(savedata, "%s\t %s\t\t %s\t\t %s\n", "misplacement_X", "SE", "SRE SE", "SE Difference");  
	fprintf(savedata,"\n");
	std::cout << "Now Computing..." << endl; 

	 //---> major for loop (begin) <---//
	for (int step = 1; step <= 100; step++) {

	double misplacement_X = ((double)1.0 - ( cos( (double)theta + STEP_A) + sin( (double)theta + STEP_A) ) + (double)x_misplacement_X + (double)STEP_MX);
    
	misplacement_X = ((double)misplacement_X/XPixelSize);

	
	  //////////////////***********//////////////////////
	  // Above formula scales the misplacement to the  //
	  // pixel size the same way the following formula //
	  // would do: (min - misplacement)/(min - max)    //  
	  //////////////////***********//////////////////////
	
	// compute thetas (begin)
	double add1, add2;

   for (int i1=PAD+1; i1 < n1-PAD-1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim


		if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) { 


		if ( i1 == (int)firstLine || i1 == ((int)firstLine + 1) || i1 == ((int)firstLine + 2) )
		{ // calculate the first three horizonatal lines of the signal (begin)


		add2 = ((double) ( SRE.pointer->the_fMRI[i1+2+XNEI][i2] + SRE.pointer->the_fMRI[i1+1+XNEI][i2] + 
			               SRE.pointer->the_fMRI[i1-1+XNEI][i2] + SRE.pointer->the_fMRI[i1-2+XNEI][i2] ));       
        
		add1 = ((double)   SRE.pointer->the_fMRI[i1+1+XNEI][i2] + SRE.pointer->the_fMRI[i1-1+XNEI][i2] );       
            		       

		// calculate the first three horizonatal lines of the signal (end)
		} else 

		if ( i1 == ((int)lastLine) ||  i1 == ((int)lastLine - 1) || i1 == ((int)lastLine - 2) )
		{ // calculate the last three horizonatal lines of the signal (begin)

	
		add2 = ((double) ( SRE.pointer->the_fMRI[i1+2-XNEI][i2] + SRE.pointer->the_fMRI[i1+1-XNEI][i2] + 
			               SRE.pointer->the_fMRI[i1-1-XNEI][i2] + SRE.pointer->the_fMRI[i1-2-XNEI][i2] ));       
        
		add1 = ((double)   SRE.pointer->the_fMRI[i1+1-XNEI][i2] + SRE.pointer->the_fMRI[i1-1-XNEI][i2] );       
    
			
		} else { // calculate the rest of the lines of the signal (begin)

	
		add2 = ((double) ( SRE.pointer->the_fMRI[i1+2][i2] + SRE.pointer->the_fMRI[i1+1][i2] + 
			               SRE.pointer->the_fMRI[i1-1][i2] + SRE.pointer->the_fMRI[i1-2][i2] ));       
        
		add1 = ((double)   SRE.pointer->the_fMRI[i1+1][i2] + SRE.pointer->the_fMRI[i1-1][i2] );       
    

		}// calculate the rest of the lines of the signal (end)


		SRE.pointer->theta[i1][i2] =  ((double) ( -2.0 * add1 ) + ( 2.0 * add2 ) );

		SRE.pointer->theta1[i1][i2] = ((double) ( -2.0 * add1 * SRE.pointer->fMRI[i1][i2] ) + ( 2.0 * add2 * SRE.pointer->fMRI[i1][i2] ) ); 
		
		SRE.pointer->theta2[i1][i2] = ((double) ( (3.0/2.0) * add1 * SRE.pointer->fMRI[i1][i2] ) - ( 0.5 * add2 * SRE.pointer->fMRI[i1][i2] ) + ( 2.0 * add1 * add1 ) - ( add1 * add2 ) - ( (20.0/9.0) * add1 * add2 ) - ( (7.0/3.0) * add2 * add2 ) );

		SRE.pointer->theta3[i1][i2] = ((double) ( (1.0/6.0) * add1 * add1 ) - ( (1.0/2.0) * add1 * add2 ) - ( (15.0/6.0) * add1 * add2 ) + ( (23.0/18.0) * add2 * add2 ) );

		SRE.pointer->theta4[i1][i2] = ((double) ( - add1 * add1 ) + ( (1.0/2.0) * add1 * add2 ) + ( (5.0/6.0) * add1 * add2 ) - ( (1.0/3.0) * add2 * add2 ) );

        SRE.pointer->theta5[i1][i2] = ((double) ( (3.0/10.0) * add1 * add1 ) - ( (1.0/5.0) * add1 * add2 ) + ( (1.0/30.0) * add2 * add2 ) );
		

		} else if ( (double)SRE.pointer->fMRI[i1][i2] <= TH ) {

		
		SRE.pointer->theta[i1][i2] = (double)0.0;

		SRE.pointer->theta1[i1][i2] = (double)0.0; 
		
		SRE.pointer->theta2[i1][i2] = (double)0.0;

		SRE.pointer->theta3[i1][i2] = (double)0.0;

		SRE.pointer->theta4[i1][i2] = (double)0.0;

        SRE.pointer->theta5[i1][i2] = (double)0.0;
	

		}

		
		} // y dim
        
	}  // x dim
	// compute thetas (end)

	// compute the h functions (begin)
	double h_LGR1 = (double) SCALE * (double) exp ( - (double)( 0.5 * misplacement_X * misplacement_X * misplacement_X - misplacement_X * misplacement_X - 0.5 * misplacement_X + 1.0 ) );

	double h_LGR2 = (double) SCALE * (double) exp ( - (double)( -(1.0/6.0) * misplacement_X * misplacement_X * misplacement_X + misplacement_X * misplacement_X - (11.0/6.0) * misplacement_X + 1.0 ) ) ;
	// compute the h functions (end)

	std::cout << "Thetas Calculated" << endl;

	
	// compute XSRE (begin)
	double q, p, z1, z, y;
	double fi0, fi1, fi2, fi3, cub;

	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
  
		
		 if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) { 


			fi0 = ((double) SRE.pointer->theta2[i1][i2]);

			fi1 = ((double) 2.0 * SRE.pointer->theta3[i1][i2]);

			fi2 = ((double) 3.0 * SRE.pointer->theta4[i1][i2]);

			fi3 = ((double) 4.0 * SRE.pointer->theta5[i1][i2]);


			if ( (double)fi3 != 0.0 )
			{

				q = ((double) ( (2.0/27.0) * ( fi2 / fi3 ) * ( fi2 / fi3 ) ) - (  ( fi2 / ( 3.0 * fi3 ) ) * ( fi1 / fi3 ) ) + ( fi0 / fi3 ));

				p = ((double) (  - ( fi2 / ( 9.0 * fi3 ) ) * ( fi2 / ( 9.0 * fi3 ) )  ) + ( fi1 / fi3 ));

			} else {


				q = (double)0.0;

				p = (double)0.0;

			}


			 if (SREPLUS == 1)
			 { //if

			 if ( ( ( q / 2.0 ) * ( q / 2.0 ) + ( ( p * p * p ) / 27.0 ) ) >= 0 )

				   z1 = ((double) ( q / 2.0 ) + (double) sqrt( (double) ( q / 2.0 ) * ( q / 2.0 ) + ( ( p * p * p ) / 27.0 )  ) );
			 
			 else  z1 = ((double) q / 2.0 );

			 } else if ( SREPLUS == 0 ) { //if

			 if ( ( ( q / 2.0 ) * ( q / 2.0 ) + ( ( p * p * p ) / 27.0 ) ) >= 0 )

				   z1 = ((double) ( q / 2.0 ) - (double) sqrt( (double) ( q / 2.0 ) * ( q / 2.0 ) + ( ( p * p * p ) / 27.0 )  ) );

			 else  z1 = ((double) q / 2.0 );

			 } // if


			 cub  = ((double)1.0 / 3.0 );

			 z  = (double) pow( (double)z1, (double)cub );
    
	
			 if ( (double)p == 0.0 && (double)z == 0.0 )		y = ((double) -1.0/3.0); // de l'Hopital

			 else if ( (double)p != 0.0 && (double)z == 0.0 )	y = (double)0.0;

			 else if ( (double)p == 0.0 && (double)z != 0.0 )	y = (double)z;
			 
			 else if ( (double)p != 0.0 && (double)z != 0.0 )	y = ((double) ( z - ( p / ( 3.0 * z ) ) ));


			 if ( (double)y != 0.0 && (double)fi2 != 0.0 && (double)fi3 != 0.0 ) 
				 
				 SRE.pointer->Xsre[i1][i2] = ((double) ( y - ( fi2 / ( 3.0 * fi3 ) ) ));

			 else if ( (double)y == 0.0 && (double)fi2 == 0.0 && (double)fi3 == 0.0 ) 

			     SRE.pointer->Xsre[i1][i2] = ((double) -1.0/3.0); // de l'Hopital
	
			 else if ( (double)y != 0.0 && (double)fi2 == 0.0 && (double)fi3 == 0.0 ) 

				 SRE.pointer->Xsre[i1][i2] = ((double) y - (1.0/3.0));

			 else if ( (double)y == 0.0 && (double)fi2 != 0.0 && (double)fi3 == 0.0 ) 
				 
				 SRE.pointer->Xsre[i1][i2] = ((double) 0.0);

			 else if ( (double)y == 0.0 && (double)fi2 == 0.0 && (double)fi3 != 0.0 ) 

				 SRE.pointer->Xsre[i1][i2] = ((double)0.0);


			 else if ( (double)y != 0.0 && (double)fi2 != 0.0 && (double)fi3 == 0.0 ) 

				  SRE.pointer->Xsre[i1][i2] = ((double) 0.0);

			 else if ( (double)y != 0.0 && (double)fi2 == 0.0 && (double)fi3 != 0.0 ) 

				  SRE.pointer->Xsre[i1][i2] = ((double) y);

			 else if ( (double)y == 0.0 && (double)fi2 != 0.0 && (double)fi3 != 0.0 ) 

				  SRE.pointer->Xsre[i1][i2] = ((double) - ( fi2 / ( 3.0 * fi3 ) ) );




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

           if ( max == min ) SRE.pointer->Xsre[i1][i2] = (double)0.0;

           else SRE.pointer->Xsre[i1][i2] = ( (double) (misplacement_X) * 
			                                  (double) fabs ( (min - SRE.pointer->Xsre[i1][i2]) / (min - max) ) );
              
		   } else if ( (double)SRE.pointer->fMRI[i1][i2] <= TH ) {

			    SRE.pointer->Xsre[i1][i2] = (double)0.0;
		   }

		} // y dim
        
	}  // x dim
	

	std::cout << "XSRE Calculated" << endl;
	// scale the XSRE to the misplacement (end)
	

	// compute DELTA (begin)
	double den_DELTA, num_DELTA;
	double powr5, powr4, powr3, powr2;


	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) {  


		powr5 = (double) pow( (double)SRE.pointer->Xsre[i1][i2] , (double)5.0 );
       
        powr4 = (double) pow( (double)SRE.pointer->Xsre[i1][i2] , (double)4.0 );

        powr3 = (double) pow( (double)SRE.pointer->Xsre[i1][i2] , (double)3.0 );

		powr2 = (double) pow( (double)SRE.pointer->Xsre[i1][i2] , (double)2.0 );


		den_DELTA = ((double) ( powr5 * SRE.pointer->theta5[i1][i2] ) + 
			                  ( powr4 * SRE.pointer->theta4[i1][i2] ) + 
							  ( powr3 * SRE.pointer->theta3[i1][i2] ) + 
							  ( powr2 * SRE.pointer->theta2[i1][i2] ) + 
					          ( SRE.pointer->Xsre[i1][i2] * SRE.pointer->theta1[i1][i2] ));
	

		powr5 = (double) pow( ((double)SRE.pointer->Xsre[i1][i2] - misplacement_X), (double)5.0 );

		powr4 = (double) pow( ((double)SRE.pointer->Xsre[i1][i2] - misplacement_X), (double)4.0 );

		powr3 = (double) pow( ((double)SRE.pointer->Xsre[i1][i2] - misplacement_X), (double)3.0 );

		powr2 = (double) pow( ((double)SRE.pointer->Xsre[i1][i2] - misplacement_X), (double)2.0 );


		num_DELTA = ((double) ( powr5 * SRE.pointer->theta5[i1][i2] ) + 
			                  ( powr4 * SRE.pointer->theta4[i1][i2] ) + 
							  ( powr3 * SRE.pointer->theta3[i1][i2] ) + 
							  ( powr2 * SRE.pointer->theta2[i1][i2] ) + 
					          ( (SRE.pointer->Xsre[i1][i2] - misplacement_X) * SRE.pointer->theta1[i1][i2] ));

		if ( (double)den_DELTA == 0.0 )							 
			
			SRE.pointer->DELTA[i1][i2] = (double)0.0;
		
		else if ( (double)num_DELTA == 0.0 && (double)den_DELTA == 0.0 ) 
			
			SRE.pointer->DELTA[i1][i2] = (double)1.0; // de L'Hopital

		else if ( (double)num_DELTA != 0.0 && (double)den_DELTA == 0.0 ) 
			
			SRE.pointer->DELTA[i1][i2] = (double)0.0; 
		else 												 
			
			SRE.pointer->DELTA[i1][i2] = ((double) num_DELTA / den_DELTA);
	

		} else if ( (double)SRE.pointer->fMRI[i1][i2] <= TH ) {

			   SRE.pointer->DELTA[i1][i2] = (double)0.0;

		}

		} // y dim
        
	}  // x dim
	// compute DELTA (end)

	std::cout << "DELTA Calculated" << endl;	

		/// Calculation of the Intensity-Curvature Functional (begin) ///
	double num, den;
	for (int i1=0; i1 < n1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) {  


		powr5 = (double) pow( (double)misplacement_X , (double)5.0 );
       
        powr4 = (double) pow( (double)misplacement_X , (double)4.0 );

        powr3 = (double) pow( (double)misplacement_X , (double)3.0 );

		powr2 = (double) pow( (double)misplacement_X , (double)2.0 );


		den       = ((double) ( powr5 * SRE.pointer->theta5[i1][i2] ) + 
			                  ( powr4 * SRE.pointer->theta4[i1][i2] ) + 
							  ( powr3 * SRE.pointer->theta3[i1][i2] ) + 
							  ( powr2 * SRE.pointer->theta2[i1][i2] ) + 
					          ( misplacement_X * SRE.pointer->theta1[i1][i2] ));

		num = (double)SRE.pointer->fMRI[i1][i2] * misplacement_X * SRE.pointer->theta[i1][i2];
	

		if ( (double)den == 0.0 )							 
			
			SRE.pointer->ICF[i1][i2] = (double)0.0;
		
		else if ( (double)num == 0.0 && (double)den == 0.0 ) 
			
			SRE.pointer->ICF[i1][i2] = (double)1.0; // de L'Hopital

		else if ( (double)num != 0.0 && (double)den == 0.0 ) 
			
			SRE.pointer->ICF[i1][i2] = (double)0.0; 
		else 												 
			
			SRE.pointer->ICF[i1][i2] = ((double) num / den);
	

		} else if ( (double)SRE.pointer->fMRI[i1][i2] <= TH ) {

			SRE.pointer->ICF[i1][i2] = (double)0.0;

		}

		} // y dim
        
	}  // x dim
	/// Calculation of the Intensity-Curvature Functional (end) ///

	// calculate the re-computed locations (begin)
	double n, s, r, v, m;
	double alpha_1, beta_1, gamma_1;
	double A, B, C, D;
	double o1, o2, df, f_df_sre;
	double PSI, f_1, f_2;
	double E, F, G, H, SQRTTERM, fi4, f2;

	for (int i1=PAD+1; i1 < n1-PAD-1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

		if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) { 


		if ( i1 == (int)firstLine || i1 == ((int)firstLine + 1) || i1 == ((int)firstLine + 2) )
		{ // calculate the first three horizonatal lines of the signal (begin)

		f_2 = ((double) ( SRE.pointer->the_fMRI[i1+2+XNEI][i2] + SRE.pointer->the_fMRI[i1+1+XNEI][i2] + 
			              SRE.pointer->the_fMRI[i1-1+XNEI][i2] + SRE.pointer->the_fMRI[i1-2+XNEI][i2]  ));       
   
            
		f_1 = ((double)  SRE.pointer->the_fMRI[i1+1+XNEI][i2] + SRE.pointer->the_fMRI[i1-1+XNEI][i2]  );  
		
		// calculate the first three horizonatal lines of the signal (end)
		} else 

		if ( i1 == ((int)lastLine) ||  i1 == ((int)lastLine - 1) || i1 == ((int)lastLine - 2) )
		{ // calculate the last three horizonatal lines of the signal (begin)


		f_2 = ((double) ( SRE.pointer->the_fMRI[i1+2-XNEI][i2] + SRE.pointer->the_fMRI[i1+1-XNEI][i2] + 
			              SRE.pointer->the_fMRI[i1-1-XNEI][i2] + SRE.pointer->the_fMRI[i1-2-XNEI][i2]  ));       
   
            
		f_1 = ((double)  SRE.pointer->the_fMRI[i1+1-XNEI][i2] + SRE.pointer->the_fMRI[i1-1-XNEI][i2]  );  
	
		} else { // calculate the rest of the lines of the signal (begin)

		f_2 = ((double) ( SRE.pointer->the_fMRI[i1+2][i2] + SRE.pointer->the_fMRI[i1+1][i2] + 
			              SRE.pointer->the_fMRI[i1-1][i2] + SRE.pointer->the_fMRI[i1-2][i2]  ));       
   
            
		f_1 = ((double)  SRE.pointer->the_fMRI[i1+1][i2] + SRE.pointer->the_fMRI[i1-1][i2]  );  
	

		}// calculate the rest of the lines of the signal (end)



		fi4 =  ((double) ( (3.0/2.0) * f_1 * f_1 ) - ( 0.5 * f_1 * f_2 ) - ( 0.5 * f_1 * f_2 ) + ( (1.0/6.0) * f_2 * f_2 ) );

      
		fi3 =  ((double) ( -4.0 * f_1 * f_1 ) + ( 2.0 * f_1 * f_2 ) + ( (10.0/3.0) * f_1 * f_2 ) - ( (4.0/3.0) * f_2 * f_2 ) );


		fi2 =  ((double) ( -(7.0/2.0) * f_1 * f_1 ) - ( (3.0/2.0) * f_1 * f_2 ) - ( (15.0/2.0) * f_1 * f_2 ) + ( (23.0/6.0) * f_2 * f_2 ) );


		fi1 =  ((double) ( 3.0 * f_1 * SRE.pointer->the_fMRI[i1][i2] ) - ( SRE.pointer->the_fMRI[i1][i2] * f_2 ) + ( 4.0 * f_1 * f_1 ) - ( 2.0 * f_1 * f_2 ) + ( (20.0/3.0) * f_1 * f_2 ) - ( (14.0/3.0) * f_2 * f_2 ) );




		o1 = ((double) ( 0.5 * SRE.pointer->Xsre[i1][i2] * SRE.pointer->Xsre[i1][i2] * SRE.pointer->Xsre[i1][i2] ) - 
			
			                 ( SRE.pointer->Xsre[i1][i2] * SRE.pointer->Xsre[i1][i2] ) -
							 
							 ( 0.5 * SRE.pointer->Xsre[i1][i2] ) + 1.0 );

     	o2 = ((double) ( -(1.0/6.0) * SRE.pointer->Xsre[i1][i2] * SRE.pointer->Xsre[i1][i2] * SRE.pointer->Xsre[i1][i2] ) + 
			
			                    ( SRE.pointer->Xsre[i1][i2] * SRE.pointer->Xsre[i1][i2] ) - 
								
								( (11.0/6.0) * SRE.pointer->Xsre[i1][i2] ) + 1.0 );


     	f2 = ((double) SRE.pointer->the_fMRI[i1][i2] + ( f_1 * o1 ) + ( f_2 * o2 )); 


     	df = ((double) ( f_1 * ( ( 3.0 * SRE.pointer->Xsre[i1][i2] ) - 2.0 ) ) + 
			           ( f_2 * ( - SRE.pointer->Xsre[i1][i2] + 2.0 ) ));

		
		f_df_sre = ((double)f2 * df);

       	r = ((double)SRE.pointer->DELTA[i1][i2]);

		fi0 =  ((double) ( -2.0 * f_1 * SRE.pointer->the_fMRI[i1][i2] ) + ( 2.0 * SRE.pointer->the_fMRI[i1][i2] * f_2 ) - 
			             ( 2.0 * f_1 * f_1 ) + ( 2.0 * f_2 * f_2 ) - ( r * f_df_sre ));
 

		//if ( fi4 != 0.0 )
		if ( (double)fi3 != 0.0 && (double)fi4 != 0.0 )
		{ // if

			A = ((double) ( fi3 / ( 4.0 * fi4 ) ) * ( fi3 / ( 4.0 * fi4 ) ) * ( fi3 / ( 4.0 * fi4 ) ) * ( fi3 / ( 4.0 * fi4 ) ) );

			B = ((double) ( fi3 / ( 4.0 * fi4 ) ) * ( fi3 / ( 4.0 * fi4 ) ) * ( fi3 / ( 4.0 * fi4 ) ) );

			C = ((double) ( fi3 / ( 4.0 * fi4 ) ) * ( fi3 / ( 4.0 * fi4 ) ) );
	
			D = ((double) ( fi3 / ( 4.0 * fi4 ) ) );

		} else if ( (double)fi3 == 0.0 && (double)fi4 == 0.0 )
		{ // if

			A = (double)1.0; // de L'Hopital

			B = (double)1.0;

			C = (double)1.0;
	
			D = (double)1.0;

		} else if ( (double)fi3 != 0.0 && (double)fi4 == 0.0 )
		{ // if

			A = (double)0.0; 

			B = (double)0.0;

			C = (double)0.0;
	
			D = (double)0.0;

		} else if ( (double)fi3 == 0.0 && (double)fi4 != 0.0 )
		{ // if

			A = (double)0.0; 

			B = (double)0.0;

			C = (double)0.0;
	
			D = (double)0.0;

		}


		if ( (double)fi4 != 0.0 )
		{ // if

			alpha_1 = ((double) ( 6.0 * C - 3.0 * D * ( fi3 / fi4 ) + ( fi2 / fi4 ) ));

			beta_1 =  ((double) ( -4.0 * B + 3.0 * C * ( fi3 / fi4 ) - 2.0 * D * ( fi2 / fi4 ) + ( fi1 / fi4 ) ));

			gamma_1 = ((double) ( A - B * ( fi3 / fi4 ) + C * ( fi2 / fi4 ) - D * ( fi1 / fi4 ) + ( fi0 / fi4 ) )); 

		} else if ( fi4 == 0.0 ) {  // if

			alpha_1 = ((double) 0.0);

			beta_1 =  ((double) 0.0);

			gamma_1 = ((double) 0.0); 

		} // if


		E = ((double) alpha_1 * alpha_1);

		F = ((double) alpha_1 * alpha_1 * alpha_1);

		G = ((double) beta_1 * beta_1);


		p = ((double) (  ( (5.0/6.0) * (5.0/3.0) * E ) + ( (5.0/6.0) * (5.0/6.0) * E ) - ( (5.0/2.0) * (5.0/3.0) * E )  + ( 2.0 * E - gamma_1 )  )  );

		H = ((double) (  F - ( alpha_1 * gamma_1 ) - ( G / 4.0 )  ) / 2.0);

		q = ((double) - ( (5.0/6.0) * (5.0/6.0) * (5.0/6.0) * F ) + ( (5.0/6.0) * (5.0/6.0) * (5.0/2.0) * F ) - ( 2.0 * E - gamma_1 ) * (5.0/6.0) * alpha_1 + H);


		if (NOVELPLUS == 1) {  // if

			if ( ( ( q / 2.0 ) * ( q / 2.0 ) + ( ( p * p * p ) / 27.0 ) ) >= 0 )
	
				   z1 = ((double) ( q / 2.0 ) + (double) sqrt( (double) ( q / 2.0 ) * ( q / 2.0 ) + ( ( p * p * p ) / 27.0 ) ) );
	
			else   z1 = ((double) q / 2.0 );

		} else if (NOVELPLUS == 0) {// if


			if ( ( ( q / 2.0 ) * ( q / 2.0 ) + ( ( p * p * p ) / 27.0 ) ) >= 0 )
	
				   z1 =  ((double) ( q / 2.0 ) - (double) sqrt( (double) ( q / 2.0 ) * ( q / 2.0 ) + ( ( p * p * p ) / 27.0 ) ) ) ;

			else   z1 = ((double) q / 2.0 );

		} // if


		cub  = ((double) 1.0 / 3.0);

		z = (double) pow( (double)z1, (double)cub );
    

		if (  (double)z != 0.0 ) { // if


			v = ((double) (  ( p / ( 3.0 * z ) ) - z ));

			y = ((double) ( v - ( (5/6) * alpha_1 ) ));

			m = ((double) ( 2.0 * (alpha_1 + 2 * y ) ));

		} 
		else if (  (double)z == 0.0 && (double)p == 0.0 ) { // if
		
		
			v = (double) 1.0; // de l'Hopital

			y = ((double) ( v - ( (5/6) * alpha_1 ) ));

			m = ((double) ( 2.0 * (alpha_1 + 2 * y ) ));

		} // if
		else if (  (double)z != 0.0 && (double)p == 0.0 ) { // if

			
			v = ((double) (  ( p / ( 3.0 * z ) ) - z ));

			y = ((double) ( v - ( (5/6) * alpha_1 ) ));

			m = ((double) ( 2.0 * (alpha_1 + 2 * y ) ));

		} // if
		else if (  (double)z == 0.0 && (double)p != 0.0 ) { // if

			
			v = (double) 0.0;

			y = ((double) ( v - ( (5/6) * alpha_1 ) ));

			m = ((double) ( 2.0 * (alpha_1 + 2 * y ) ));

		} // if


		if ( (alpha_1 + 2.0 * y ) > 0.0 ) { // if
	
			n = ((double) beta_1 / ( 2.0 * (double) sqrt( (double)alpha_1 + 2.0 * y ) ));
	
			s = ((double) sqrt( (double)alpha_1 + 2.0 * y ));

		} else { // if

			n = (double)0.0;
			s = (double)0.0;

		} // if

		
		SQRTTERM = ((double) ( ( s * s ) - 4.0 * ( alpha_1 + y + n ) ));

		if ( (SQRTTERM >= 0.0) && (fi4 != 0.0) )
		{ // if

			if ( SQRTNOVELPLUS == 1 ) 
			
				PSI  =  ((double) ( ( s + (double)sqrt( (double)SQRTTERM ) ) / 2.0 ) - ( fi3 / ( 4.0 * fi4 ) ) ); 
		
			else if ( SQRTNOVELPLUS == 0 ) 

				PSI  = -((double) ( ( s + (double)sqrt( (double)SQRTTERM ) ) / 2.0 ) - ( fi3 / ( 4.0 * fi4 ) ) );  
		

				SRE.pointer->X_recomputed[i1][i2] = (double)PSI;


		} else if ( (SQRTTERM > 0.0) && (fi4 == 0.0) ) { // else if
			
			    if ( SQRTNOVELPLUS == 1 ) 
			
				PSI  = (double) 0.0;  
		
			    else if ( SQRTNOVELPLUS == 0 ) 

				PSI  = (double) 0.0;  
		

				SRE.pointer->X_recomputed[i1][i2] = (double)PSI;	
		
		} else if ( (SQRTTERM > 0.0) && (fi4 == 0.0) && (fi3 == 0.0) ) { // else if
			
			    if ( SQRTNOVELPLUS == 1 ) 
			
				PSI  =  ((double) ( ( s + (double)sqrt( (double)SQRTTERM ) ) / 2.0 ) - ( 1.0 / 4.0 ) ); 
		
			else if ( SQRTNOVELPLUS == 0 ) 

				PSI  = -((double) ( ( s + (double)sqrt( (double)SQRTTERM ) ) / 2.0 ) - ( 1.0 / 4.0 ) ); 
		

				SRE.pointer->X_recomputed[i1][i2] = (double)PSI;	
		
		} else if ( (SQRTTERM < 0.0) && (fi4 != 0.0) ) { // else if

				if ( SQRTNOVELPLUS == 1 ) 
			
				PSI  =  ((double) (s / 2.0) - ( fi3 / ( 4.0 * fi4 ) ) ); 
		
			    else if ( SQRTNOVELPLUS == 0 ) 

				PSI  = -((double) (s / 2.0) - ( fi3 / ( 4.0 * fi4 ) ) );  
		

				SRE.pointer->X_recomputed[i1][i2] = (double)PSI;

		} else if ( (SQRTTERM < 0.0) && (fi4 == 0.0) ) { // else if

				if ( SQRTNOVELPLUS == 1 ) 
			
				PSI  =  ((double) (s / 2.0) ); 
		
			    else if ( SQRTNOVELPLUS == 0 ) 

				PSI  = -((double) (s / 2.0) );  
		

				SRE.pointer->X_recomputed[i1][i2] = (double)PSI;

		} // else if

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

		if ( (double)SRE.pointer->the_fMRI[i1][i2] > TH ) { 
	

		if( SRE.pointer->X_recomputed[i1][i2] > (double)max ) max = (double)SRE.pointer->X_recomputed[i1][i2];
              
		if( SRE.pointer->X_recomputed[i1][i2] < (double)min ) min = (double)SRE.pointer->X_recomputed[i1][i2];
	
		
		}

		} // y dim
        
	}  // x dim


	for (int i1=0; i1 < n1; i1++) {// x dim

		for (int i2=0; i2 < n2; i2++) { // y dim

		   if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) {
	
           if ( max == min ) SRE.pointer->X_recomputed_scaled[i1][i2] = (double)0.0;

           else  SRE.pointer->X_recomputed_scaled[i1][i2] = ((double) ( XPixelSize ) * 
			                                                 (double) fabs ( ( min - SRE.pointer->X_recomputed[i1][i2] ) / ( min - max )) );
                    
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
			
			SRE.pointer->h_LGR1_recomp[i1][i2] = (double) SCALE * (double) exp ( (double) ( - ( 0.5 * SRE.pointer->X_recomputed_scaled[i1][i2] * SRE.pointer->X_recomputed_scaled[i1][i2] * SRE.pointer->X_recomputed_scaled[i1][i2] - SRE.pointer->X_recomputed_scaled[i1][i2] * SRE.pointer->X_recomputed_scaled[i1][i2] - 0.5 * SRE.pointer->X_recomputed_scaled[i1][i2] + 1.0 ) ));


			SRE.pointer->h_LGR2_recomp[i1][i2] = (double) SCALE * (double) exp ( -(double) ( -(1.0/6.0) * SRE.pointer->X_recomputed_scaled[i1][i2] * SRE.pointer->X_recomputed_scaled[i1][i2] * SRE.pointer->X_recomputed_scaled[i1][i2] + SRE.pointer->X_recomputed_scaled[i1][i2] * SRE.pointer->X_recomputed_scaled[i1][i2] - (11.0/6.0) * SRE.pointer->X_recomputed_scaled[i1][i2] + 1.0 ) );


			} else if ( (double)SRE.pointer->fMRI[i1][i2] <= TH ) {


			SRE.pointer->h_LGR1_recomp[i1][i2] = (double) 0.0;


			SRE.pointer->h_LGR2_recomp[i1][i2] = (double) 0.0;


			}

		} // y dim
        
	}  // x dim
	// compute the h functions (end)

	std::cout << "H functions Calculated" << endl;

	// shift the image of the given misplacement (begin)
	double add_up, sum;

	for (int i1=PAD+1; i1 < n1-PAD-1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim

			if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) {
           
			if ( i1 == (int)firstLine || i1 == ((int)firstLine + 1) || i1 == ((int)firstLine + 2) )
			{ // calculate the first three horizonatal lines of the signal (begin)

		    add_up = ((double) ( SRE.pointer->the_fMRI[i1+2+XNEI][i2] + SRE.pointer->the_fMRI[i1+1+XNEI][i2] + 
				                 SRE.pointer->the_fMRI[i1-1+XNEI][i2] + SRE.pointer->the_fMRI[i1-2+XNEI][i2]  ));       
   
            
			sum = ((double) SRE.pointer->the_fMRI[i1+1+XNEI][i2] + SRE.pointer->the_fMRI[i1-1+XNEI][i2]  );   
			
			// calculate the first three horizonatal lines of the signal (end)
			} else 

			if ( i1 == ((int)lastLine) ||  i1 == ((int)lastLine - 1) || i1 == ((int)lastLine - 2) )
			{ // calculate the last three horizonatal lines of the signal (begin)

			add_up = ((double) ( SRE.pointer->the_fMRI[i1+2-XNEI][i2] + SRE.pointer->the_fMRI[i1+1-XNEI][i2] + 
				                 SRE.pointer->the_fMRI[i1-1-XNEI][i2] + SRE.pointer->the_fMRI[i1-2-XNEI][i2]  ));       
   
            
			sum = ((double) SRE.pointer->the_fMRI[i1+1-XNEI][i2] + SRE.pointer->the_fMRI[i1-1-XNEI][i2]  );   
			
			} else { // calculate the rest of the lines of the signal (begin)

			
			add_up = ((double) ( SRE.pointer->the_fMRI[i1+2][i2] + SRE.pointer->the_fMRI[i1+1][i2] + 
				                 SRE.pointer->the_fMRI[i1-1][i2] + SRE.pointer->the_fMRI[i1-2][i2]  ));       
   
            
			sum = ((double) SRE.pointer->the_fMRI[i1+1][i2] + SRE.pointer->the_fMRI[i1-1][i2]  );   
			

			}// calculate the rest of the lines of the signal (end)

			SRE.pointer->shifted_fMRI[i1][i2] = ((double) SRE.pointer->fMRI[i1][i2] + ( sum * h_LGR1 ) + ( add_up * h_LGR2 ) );

			} else if ( (double)SRE.pointer->fMRI[i1][i2] <= TH ) {

				   SRE.pointer->shifted_fMRI[i1][i2] = (double)0.0;

			}

 		} // y dim
        
	}  // x dim
	// shift the image of the given misplacement (end)

	std::cout << "The Image is Shifted" << endl;

	// motion correct the shifted image of the given misplacement,
	// and motion correct the shifted image of the recomputed location,
	// and compute SE (square error) -> (begin)

	double sum_SE = 0.0;
	double sum_SE_SRE = 0.0;

	for (int i1=PAD+1; i1 < n1-PAD-1; i1++) {// x dim
       	
		for (int i2=0; i2 < n2; i2++) { // y dim
    
			if ( (double)SRE.pointer->fMRI[i1][i2] > TH ) {

			if ( i1 == (int)firstLine || i1 == ((int)firstLine + 1) || i1 == ((int)firstLine + 2) )
			{ // calculate the first three horizonatal lines of the signal (begin)

		    add_up = ((double) ( SRE.pointer->the_fMRI[i1+2+XNEI][i2] + SRE.pointer->the_fMRI[i1+1+XNEI][i2] + 
				                 SRE.pointer->the_fMRI[i1-1+XNEI][i2] + SRE.pointer->the_fMRI[i1-2+XNEI][i2]  ));       
   
            
			sum = ((double) SRE.pointer->the_fMRI[i1+1+XNEI][i2] + SRE.pointer->the_fMRI[i1-1+XNEI][i2]  );   
			
			// calculate the first three horizonatal lines of the signal (end)
			} else 

			if ( i1 == ((int)lastLine) ||  i1 == ((int)lastLine - 1) || i1 == ((int)lastLine - 2) )
			{ // calculate the last three horizonatal lines of the signal (begin)

			add_up = ((double) ( SRE.pointer->the_fMRI[i1+2-XNEI][i2] + SRE.pointer->the_fMRI[i1+1-XNEI][i2] + 
				                 SRE.pointer->the_fMRI[i1-1-XNEI][i2] + SRE.pointer->the_fMRI[i1-2-XNEI][i2]  ));       
   
            
			sum = ((double) SRE.pointer->the_fMRI[i1+1-XNEI][i2] + SRE.pointer->the_fMRI[i1-1-XNEI][i2]  );   
			
			} else { // calculate the rest of the lines of the signal (begin)

			
			add_up = ((double) ( SRE.pointer->the_fMRI[i1+2][i2] + SRE.pointer->the_fMRI[i1+1][i2] + 
				                 SRE.pointer->the_fMRI[i1-1][i2] + SRE.pointer->the_fMRI[i1-2][i2]  ));       
   
            
			sum = ((double) SRE.pointer->the_fMRI[i1+1][i2] + SRE.pointer->the_fMRI[i1-1][i2]  );   
			

			}// calculate the rest of the lines of the signal (end)

			SRE.pointer->the_m_corr_recomp_fMRI[i1][i2] = ((double) SRE.pointer->fMRI[i1][i2] + ( sum * SRE.pointer->h_LGR1_recomp[i1][i2] ) + ( add_up * SRE.pointer->h_LGR2_recomp[i1][i2] ) );

			} else if ( (double)SRE.pointer->fMRI[i1][i2] <= TH ) {

					SRE.pointer->the_m_corr_recomp_fMRI[i1][i2] = (double)0.0;
		
			}
		
		} // y dim
        
	}  // x dim

	std::cout << "Image SRE Shifted" << endl;
	std::cout << endl;

	// calculate Square Errors (begin)
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
	std::cout << "SE without SRE is: " << ((double)sum_SE/( counter ))  << endl;
	std::cout << endl;
	std::cout << "SE with SRE is: " << ((double)sum_SE_SRE/( counter ))  << endl;
	std::cout << endl;
	//----->>> print and save Square Errors (end)
	// calculate Square Errors (end)


	// motion correct the shifted image of the given misplacement,
	// and motion correct the shifted image of the recomputed location,
	// and compute SE (square error) -> (end)

	double DIFF = ((double) ((double)sum_SE_SRE/( counter )) - ((double)sum_SE/( counter )) );

	fprintf(savedata, "%e\t %e\t %e\t %e\n", misplacement_X, ((double)sum_SE/( counter )), ((double)sum_SE_SRE/( counter )), DIFF);

    STEP_A = ((double)STEP_A + 0.005);
             
    STEP_MX = ((double)STEP_MX + 0.005);
            
	}  //---> major for loop (end) <---//

	fprintf(savedata,"\n");
	fprintf(savedata,"%s\n", "End of Computation...");
	std::cout << "End of Computation..." << endl;
	std::cout << endl;

	
	fclose(savedata);
	delete SRE.pointer;
	SRE.~LGRSRE2013();
	} // processing (end)

	} // run the program (end)

	system("java processCompleted");

	return 0;
} // end of main 