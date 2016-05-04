#include <stdio.h>
#include <string.h>
#include <iostream>
#include <iomanip>
#include <array>
#include <fstream>
#include <sstream>
using namespace std;

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include <unistd.h>
#define BUFFSIZE 8192


////  Definitions
string field_name = "sg0023+0423";
string version = "v0.0.1";

string base_dir = "/Volumes/PHOENIX/atomczak/DATA/ORELSE/catalogs/v001/Cl_0023+0423/";




//// Count the number of readable rows and
//// columns and set to global variables.
long nobj = 0;
long ncol = 0;
void count_nobj_ncol(string CATALOG_FILE) {
    FILE *fp;
    char buff[BUFFSIZE];

    if (!(fp = fopen(CATALOG_FILE.c_str(), "r"))) {
        fprintf(stderr,"Hey, I can't read from %s!\n", CATALOG_FILE.c_str());
        exit(1);
    }

    while (fgets(buff, BUFFSIZE, fp)) {
    	if (buff[0] == '#') ncol++;
        if (buff[0] != '#') ++nobj;
    }
    fclose(fp);

    printf("Catalog file contains %ld objects and %ld columns.\n", nobj, ncol);
}

void readdata(string CATALOG_FILE, double *catalog[nobj]) {
    FILE *fp;
    long i,j,k;
    char buff[BUFFSIZE], *arg;
    char  delim[]=" \n\t,\r";

    if (!(fp = fopen(CATALOG_FILE.c_str(), "r"))) {
        fprintf(stderr,"Hey, I can't read from %s!\n", CATALOG_FILE.c_str());
        exit(1);
    }
    fclose(fp);



    // ** Emulate sm -- assume lines starting with "#" are comments to be ignored!
    if (nobj<1) {
        fprintf(stderr,
            "Hm... a null catalog file is not very interesting. Quitting.\n");
        exit(1);
    }



    // //// Reopen file, now that we know whether or not we've read in readshifts
    // fp = fopen(CATALOG_FILE,"r");
    // fgets(buff,BUFFSIZE,fp);
    // get_column_defs(buff);
    // fclose(fp);

    // Now to actually read in the data
    if (!(fp = fopen(CATALOG_FILE.c_str(),"r"))) {
        fprintf(stderr, "Hey, I can't read from %s\n",CATALOG_FILE.c_str());
        exit(1);
    }  

    j = 0;
    while (fgets(buff, BUFFSIZE, fp))
        if (buff[0] != '#') {   // ignore comment "#" lines
        for (i=0; i<ncol; ++i) {
            if (i)    
                arg=strtok(NULL,delim);
            else
                arg=strtok(buff,delim);   
            if (!arg) {
                fprintf(stderr,
                  "Something is wrong with the data in line %ld: \n",
                  j+1);
                fprintf(stderr, "%s\n", buff);
                fprintf(stderr, "Maybe this line doesn't have the number of columns expected from the header?\n");
                exit(1);
            }
	        catalog[j][i] = atof(arg);
        }
        ++j;
    }
}





int main()
{
    ////  counting rows and columns and initializing 2d array
    cout <<  endl << "reading in catalogs..." << endl;
    count_nobj_ncol(base_dir + "images/sextractor/catalog_B.cat");

	cout << endl << field_name << '_' << version << endl << endl;


	string asdf = field_name + "_" + version;
	cout << asdf << endl << endl;



	// an array containing nobj pointers
    double *catalog_B[nobj];
    double *catalog_V[nobj];
    double *catalog_Rsup[nobj];
    double *catalog_Isup[nobj];
    double *catalog_r[nobj];
    double *catalog_i[nobj];
    double *catalog_z[nobj];
    double *catalog_J[nobj];
    double *catalog_K[nobj];
    for(int i = 0; i < nobj; i++){
        catalog_B[i] = new double[ncol];
        catalog_V[i] = new double[ncol];
        catalog_Rsup[i] = new double[ncol];
        catalog_Isup[i] = new double[ncol];
        catalog_r[i] = new double[ncol];
        catalog_i[i] = new double[ncol];
        catalog_z[i] = new double[ncol];
        catalog_J[i] = new double[ncol];
        catalog_K[i] = new double[ncol];
    }

    double *catalogs[9] = {*catalog_B, 
    	                   *catalog_V, 
    	                   *catalog_Rsup, 
    	                   *catalog_Isup, 
    	                   *catalog_r, 
    	                   *catalog_i, 
    	                   *catalog_z, 
    	                   *catalog_J, 
    	                   *catalog_K};

    readdata(base_dir + "images/sextractor/catalog_B.cat", catalog_B);

    cout << catalogs << endl;


}










