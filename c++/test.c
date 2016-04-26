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

	double MW_extinction_mag[9] = {0.078, 0.056, 0.049, 0.036, 0.049, 0.036, 0.027, 0.015, 0.006};


	double MW_extinction_scale[9];
	for(int i=0; i<9; i++) {
		MW_extinction_scale[i] = pow(10, MW_extinction_mag[i] / 2.5);
	}

	string a = "hello";
	string b = "again";
	fprintf(stderr, "%s world %s!\n", a.c_str(), b.c_str());














}










