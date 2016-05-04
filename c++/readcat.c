#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <array>
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



////  Count the number of readable rows and
////  columns and set to global variables.
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

    readdata(base_dir + "images/sextractor/catalog_B.cat", catalog_B);
    readdata(base_dir + "images/sextractor/catalog_V.cat", catalog_V);
    readdata(base_dir + "images/sextractor/catalog_Rsup.cat", catalog_Rsup);
    readdata(base_dir + "images/sextractor/catalog_Isup.cat", catalog_Isup);
    readdata(base_dir + "images/sextractor/catalog_r.cat", catalog_r);
    readdata(base_dir + "images/sextractor/catalog_i.cat", catalog_i);
    readdata(base_dir + "images/sextractor/catalog_z.cat", catalog_z);
    readdata(base_dir + "images/sextractor/catalog_J.cat", catalog_J);
    readdata(base_dir + "images/sextractor/catalog_K.cat", catalog_K);
    cout << "done!" << endl << endl;



    float zeropoints[9] = {25., 25., 27.44, 27.29, 26.50, 26.20, 25.10, 28.17 + 0.922, 27.55 + 1.900};
    float flux_errors[9] = {0.0414, 0.0625, 1.0644, 1.5525, 0.5309, 0.6829, 0.9143, 108.63, 178.05};

    float MW_extinction_mag[9] = {0.078, 0.056, 0.049, 0.036, 0.049, 0.036, 0.027, 0.015, 0.006};
    double MW_extinction_scale[9];
    for(int i=0; i<9; i++) {
        MW_extinction_scale[i] = pow(10, MW_extinction_mag[i] / 2.5);
    }



    ////  writing .translate file
    ofstream translateFile;
    translateFile.open ("test.translate");

    translateFile << "fluxaper_B   F8\nerraper_B   E8\n";
    translateFile << "fluxaper_V   F9\nerraper_V   E9\n";
    translateFile << "fluxaper_Rplus  F14\nerraper_Rplus  E14\n";
    translateFile << "fluxauto_Rplus  TOT14\n";
    translateFile << "fluxaper_Iplus  F15\nerraper_Iplus  E15\n";
    translateFile << "fluxaper_r   F5\nerraper_r   E5\n";
    translateFile << "fluxaper_i   F6\nerraper_i   E6\n";
    translateFile << "fluxaper_z   F7\nerraper_z   E7\n";
    translateFile << "fluxaper_J   F17\nerraper_J   E17\n";
    translateFile << "fluxaper_K   F19\nerraper_K   E19\n";
    translateFile << "fluxcolor_ch1  F20\nerrcolor_ch1  E20\n";
    translateFile << "fluxcolor_ch2  F21\nerrcolor_ch2  E21\n";
    translateFile.close();



    ////  setting up output catalogs
    ofstream fluxCatalogFile;
    ofstream zspecFluxCatalogFile;
    ofstream magCatalogFile;

    fluxCatalogFile.open ("catalog.cat");
    magCatalogFile.open ("catalog.mag");
    zspecFluxCatalogFile.open ("catalog.zspec.cat");

    fluxCatalogFile << "# id z_spec x y ra dec ";
    fluxCatalogFile << "fluxaper_B erraper_B ";
    fluxCatalogFile << "fluxaper_V erraper_V ";
    fluxCatalogFile << "fluxaper_Rsup erraper_Rsup ";
    fluxCatalogFile << "fluxaper_Isup erraper_Isup ";
    fluxCatalogFile << "fluxaper_r erraper_r ";
    fluxCatalogFile << "fluxaper_i erraper_i ";
    fluxCatalogFile << "fluxaper_z erraper_z ";
    fluxCatalogFile << "fluxaper_J erraper_J ";
    fluxCatalogFile << "fluxaper_K erraper_K ";
    fluxCatalogFile << "\n";

    magCatalogFile << "# id z_spec x y ra dec ";
    magCatalogFile << "magaper_B erraper_B ";
    magCatalogFile << "magaper_V erraper_V ";
    magCatalogFile << "magaper_Rsup erraper_Rsup ";
    magCatalogFile << "magaper_Isup erraper_Isup ";
    magCatalogFile << "magaper_r erraper_r ";
    magCatalogFile << "magaper_i erraper_i ";
    magCatalogFile << "magaper_z erraper_z ";
    magCatalogFile << "magaper_J erraper_J ";
    magCatalogFile << "magaper_K erraper_K ";
    magCatalogFile << "\n";

    zspecFluxCatalogFile << "# id z_spec x y ra dec ";
    zspecFluxCatalogFile << "fluxaper_B erraper_B ";
    zspecFluxCatalogFile << "fluxaper_V erraper_V ";
    zspecFluxCatalogFile << "fluxaper_Rsup erraper_Rsup ";
    zspecFluxCatalogFile << "fluxaper_Isup erraper_Isup ";
    zspecFluxCatalogFile << "fluxaper_r erraper_r ";
    zspecFluxCatalogFile << "fluxaper_i erraper_i ";
    zspecFluxCatalogFile << "fluxaper_z erraper_z ";
    zspecFluxCatalogFile << "fluxaper_J erraper_J ";
    zspecFluxCatalogFile << "fluxaper_K erraper_K ";
    zspecFluxCatalogFile << "\n";


    for(int i=0; i<nobj; i++) {

    	for(int filter_i=0; filter_i<n_filters; filter_i++) {
    		pass;
    	}

    }









    fluxCatalogFile.close();
    zspecFluxCatalogFile.close();
    magCatalogFile.close();





}










