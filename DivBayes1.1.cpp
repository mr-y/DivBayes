/********************************************************************
Copyright (C) 2010 Martin Ryberg

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

contact: kryberg@utk.edu
*********************************************************************/

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>

using namespace std;

/*** Function declarations ***/
//Probabilities
double LikeValue (int *nspecies, double *ages, double *netdivvalues, double *extinctpropvalues, int ntaxa );
double PriorProb (char distribution[10], double *x, int a, double b, int ntaxa);
double PriorProb (char distribution[10], double *x, int *a, double *b, int ntaxa);
double AgePriorProb (char distribution[10], double *x, double *mean, double *sd, int ntaxa);

//Distributions
double gamma(double x, int k, double b);
double uniform(double x, double a, double b);
double normal(double a);
double normalprob (double x, double mean, double sd);

//Math
int factorial(int x);

//Simulations
int MCMC(double *netdivvalues, int *equalNDs, double *extinctpropvalues, int *equalEPs, int *nspecies, double *ages, double *agessd, int ngen, int sampfreq, char model, double multiplier, double *meanGamma, int *shapeGamma, int ntaxa);

//Misc
int WhichParameter (int value, int *equalvalues, int ntaxa);
void SetSameRates (double *values, int *equalvalues, int n);
void rearrange (double *values, int *rearrangearray, int n);
int nparameters (int *model, int n);
void MCMC_summarize (char filename[], int burnin);
void help();
/***/

/*** Classes ***/
/*** Class to handel data as a binary tree ***/
class data_tree {
    public:
        struct node { //store the information
            double value; //store a floating point value
            node *left; //points to lower values
            node *right; //point to higher values
        };
        data_tree () { root = 0; total = 0; number_of_values = 0; } //initilize the values
        ~data_tree () { destroy_tree(root); } //deletes the tree starting from root
        void insert_data (double x); //function ti insert data
        double mean() { return total/number_of_values; } //calculate mean
        double median() { //get the median value
            if ( n_nodes() % 2 == 0 ) return (get_value((n_nodes()/2)) + get_value((n_nodes()/2)+1))/2; //if even number of nodes get the mean of the equaly middle values
            else return get_value((n_nodes()/2)+1); //get middle vale
        }
        double sum() { //sum up values in tree
            if (root == 0) return 0; //if the root is empty there is no sum
            return sum(root); //call sum function
        }
        int n_nodes() { //get number of nodes
            if (root == 0) return 0; //if the root is empty there is no values to sum up
            return n_descendants(root); //calculate the sum from the root up
        }
        double sd() { //calculate standard deviation
            if ( root == 0) return 0; //if root is empty there is no values
            return sqrt(sum_square(root)/n_nodes()); //take the square root of the mean of sum of the square differences between values and the mean
        }
        double max() { //get the highest value in the tree
            node *present = root; //start at root
            if (root == 0) return 0; //if root emty there are no values
            while (present->right != 0) present=present->right; //while not at a tip go to the right node
            return present->value; //the highest value is in the right most node
        }
        double min() { //get the lowest value in the tree
            if ( root == 0 ) return 0; //if the root is empty there are no values
            node *present = root; //start at root
            while (present->left != 0) present=present->left; //while not at a tip go to the left
            return present->value; //the left most tip has the lowest value
        }
        double get_value( int i ); //function to get a value counted from the lowest (value 1 beeing the lowest)

    private:
        data_tree (const data_tree&, data_tree vector_copy); //objects of the class can not be copied
        void destroy_tree ( node *leaf ) { //delete the tree from the given node
            if(leaf!=0) { //while not at a tip delete next nodes
                destroy_tree ( leaf->left ); //delete node to the left
                destroy_tree ( leaf->right ); //delete node to the right
                delete leaf; //delete present node
            }
        }
        void insert_data ( double x, node *leaf ); //inser data somwhere after given node
        double sum(node *leaf) { //sum the values in the tree from the given node
            if ( leaf == 0 ) return 0; //if leaf is empty return 0
            return leaf->value + sum(leaf->left) + sum(leaf->right); //return the sum of the value of the leaf and all values to the left and right
        }
        double sum_square ( node *leaf ) { //calculate the sum of the square of the difference between the mean and the value for all nodes higher in the tree
            if ( leaf == 0 ) return 0; //if the leaf is empty return 0
            return (( leaf->value-mean() ) * ( leaf->value-mean() )) + sum_square(leaf->left) + sum_square(leaf->right); //return the sum of the square difference of this leaf and all nodes to the left and right
        }
        int n_descendants ( node *leaf ) { //returns the number of daughter nodes for a node +1
            if (leaf == 0) return 0; //if a empty leaf return 0
            return 1 + n_descendants(leaf->left) + n_descendants(leaf->right); //return 1 plus the number of descendants to the left and right
        }

        node *root; //store the root location of the tree
        double total; //store the total off all values in the tree, needed to make mean() fast so sd() is resonable fast
        int number_of_values; //store the total number of nodes in the tree for the same reason total is stored
};

/** Functions for sorted vector that ***
*** should not be inlined            **/
void data_tree::insert_data ( double x ) { //function to insert data
    if ( root == 0 ) { //if the node has not been initiated do so
        root = new node; //allocate memory for the root
        root->value = x; //store the value
        root->left = 0; //set the left leaf to a null pointer
        root->right = 0; //set the right leaf to a null pointer
    }
    else { insert_data( x, root ); } //otherwise insert data, starting to look for where at the root

}

void data_tree::insert_data ( double x, node *leaf ) { //insert data somwhere after the given node
    if ( x < leaf->value ) { //if the value is less than the value in the present leaf it belonges to the left
        if ( leaf->left != 0 ) insert_data ( x, leaf->left); //if not at a tip proceed (recursive function)
        else { //if at tip create node
            leaf->left = new node; //allocate memory for the left node
            leaf->left->value = x; //insert value
            leaf->left->left = 0; //set left leaf to null pointer
            leaf->left->right = 0; //set right leaf to null pointer
        }
    }
    else if ( x >= leaf->value ) { //if value higher (or equal) to present node insert it to the right
        if ( leaf->right !=0 ) insert_data ( x, leaf->right ); //if not at a tip proceed (recursive finction)
        else { //if at tip insert value
            leaf->right = new node; //allocate memory for right node
            leaf->right->value = x; //insert value
            leaf->right->left = 0; //set left leaf to null pointer
            leaf->right->right = 0; //set right leaf to null pointer
        }
    }
    total += x; //add value to the total
    ++number_of_values; //add one to the number of nodes
}

double data_tree::get_value( int i ) { //get the i:th lowest value
    if ( root == 0 ) { cerr << "Trying to get value from uninitiated data_tree; returning 0." << endl; return 0; } //if empty tree return 0, this is not strictly correct as it should return nothing
    if ( i < 0 || i > n_nodes() ) { cerr << "get_value function in data_tree (class) called with value out of range; returning 0." << endl; return 0; } //if the node is out of range return 0, should return nothing
    node *present = root; //start at root
    while ( i > 0 ) { //loop, should stop when reaching right value but if i is < 0 something is wrong
        if (i <= n_descendants(present->left) ) present = present->left; //if there is more nodes than i nodes to the left the value we are looking for is to the left
        else if ( i == n_descendants(present->left) + 1 ) return present->value; //when we have i-1 nodes to the left we are at the value we are looking for
        else if ( i > n_descendants(present->left) ) { //if i is higher than the number of nodes to the left we are looking for a higher value
            i -= n_descendants(present->left) + 1; //after we moved to the right we will have excluded the number of nodes to the left plus one from our search (as lower values)
            present = present->right; //the value we are looking for is higher than the present node and nodes to the left so ove to the right
        }
    }
    return present->value; //should not get here so shuld probably give error message and return 0 instead, but the value we are at should be closer to the real value we are looking for
}

/***********************************
************************************
*** Main function of the program ***
************************************
***********************************/
main (int argc, char *argv []) {

/*** Declarations of variables ***/
   int ngen=100;            //Number of generations
   int sampfreq=10;         //Sample frequency
   int i=0;                 //Counter
   double LHvalue=1;        //log likelihood value   
   int acceptrate=0;        //Variable to store the number of accepted state changes in the MCMC
   int BayesOrLH=0;         //Decides if to run Bayesia analysis of just calculate likelihood   
   char model= 'y';         //Variable to chose model d=Birth-Death model, y=Pure-Birth (Yule)
   double windowwidth=0.2;  //This number adjust the window width of the proposal dencity
   char temp[100];          //Variable to store words input from file
   int ntaxa=0;             //Variable to store number of clades
   char *filename;          //Pointer to stored file name
   int bracket_balance=0; //Flag to signal that we are in a comment that should not be read by program
   srand(time(NULL));       //Set the random seed by clock
/***/

/*** This will print the line that was used to call the program ***/
   cout<< "The program was called using the following command: " << endl;
   while (argv[i]) cout << argv[i++] << " ";
   cout << endl<< endl;
/***************************************************************/

/*** Parsing the arguments given to the program ***/
//   if (argc>2) {cout << "Too many arguments, just give name of the infile or --help (-h).\n"; return 1;}
   if (argc>2) {
      int burnin = 0;
      for (i=1; i < argc; ++i) {
         if (!strcmp(argv[i],"-s") || !strcmp(argv[i],"--summerize")) { filename = argv[++i]; }
         else if (!strcmp(argv[i],"-b") || !strcmp(argv[i],"--burnin")) { burnin = atoi(argv[++i]); }
         else {cout << "Do not recognize arguments, just give name of the infile, -s and name of MCMC file, or --help (-h).\n"; return 1; }
      }
      MCMC_summarize(filename, burnin); //start summarizing
      return 0;
   }

   if (argc==1) {cout << "Too few arguments, give name of the infile or --help (-h).\n"; return 1;}

   if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) { help(); return 1;}

   if (argc==2) { filename = argv[1];}

   ifstream infile (filename);

   if (!infile.is_open()) { cout << "Could not open \"" << argv[1] << "\". Quitting.\n"; return 1;}

// Find the number of taxa
   

   while (infile) {
      infile >> temp;
      /*** Remove comments in square brackets ***/
      i=0; //Set counter to 0
      while ( temp[i] != '\0') { //read each character till end of string
         if (bracket_balance < 0) { cout << "Found an unmatched ]. Square brackets need to be balanced. Quitting." << endl; return 1; }
         else if ( temp[i] == '[' ) { //Start of comment to be removed
            bracket_balance++;   //Note one [ for balance
            int j=i;
            while (temp[j] !='\0') temp[j++]=temp[j+1]; //Shift string one character left
            continue; //check next character
         }
         else if ( temp[i] == ']' ) { //End of comment
            bracket_balance--;        //Note one ] for balance
            int j=i;
            while (temp[j] !='\0') temp[j++]=temp[j+1]; //Shift string one character left
            continue; //check next character
         }
         else if (bracket_balance > 0) { //Means we are still in a comment
            int j=i;
            while (temp[j] !='\0') temp[j++]=temp[j+1]; //Shift string one character left
            continue; //check next character
         }
         i++; //If not square bracket or between square brackets check next character
      }
      /******************************************/
      if (temp[0]=='\0') continue; //If nothing in string read next word
      if (!strcmp(temp,"ntaxa")) {
         infile >> temp;
         ntaxa=atoi(temp);
         break;
         }
      //infile >> temp;
   }
   infile.close();

  //Check if brackets are balanced
  if (bracket_balance > 0) { cout << bracket_balance << " too many [. Square brackets need to be balanced. Quitting." << endl; return 1; }
  else if (bracket_balance < 0) { cout << bracket_balance * -1 << " too many ]. Square brackets need to be balanced. Quitting." << endl; return 1; }

  //Check if ntaxa could be read
  if (ntaxa == 0) { cout << "Could not find number of taxa (ntaxa). ntaxa need to be > 0. Quitting." << endl; return 1; }

/*** Asign memory to arrays depending on the number of taxa ***********/
   double *extinctpropvalues= new double[ntaxa];  //pointer to memory to store extinction proportions
   double *netdivvalues= new double[ntaxa];       //pointer to memory to store netdiversification rates
   int *nspecies= new int[ntaxa];                 //pointer to memory to store the number of species for the different taxa
   double *ages= new double[ntaxa];               //pointer to memory to store the ages of different taxa
   double *agessd= new double[ntaxa];             //pointer to memory to store the variance in ages of different taxa
   int *shapeGamma= new int[ntaxa];               //pointer to memory to store the shape parameter for the prior for different tax, not used yet
   double *meanGamma= new double[ntaxa];          //pointer to memory to store the mean for the prior for different tax, not used yet
   int *equalNDs= new int[ntaxa];                 //pointer to memory to store array of which taxa have same net diversification rate
   int *equalEPs= new int[ntaxa];                 //pointer to memory to store array of which taxa have same extinction proportion
/**********************************************************************/             
// Set all variables to default
   for (i=0; i<ntaxa; i++) { shapeGamma[i]=3;}
   for (i=0; i<ntaxa; i++) { meanGamma[i]=0.1;}
   for (i=0; i<ntaxa; i++) { netdivvalues[i]=0.1; }
   for (i=0; i<ntaxa; i++) { nspecies[i]=0;}
   for (i=0; i<ntaxa; i++) { ages[i]=0;}
   for (i=0; i<ntaxa; i++) { agessd[i]=0;}
   for (i=0; i<ntaxa; i++) { extinctpropvalues[i]=0.1; }
   for (i=0; i<ntaxa; i++) { equalNDs[i]=1;}
   for (i=0; i<ntaxa; i++) { equalEPs[i]=1;}

// Begin parsing file for commands and data
// for the function of each command se software manual
   if (ntaxa>0) { 
      ifstream infile (filename);
      if (!infile.is_open()) { cout << "Could not open \"" << argv[1] << "\". Quitting.\n"; return 1;}
      while(infile) {
         infile >> temp;
         //if (!infile) {break;}
         /*** Remove comments in square brackets ***/
         i=0; //Set counter to 0
         while ( temp[i] != '\0') { //read each character till end of string
            if (bracket_balance < 0) { cout << "Found an unmatched ]. Square brackets need to be balanced. Quitting." << endl; return 1; }
            else if ( temp[i] == '[' ) { //Start of comment to be removed
               bracket_balance++;   //Note one [ for balance
               int j=i;
               while (temp[j] !='\0') temp[j++]=temp[j+1]; //Shift string one character left
               continue; //check next character
            }
            else if ( temp[i] == ']' ) { //End of comment
               bracket_balance--;        //Note one ] for balance
               int j=i;
               while (temp[j] !='\0') temp[j++]=temp[j+1]; //Shift string one character left
               continue; //check next character
            }
            else if (bracket_balance > 0) { //Means we are still in a comment
               int j=i;
               while (temp[j] !='\0') temp[j++]=temp[j+1]; //Shift string one character left
               continue; //check next character
            }
            i++; //If not square bracket or between square brackets check next character
         }
         /******************************************/
         if (temp[0]=='\0') {continue;} //If nothing in string read next word

         else if (!strcmp(temp,"ntaxa")) { infile >> temp; continue; } //have already read ntaxa but need to state that this is a valid option

         else if (!strcmp(temp,"extinctpropvalues") || !strcmp(temp,"rel_extinct_rates")) {
            for (i=0; i<ntaxa; i++) {
               infile >> temp;
               extinctpropvalues[i] = atof(temp);
            }
         }
         else if (!strcmp(temp,"netdivvalues")) {
            for (i=0; i<ntaxa; i++) {
               infile >> temp;
               netdivvalues[i]=atof(temp);
            }
         }
         else if (!strcmp(temp,"ages")) {
            for (i=0; i<ntaxa; i++) {
               infile>>temp;
               ages[i]=atof(temp);
            }
         }
         else if (!strcmp(temp,"agessd")) {
            for (i=0; i<ntaxa; i++) {
               infile>>temp;
               agessd[i]=atof(temp);
            }
         }
         else if (!strcmp(temp,"nspecies")) {
            for (i=0; i<ntaxa; i++) {
               infile>>temp;
               nspecies[i]=atoi(temp);
            }
         }
         else if (!strcmp(temp,"shapeGamma")) {
            for (i=0; i<ntaxa; i++) {
               infile>>temp;
               shapeGamma[i]=atoi(temp);
            }
         }
         else if (!strcmp(temp,"meanGamma")) {
            for (i=0; i<ntaxa; i++) {
               infile>>temp;
               meanGamma[i]=atof(temp);
            }
         }
         else if (!strcmp(temp, "equalEP") || !strcmp(temp, "equalRER")) {
            infile>>temp;
            if (!strcmp(temp, "no")||!strcmp(temp, "No")) {
               for (i=0; i<ntaxa; i++) { equalEPs[i]=i+1;}  
               cout << "Extinction proportions set to all independent.\n";
            }
            else if (!strcmp(temp, "yes")||!strcmp(temp, "Yes")) {
               for (i=0; i<ntaxa; i++) { equalEPs[i]=1;}
               cout << "Extinction proportions set to all equal.\n";
            }
            else if (!strcmp(temp, "string")||!strcmp(temp, "String")) {
               for (i=0; i<ntaxa; i++) { 
                  infile>>temp;
                  equalEPs[i]=atoi(temp);
               }  
            }
            else { cout << "The command equalEP must be followed by: no, yes or string (see manual). Quitting." << endl; return 1; }
         }
         else if (!strcmp(temp, "equalND")) {
            infile>>temp;
            if (!strcmp(temp, "no")||!strcmp(temp, "No")) {
               for (i=0; i<ntaxa; i++) { equalNDs[i]=i+1;}  
               cout << "Net diversification rates set to all independent.\n";
               }
            else if (!strcmp(temp, "yes")||!strcmp(temp, "Yes")) {
               for (i=0; i<ntaxa; i++) { equalNDs[i]=1;}
               cout << "Net diversification rates set to all equal.\n";
               }

            else if (!strcmp(temp, "string")||!strcmp(temp, "String")) {
               for (i=0; i<ntaxa; i++) { 
                  infile>>temp;
                  equalNDs[i]=atoi(temp);
               }
            }
            else { cout << "The command equalND must be followed by: no, yes or string (see manual). Quitting." << endl; return 1; }
         }
         else if (!strcmp(temp, "ngen"))       {
            infile>>temp; 
            ngen=atoi(temp);
         }          
         else if (!strcmp(temp, "sampfreq"))   {
            infile>>temp;
            sampfreq=atoi(temp);
         }      
         else if (!strcmp(temp, "Bayes"))      {BayesOrLH=0;} //This is an old and no longer used command
         else if (!strcmp(temp, "model"))      {
            infile>>temp;
            model=temp[0];
         }         
         else if (!strcmp(temp, "windowidth") || !strcmp(temp, "windowwidth")) {
            infile>>temp; 
            windowwidth=atof(temp);
         }  
         else if (infile) { cout << "Could not understand \"" << temp << "\" in the input file. Please check the manual for available commands." << endl;
          return 1;
         }
      }
      infile.close();
   }
   else { cout << "Need to know how many different groups (ntaxa) there are.\n"; return 1;}   

   //Check if brackets are balanced
   if (bracket_balance > 0) { cout << bracket_balance << " too many [. Square brackets need to be balanced. Quitting." << endl; return 1; }
   else if (bracket_balance < 0) { cout << bracket_balance * -1 << " too many ]. Square brackets need to be balanced. Quitting." << endl; return 1; }

   /*** Test the different variables ***/
   for (i=0; i<ntaxa; i++) {
      if (nspecies[i]<1) {
         cout << "Taxa " << i+1 << " has no species. Each taxon must have at least 1 species (given by nspecies). Quitting." << endl << endl;
         return 1;
      }
   }
   for (i=0; i<ntaxa; i++) {
      if (ages[i]<=0) {
         cout << "Taxa " << i+1 << " has age zero (or less). Each taxon must have the age > 0 (given by ages). Quitting." << endl << endl;
         return 1;
      }
   }
   if (ngen<=0) {cout << "The number of generations (ngen) to run the MCMC is set to zero (or less). Must be a positive integer. Quitting." << endl; return 1;}
   //Test so that the sampling is more frequent than the number of generations, else quit
   if (sampfreq>ngen || sampfreq <= 0) {cout << "The sample frequency (sampfreq) needs to be less than the number of generations (ngen) but larger than zero to get any estimation of the posterior distribution.\n"; return 1;}

   //Test so a viable model is given, else quit
   if (!(model=='d' || model=='y')) { cout << "You set model: " << model << ". y (Yule; pure-birth) or d (birth-death) are available.\n"; return 1;}

   for (i=0; i<ntaxa; i++) {
      if (agessd[i]<0) {
         agessd[i]=0;
         cout << endl << "WARNING!!! The standard deviation of taxon " << i+1 << " is set to less than zero. Changed to default value (0)." << endl << endl;
      }
   }
   for (i=0; i<ntaxa; i++) {
      if (equalNDs[i]<=0) {
         equalNDs[i]=1;
         cout << endl << "WARNING!!! Net diversification rate category (equalND) of taxon " << i+1 << " set to zero (or less). Must be posetive integer. Changed to default value (1)." << endl << endl;
      }
   }
   for (i=0; i<ntaxa; i++) {
      if (model != 'y' && equalEPs[i]<=0) {
         equalEPs[i]=1;
         cout << endl << "WARNING!!! Extinction proportion category (equalEP) of taxon " << i+1 << " set to zero (or less). Must be posetive integer. Changed to default value (1)." << endl << endl;
      }
   }
   for (i=0; i<ntaxa; i++) {
      if (meanGamma[i]<=0) {
         meanGamma[i]=0.1;
         cout << endl << "WARNING!!! Mean for net diversification prior (meanGamma) of taxon " << i+1 << " set to zero (or less). Changed to default value (0.1)." << endl << endl;
      }
   }
   for (i=0; i<ntaxa; i++) { 
      if (shapeGamma[i]<=0) { 
         shapeGamma[i]=3;
         cout << endl << "WARNING!!! Shape for net diversification (shapeGamma) prior of taxon " << i+1 << " set to zero (or less). Changed to default value (3)." << endl << endl;
      }
   }
   for (i=0; i<ntaxa; i++) {
      if (netdivvalues[i]<=0) {
         netdivvalues[i]=0.1;
         cout << endl << "WARNING!!! Starting value for net diversification rate (netdivvalues) of taxon " << i+1 << " set to zero (or less). Changed to default value (0.1)." << endl << endl;
      }
   }
   for (i=0; i<ntaxa; i++) {
      if (model != 'y' && (extinctpropvalues[i] <0 || extinctpropvalues[i] >1)) {
         extinctpropvalues[i]=0.01;
         cout << endl << "WARNING!!! Starting value for extinction proportion (extinctpropvalues) of taxon " << i+1 << " set to less than zero or more than 1. Changed to default value (0.01)." << endl << endl;
      }
   }
   if (windowwidth <=0) {
      windowwidth=0.2;
      cout << endl << "WARNING!!! The window width (windowidth) is set to less than zero. Changed to default (0.2)." << endl << endl; 
   }
/**************************************************/

// Print the settings of the run
   cout << "Number of taxa/groups (ntaxa): " << ntaxa << endl;
   cout << "Number of generations (ngen): " << ngen << endl;
   cout << "Sample frequency (sampfreq): " << sampfreq << endl;
   cout << "The model used (model): ";
   if (model == 'y') cout << "Yule" << endl;
   else if (model == 'd') cout << "birth-death" << endl;
   else cout << "?" << endl;
   cout << "The window width multiplier (windowidth): " << windowwidth << endl << endl;
   cout << "Taxa\tNumber_of_species\tAge_mean\tAge_sd\tNet_div_rate_cat\t";
   if (model == 'd') cout << "Rel_ext_rate_cat\t";
   cout << "Net_div_prior_mean\tNet_div_prior_shape\tNet_div_start";
   if (model == 'd') cout << "\tRel_ext_rate_start";
   cout << endl;
   for (i=0; i<ntaxa; i++) {
      cout << "Taxa_" << i+1 << "\t" << nspecies[i] << "\t" << ages[i] << "\t" << agessd[i] << "\t" << equalNDs[i] << "\t";
      if (model == 'd') cout << equalEPs[i] << "\t";
      cout << meanGamma[i] << "\t" << shapeGamma[i] << "\t" << netdivvalues[i];
      if (model == 'd') cout << "\t" << extinctpropvalues[i];
      cout << endl;
   }
   cout << endl;
   cout << "Starting MCMC" << endl << "*******************************************************" << endl << endl;

/*** Start the MCMC ***/
   if (BayesOrLH==0) { acceptrate= MCMC(netdivvalues, equalNDs, extinctpropvalues, equalEPs, nspecies, ages, agessd, ngen, sampfreq, model, windowwidth, meanGamma, shapeGamma, ntaxa); }

   if (acceptrate<0) { cout << " Hey, something went wrong in the MCMC function... Quitting... Try to have a nice day despite of this... Sorry... Bye...\n"; return 1; }

/***********************************************************/
   
/**************************************************************/
// Clean up
   delete [] extinctpropvalues;
   delete [] netdivvalues;
   delete [] nspecies;
   delete [] ages;
   delete [] agessd;
   delete [] shapeGamma;
   delete [] meanGamma;
   delete [] equalNDs;
   delete [] equalEPs;
}

/*********************************************************/
/*** Functions *******************************************/
/*********************************************************/

void help () {  //Print the help
   cout << "This program takes a text (ASCII) file as the only argument and reads the commands in that file.\n";
   cout << "\nThe file needs to give the number of taxa as ntaxa followed by space and the number of taxa (e.g. ntaxa 5).\n";
   cout << "It also needs the number of species of each taxon given as nspecies followed by a space separated vector of the\n   number of species in each taxon (e.g. nspecies 10 20 30 40 50).\n";
   cout << "In addition you should give the age of each taxon (in the same order as you gave the number of species) as ages\n   followed by a space separated vector (e.g. ages 50 40 30 20 10).\n";
   cout << "If you want to use a prior distributions instead of fixed ages, you should give the standard deviation of the\n   distributions using agessd (e.g. agessd 5 4 3 2 1). At present the normal distribution is the only available distribution.\n"; 
   cout << "\nOther arguments are:\n";
   cout << "ngen - number of generations for the Markov chain (e.g. ngen 100000)\n";
   cout << "sampfreq - the sample frequency (e.g. sampfreq 1000).\n";
   cout << "windowidth - the width of the proposal window. Should be adjusted to get an as efficient chain as possible\n   (acceptance rates 0.2-0.4; e.g. windowidth 0.3).\n";
   cout << "equalND - Equal net diversification in all taxa, Yes or No (e.g. equalND Yes). You may also give string\n";
   cout << "   followed by a space-separated vector with the same integer numbers for taxa with the same net\n";
   cout << "   diversification rate (equalND string 1 1 2 2 3).\n";
   cout << "equalEP - Equal extinction proportion, Yes or No (e.g. equalEP Yes). You may also give string followed\n";
   cout << "   by a space-separated vector with the same integer numbers for taxa with the same extinction proportion\n";
   cout << "   (equalEP string 1 2 1 3 4).\n"; 
   cout << "model d - Allows the extinction proportion to be >0, i.e. use the Birth-Death model.\n";
   cout << "netdivvalues - Starting values for the net diversification rates, given for each taxon as a space separated\n   string.\n";
   cout << "extinctpropvalues - Starting values for the extinction proportions, given for each taxon as a space separated\n   string.\n";
   cout << "meanGamma - The mean of the net diversification rate prior, given for each taxon as a space separated string.\n";
   cout << "shapeGamma - The shape of the net diversification rate prior (integer number), given for each taxon as a space\n   separated string.\n";

   cout << endl;
   cout << "The output is a Markov chain as a tab separated table that can be analyzed in programs such as R or Excel (if deleting the lines before the heading of the Markov chain).\n";

   return;
}

void MCMC_summarize (char filename[], int burnin) {
    char temp[100]; //to store what is read from file
    int ntaxa; //stor number of taxa
    double sum_LH=0; //the sum of the likelihood values (LH)
    int number_LH=0; //the number of LH values read
    int n_generations = 0; //the number of generations of the chain
    char generations_err = 'F'; //if the given number of generations differes from the read number warn user
    int sample_frequency = 0; //store the sample frequency
    int iteration=0; //store which generation is being read
    char netdiv_err = 'F'; //flag for errors in reading net diversification rate
    char ext_prop_err = 'F'; //flag for errors in reading extinction proportions
    char age_err = 'F'; //Flag for suspicious age values
    double acceptance; //to store acceptance rate
    char model = 'y'; //Yule (y) or birth-death (b) model?
    char var_age = 'F'; //ages treated as priors False or True
    int ncolumns;

    ifstream infile (filename);
    if (!infile.is_open()) { cout << "Could not open \"" << filename << "\". Quitting.\n"; exit(1);}
    { //start scope for first round of parsing file, also scope for data_trees (see below)
        data_tree netdiv; //data_tree to store net diversification values
        data_tree ext_prop; //data_tree to stor extinction proportion values
        data_tree age; //data_tree to store ages

        while (infile) {
            infile >> temp;
            if ( !strcmp(temp, "(ntaxa):") ) { infile >> temp; ntaxa = atoi(temp); }
            else if ( !strcmp(temp, "(ngen):") ) { infile >> temp; n_generations = atoi(temp); }
            else if ( !strcmp(temp, "(model):") ) { 
                infile >> temp;  
                if ( !strcmp(temp, "Yule") ) model = 'y';
                else if ( !strcmp(temp, "birth-death") ) model = 'b';
            }
            else if ( !strcmp(temp, "(sampfreq):") ) { infile >> temp; sample_frequency = atoi(temp); }
            else if ( temp[0] == 'T' && temp[1] == 'a' && temp[2] == 'x' && temp[3] == 'a' && temp[4] == '_' ) {
                infile >> temp; infile >> temp; infile >> temp;
                if ( atof(temp) > 0 ) var_age = 'T';
            }
            else if ( !strcmp(temp, "Acceptance_rate") ) {
                int column = 1;
                ncolumns = 5 + ntaxa;
                if ( model == 'b' ) ncolumns += ntaxa;
                if ( var_age == 'T' ) ncolumns += ntaxa;
            
                while (infile) {
                    infile >> temp;
                    if (column > ncolumns) column = 1;
                    if (column == 1) { 
                        iteration = atoi(temp);
                        if ( iteration > n_generations ) generations_err = 'T';
                    }
                    else if ( column == 2 && iteration > burnin ) {
                        sum_LH += 1/atof(temp);
                        ++number_LH;
                    }
                    else if ( column == 5 && iteration > burnin ) {
                        netdiv.insert_data(atof(temp)); //add value to the netdiv data_tree
                        if (atof(temp) < 0) netdiv_err = 'T'; //if negative value issue warning
                        if ( model == 'b' ) {
                            infile >> temp;
                            ++column;
                            ext_prop.insert_data(atof(temp)); //add value to extinction proportion tree
                            if (atof(temp) < 0 || atof(temp) > 1) ext_prop_err = 'T'; //if extinction proportion is over 1 or negative issue warning
                        }
                        if ( var_age == 'T' ) {
                            infile >> temp;
                            ++column;
                            age.insert_data(atof(temp)); //add value to extinction proportion tree
                            if ( atof(temp) < 0 ) age_err = 'T'; //if extinction proportion is over 1 or negative issue warning
                        }
                    }
                    else if ( column == ncolumns ) acceptance = atof(temp);
                    ++column;
                }
            }
        }
        cout << "Summarizing MCMC" << endl << "******************************" << endl;
        if ( generations_err == 'T' ) cout << "WARNING!!! The number of generations read does not agree with what was given at start!!!" << endl << endl; //if warning issued for the generation number read print it
        cout << "The MCMC to summarize had " << n_generations << " generations, and were sampled every " << sample_frequency << " generations." << endl; //basic statistics on the chain
        cout << "Acceptance rate: " << acceptance << endl;
        cout << "Harmonic mean of the log likelihood: " << number_LH/sum_LH << endl << endl; //read the text
        cout << "Taxa 1:" << endl;
        if ( netdiv_err == 'T' ) cout << "WARNING!!! Net diversification rate had strange values!!!" << endl << endl; //if net diversification warning issued print it
        if ( ext_prop_err == 'T' ) cout << "WARNING!!! Relative extinction rate had strange values!!!" << endl << endl; //if extinction proportion warning issued print it
        if ( age_err == 'T' ) cout << "WARNING!!! Age estimates had strange values!!!" << endl << endl; //if age warning issued print it
        cout << "Mean net diversification rate: " << netdiv.mean() << ", median: " << netdiv.median() << ", standard deviation: " << netdiv.sd();
        cout << ", 95\% credibility interval: " << netdiv.get_value(int(netdiv.n_nodes()*0.025)) << "-" << netdiv.get_value(int(netdiv.n_nodes()*0.975));
        cout << ", min: " << netdiv.min() << ", max: " << netdiv.max() << endl;
        if ( model == 'b' ) {
            cout << "Mean relative extinction rate: " << ext_prop.mean() << ", median: " << ext_prop.median() << ", standard deviation: " << ext_prop.sd();
            cout << ", 95\% credibility interval: " << ext_prop.get_value(int(ext_prop.n_nodes()*0.025)) << "-" << ext_prop.get_value(int(ext_prop.n_nodes()*0.975));
            cout << ", min: " << ext_prop.min() << ", max: " << ext_prop.max() << endl;
        }
        if ( var_age == 'T' ) {
            cout << "Mean estimated age: " << age.mean() << ", median: " << age.median() << ", standard deviation: " << age.sd();
            cout << ", 95\% credibility interval: " << age.get_value(int(age.n_nodes()*0.025)) << "-" << age.get_value(int(age.n_nodes()*0.975));
            cout << ", min: " << age.min() << ", max: " << age.max() << endl;
        }
    } //end of first round parsing the file, also end of scope for data_trees for age, netdiv, and ext_prop
   for ( int i = 1; i < ntaxa; ++i ) {
        infile.clear();              // forget we hit the end of file
        infile.seekg(0, ios::beg);   // move to the start of the file
        data_tree netdiv;
        data_tree ext_prop;
        data_tree age;
        netdiv_err = 'F';
        ext_prop_err = 'F';
        age_err = 'F';
        while (infile) {
            infile >> temp;
            if ( !strcmp(temp, "Acceptance_rate") ) {
                int column = 1;
                int netdiv_column = 5+i;
                if ( model == 'b' ) netdiv_column += i;
                if ( var_age == 'T' ) netdiv_column += i;
                while (infile) {
                    infile >> temp;
                    if (column > ncolumns) column = 1;
                    if (column == 1) {
                        iteration = atoi(temp); 
                        if ( iteration > n_generations ) generations_err = 'T';
                    }
                    else if ( column == netdiv_column && iteration > burnin ) {
                        netdiv.insert_data(atof(temp)); //add value to the netdiv data_tree
                        if (atof(temp) < 0) netdiv_err = 'T'; //if negative value issue warning
                        if ( model == 'b' ) {
                            infile >> temp;
                            ++column;
                            ext_prop.insert_data(atof(temp)); //add value to extinction proportion tree
                            if (atof(temp) < 0 || atof(temp) > 1) ext_prop_err = 'T'; //if extinction proportion is over 1 or negative issue warning
                        }
                        if ( var_age == 'T' ) {
                            infile >> temp;
                            ++column;
                            age.insert_data(atof(temp)); //add value to extinction proportion tree
                            if ( atof(temp) < 0 ) age_err = 'T'; //if extinction proportion is over 1 or negative issue warning
                        }
                    }
                ++column;
                }
            }
        }
        cout << "Taxa " << i+1 << ":" << endl; 
        if ( netdiv_err == 'T' ) cout << "WARNING!!! Net diversification rate had strange values!!!" << endl << endl; //if net diversification warning issued print it
        if ( ext_prop_err == 'T' ) cout << "WARNING!!! Relative extinction rate had strange values!!!" << endl << endl; //if extinction proportion warning issued print it
        if ( age_err == 'T' ) cout << "WARNING!!! Age estimates had strange values!!!" << endl << endl; //if age warning issued print it
        cout << "Mean net diversification rate: " << netdiv.mean() << ", median: " << netdiv.median() << ", standard deviation: " << netdiv.sd();
        cout << ", 95\% credibility interval: " << netdiv.get_value(int(netdiv.n_nodes()*0.025)) << "-" << netdiv.get_value(int(netdiv.n_nodes()*0.975));
        cout << ", min: " << netdiv.min() << ", max: " << netdiv.max() << endl;
        if ( model == 'b' ) {
            cout << "Mean relative extinction rate: " << ext_prop.mean() << ", median: " << ext_prop.median() << ", standard deviation: " << ext_prop.sd();
            cout << ", 95\% credibility interval: " << ext_prop.get_value(int(ext_prop.n_nodes()*0.025)) << "-" << ext_prop.get_value(int(ext_prop.n_nodes()*0.975));
            cout << ", min: " << ext_prop.min() << ", max: " << ext_prop.max() << endl;
        }
        if ( var_age == 'T' ) {
            cout << "Mean estimated age: " << age.mean() << ", median: " << age.median() << ", standard deviation: " << age.sd();
            cout << ", 95\% credibility interval: " << age.get_value(int(age.n_nodes()*0.025)) << "-" << age.get_value(int(age.n_nodes()*0.975));
            cout << ", min: " << age.min() << ", max: " << age.max() << endl;
        }
    }
    infile.close();
}

/*** Probability ***/
// Calculate the log likelihood
double LikeValue (int *nspecies, double *ages, double *netdivvalues, double *extinctpropvalues, int ntaxa ) {
   double Lh=0;
   for (int i=0; i<ntaxa; i++) {
   double beta=((exp(netdivvalues[i]*ages[i])-1)/(exp(netdivvalues[i]*ages[i])-(1/(((1-extinctpropvalues[i])/extinctpropvalues[i])+1))));
   Lh += log((1-beta))+(log(beta)*(nspecies[i]-1));
   }
   return Lh;
}

// Calculate the prior probabilities
double PriorProb (char distribution[10], double *x, int a, double b, int ntaxa) {
   double pb=0;
   for (int i=0; i<ntaxa; i++) {
      if (x[i]<0) return log(0);
      if (!strcmp(distribution,"gamma")) pb+=log(gamma(x[i],a,1/b));            //a represents shape and b scale
      if (!strcmp(distribution,"uniform")) pb+=log(uniform(x[i],a,b));          //a represents lower bound and b higher bound
   }
   return pb;
}

double PriorProb (char distribution[10], double *x, int *a, double *b, int ntaxa) {
   double pb=0;
   for (int i=0; i<ntaxa; i++) {
      if (x[i]<0) return log(0);
      if (!strcmp(distribution,"gamma")) pb+=log(gamma(x[i],a[i],1/b[i]));                 //a represents shape and b scale
      else if (!strcmp(distribution,"uniform")) pb+=log(uniform(x[i],a[i],b[i]));          //a represents lower bound and b higher bound
   }
   return pb;
}

// Calculates the prior probability of the ages
double AgePriorProb (char distribution[10], double *x, double *mean, double *sd, int ntaxa) {
   double pb=1;
   for (int i=0; i<ntaxa; i++) {
      if (!strcmp(distribution,"normal")) pb+=log(normalprob(x[i],mean[i],sd[i]));
   }
   return pb;
}
/***/

/*** Distributions ***/
/* Probability distributions */
double gamma(double x, int k, double b) { //Takes value for which to calculate probability (x), shape (k), and inv. scale (b).

   return pow(x,k-1)*(exp(-x/(1/b))/(pow(1/b,k)*factorial(k-1)));

}

double uniform (double x, double a, double b) {
   if (a > b) return 0;
   if (x < a) return 0;
   if (x > b) return 0;
   else return 1/(b-a);
}

double normalprob (double x, double mean, double sd) {
   if (sd==0 && x==mean) return 1;
   else if (sd==0 && x!=mean) return 0; 
   else return (1/sqrt(2*3.14159265*pow(sd,2)))*exp(-(pow(x-mean,2)/(2*pow(sd,2))));
}

/***/

/* Random distributions */
//Generate normal distributed random value using the Marsaglia polar method
double normal(double a) {
   double s=2;
   double x;
   double y;

   while (!(s<1)) {
      x = (rand()/((double)RAND_MAX/2))-1;
      y = (rand()/((double)RAND_MAX/2))-1;
      s = (pow(x,2)+pow(y,2));
   }
   return x*sqrt(-2*log(s)/s)*a;

}
/***/
/***/

/*** Math ***/
int factorial(int x){
 int i=x;
 x=1;
 while (i) x=x*i--;
 return x; 

}
/***/

/**********************/
/*** MCMC functions ***/
/**********************/

int MCMC(double *netdivvalues, int *equalNDs, double *extinctpropvalues, int *equalEPs, int *nspecies, double *ages, double *agessd, int ngen, int sampfreq, char model, double multiplier, double *meanGamma, int *shapeGamma, int ntaxa) {
/** Variables *******************************************************/
   double stateLH;                     //Variables to store the probabilities of state i-1
   double statePr=1; 
   double propLH;                      //Variables to store the probabilities of the proposed state i
   double propPr;
   double *propb= new double[ntaxa];   //Variables to store the proposed parameters
   double *propd= new double[ntaxa];
   double *propAges= new double[ntaxa];
   double *Agesmean= new double[ntaxa];
   int accept=0;                       //Counter of number of acceptet proposals
   int i=0;                            //Counter
   int updateturn=0;                   //Counter to se what value to uppdate
   int varages=0;                      //Flag to see if ages are given exactly or as priors
   char gamm[] = "gamma";              //Strings to path to function
   char norm[] = "normal";
   char unif[] = "uniform";
/*******************************************************************/
   // Check if prior distributions are used for ages
   for (i=0; i<ntaxa; i++) {
      propAges[i]=ages[i];
      if (agessd[i]>0) varages=1;
      Agesmean[i]=ages[i]; 
   }

   if (model=='y') for (i=0; i<ntaxa; i++) extinctpropvalues[i]=0; //set d values to 0 if yule model

   //Set the rates to the same for the taxa with same rate category
   SetSameRates(netdivvalues,equalNDs,ntaxa);
   SetSameRates(extinctpropvalues,equalEPs,ntaxa);
   
   //Calculate likelihood for start values
   stateLH=LikeValue(nspecies, ages, netdivvalues, extinctpropvalues, ntaxa);

   //Calculate the prior probability for the starting values
   if (model=='y') for (i=0; i<ntaxa; i++) statePr=PriorProb(gamm,netdivvalues,3,0.1,ntaxa); //If Yule model

   else statePr=PriorProb(gamm,netdivvalues,shapeGamma,meanGamma,ntaxa)+PriorProb(unif,extinctpropvalues,0,1,ntaxa); //If birth-death model  

   if (varages) statePr += AgePriorProb(norm, ages, Agesmean, agessd, ntaxa); //If priors for ages

/*** Print header for MCMC *********/
   cout << "Generation\t"; 
   cout << "Log_Lh\tPrior_P\tPost_P";
   for (i=0; i<ntaxa; i++) {cout << "\tNetdiv_" << i+1; if (model != 'y') cout << "\tRel_ext_rate_" << i+1; if (varages) cout << "\tAge_" << i+1;}
   cout << "\tAcceptance_rate";
   cout << endl;
/***********************************/

/*** Print values for start point in chain **************************/
   cout << "0\t";
   cout << stateLH << "\t" << statePr << "\t" << stateLH+statePr;
   for (i=0; i<ntaxa; i++) {cout << "\t" << netdivvalues[i]; if (model != 'y') cout  << "\t" << extinctpropvalues[i]; if (varages) cout << "\t" << ages[i];}
   cout << "\t0";
   cout << endl;
/********************************************************************/



/*** Start the chain ****/
   for (i=1; i<=ngen; i++) {

      /*** null the proposed state **************/
      for (int j=0; j<ntaxa; j++) {propb[j]=netdivvalues[j];}                       //Set b(irth) to new state
      if (model!='y') for (int j=0; j<ntaxa; j++) {propd[j]=extinctpropvalues[j];}  //If not Yule model set d(eth) to new state
      if (varages) for (int j=0; j<ntaxa; j++) {propAges[j]=ages[j];}               //If ages are a distribution set ages to new state
      /******************************************/

      /*** Proposal of new states ****************************/
      int nNDparameters=nparameters(equalNDs,ntaxa);
      for (int j=0; j<nNDparameters; j++) {
         propb[WhichParameter(j+1,equalNDs,ntaxa)]=netdivvalues[WhichParameter(j+1,equalNDs,ntaxa)]-normal(multiplier); 
      }
      if (model!='y') { 
         int nEPparameters=nparameters(equalEPs,ntaxa);
         for (int j=0; j<nEPparameters; j++) {
            propd[WhichParameter(j+1,equalEPs,ntaxa)]=extinctpropvalues[WhichParameter(j+1,equalEPs,ntaxa)]-normal(multiplier);
         } 
      }
      
      SetSameRates(propb,equalNDs,ntaxa);
      SetSameRates(propd,equalEPs,ntaxa);
      /*******************************************************/

      /*** If ages are distributions propose new ages ********/
      if (varages) for (int j=0; j<ntaxa; j++) { propAges[j] = ages[j]-normal(agessd[j]*multiplier); }
      /*******************************************************/

      /*** Calculating probability of proposed state ************/
      propLH = LikeValue(nspecies, propAges, propb, propd, ntaxa);
      if (model=='y') propPr=PriorProb(gamm,propb,shapeGamma,meanGamma,ntaxa);
      else propPr=PriorProb(gamm,propb,shapeGamma,meanGamma,ntaxa)+PriorProb(unif,propd,0,1,ntaxa);
      if (varages) propPr += AgePriorProb(norm, propAges, Agesmean, agessd, ntaxa);
      /**********************************************************/

      /*** Check if to accept the proposed state and if so uppdate it ***/
      if (exp((propLH+propPr)-(stateLH+statePr)) > (rand()/(double)RAND_MAX)) { 
         for (int j=0; j<ntaxa; j++) { netdivvalues[j]=propb[j];}                         //Set b(irth) to new state
         if (model!='y') {for (int j=0; j<ntaxa; j++) { extinctpropvalues[j]=propd[j];}}  //If not Yule model set d(eth) to new state
         if (varages) {for (int j=0; j<ntaxa; j++) { ages[j]=propAges[j]; }}              //If ages are distributions set ages to new state
         stateLH=propLH;           //Set probabilities for printing and
         statePr=propPr;           //not having to calculate again for next iteration
         accept++;                 //count # acceptens
      }
      /******************************************************************/

      /*** Print the values of the chain with sampfreq intervals ******/
      if (fmod(i,sampfreq)==0) {
         cout << i << "\t"; 
         cout << stateLH << "\t" << statePr << "\t" << stateLH+statePr; 
         for (int j=0; j<ntaxa; j++) {cout << "\t" << netdivvalues[j]; if (model != 'y') cout  << "\t" << extinctpropvalues[j]; if (varages) cout << "\t" << ages[j];} 
         cout << "\t" << accept/(double)i;
         cout << endl;
      }
      /*****************************************************************/
   }
/************************/
// Cleaning up
   delete [] propb;
   delete [] propd;
   delete [] propAges;
   delete [] Agesmean;


   return accept;   //Return the acceptence rate

}

/***********************************************************************************************
************************************************************************************************
***********************************************************************************************/

// Gives where in the string the n:th parameter is
int WhichParameter (int value, int *equalvalues, int ntaxa) {
   int k=0;
   for (int i=0; i<ntaxa; i++) {
      for (int j=0; j<=i; j++) {
         if (i==j) {k++; if (k==value) {return i;}}
         else if (equalvalues[i] == equalvalues[j]) {break;}
      }
   }
}

//Sets the values with the same category to the first value that apear in the category
void SetSameRates (double *values, int *equalvalues, int n) {
   int *identity = new int[n];
   for (int i=0; i<n; i++) {
      for (int j=0; j<=i; j++) {
         if (equalvalues[i] == equalvalues[j]) {identity[i]=j; break;}
      }
   }
   for (int i=0; i<n; i++) {values[i]=values[identity[i]]; }

   delete [] identity;
}

//Rearanges the values in an array according to given array
void rearrange (double *values, int *rearrangearray, int n) {
   double *copy = new double[n];
                               
   for (int i=0; i<n; i++) { copy[i]=values[i]; }
   for (int i=0; i<n; i++) { values[i]=copy[rearrangearray[i]-1];}
   delete [] copy;
   return;
}
/***/

// Calculate number of parameters
int nparameters (int *model, int n) {
   int npar=0;

   for (int i=0; i<n; i++) {
      for (int j=0; j<=i; j++) {
         if (j==i) { npar++; break;}
         if (model[i]==model[j]) {break;}
      }
   }
   return npar;
}
/***/
