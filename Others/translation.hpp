#include <condefs.h>
#include <iostream.h>
#include <stdio.h>
#include <conio.h>
#include <fstream.h>
#include <math.h>
#include <iomanip.h>

#define MAXLINE 512
#define MAXSEQ  50000

#pragma hdrstop

char dbname[MAXLINE];
bool REVERSE;
/*
    DSEQ--struct hold the DNA sequence
*/
struct DSEQ
{ char name[MAXLINE], seq [MAXSEQ], r_seq[MAXSEQ];
  int  len;
  bool reverse;
};

void Reverse (DSEQ* dseq);
void Translation (char name[], char seq[]);
DSEQ * Get_Lib_Seq ();
char Translate(char s[]); 		//translate one aa per 3 bp.
int Convert (char c);				//convert (A,C,G,T) to (0,1,2,3),
									//called from Translate


//ifstream infile ("D:\\GeneWork\\Search Result\\card\\SEQ.EST08.txt");
ifstream infile ("DNA.txt");
//ifstream infile ("E:\\Flysequence\\geno\\transf3k.txt");
//ifstream infile ("E:\\Human\\Unigene\\HsUni.txt");
//ifstream infile ("E:\\Flysequence\\Est\\ESTtest.txt");
//ofstream outfile ("S:\\human\\gd\\hs.gd.pep");
ofstream outfile ("pep.txt");
