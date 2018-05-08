/**********                Transform.hpp  			*************
 *********               copy right Lei Zhou            *************/
/* this programm is intended to trnasform long (>100K) genomic sequences
into realitively shorter sequences (eg. 3k) to faciliate searching
                                                                    */

#include <condefs.h>
#include <iostream.h>
#include <stdio.h>
#include <conio.h>
#include <fstream.h>
#include <iomanip.h>

#define FALSE 0
#define TRUE 1
#define OVERLAP 180         //To compensate breaking of sequence

#define ARRSIZ 256
#define MAXLINE 512

#define MAXSEQ  100000000      // max len of input seq
#define OUTSEQL 10000        //len of out put seq

#pragma hdrstop

//------ Globals -------------
int num_gnt;
char dbname[MAXLINE];

struct GSEQ
{ char name[MAXLINE], seq [MAXSEQ], rseq[MAXSEQ];
  int  len, hit_p;
  bool reverstrand;
};
//----------FUNCTIONS & PROCEDURES--------


GSEQ* Get_Lib_Seq ();        //read from mutiple FASTA, buid GSEQ
void Reverse ( GSEQ * Gseq);
void Transform (GSEQ* gseq);
void Output_seq (int start, char name[], bool reverse, char seq[]);

//-------FILE HANDLER-------------------

//ifstream lfile ("E:\\Flysequence\\geno\\P1test.txt");
ifstream lfile ("S:\\ins\\beetle\\TCgenome2.0");
//ifstream lfile ("E:\\Worm\\allcmid.txt");

ofstream seq_file ("S:\\ins\\beetle\\Tcg2.10k");

