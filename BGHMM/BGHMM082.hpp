/**********                   BGHMM.hpp    			*************
 *********               copy right Lei Zhou            **************/


#include <condefs.h>
#include <iostream.h>
#include <stdio.h>
#include <conio.h>
#include <fstream.h>
#include <math.h>
#include <iomanip.h>
#include <time.h>

#define FALSE 0
#define TRUE 1
#define MAXMATRIX 60  //length of motifs
#define ALPHEBIT 20
#define ARRSIZE 256
#define MAXMOTIFS 30
#define MAXLINE 512
#define MAXHITS 50
#define MAXSEQ  41000      //greatly affacts searching time
#define MAXOUTSEQ 15000

/*searh parameter, could changed to input driven
																	*/
#define CHUNKS 2000
#define MARGIN 30                    //for report & Check_Redundancy

#define NAMESPACE  20             //for extra entry of same sequence
#pragma hdrstop

//------ Globals -------------

bool generate_seq_file, genomic;
char dbname[MAXLINE];
ifstream lfile, m_file;
ofstream seq_file, rfile;

/* MOTIF--holds matrix for one motif;
    MOTIFLIST--all motifs
                                                                    */
struct MOTIF
{ char name[MAXLINE];
  float matrix[MAXMATRIX][ARRSIZE];
  float cut_off;
  int len;
  int Check_up;  // for checking stop codons
  int Check_down;
};

struct MOTIFLIST
{ MOTIF* motifs[MAXMOTIFS];
  int num_motifs;
};
MOTIFLIST* mlp;

/* PSEQ--struct hold the protein sequence
                                                                    */
struct PSEQ
{ char name[MAXLINE], seq [MAXSEQ];
  int  len;
};

/* SCORE--holding scores as parsing the seq against the HMM
                                                                    */
struct SCORES
{float seq_scores[MAXMOTIFS][MAXSEQ];
};


/* HIT--holds one potential hit;
    HITLIST--hold all qualified hits for one motif
    HITTABLE--holds all data for output
                                                                    */
struct HIT
{   int motif_n;
    PSEQ *pseq;
    float score;
    int hit_p;      /*postion of first match*/
    int more_entries;
    char ext_names[MAXLINE];
};

struct HITLIST
{   int motif_n;
    char motif_name;
    HIT* hitlist[MAXHITS];
    int n_hits;
};

struct HITTABLE
{   HITLIST* lp[MAXMOTIFS];
    int num_motifs;
};
HITTABLE* htp;

//----------FUNCTIONS & PROCEDURES--------

MOTIFLIST* Input_Matrix ();     //generate matrixlis from input file
PSEQ* Get_Lib_Seq ();        //read from mutiple FASTA, buid PSEQ

SCORES* Score_it
(PSEQ* pseq, MOTIFLIST* mlp);

/*chech stopcodon before scoreing procedure */
int Check_Stopcodon (char seq[]);

/*sort score to get the highest score, with consideration of ORF ect.*/
void Sort_Scores
(  SCORES* scores,
   PSEQ* pseq,
   MOTIFLIST* mlp);

/* check to see if needs to add to hitlist */
bool Check_Hit (HIT* hit);
bool Check_Redundancy (HIT* new_hit, HIT* old_hit);
/* add hit to hitlist[motif_n]*/
void Insert_Hit(HIT* hit, int insertp);
void Output
(char* search_name, HITTABLE* htp, MOTIFLIST* mlp);
void Output_seq
(HITTABLE* htp, MOTIFLIST* mlp);


//-------FILE HANDLER-------------------


//ifstream m_file ("E:\\Bagua\\matrix\\celldeath0409.hmm.txt");
//ifstream m_file ("E:\\Bagua\\matrix\\RGH\\IBM.hmm.txt");

//ifstream lfile ("E:\\Flysequence\\Est\\EstP.txt");
//ifstream lfile ("S:\\fish\\tet\\pep.txt");
//ifstream lfile ("S:\\protein\\swp.txt");
//ifstream lfile ("S:\\ins\\silk\\g10k.pep.txt");

//ofstream rfile ("E:\\Bagua\\output\\IBM_swp.txt");
//ofstream seq_file ("E:\\Bagua\\output\\IBM_swp.seq.txt");

//ofstream rfile ("S:\\ins\\silk\\IBM_SWg.txt");
//ofstream seq_file ("S:\\ins\\silk\\IBM_SWg.seq.txt");
