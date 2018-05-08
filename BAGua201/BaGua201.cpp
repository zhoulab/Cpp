/* FROM BaGua20: in order to get single file input for motifs as well as
dbfiles insteadof on screen inpout, BaGua 2.01 using Search
( char*dfname) to replace part of the main code and the
Get_Lib_Seq () function; also changed are the various file handling
procedures. */

#include <vcl\condefs.h>
#include <iostream.h>
#include <stdio.h>
#include <conio.h>
#include <fstream.h>
#include <math.h>
#include <iomanip.h>

#define FALSE 0
#define TRUE 1
#define MAXMATRIX 40
#define MAXQ 100
#define ARRSIZ 256
#define LINE 80
#define MAXLINE 512
#define MAXMOTIFS 6
#define MAXSEQ  120000
#define FILENAMELEN 80

#pragma hdrstop

int dlength, qlength[MAXMOTIFS], nom,matches[MAXMOTIFS];
bool option,reverse;
char matrix[MAXMOTIFS][ARRSIZ][MAXMATRIX];
char qseq[MAXMOTIFS][MAXQ], dseq[MAXSEQ],dseq_r[MAXSEQ];
char efname[FILENAMELEN], rfname[FILENAMELEN], dfname[FILENAMELEN];
char dseqname[71];


//----------FUNCTIONS & PROCEDURES--------
void Search (char dfname[FILENAMELEN]);
char* Logo (int n);
void Make_Matrix (int n, char q[]);
int Get_Lib_Seq (char string[]);
void Compare();
void PcompareN(char seq[]); //search protein pattern against neucleotide
                                	//sequences.
char Translate (char s[]);  //translate one aa per 3 bp.
int Convert (char c);//convert (A,C,G,T) to (0,1,2,3), called from Translate
void Reverse (char seq[]);
void Output (int s, int nma, int n, int f, int i, char seq[]); //s:swutch; f:frame;
void Output2 (int tm);                                         //nma:num. of Matches.
void PrintHit (char hit[]); //for print out matched seq. called from Output.

//-------FILE HANDLER-------------------

//ifstream lfile ("D:\\GeneWork\\test.seq.txt");
//ifstream lfile ("E:\\Flysequence\\p1\\P1.txt");
//ifstream lfile ("E:\\Flysequence\\Est\\EST0331.txt");
//ifstream lfile ("E:\\Flysequence\\STS\\STS.txt");
ofstream rfile ("D:\\GeneWork\\testout.txt");
ofstream rfile2 ("D:\\GeneWork\\testout2.txt");


int main(int argc, char **argv)
{	char temp[MAXQ],qfname[FILENAMELEN];

	ifstream qfile ;
    option=1;
//======TITLE======
   	cout<<"------BAGUA20-------"<<endl;
   	cout<<"-Copy Right Lei Zhou-"<<endl;


//------get input------
	nom=0;
    strcpy(qfname,"D:\\GeneWork\\mysearch_bagua\\query.txt");
    cout<<qfname<<endl;
    qfile.open(qfname);
   while (qfile.getline (temp, MAXQ))
   {	if (temp[0]=='@')
   		{	strcpy(dfname, &temp[1]);
            Search (dfname);
            break;
        }
        if (temp[0]=='/')
        {	rfile<<"search profile:"<<temp<<endl;
            cout<<"search profile:"<<temp<<endl;
        }

        else
        {	strcpy( qseq[nom],&temp[0]);
        	rfile<<qseq[nom]<<endl;
            Make_Matrix (nom, qseq[nom]);
        	nom=nom+1;
        }


   		/****consol input ************
        for (nom=0; nom<MAXMOTIFS;nom++)
   		cout<<"motif "<<(nom+1)<<" ?:";
   		gets (qseq[nom]);
        if (qseq[nom][0]=='0') break;
        Make_Matrix (nom,qseq [nom]);*/
    }

   rfile<<"total motif:"<<nom<<endl<<endl;
   for (int i=0; i<nom; i++)
   		cout<<qseq[i]<<endl;

   /*cout<<"library file:";
   gets(lname); */

   /*cout<<"Search Option(1=P->N):";
   cin>>option;****for p->p ot n->n option****/

   cout<<"search finished"<<endl;
   getchar();
   qfile.close();

  	return 0;
}

void Search (char dfname[FILENAMELEN])
{   char string [MAXLINE];
	int n=0;   /*for desq*/
    ifstream lfile;
	lfile.open(dfname);
    if (lfile.fail())
    	cout<<"failed to open lfile"<<endl;
	lfile.getline (string,80);
    if (string[0]=='>') strncpy (dseqname,string,LINE);
    cout<<string<<endl;
    while (!lfile.eof())
	{
    	 //----------get seq from database file--------

    	while (lfile.getline (string, LINE))
    	{	if (string[0]=='>')
        	{   cout<<"dseq: "<<dseq<<endl;
            	rfile<<"dseq:"<<dseq<<endl;
            	PcompareN(dseq);
				Reverse(dseq);
        		PcompareN (dseq_r);
                reverse=0;
            	strcpy(dseqname, string);

             	break;
            }
            else
            {   int len=strlen(string);
            	for (int i=0;i<len; i++)
                {	if (n>=MAXSEQ) break;
                	dseq[n++]=string[i];
                }
            }
        }
        //-----------------compare---------------------

    }
    lfile.close();
}


void Make_Matrix (int n, char q[])
{	int i, j, bracket, k,count,q_length;
	count = bracket = 0;
    q_length=strlen (q);
    for(i = 0;i<MAXMATRIX;i++)
         for(j=0;j<ARRSIZ;j++)
              matrix[n][j][i] = 0;        /* set whole setmatrix to zero */
    cout<<"matrixzeroed"<<endl;

    for(i=0; i<q_length; i++)
    {	if (qseq[n][i] == '[')
    	{	bracket = 1;           /* increment the counter */
        	continue;
        }
        if (qseq[n][i] ==']')
        {   bracket = 0;
        	count++;
            continue;
        }
        if((qseq[n][i] == 'X') || (qseq[n][i] =='-'))
        {
        	for(k=0;k<ARRSIZ;k++)
            {	matrix[n][k][count] = 1;
            }
        }
         j = (int)qseq[n][i];//get the int value of qseq[i]
         if (bracket == 0)
         {	matrix[n][j][count] = 1;
         	count++;
         }
         if (bracket  == 1)
         {	 matrix[n][j][count] = 1;
         	continue;
         }
    }                                  /* invert the column of the*/
         qlength[n] = count;                   /* setmatrix */

     /* adjust qlength to account for brackets */

     if(bracket == 1)
     {
         cout<<"ERROR: bracket not closed in query motif"<<endl;
        //needs exit(0) eqivalent here
     }

    cout<<"qlenth ["<<(n+1)<<"]"<<qlength[n]<<endl;
    getchar();

}

/*void Compare()
{
    int i, j, score, k, len;
    int matches[MAXMOTIFS];
    int matches[MAXMOTIFS];
    matches[n] = 0;
    dlength = dlength + 1;

    len = dlength-qlength[n];
    for(i=0;i<len;i++)
    {   for(j=0,score=0;j<qlength[n];j++)
        {  score = score + matrix[n][(int)dseq[i+j]][j];
            if (score <= j) break;
        }

	   	if (score == qlength[n])
    	{	matches[n] = matches[n] + 1;
         	if (matches[n] == 1)
         	{	rfile<<dseqname;
         		cout<<dseqname;
         		rfile<<"matches at : "<<(i+1)<<endl;
            	for (k=i;k<i+qlength[n];k++)
                   rfile<<dseq[k];
         		cout<<"matches at:"<<(i+1)<<endl;
         	}

         	if (matches[n] > 1)
         	{	cout<<endl;
         		rfile<<endl;
            	rfile<<"match at"<<(i+1);
            	for (k=i;k<i+qlength[n];k++)
                   rfile<<dseq[k];
            	cout<<"matech at:"<<(i+1);
         	}
    	}
    	if(matches[n] > 0)
    	{   rfile<<endl;
         	rfile<<"total matches for motif ["<<(n+1)<<"]:"
            <<matches<<endl;
            totalmatch=totalmatch+1
	   	}
	}
} */

void PcompareN(char seq[])
{   int j,n,i,f,len, score, totalmatches, matches[MAXMOTIFS];
    int s=0;
    totalmatches=0;

    for (n=0; n<nom; n++)
    {	matches[n] = 0;
    	len = strlen(seq)-qlength[n];
        for ( f=0; f<3; f++)
		{	for( i=f;i<len;i=i+3)
    		{	for(j=0,score=0;j<qlength[n];j++)
        		{	score = score + matrix[n][(int)Translate(&seq[i+3*j])][j];
            		if (score <= j) break;
        		}
    			if (score == qlength[n])
            	{	matches[n]=matches[n]+1;
    				Output (s,matches[n],n,f,i,seq);
                    s=+1;
				}
			}
        }
     	if (matches[n]>0) totalmatches=totalmatches+1;
    }
    if (totalmatches==(nom+1)) Output2(totalmatches);
}

char Translate(char s[])  //take 3bp return 1aa
{   int c1,c2,c3;
	char P, code[3];

        //translation table, A(0),C(1), G(2), T(3)-----------
    char table [4][4][4]=
    {{{'K','N','K','N'},{'T','T','T','T'},{'R','S','R','S'},{'I','I','M','I'}},
     {{'Q','H','Q','H'},{'P','P','P','P'},{'R','R','R','R'},{'L','L','L','L'}},
     {{'E','D','E','D'},{'A','A','A','A'},{'G','G','G','G'},{'V','V','V','V'}},
     {{'*','Y','*','Y'},{'S','S','S','S'},{'*','C','W','C'},{'L','F','L','F'}}};

    strncpy (code, s, 3);
    c1=Convert(code[0]);
    c2=Convert(code[1]);
    c3=Convert(code[2]);
    if (c1>=4||c2>=4||c3>=4)
    	P='A';  //can be Optimized further here by considering....
    else
     	P=table[c1][c2][c3];
    //P=table[Convert(code[0])][Convert(code[1])][Convert(code[2])];

    return (P);
}

int Convert (char c) //convert ACGT to Int for Translate
{ 	char s=c;
	if (s=='A'||s=='a') return (0);
    if (s=='C'||s=='c') return (1);
    if (s=='G'||s=='g') return (2);
    if (s=='T'||s=='t'||s=='U'||s=='u') return (3);
    if (s=='N'||s=='n') return (4);
    else return (5);
}

void Reverse (char seq[]) //Reverse dseq

{	int i,j;

	j=0;

    for (i=(strlen(dseq)-1);i>0;i--)
    {   if (seq[i]=='A'||seq[i]=='a')
    		dseq_r[j++]='T';
        if (seq[i]=='C'||seq[i]=='c')
        	dseq_r[j++]='G';
        if (seq[i]=='G'||seq[i]=='g')
             dseq_r[j++]='C';
        if (seq[i]=='T'||seq[i]=='t')
        	dseq_r[j++]='A';
    }
    reverse=1;
}
//--------for every match----------------
// s:number of matches--as output switch;
//n:motif No.; f:fram; i:position
void Output (int s, int nma, int n, int f, int i, char seq[])
{   int k,j,m;
	char hit[600], temp[66];
    j=m=0;
    if (s==0)
    {	cout<<"dseqname:"<<dseqname<<endl;
    	rfile<<"dseqname:"<<dseqname<<endl;
    	if (reverse==1)
    	{	cout<<"---Reverse strand"<<endl;
    		rfile<<"---Reverse strand"<<endl;
    	}
    }
    if (nma==1)
    {	cout<<setw(5)<<"Motif ["<<(n+1)<<"]"<<endl;
    	rfile<<setw(5)<<"Motif ["<<(n+1)<<"]"<<endl;
    }
	cout<<setw(5)<<"Matches at:"<<i;
    rfile<<setw(5)<< "Matches at:"<<i;
    cout<<setw(10)<<"Frame:"<<f<<endl;
    rfile<<setw(10)<<"Frame:"<<f<<endl;
    for (k=(i-15);k<(i+3*qlength[n]+15);k++)
    	hit[j++]=seq[k];
    hit [j++]='\0';                   //end of array;
    if (strlen(hit)<=66)
    	PrintHit (&hit[0]);
    if (strlen(hit)>66)
    {	for (j=0;j<=65;j++)
    		temp[m++]=hit[j];
    	PrintHit (&temp[0]);
        PrintHit (&hit[66]);
    }
    rfile<<endl;
}
//---------output2 for multiple hits------------
void Output2 (int tm)
{	cout<<dseqname<<endl;
	rfile2<<dseqname<<endl;
    if (reverse==1)
    {	cout<<"--reverse strand--"<<endl;
    	rfile2<<"--reverse strand--"<<endl;
    }
    cout<<"totalmatches :"<<tm<<endl;
    rfile2<<"totalmatches :"<<tm<<endl;
}

void PrintHit (char hit[])
{
    char temp[3];
    int len=strlen(hit);
    for (int i=0;i<len;i++)
    {	cout<<hit[i];
    	rfile<<hit[i];
    }
    cout<<endl;
    rfile<<endl;
    for (int i=0; i<len;i=i+3)
    { 	strncpy (temp, &hit[i],3);
    	cout<<Translate(temp)<<"  ";
        rfile<<Translate(temp)<<"  ";
    }
    cout<<endl;
    rfile<<endl;
}

//----------------------------------------------------------



