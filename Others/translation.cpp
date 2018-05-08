

/* TRANSLATION: 3 or 6 frame translate cDNA sequences				*/

//---------------------------------------------------------------------------
#include "translation.hpp"

int main(int argc, char **argv)
{   int num_seq=0;
	char string[MAXLINE];
    DSEQ * dseq;

    cout<<"translate reverse strand? (0=No, 1=Yes): ";
    cin>>REVERSE;
   
	infile.getline (string,MAXLINE);
    if (string[0]=='>') strncpy (dbname,string,MAXLINE);
    while (!infile.eof())
    {   dseq=Get_Lib_Seq ();
        Translation (&dseq->name[1], dseq->seq);
        if (REVERSE==1)
            Translation (&dseq->name[1], dseq->r_seq);
    	num_seq++;
        if (num_seq%1000==0)
        {   cout<<num_seq<<endl;
            cout<<dseq->name<<endl;
        }
        delete dseq;
    }

    infile.close();
    outfile.close();
    cout<<num_seq<<" translated"<<endl;
    getch();
	return 0;
}

DSEQ* Get_Lib_Seq()
{   int i,n;
    char str[MAXLINE];
    DSEQ* dseq;
    n = 0;
    dseq=new DSEQ;
    strcpy (dseq->name, dbname);

    while(infile.getline(str,MAXLINE))
    {	if (str[0] == '>')
        {     strcpy( dbname, str);
              break;
        }
	   	for(i=0;i<strlen(str);i++)
        {	if(n==MAXSEQ) break;
        	dseq->seq[n++] = str[i];
        }
    }
    dseq->seq[n]='\0';

    if(n==MAXSEQ)
       cout<<"WARNING: sequence"<<dbname<<"too long!"<<endl;
    dseq->len=n;
    if (REVERSE==1)
        Reverse (dseq);
    else
        dseq->reverse=0;
    return dseq;
}

void Reverse (DSEQ* dseq) //Reverse dseq
{	int i,j;
    j=0;
    for (i=(dseq->len-1);i>0;i--)
    {   if (dseq->seq[i]=='A'||dseq->seq[i]=='a')
    		dseq->r_seq[j++]='T';
        if (dseq->seq[i]=='C'||dseq->seq[i]=='c')
        	dseq->r_seq[j++]='G';
        if (dseq->seq[i]=='G'||dseq->seq[i]=='g')
             dseq->r_seq[j++]='C';
        if (dseq->seq[i]=='T'||dseq->seq[i]=='t')
        	dseq->r_seq[j++]='A';
        if (dseq->seq[i]=='N'||dseq->seq[i]=='n')
        	dseq->r_seq[j++]='N';
    }
    dseq->r_seq[j++]='\0';
    dseq->reverse=1;

}
void Translation (char name[], char seq[])
{   char ppseq[MAXSEQ/3];

    for (int f=0; f<3; f++)
    {   outfile<<">"<<"F_"<<f<<name<<endl;
        int j=0;
    	int len=strlen(seq);
        for( int i=f; i<len; i=i+3)
        	ppseq[j++]=Translate(&seq[i]);
        ppseq[j++]='\0';
        int m=strlen(ppseq)/50;           	// output 50 aa per line
        for (int n=0; n<=m; n++)
    	{	for (int i=n*50; i<50*(n+1); i++)
        	{	outfile<<ppseq[i];
            	if (ppseq[i]=='\0') break;
            }
        	outfile<<endl;
    	}
    }
}


char Translate(char s[])
{   int c1,c2,c3;
	char P, code[3];

//***standard translation table, A(0),C(1), G(2), T(3)*****

    char table [4][4][4]=
    {{{'K','N','K','N'},{'T','T','T','T'},{'R','S','R','S'},{'I','I','M','I'}},
     {{'Q','H','Q','H'},{'P','P','P','P'},{'R','R','R','R'},{'L','L','L','L'}},
     {{'E','D','E','D'},{'A','A','A','A'},{'G','G','G','G'},{'V','V','V','V'}},
     {{'*','Y','*','Y'},{'S','S','S','S'},{'*','C','W','C'},{'L','F','L','F'}}};

//*********** table2 for n at 3rd position********************
char table2 [4][4]={{'X','T','X','X'},{'X','P','R','L'},
    					{'X','A','G','V'},{'X','S','X','X'}};
    strncpy (code, s, 3);
    c1=Convert(code[0]);
    c2=Convert(code[1]);
    c3=Convert(code[2]);
    if (c1>=4 || c2>=4)
    P='X';  //can be Optimized further here by considering....
	else
    {   if (c3>=4)
    		P=table2[c1][c2];
    	else
    		P=table[c1][c2][c3];
		//P=table[Convert(code[0])][Convert(code[1])][Convert(code[2])];
     }
    return (P);
}

int Convert (char c)
{ 	char s=c;
	if (s=='A'||s=='a') return (0);
    if (s=='C'||s=='c') return (1);
    if (s=='G'||s=='g') return (2);
    if (s=='T'||s=='t'||s=='U'||s=='u') return (3);
    if (s=='N'||s=='n') return (4);
    else return (5);
}


