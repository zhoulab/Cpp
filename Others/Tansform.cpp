#include "Transform.hpp"

int main(int argc, char **argv)
{   int num_tf=0;
    char string[MAXLINE];
    GSEQ* Gseq;
    num_gnt=0;
//-----get seq from database file--------
	lfile.getline (string,MAXLINE);
    if (string[0]=='>') strncpy (dbname,string,MAXLINE);

    while (! lfile.eof())
    {   Gseq=Get_Lib_Seq ();
        Transform (Gseq);
        delete Gseq;
        num_tf+=1;
        if (num_tf%10==0)
        cout<<num_tf<<setw(2)<<Gseq->name<<endl;
     }
   cout<<num_tf<<" entries transformed"<<endl;
   cout<<num_gnt<<" entries generated"<<endl;
   cout<<"Transformation finished"<<endl;
   lfile.close();

   seq_file.close();
   getchar();
   return 0;
}

/* read from mutiple fatsa file and generate Gseq
                                                                */
GSEQ* Get_Lib_Seq()
{   int i,n;
    char str[MAXLINE];
    GSEQ* Gseq;
    n = 0;
    Gseq=new GSEQ;
    strcpy (Gseq->name, dbname);

    while(lfile.getline(str,MAXLINE))
    {	if (str[0] == '>')
        {     strcpy( dbname, str);
              break;
        }
	   	for(i=0;i<strlen(str);i++)
        {	if(n==MAXSEQ) break;
        	Gseq->seq[n++] = str[i];
        }
    }
    Gseq->seq[n++]='\0';

    if(n==MAXSEQ)
       cout<<"WARNING: sequence"<<dbname<<"too long!";
    Gseq->len=n;
    Reverse(Gseq);

    return Gseq;

}
void Reverse (GSEQ* Gseq) //Reverse dseq

{	int i,j;

	j=0;

    for (i=(Gseq->len-1);i>0;i--)
    {   if (Gseq->seq[i]=='A'||Gseq->seq[i]=='a')
    		Gseq->rseq[j++]='T';
        if (Gseq->seq[i]=='C'||Gseq->seq[i]=='c')
        	Gseq->rseq[j++]='G';
        if (Gseq->seq[i]=='G'||Gseq->seq[i]=='g')
             Gseq->rseq[j++]='C';
        if (Gseq->seq[i]=='T'||Gseq->seq[i]=='t')
        	Gseq->rseq[j++]='A';
        if (Gseq->seq[i]=='N'||Gseq->seq[i]=='n')
        	Gseq->rseq[j++]='N';
    }
    Gseq->rseq[j++]='\0';
    Gseq->reverstrand=1;
}
void Transform (GSEQ* gseq)
{   char temp[OUTSEQL];
    int n=gseq->len/OUTSEQL;
    for (int i=0; i<=(n+1); i++)
    {   int start=i*(OUTSEQL-OVERLAP);
        strncpy (temp, &gseq->seq[start], OUTSEQL);
        Output_seq (start, gseq->name, 0, temp);
    }
     for (int i=0; i<=(n+1); i++)
    {   int start=i*(OUTSEQL-OVERLAP);
        strncpy (temp, &gseq->rseq[start], OUTSEQL);
        Output_seq (start, gseq->name, 1, temp);
    }
}

void Output_seq (int start, char name[], bool reverse, char seq[])
{  /* output seq name,(reverse), start-stop, annotation
                                                                 */
    int len=strlen(seq);
    for(int i=0; i<10; i++)
    {   seq_file<<name[i];
        if (name[i]==' ') break;
    }
	if (reverse==1)
    	seq_file<<" -REV. STD"<<setw(1);
    seq_file<<" ("<<start<<"-"<<(start+len)<<")";

    for(int i=10; i<MAXLINE; i++)
    {   if( name[i]=='\0') break;
        seq_file<<name[i];
    }
    seq_file<<endl;

    //-------out put sequence--------------
    int m=len/150;           // output 150 bp per line
    for (int n=0; n<=m; n++)
    {	for (int i=n*150; i<150*(n+1); i++)
        {    if (seq[i]=='\0') break;
    		seq_file<<seq[i];
        }
    	seq_file<<endl;
    }
    num_gnt+=1;
}
