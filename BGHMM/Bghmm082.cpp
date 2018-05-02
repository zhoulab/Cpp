
/**********************************************************************/

/*    BaGua-HMM
      copy right Lei Zhou
*/

/*********************************************************************/

/* this programm was developed on the base of BaGuaG to search motifs
using the HMM model scoring matrix
*/

/*  modified 99/2-3 to search mutiple motifs against the database
*/

/*  modified 04/12 to take file name at run time
*/

/* 05/04 - put former golbal variables Check_up, Check_down into motifs{}
*/

/* 07/08 - minor changes, Namespace = 20. Check redundancy from p instead of p-MARGIN
*/

/* 08/10 - minor changes, MaSeq=21,000
*/


#include "BGHMM082.hpp"

using namespace std;

int main(int argc, char **argv)
{   int Num_searched=0;
	char string[MAXLINE], l_file[MAXLINE], out_file[MAXLINE], out_seq[MAXLINE];
    char* search_name;
    PSEQ* pseq;
    MOTIFLIST *mlp;
    SCORES* scores;

//======TITLE======

   	cout<<"------BGHMM 0.8-------"<<endl;
   	cout<<"-Copy Right Lei Zhou-"<<endl;

    // get matrix


     //------get input------

    cout<<"search name: ";
    cin>>search_name;

    cout<<"sequence library to search: ";
    cin>>l_file;
    //l_file = "S:\\fly\\Dmel.allp";
    lfile.open(l_file);  // seq file handler
    if (!lfile) cout<<"can not open sequence file";

    cout<<"file for search result: ";
    cin>>out_file;

    rfile.open(out_file);  // result file handler
    if (!rfile) cout<<"can not open result file";

    cout<<"Generate Sequence File?(1=yes)";
    cin>>generate_seq_file;
    if (generate_seq_file==1)
    {   cout<<"file for holding sequences: ";
        cin>>out_seq;
        seq_file.open(out_seq);  // out_seq file handler
    }

     mlp=Input_Matrix ();
    cout<<"total "<<mlp->num_motifs<<" matrix built"<<endl;

//---Initialize the list-----
    htp=new HITTABLE;
    htp->num_motifs=mlp->num_motifs;
    for(int i=0; i<htp->num_motifs; i++)
    {   htp->lp[i]=new HITLIST;
        htp->lp[i]->n_hits=0;
    }
//-----get seq from database file--------
	lfile.getline (string,MAXLINE);
    if (string[0]=='>') strncpy (dbname,string,MAXLINE);

    while (! lfile.eof())
    {   pseq=Get_Lib_Seq ();

        if (pseq->len>=MAXMATRIX)

    	{   Num_searched++;
            if ((Num_searched%1000)==0)
            {   cout<<Num_searched<<" ";
            //<<setw(5)<<pseq->name<<" length: "<<pseq->len<<endl;
            }
            scores=Score_it(pseq,mlp);

            Sort_Scores (scores, pseq, mlp);
        }

    }

   Output(search_name, htp, mlp);
   if (generate_seq_file==1)
        Output_seq (htp, mlp);

   cout<<Num_searched<<" entries searched"<<endl;
   cout<<"DONE"<<endl;
   rfile<<Num_searched<<" entries searched"<<endl;
   lfile.close();

   rfile.close();
   seq_file.close();
   getchar();
   return 0;
}

/* read from mutiple fatsa file and generate pseq
                                                                */

PSEQ* Get_Lib_Seq()
{   int i,n;
    char str[MAXLINE];
    PSEQ* pseq;
    n = 0;
    pseq=new PSEQ;                      //delete after Score in Main
    strcpy (pseq->name, dbname);

    while(lfile.getline(str,MAXLINE))
    {	if (str[0] == '>')
        {     strcpy( dbname, str);
              break;
        }
	   	for(i=0;i<strlen(str);i++)
        {	if(n==MAXSEQ) break;
        	pseq->seq[n++] = str[i];
        }
    }
    pseq->seq[n]='\0';

    if(n==MAXSEQ)
       cout<<"WARNING: sequence"<<dbname<<"too long!"<<endl;
    pseq->len=n;
    return pseq;

}

/* read from matrix.txt file buid matrix[][]
                                                                    */
MOTIFLIST* Input_Matrix ()
{   MOTIFLIST* mlp;
    int n=0;
    char alphebit[ALPHEBIT]="ACDEFGHIKLMNPQRSTVWY";
    char str[MAXLINE], mfile [MAXLINE];

    mlp=new MOTIFLIST;
    mlp->num_motifs=0;

    cout<<"matrix file name: ";
    cin>>mfile;
    //mfile="E:\\Bagua\\matrix\\corp\\corp.m0502.txt";
    ifstream m_file(mfile);
    if (!m_file) cout<<"can not open matrix file";

    while (!m_file.eof())
    {   m_file.getline(str, MAXLINE);
        if (str[0]=='#') break;
        if (str[0]=='M')
        {   mlp->motifs[n]=new MOTIF;
            strcpy(mlp->motifs[n]->name, &str[6]);
            m_file>>mlp->motifs[n]->len;
            m_file>>mlp->motifs[n]->cut_off;
            m_file>>mlp->motifs[n]->Check_up;
            m_file>>mlp->motifs[n]->Check_down;

            for(int i = 0;i<MAXMATRIX;i++)
                for(int j=0;j<ARRSIZE;j++)
                    mlp->motifs[n]->matrix[i][j] = 0;        	/* set whole matrix to zero */

            for(int i=0; i<mlp->motifs[n]->len; i++)
            {   for (int j=0; j<ALPHEBIT; j++)
                m_file>>mlp->motifs[n]->matrix[i][int (alphebit[j])];
            }
            n+=1;
            cout<<n<<" matrix build "<<&str[6]<<endl;
        }
    }
    mlp->num_motifs=n;
    m_file.close();
    getchar();
    return (mlp);
}

/* score one sequence against the models, generate score[] in seq_score
                                                                    */
SCORES* Score_it(PSEQ* pseq, MOTIFLIST* mlp)
{   int i, j, len;
    char temp[MAXMATRIX];

    SCORES* scores=new SCORES;

    /*Zero the 2d array*/
    for (int i=0;i<MAXMOTIFS; i++)
        for (int j=0; j<MAXSEQ; j++)
           scores->seq_scores[i][j]=0;

    for (int n=0; n<mlp->num_motifs; n++)
    {   len = pseq->len - mlp->motifs[n]->len;

        for (i=0; i<len;i++)
        {   if (i > mlp->motifs[n]->Check_up && (len-i)> mlp->motifs[n]->Check_down)        //checking for stop codon
            strncpy(temp, &pseq->seq[i-mlp->motifs[n]->Check_up], (mlp->motifs[n]->len+mlp->motifs[n]->Check_up
                + mlp->motifs[n]->Check_down) );
            else strncpy (temp, &pseq->seq[i], mlp->motifs[n]->len);
            if (Check_Stopcodon(temp)==0)      //score=0 if * found
                scores->seq_scores[n][i]=0;
            else
            {   float s=0;
                for(j=0;j<mlp->motifs[n]->len;j++)
                    s = s+mlp->motifs[n]->matrix[j][(int)temp[j]];
                scores->seq_scores[n][i]=s;
            }
        }
     }
     return (scores);
}
int Check_Stopcodon (char seq[])
{   for (int i=0; i<strlen(seq); i++)
        if (seq[i]=='*'||seq[i]=='X') return (FALSE);

    return (TRUE);
}

void Sort_Scores (SCORES* scores,PSEQ*pseq, MOTIFLIST* mlp)
{   int removepseq=0;
    for (int n=0; n<mlp->num_motifs; n++)
    {   float hscore=0;
        int p=0;
        HIT* hit=0;

        /* find highest score for each motif-seq pair*/
        for (int i=0; i<(pseq->len - mlp->motifs[n]->len); i++)
            if (scores->seq_scores[n][i]>hscore)
            {   hscore=scores->seq_scores[n][i];
                p=i;
            }

        if (hscore>mlp->motifs[n]->cut_off)
        {   hit= new HIT;         //delete in Check_hit if not chosen
            hit->motif_n=n;
            hit->score=hscore;
            hit->pseq=pseq;
            hit->hit_p=p;
            hit->more_entries=0;
            hit->ext_names[0]=NULL;

            int n=Check_Hit (hit);
            removepseq+=n;
        }
        else
        {    delete hit;
        }
      }
     /* clearing */
     delete scores;
     if (removepseq==0) //pseq is not used
        delete pseq;
}

/* check to see if the score is high enough for output return TRUE
if pseq is needed for output
                                                                    */
bool Check_Hit (HIT* hit)
{   int n=hit->motif_n;
    if (htp->lp[n]->n_hits==0)
    {   Insert_Hit (hit, 0);
        return (TRUE);
    }
    else
    {   for (int i=0; i<htp->lp[n]->n_hits; i++)
        {   /**for identical score, check redundancy**/
            //cout<<"hitscore="<<hit->score<<" list="<<htp->lp[n]->hitlist[i]->score<<endl;
            if (hit->score == htp->lp[n]->hitlist[i]->score )
            {   if (Check_Redundancy (hit, htp->lp[n]->hitlist[i])==TRUE)
                   return (FALSE);
                else
                {   Insert_Hit (hit, i);
                    return (TRUE);
                }
            }
            /** if larger then an earlier one, insert*/
            if (hit->score > htp->lp[n]->hitlist[i]->score)
            {   Insert_Hit (hit, i);
                return (TRUE);
            }
        }
        /** append if n_hits<MAXHITS*/
        if (htp->lp[n]->n_hits<MAXHITS)
        {   htp->lp[n]->hitlist[htp->lp[n]->n_hits]=hit;
            htp->lp[n]->n_hits++;
            return (TRUE);
        }

        else  /** delete hit if failed to make into the list*/
        {   
            delete hit;
            return (FALSE);
        }
    }
}

/* check reduncy, if fund, just add new name to old hit
*/
bool Check_Redundancy (HIT* new_hit, HIT* old_hit)
{   cout<<"C-R called-- ";
    int n= new_hit->motif_n;
    int len= MARGIN;
    int p1=new_hit->hit_p;
    int p2=old_hit->hit_p;
    
    if (memcmp(&new_hit->pseq->seq[p1], &old_hit->pseq->seq[p2],
            len)==0)
    {   cout<<"compared ";
        strncat(old_hit->ext_names, new_hit->pseq->name, NAMESPACE);
        old_hit->more_entries ++;
        cout<<"added ";
        delete new_hit;
        cout<<"deleted "<<endl;
        return (TRUE);
    }
    else
        return (FALSE);
}

/* add hit to out put list
                                                                    */
void Insert_Hit (HIT* hit, int insertp)
{   int n=hit->motif_n;
    if (htp->lp[n]->n_hits==MAXHITS)    //full list, remove the last one
    {   delete htp->lp[n]->hitlist[MAXHITS-1];
        for (int j=MAXHITS-1; j>insertp; j--)
            htp->lp[n]->hitlist[j] = htp->lp[n]->hitlist[j-1];

        htp->lp[n]->hitlist[insertp]=hit;
        return;
    }

    /*** move everything after insertp one postion towards end  ****/

    for (int j=htp->lp[n]->n_hits; j>insertp; j--)
            htp->lp[n]->hitlist[j] = htp->lp[n]->hitlist[j-1];

    htp->lp[n]->hitlist[insertp]=hit;

    if (htp->lp[n]->n_hits<MAXHITS)
        htp->lp[n]->n_hits++;
    cout<<"motif "<<(n+1)<<" has "<<htp->lp[n]->n_hits<<" hits "
        <<htp->lp[n]->hitlist[insertp]->score <<endl;
}

//--------for every match----------------
void Output (char* search_name,HITTABLE* htp,MOTIFLIST* mlp)
{   time_t* dt=new time_t;
    time (dt);
    rfile<<setw(45)<<"------BGHMM 0.8-------"<<endl;
   	rfile<<setw(45)<<"-Copy Right Lei Zhou-"<<endl<<endl;
    rfile<<setw(15)<<search_name<<" done at-- "<<ctime(dt)<<endl;
    rfile<<endl;


    rfile<<setw(15)<<"MOTIF"<<setw(10)<<"CUT_OFF"<<setw(10)
    <<"Check_up"<<setw(10)<<"Check_down"<<setw(10)
                                            <<"No. HITS"<<endl;
    rfile<<setw(15)<<"---------"<<setw(10)<<"--------"
    <<setw(10)<<"--------"<<setw(10)<<"--------"
                                    <<setw(10)<<"--------"<<endl;
  /* generate summary of search result */
    for (int n=0; n<htp->num_motifs; n++)
    {   rfile<<setw(15)<<mlp->motifs[n]->name;
        rfile<<setw(10)<<mlp->motifs[n]->cut_off;
        rfile<<setw(10)<<mlp->motifs[n]->Check_up;
        rfile<<setw(10)<<mlp->motifs[n]->Check_down;
        rfile<<setw(10)<<htp->lp[n]->n_hits;
        rfile<<endl;
    }
   /* generate reports for individual motif*/
    for (int n=0; n<htp->num_motifs; n++)
    {   rfile<<endl<<setw(15)<<"##   "<<mlp->motifs[n]->name<<"   ##"
            <<endl<<setw(20)<< htp->lp[n]->n_hits<<" HITS"
                <<endl<<endl;

        for(int i=0; i<htp->lp[n]->n_hits; i++)
        {   int p= htp->lp[n]->hitlist[i]->hit_p;
            int len=mlp->motifs[n]->len;
            rfile<<"*****************************"<<endl;
            rfile<<"dbname:"<<htp->lp[n]->hitlist[i]->pseq->name<<setw(5)
                    <<"len "<<htp->lp[n]->hitlist[i]->pseq->len<<endl;
            if ( htp->lp[n]->hitlist[i]-> more_entries > 0)
                rfile<<setw(10)<<"RD:"<<htp->lp[n]->hitlist[i]->ext_names
                        <<endl ;
   	        cout<<setw(5)<<"Matches at:"<<p<<endl;
            cout<<"score "<<htp->lp[n]->hitlist[i]->score<<endl;
            rfile<<setw(5)<< "Matches at:"<<p<<endl;
            rfile<<"score "<<htp->lp[n]->hitlist[i]->score<<endl;
            //--output sequence
            for (int j=p-MARGIN; j<p; j++)
                rfile<< htp->lp[n]->hitlist[i]->pseq->seq[j];
            rfile<<"  ";
            for (int j=p; j<(p+ len); j++)
                rfile<<htp->lp[n]->hitlist[i]->pseq->seq[j];
            rfile<<"  ";
            for (int j=(p+ len);j<(p+len+MARGIN); j++)
                rfile<<htp->lp[n]->hitlist[i]->pseq->seq[j];

            rfile<<endl;
        }
    }
}

void Output_seq (HITTABLE* htp, MOTIFLIST* mlp)
{   for (int n=0; n<htp->num_motifs; n++)
    {   for(int i=0; i<htp->lp[n]->n_hits; i++)
        {   /* output seq name. (reverse) */
            seq_file<<endl;
	        seq_file<<htp->lp[n]->hitlist[i]->pseq->name<<endl;

        /*********  ouput sequence   *************/

            char temp[MAXSEQ];
            strcpy (temp,htp->lp[n]->hitlist[i]->pseq->seq);
            int m=strlen(temp)/50;           // output 50 aa per line
            for (int n=0; n<=m; n++)
            {    for (int i=n*50; i<50*(n+1); i++)
    	        {   if (temp[i]=='\0') break;
                seq_file<<temp[i];
                }
    	        seq_file<<endl;
            }
         }
    }
}

//----------------------------------------------------------









