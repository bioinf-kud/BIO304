#include <stdio.h>
#include <stdlib.h>
#include <math.h>
struct score_matrix{
    char a;
    char b;
    int score;
};
struct node{
    int father[2]; 
    int gapn;
    int score;
};
struct result{
    char seq1[1000];
    char seq2[1000];
    char seq3[1000];
};
struct treenode{
    int id;
    int son_id[2];
};
struct node_dist{
    int id1;
    int id2;
    int dist;
};
struct max_dist{
    int id1;
    int id2;
    int dist;
};
int length(char*a);
void copystr(char*a,char*b);
int alignDNA(char*seq1,char*seq2,struct result*result);
int alignProtein(char*seq1,char*seq2,struct score_matrix*protmatrix[400],struct result*result);
int main(int argc, char **argv){
    char seq[5][100]={"PEEKSAVTALWGKVNVDEYGG","GEEKAAVLALWDKVNEEEYGG","PADKTNVKAAWGKVGAHAGEYGA","AADKTNVKAAWSKVGGHAGEYGA","AATNVKTAWSSKVGGHAPAA"};
    struct score_matrix*protmatrix[400];
    FILE*fp=fopen("/Users/sunkai/Desktop/bio304/BLOSUM62.txt","r");
    FILE*seq1f=fopen("/Users/sunkai/Desktop/bio304/seq1.fasta","r");
    FILE*seq2f=fopen("/Users/sunkai/Desktop/bio304/seq2.fasta","r");
    for(int i=0;i<400;i++){
        protmatrix[i]=(struct score_matrix*)malloc(sizeof(struct score_matrix));
        fscanf(fp,"%c %c %d\n",&(protmatrix[i]->a),&(protmatrix[i]->b),&(protmatrix[i]->score));
    }
    fclose(fp);
    struct result result;
    struct treenode*tree=(struct treenode*)malloc(sizeof(struct treenode)*9);
    for(int i=0;i<9;i++){
        tree[i].id=i;
        tree[i].son_id[0]=-1;
        tree[i].son_id[1]=-1;
    }
    struct node_dist*dist=(struct node_dist*)malloc(sizeof(struct node_dist)*81);
    for(int i=0;i<9;i++){
        for(int j=i+1;j<9;j++){
            dist[i*9+j].id1=i;
            dist[i*9+j].id2=j;
            dist[i*9+j].dist=-999999;
        }
    }
    for(int i=0;i<5;i++){
        for(int j=i+1;j<5;j++){
            if(i==j)
                continue;
            int score=alignProtein(seq[i],seq[j],protmatrix,&result);
            printf("%d\n",score);
            dist[i*9+j].dist=score;
            dist[j*9+i].dist=score;
        }
    }
    int coldata[5]={0,1,2,3,4};
    int newcoldata[5]={0,1,2,3,4};
    for(int i=5;i<9;i++){
        for(int j=0;j<5;j++){
            coldata[j]=newcoldata[j];
        }
        struct max_dist maxdist;
        maxdist.dist=-999999;
        int score[10-i][10-i];
        for(int j=0;j<10-i;j++)
            for(int k=j+1;k<10-i;k++)
                for(int l=0;l<81;l++)
                    if(coldata[j]==dist[l].id1&&coldata[k]==dist[l].id2){
                        score[j][k]=dist[l].dist;
                        score[k][j]=dist[l].dist;
                    }
        for(int j=0;j<10-i;j++){
            for(int k=j+1;k<10-i;k++){
                if(score[j][k]>maxdist.dist){
                    maxdist.id1=j;
                    maxdist.id2=k;
                    maxdist.dist=score[j][k];
                }
            }
        }
        tree[i].son_id[0]=coldata[maxdist.id1];
        tree[i].son_id[1]=coldata[maxdist.id2];
        for(int j=0;j<10-i;j++){
            if(j!=maxdist.id1&&j!=maxdist.id2){
                dist[(10-i)*9+j].id1=coldata[maxdist.id1];
                dist[(10-i)*9+j].id2=coldata[j];
                dist[(10-i)*9+j].dist=(score[maxdist.id1][j]+score[maxdist.id2][j])/2;
                dist[j*9+(10-i)].id1=coldata[maxdist.id1];
                dist[j*9+(10-i)].id2=coldata[j];
                dist[j*9+(10-i)].dist=(score[maxdist.id1][j]+score[maxdist.id2][j])/2;
            }
        }
        int cnt=0;
        for(int j=0;j<10-i;j++){
            if(j!=maxdist.id1&&j!=maxdist.id2){
                newcoldata[cnt]=coldata[j];
                cnt++;
            }
            newcoldata[cnt]=i;
        }
    }
    for(int i=0;i<9;i++){
        printf("%d %d %d\n",tree[i].id,tree[i].son_id[0],tree[i].son_id[1]);
    }
}
int length(char*a){
    int cnt=0;
    while(a[cnt]!='\0'){
        if(a[cnt]=='\n'){
            a[cnt]='\0';
            break;
        }
        cnt++;
    }
    return cnt;
}
void copystr(char*a,char*b){//b--->a
    int la=length(a),lb=length(b);
    if(la>lb){
        for(int i=0;i<la;i++){
            if(i<lb)
                *(a+i)=*(b+i);
            else
                *(a+i)='\0';
        }
    }
    else{
        for(int i=0;i<lb;i++){
            *(a+i)=*(b+i);
        *(a+lb)='\0';
    }
    }
}
int alignDNA(char*seq1,char*seq2,struct result*result){
    int l1=length(seq1);
    int l2=length(seq2);
    struct node score_matrix[l1+1][l2+1];
    for(int i=0;i<l1+1;i++){
        for(int j=0;j<l2+1;j++){
            score_matrix[i][j].father[0]=-1;
            score_matrix[i][j].father[1]=-1;
            score_matrix[i][j].gapn=0;
            score_matrix[i][j].score=0;
        }
    }
    for(int i=1;i<l1+1;i++){
        score_matrix[i][0].score=score_matrix[i-1][0].score-2;
        score_matrix[i][0].father[0]=i-1;
        score_matrix[i][0].father[1]=0;
        score_matrix[i][0].gapn=score_matrix[i-1][0].gapn+1;
    }
    for(int i=1;i<l2+1;i++){
        score_matrix[0][i].score=score_matrix[0][i-1].score-2;
        score_matrix[0][i].father[0]=0;
        score_matrix[0][i].father[1]=i-1;
        score_matrix[0][i].gapn=score_matrix[0][i-1].gapn+1;
    }
    for(int i=1;i<l1+1;i++){
        for(int j=1;j<l2+1;j++){
            int match=-1;
            if (seq1[i-1]==seq2[j-1])
                match=1;
            int score1=score_matrix[i-1][j-1].score+match;
            int score2=score_matrix[i-1][j].score-2;
            int score3=score_matrix[i][j-1].score-2;
            if(score1>=score2&&score1>=score3){
                score_matrix[i][j].score=score1;
                score_matrix[i][j].father[0]=i-1;
                score_matrix[i][j].father[1]=j-1;
                score_matrix[i][j].gapn=0;
            }
            else if(score2>=score1&&score2>=score3){
                score_matrix[i][j].score=score2;
                score_matrix[i][j].father[0]=i-1;
                score_matrix[i][j].father[1]=j;
                score_matrix[i][j].gapn=score_matrix[i-1][j].gapn+1;
            }
            else{
                score_matrix[i][j].score=score3;
                score_matrix[i][j].father[0]=i;
                score_matrix[i][j].father[1]=j-1;
                score_matrix[i][j].gapn=score_matrix[i][j-1].gapn+1;
            }
        }
    }
    int i=l1,j=l2;
    int flag=0;
    while(i!=0&&j!=0){
        if(score_matrix[i][j].father[0]==i-1&&score_matrix[i][j].father[1]==j-1){
            result->seq1[flag]=seq1[i-1];
            result->seq3[flag]=seq2[j-1];
            if (seq1[i-1]!=seq2[j-1])
                result->seq2[flag]=' ';
            else
                result->seq2[flag]='|';
            i--;
            j--;
        }
        else if(score_matrix[i][j].father[0]==i-1){
            result->seq1[flag]=seq1[i-1];
            result->seq2[flag]=' ';
            result->seq3[flag]='-';
            i--;
        }
        else{
            result->seq1[flag]='-';
            result->seq2[flag]=' ';
            result->seq3[flag]=seq2[j-1];
            j--;
        }
        flag++;
    }
    result->seq1[flag]='\0';
    result->seq2[flag]='\0';
    result->seq3[flag]='\0';
    return score_matrix[l1][l2].score;
}
int alignProtein(char*seq1,char*seq2,struct score_matrix*protmatrix[400],struct result*result){
    int gappenalty=-5;
    int l1=length(seq1);
    int l2=length(seq2);
    struct node score_matrix[l1+1][l2+1];
    for(int i=0;i<l1+1;i++){
        for(int j=0;j<l2+1;j++){
            score_matrix[i][j].father[0]=-1;
            score_matrix[i][j].father[1]=-1;
            score_matrix[i][j].gapn=0;
            score_matrix[i][j].score=0;
        }
    }
    for(int i=1;i<l1+1;i++){
        score_matrix[i][0].score=score_matrix[i-1][0].score+gappenalty;
        score_matrix[i][0].father[0]=i-1;
        score_matrix[i][0].father[1]=0;
        score_matrix[i][0].gapn=score_matrix[i-1][0].gapn+1;
    }
    for(int i=1;i<l2+1;i++){
        score_matrix[0][i].score=score_matrix[0][i-1].score+gappenalty;
        score_matrix[0][i].father[0]=0;
        score_matrix[0][i].father[1]=i-1;
        score_matrix[0][i].gapn=score_matrix[0][i-1].gapn+1;
    }
    for(int i=1;i<l1+1;i++){
        for(int j=1;j<l2+1;j++){
            int match=-1;
            for(int k=0;k<400;k++){
                if((seq1[i-1]==protmatrix[k]->a&&seq2[j-1]==protmatrix[k]->b)||(seq1[i-1]==protmatrix[k]->b&&seq2[j-1]==protmatrix[k]->a)){
                    match=protmatrix[k]->score;
                    break;
                }
            }
            int score1=score_matrix[i-1][j-1].score+match;
            int score2=score_matrix[i-1][j].score+gappenalty;
            int score3=score_matrix[i][j-1].score+gappenalty;
            if(score1>=score2&&score1>=score3){
                score_matrix[i][j].father[0]=i-1;
                score_matrix[i][j].father[1]=j-1;
                score_matrix[i][j].gapn=0;
                score_matrix[i][j].score=score1;
            }
            else if(score2>=score1&&score2>=score3){
                score_matrix[i][j].father[0]=i-1;
                score_matrix[i][j].father[1]=j;
                score_matrix[i][j].gapn=score_matrix[i-1][j].gapn+1;
                score_matrix[i][j].score=score2;
            }
            else{
                score_matrix[i][j].father[0]=i;
                score_matrix[i][j].father[1]=j-1;
                score_matrix[i][j].gapn=score_matrix[i][j-1].gapn+1;
                score_matrix[i][j].score=score3;
            }
        }
    }
    int i=l1,j=l2;
    int flag=0;
    while(i!=0&&j!=0){
        if(score_matrix[i][j].father[0]==i-1&&score_matrix[i][j].father[1]==j-1){
            result->seq1[flag]=seq1[i-1];
            result->seq3[flag]=seq2[j-1];
            if (seq1[i-1]!=seq2[j-1]){
                if(score_matrix[i][j].score-score_matrix[i-1][j-1].score>=0)
                    result->seq2[flag]='+';
                else
                    result->seq2[flag]=' ';
            }
            else
                result->seq2[flag]=seq1[i-1];
            i--;
            j--;
        }
        else if(score_matrix[i][j].father[0]==i-1){
            result->seq1[flag]=seq1[i-1];
            result->seq2[flag]=' ';
            result->seq3[flag]='-';
            i--;
        }
        else{
            result->seq1[flag]='-';
            result->seq2[flag]=' ';
            result->seq3[flag]=seq2[j-1];
            j--;
        }
        flag++;
    }
    result->seq1[flag]='\0';
    result->seq2[flag]='\0';
    result->seq3[flag]='\0';
    return score_matrix[l1][l2].score;
}