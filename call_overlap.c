#include <stdio.h>
#include <stdlib.h>
#include <string.h>
struct attribute{//used to record separated attributes
    char attribute_name[50];
    char attribute_value[50];
};
struct gtfinfo{
    char seqname[15];
    char source[20];
    char feature[12];
    int start;
    int end;
    char score;
    char strand;
    char frame;
    int attribute_number;
    struct attribute attr[30];
    };
struct gene{
    char gene_id[20];
    char gene_name[20];
    char gene_type[25];
    char chromosome[15];
    int gene_start;//position on chromosome
    int gene_end;//position on chromosome
    int gene_transcriptf_start;
    int gene_transcriptf_end;
    int selected_transcript;//use to selected one transcript and merge features
    };
struct transcript{
    int gene_flag;
    char transcript_id[20];
    char transcript_type[50];
    char transcript_strand;//+ or -
    int transcript_start;//position on chromosome
    int transcript_end;//position on chromosome
    int EXON_length;//length of all exons
    int CDS_length;//length of all CDS
    int transcript_exonf_start;
    int transcript_exonf_end;
    int transcript_CDSf_start;
    int transcript_CDSf_end;
    int transcript_UTRf_start;
    int transcript_UTRf_end;
    int transcript_otherf_start;
    int transcript_otherf_end;
    int is_canonical;//boolean 1(canonical) or 0(non-canonical)
    };
struct exon{
    int transcript_flag;
    char exon_id[20];
    int exon_number;
    int exon_start;
    int exon_end;
    int exon_length;
    };
struct CDS{
    int transcript_flag;
    int exon_flag;
    int CDS_start;
    int CDS_end;
    int CDS_length;
    int CDS_frame;
    };
struct UTR{
    int transcript_flag;
    int exon_flag;
    int UTR_start;
    int UTR_end;
    int UTR_length;
    };
struct other_feature{
    int transcript_flag;
    int exon_flag;
    char feature_type[20];//start_codon,stop_codon
    int feature_start;//position on chromosome
    int feature_end;//position on chromosome
    int feature_length;
    };
struct startpos{
    char transcript_id[20];
    char chromosome[5];
    int start_pos;
};
void fgetstre(char*a,FILE*f);//scan attribute from file
int fgetattr(char*attrn,char*attrv,FILE*f);//get one attribute from file
int scan_row(FILE*gtf,struct gtfinfo*gtfinfo);//scan a row from gtf file
void read_attribute(struct attribute *attr,FILE *f,int attrn);
int read_row(FILE*gtf,struct gtfinfo*gtfinfo);
int length(char*a);//字符串长度
void swapstr(char*a,char*b);//字符串交换
void copystr(char*a,char*b);//字符串复制
void write_gene(struct gene*gene_list,struct gtfinfo*gtf,int tempg,int tempt);
void write_transcript(struct transcript*transcript_list,struct gtfinfo*gtf,int tempg,int tempt,int tempe,int tempc,int tempu,int tempo);
void write_exon(struct exon*exon_list,struct gtfinfo*gtf,int tempt,int tempe);
void write_CDS(struct CDS*CDS_list,struct gtfinfo*gtf,int tempt,int tempe,int tempc);
void write_UTR(struct UTR*UTR_list,struct gtfinfo*gtf,int tempt,int tempe,int tempu);
void write_other(struct other_feature*other_list,struct gtfinfo*gtf,int tempt,int tempe,int tempo);
void write_close(struct gene*gene_list,struct transcript*transcript_list,int tempg,int tempt,int tempe,int tempc,int tempu,int tempo);
void read_start_pos(struct startpos*pos,struct gene*genelist,struct transcript*translist,char*chr,int trnasnum);
void sort_start_pos(struct startpos*pos,int num);
int main(){
    int genenum=0,transnum=0,exonnum=0,CDSnum=0,UTRnum=0,othernum=0;//calculate number of each feature
    FILE*gtf,*out1,*out2,*out3,*out4;//gtf file
    gtf=fopen("/Users/sunkai/Desktop/gencode.v46.chr_patch_hapl_scaff.annotation.gtf","r");//open gtf file for scan
    out1=fopen("/Users/sunkai/Desktop/Gene_transcript information.tsv","w+");//open output file
    //out2=fopen("/Users/sunkai/Downloads/longest_trans_list.v44.longest_PC.tsv","w+");
    struct gtfinfo*gtfrow;
    struct gene*genelist;
    struct transcript*translist;
    struct exon*exonlist;
    struct CDS*CDSlist;
    struct UTR*UTRlist;
    struct other_feature*otherlist;
    struct trans_structure*transstructurelist;
    struct startpos*chr1,*chr2,*chr3,*chr4,*chr5,*chr6,*chr7,*chr8,*chr9,*chr10,*chr11,*chr12,*chr13,*chr14,*chr15,*chr16,*chr17,*chr18,*chr19,*chr20,*chr21,*chr22,*chrX,*chrY;
    gtfrow=(struct gtfinfo*)malloc(sizeof(struct gtfinfo)*1);//record a row of gtf
    char chr[24][5]={"chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"};
    int chrnum[24]={0};
    //scan to count number of each feature
    printf("Start scanning the gtf file.\n");
    int sflag=scan_row(gtf,gtfrow);
    while(sflag!=EOF){
        if(strcmp(gtfrow->feature,"gene")==0)
            genenum++;
        else if(strcmp(gtfrow->feature,"transcript")==0)
            transnum++;
        else if(strcmp(gtfrow->feature,"exon")==0)
            exonnum++;
        else if(strcmp(gtfrow->feature,"CDS")==0)
            CDSnum++;    
        else if(strcmp(gtfrow->feature,"UTR")==0)
            UTRnum++;    
        else
            othernum++;
        sflag=scan_row(gtf,gtfrow);
        
    }
    fclose(gtf);
    //allocate memory for each list
    genelist=(struct gene*)malloc(sizeof(struct gene)*genenum);
    translist=(struct transcript*)malloc(sizeof(struct transcript)*transnum);
    exonlist=(struct exon*)malloc(sizeof(struct exon)*exonnum);
    CDSlist=(struct CDS*)malloc(sizeof(struct CDS)*CDSnum);
    UTRlist=(struct UTR*)malloc(sizeof(struct UTR)*UTRnum);
    otherlist=(struct other_feature*)malloc(sizeof(struct other_feature)*othernum);
    printf("Successfully scanned the gtf file!\n");
    printf("gene number:%d\ntranscript number:%d\nexon number:%d\nCDS number:%d\nUTR number:%d\nother number:%d\n",genenum,transnum,exonnum,CDSnum,UTRnum,othernum);

    //read gtf file and record information
    printf("Start reading the gtf file.\n");
    gtf=fopen("/Users/sunkai/Desktop/gencode.v46.chr_patch_hapl_scaff.annotation.gtf","r");//open gtf file
    int rflag=read_row(gtf,gtfrow);
    int gflag=-1,tflag=-1,eflag=-1,cflag=-1,uflag=-1,oflag=-1;
    while(1){
        if(strcmp(gtfrow->feature,"gene")==0){
            if(gflag>-1)
                write_close(genelist,translist,gflag,tflag,eflag,cflag,tflag,oflag);
            gflag++;
            write_gene(genelist,gtfrow,gflag,tflag);
            copystr((genelist+gflag)->chromosome,gtfrow->seqname);

        }
        else if(strcmp(gtfrow->feature,"transcript")==0){
            tflag++;
            write_transcript(translist,gtfrow,gflag,tflag,eflag,cflag,uflag,oflag);
        }
        else if(strcmp(gtfrow->feature,"exon")==0){
            eflag++;
            write_exon(exonlist,gtfrow,tflag,eflag);
        }
        else if(strcmp(gtfrow->feature,"CDS")==0){
            cflag++;
            write_CDS(CDSlist,gtfrow,tflag,eflag,cflag);
        }
        else if(strcmp(gtfrow->feature,"UTR")==0){
            uflag++;
            write_UTR(UTRlist,gtfrow,tflag,eflag,uflag);
        }
        else{
            oflag++;
            write_other(otherlist,gtfrow,tflag,eflag,oflag);
        }
        rflag=read_row(gtf,gtfrow);
        if(rflag==EOF){
            write_close(genelist,translist,gflag,tflag,eflag,cflag,tflag,oflag);
            break;
        }
    }
    
    for(int i=0;i<transnum;i++){
        translist[i].EXON_length=0;
        translist[i].CDS_length=0;
        for(int j=translist[i].transcript_exonf_start;j<=translist[i].transcript_exonf_end;j++){
            translist[i].EXON_length+=exonlist[j].exon_length;
        }
        for(int j=translist[i].transcript_CDSf_start;j<=translist[i].transcript_CDSf_end;j++){
            translist[i].CDS_length+=CDSlist[j].CDS_length;
        }
    }
    printf("successfully read the gtf file!\nassigned features:\n");
    printf("gene:%d\ntranscript:%d\nexon:%d\nCDS:%d\nUTR:%d\nother:%d\n",gflag+1,tflag+1,eflag+1,cflag+1,uflag+1,oflag+1);
    fclose(gtf);

    for(int i=0;i<24;i++){
        for(int j=0;j<transnum;j++){
            if(strcmp(genelist[translist[j].gene_flag].chromosome,chr[i])==0){
                chrnum[i]++;
            }
        }
    }
    printf("%d",chrnum[0]);
    chr1=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[0]);
    /*chr2=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[1]);
    chr3=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[2]);
    chr4=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[3]);
    chr5=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[4]);
    chr6=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[5]);
    chr7=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[6]);
    chr8=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[7]);
    chr9=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[8]);
    chr10=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[9]);
    chr11=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[10]);
    chr12=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[11]);
    chr13=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[12]);
    chr14=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[13]);
    chr15=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[14]);
    chr16=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[15]);
    chr17=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[16]);
    chr18=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[17]);
    chr19=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[18]);
    chr20=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[19]);
    chr21=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[20]);
    chr22=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[21]);
    chrX=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[22]);
    chrY=(struct startpos*)malloc(sizeof(struct startpos)*chrnum[23]);*/
    read_start_pos(chr1,genelist,translist,chr[0],transnum);
    /*read_start_pos(chr2,genelist,translist,chr[1],transnum);
    read_start_pos(chr3,genelist,translist,chr[2],transnum);
    read_start_pos(chr4,genelist,translist,chr[3],transnum);
    read_start_pos(chr5,genelist,translist,chr[4],transnum);
    read_start_pos(chr6,genelist,translist,chr[5],transnum);
    read_start_pos(chr7,genelist,translist,chr[6],transnum);
    read_start_pos(chr8,genelist,translist,chr[7],transnum);
    read_start_pos(chr9,genelist,translist,chr[8],transnum);
    read_start_pos(chr10,genelist,translist,chr[9],transnum);
    read_start_pos(chr11,genelist,translist,chr[10],transnum);
    read_start_pos(chr12,genelist,translist,chr[11],transnum);
    read_start_pos(chr13,genelist,translist,chr[12],transnum);
    read_start_pos(chr14,genelist,translist,chr[13],transnum);
    read_start_pos(chr15,genelist,translist,chr[14],transnum);
    read_start_pos(chr16,genelist,translist,chr[15],transnum);
    read_start_pos(chr17,genelist,translist,chr[16],transnum);
    read_start_pos(chr18,genelist,translist,chr[17],transnum);
    read_start_pos(chr19,genelist,translist,chr[18],transnum);
    read_start_pos(chr20,genelist,translist,chr[19],transnum);
    read_start_pos(chr21,genelist,translist,chr[20],transnum);
    read_start_pos(chr22,genelist,translist,chr[21],transnum);
    read_start_pos(chrX,genelist,translist,chr[22],transnum);
    read_start_pos(chrY,genelist,translist,chr[23],transnum);*/
    sort_start_pos(chr1,chrnum[0]);
    /*sort_start_pos(chr2,chrnum[1]);
    sort_start_pos(chr3,chrnum[2]);
    sort_start_pos(chr4,chrnum[3]);
    sort_start_pos(chr5,chrnum[4]);
    sort_start_pos(chr6,chrnum[5]);
    sort_start_pos(chr7,chrnum[6]);
    sort_start_pos(chr8,chrnum[7]);
    sort_start_pos(chr9,chrnum[8]);
    sort_start_pos(chr10,chrnum[9]);
    sort_start_pos(chr11,chrnum[10]);
    sort_start_pos(chr12,chrnum[11]);
    sort_start_pos(chr13,chrnum[12]);
    sort_start_pos(chr14,chrnum[13]);
    sort_start_pos(chr15,chrnum[14]);
    sort_start_pos(chr16,chrnum[15]);
    sort_start_pos(chr17,chrnum[16]);
    sort_start_pos(chr18,chrnum[17]);
    sort_start_pos(chr19,chrnum[18]);
    sort_start_pos(chr20,chrnum[19]);
    sort_start_pos(chr21,chrnum[20]);
    sort_start_pos(chr22,chrnum[21]);
    sort_start_pos(chrX,chrnum[22]);
    sort_start_pos(chrY,chrnum[23]);*/
    out2=fopen("/Users/sunkai/Desktop/chr.tsv","w+");
    for(int i=0;i<chrnum[0];i++){
        for(int j=i+1;j<chrnum[1];j++){
            if(chr1[j].start_pos<chr1[i].start_pos+4000){
                fprintf(out2,"%s\t%d\t%s\t%d\n",chr1[i].transcript_id,chr1[i].start_pos,chr1[j].transcript_id,chr1[j].start_pos);   
            }
        }
    }


    



    
    fclose(out2);

    return 0;
}
void fgetstre(char*a,FILE*f){
    int cnt=0;
    while(1){
        char c=fgetc(f);
        if(c=='\n'){
            *(a+cnt)='\0';
            break;
        }
        *(a+cnt)=c;
        cnt++;
    }
}
int fgetattr(char*attrn,char*attrv,FILE*f){
    int cnt=0;
    while(1){
        char c=fgetc(f);
        if(c==' '){
            *(attrn+cnt)='\0';
            break;
        }
        *(attrn+cnt)=c;
        cnt++;
    }
    cnt=0;
    while(1){
        char c=fgetc(f);
        if(c==';'){
            char d=fgetc(f);
            if(d=='\n'){
                *(attrv+cnt)='\0';
                return 0;
            }
            *(attrv+cnt)='\0';
            return 1;
        }
        if(c=='"'){
            continue;
        }
        *(attrv+cnt)=c;
        cnt++;
    }
}
int scan_row(FILE*gtf,struct gtfinfo*gtfinfo){
    char *temp;
    temp=(char*)malloc(sizeof(char)*1000);
    int a=fscanf(gtf,"%s\t%s\t%s\t%d\t%d\t%c\t%c\t%c\t",gtfinfo->seqname,gtfinfo->source,gtfinfo->feature,&(gtfinfo->start),&(gtfinfo->end),&(gtfinfo->score),&(gtfinfo->strand),&(gtfinfo->frame));
    if (a==EOF)
        return EOF;
    fgetstre(temp,gtf);
    return a;
}
int read_row(FILE*gtf,struct gtfinfo*gtfinfo){
    int a=0;
    a=fscanf(gtf,"%s\t%s\t%s\t%d\t%d\t%c\t%c\t%c\t",gtfinfo->seqname,gtfinfo->source,gtfinfo->feature,&(gtfinfo->start),&(gtfinfo->end),&(gtfinfo->score),&(gtfinfo->strand),&(gtfinfo->frame));
    if (a==EOF)
        return EOF;
    int count=0;
    int flag=1;
    while(flag==1){
        flag=fgetattr((gtfinfo->attr[count].attribute_name),(gtfinfo->attr[count].attribute_value),gtf);
        count++;
    }
    gtfinfo->attribute_number=count;
    return a;
}
void swapstr(char*a,char*b){
    char c[99];
    copystr(c,a);
    copystr(a,b);
    copystr(b,c);
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
void write_gene(struct gene*gene_list,struct gtfinfo*gtf,int tempg,int tempt){
    copystr((gene_list+tempg)->chromosome,gtf->seqname);
    copystr((gene_list+tempg)->gene_id,gtf->attr[0].attribute_value);
    copystr((gene_list+tempg)->gene_name,gtf->attr[2].attribute_value);
    copystr((gene_list+tempg)->gene_type,gtf->attr[1].attribute_value);
    (gene_list+tempg)->gene_start=gtf->start;
    (gene_list+tempg)->gene_end=gtf->end;
    (gene_list+tempg)->gene_transcriptf_start=tempt+1;
}
void write_transcript(struct transcript*transcript_list,struct gtfinfo*gtf,int tempg,int tempt,int tempe,int tempc,int tempu,int tempo){
    if(tempt>0){
        transcript_list[tempt-1].transcript_exonf_end=tempe;
        transcript_list[tempt-1].transcript_CDSf_end=tempc;
        transcript_list[tempt-1].transcript_UTRf_end=tempu;
        transcript_list[tempt-1].transcript_otherf_end=tempo;
    }
    transcript_list[tempt].gene_flag=tempg;
    copystr((transcript_list+tempt)->transcript_id,gtf->attr[1].attribute_value);
    copystr((transcript_list+tempt)->transcript_type,gtf->attr[4].attribute_value);
    transcript_list[tempt].transcript_strand=gtf->strand;
    transcript_list[tempt].transcript_start=gtf->start;
    transcript_list[tempt].transcript_end=gtf->end;
    transcript_list[tempt].transcript_exonf_start=tempe+1;
    transcript_list[tempt].transcript_CDSf_start=tempc+1;
    transcript_list[tempt].transcript_UTRf_start=tempu+1;
    transcript_list[tempt].transcript_otherf_start=tempo+1;
    transcript_list[tempt].is_canonical=0;
    for(int i=0;i<gtf->attribute_number;i++){
        if(strcmp(gtf->attr[i].attribute_name,"tag")==0){
            if(strcmp(gtf->attr[i].attribute_value,"Ensembl_canonical")==0)
                transcript_list[tempt].is_canonical=1;
        }
    }
}
void write_exon(struct exon*exon_list,struct gtfinfo*gtf,int tempt,int tempe){
    exon_list[tempe].transcript_flag=tempt;
    copystr((exon_list+tempe)->exon_id,gtf->attr[7].attribute_value);
    (exon_list+tempe)->exon_number=atoi(gtf->attr[6].attribute_value);
    (exon_list+tempe)->exon_start=gtf->start;
    (exon_list+tempe)->exon_end=gtf->end;
    (exon_list+tempe)->exon_length=gtf->end-gtf->start+1;
}
void write_CDS(struct CDS*CDS_list,struct gtfinfo*gtf,int tempt,int tempe,int tempc){
    CDS_list[tempc].transcript_flag=tempt;
    CDS_list[tempc].exon_flag=tempe;
    (CDS_list+tempc)->CDS_start=gtf->start;
    (CDS_list+tempc)->CDS_end=gtf->end;
    (CDS_list+tempc)->CDS_length=gtf->end-gtf->start+1;
    (CDS_list+tempc)->CDS_frame=atoi(&gtf->frame);
}
void write_UTR(struct UTR*UTR_list,struct gtfinfo*gtf,int tempt,int tempe,int tempu){
    UTR_list[tempu].transcript_flag=tempt;
    UTR_list[tempu].exon_flag=tempe;
    (UTR_list+tempu)->UTR_start=gtf->start;
    (UTR_list+tempu)->UTR_end=gtf->end;
    (UTR_list+tempu)->UTR_length=gtf->end-gtf->start+1;
}
void write_other(struct other_feature*other_list,struct gtfinfo*gtf,int tempt,int tempe,int tempo){
    other_list[tempo].transcript_flag=tempt;
    other_list[tempo].exon_flag=tempe;
    copystr((other_list+tempo)->feature_type,gtf->feature);
    (other_list+tempo)->feature_start=gtf->start;
    (other_list+tempo)->feature_end=gtf->end;
    (other_list+tempo)->feature_length=gtf->end-gtf->start+1;
}
void write_close(struct gene*gene_list,struct transcript*transcript_list,int tempg,int tempt,int tempe,int tempc,int tempu,int tempo){
    gene_list[tempg].gene_transcriptf_end=tempt;
    transcript_list[tempt].transcript_exonf_end=tempe;
    transcript_list[tempt].transcript_CDSf_end=tempc;
    transcript_list[tempt].transcript_UTRf_end=tempu;
    transcript_list[tempt].transcript_otherf_end=tempo;
}
void read_start_pos(struct startpos*pos,struct gene*genelist,struct transcript*translist,char*chr,int trnasnum){
    int cnt=0;
    for(int i=0;i<trnasnum;i++){
        if(strcmp(genelist[translist[i].gene_flag].chromosome,chr)==0){
            copystr((pos+cnt)->chromosome,chr);
            copystr((pos+cnt)->transcript_id,translist[i].transcript_id);
            if(translist[i].transcript_strand=='+')
                (pos+cnt)->start_pos=translist[i].transcript_start;
            else
                (pos+cnt)->start_pos=translist[i].transcript_end;
            cnt++;
        }
    }
}
void sort_start_pos(struct startpos*pos,int num){
    for(int i=0;i<num;i++){
        for(int j=0;j<num-i-1;j++){
            if((pos+j)->start_pos>(pos+j+1)->start_pos){
                swapstr((pos+j)->chromosome,(pos+j+1)->chromosome);
                swapstr((pos+j)->transcript_id,(pos+j+1)->transcript_id);
                int temp=(pos+j)->start_pos;
                (pos+j)->start_pos=(pos+j+1)->start_pos;
                (pos+j+1)->start_pos=temp;
            }
        }
    }
}
