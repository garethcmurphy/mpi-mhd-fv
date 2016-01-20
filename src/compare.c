/*
  WRITTEN S O'SULLIVAN JUNE3 2003 GSFC
  EXECUTE THIS ROUTINE SUCCESSIVELY TO TKDIFF .F FILES
  LISTED IN MAKEFILE FOLLOWING "ajax_F_files ="
  MAINTAINS A RECORD OF NEXT FILE TO VISIT IN .info.compare
  EDIT COMMANDNAME TO CHANGE TKDIFF DIRECTORY
*/




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main()
{
     char 
	  ch, str[10], commandname[100];
     int n, cnt;

     FILE 
	  *finitial, *finfo, *finfo2;     
     
     
     if((finfo=fopen("./.info.compare","r"))==NULL)
       {
	 printf("problems reading .info.compare!!!!!!\n");
	 exit(0);
       }
     
     while (!feof(finfo)) {
       fscanf(finfo,"%d",&cnt);
     }

     if((finfo2=fopen("./.info.compare.tmp","w"))==NULL)
       {
	 printf("problems reading .info.compare.tmp!!!!!!\n");
	 exit(0);
       }


     system("rm .compare_filenames");
     system("ls *.cpp *.h > .compare_filenames");
     
     if((finitial=fopen("./.compare_filenames","r"))==NULL)
       {
	 printf("problems reading .compare_filenames!!!!!!\n");
	 exit(0);
       }
     

/*
 *     COMMENTS MAY BE PLACED ANYWHERE BEFORE MANDATORY :
 */

     n=0;
     
     while (!feof(finitial)) {

	 fscanf(finitial,"%s",str);
	 fscanf(finitial,"%s",str);
	 while (!feof(finitial)) { 

	   fscanf(finitial,"%s",str);
	   if (strcmp(str, "OBJS")==0) {
	     fprintf(finfo2,"%d\n",1);
	     break;
	   }else { /* start reading in filenames and counting */
	     n++;
	     if(n==cnt) { /* have found file number cnt */
	       fprintf(finfo2,"%d\n",n+1);
	       break;
	     }
	   }
	 }
	 if (strcmp(str, "OBJS")==0 || n==cnt) {
	   break;
	 }
     }
     
     fclose(finitial);
     fclose(finfo);
     fclose(finfo2);
     system("rm .info.compare");
     system("mv .info.compare.tmp .info.compare");

     commandname[0]='\0';
     strcat(commandname,"tkdiff ../../ss/Shamrock/src ");
     strcat(commandname,str);
     strcat(commandname," &");
     printf("%s\n",commandname);
     system(commandname);
     

     return(1);
}


