#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <iostream> //
#include <sstream>   //
#include <tgmath.h> //
#include <getopt.h>
#include <assert.h>
#include <sys/time.h>

//#include <omp.h>


#include <sys/resource.h>
#include <unistd.h>


#include <sys/resource.h>
#include <unistd.h>

using namespace std;

int main(int argc, char **argv)
{
    FILE * R;

    FILE * W;

    char * tempr=new char[100];
    strcpy ( tempr, argv[1] );

    char * tempw=new char[100];
    strcpy ( tempw, argv[2] );
    //W.open (tmpw);

    R=fopen(tempr,"r");W=fopen(tempw,"w");
    int i=0;
    if ( R != NULL )
   {
      char line [ 128*10 ]; /* or other suitable maximum line size */
      while ( fgets ( line, sizeof line, R ) != NULL ) /* read a line */
      {
          i++;
         //fputs ( line, stdout ); /* write the line */
         if (i%4==1)
        {//W<<">"<<line<<"\n";
        fprintf(W,">%s",line);}
         if (i%4==2){
        //W<<">"<<line<<"\n";
        fprintf(W,"%s",line);}

      }
      fclose ( R );
   }
   else
   {
      perror ( tempr ); /* why didn't the file open? */
   }

   cout << i<<endl;

  // W.close();
   fclose(W);
   delete[] tempr;
   delete[] tempw;

}
