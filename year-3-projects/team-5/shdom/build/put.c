/*  PUT - puts command line arguments to stdout, for use in batch files  */

#include <stdio.h>

main(argc,argv)
int argc;
char *argv[];
{
    int i;

    for (i = 1; i < argc; i++)
        fprintf (stdout, "%s\n", argv[i]);

}


