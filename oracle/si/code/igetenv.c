/* FileName    : igetenv.c
   Objectfile  : igetenv.o
   Purpose     : This function provides FORTRAN programs with strings from
                 the shell environment.
   Calling Seq.: length = GETENV( EnvVarName, EnvVarValue )
   Input       : EnvVarName : The name of the environment variable, a
                              FORTRAN character-type expression.  It is *not*
                              necessary to terminate it with a null byte.
   Output      : EnvVarValue: The value of the environment variable, padded
                              on the right with blanks.  If the function call
                              was for some reason unsuccessful, this is a
                              blank string.  Provide a character variable of
                              sufficient length for this argument.
                 length     : If the function call was successful, `length'
                              contains the number of significant characters
                              (ignoring trailing blanks).  If the variable
                              `EnvVarValue' is of insufficient size, `length'
                              receives the negated length of the environment
                              variable's value (ie. the minimum size required
                              for `EnvVarValue').  If the environment variable
                              doesn't exist or some other error occurs, `length'
                              receives the value 0.

   Example     : C     Open a File (my.dat) on HOME-Directory

                       PROGRAM MyFortranProg
                       CHARACTER*100 env_val
                       INTEGER*4 igetenv,length
                       ......
                       length = igetenv( 'HOME', env_val )
                       if ( length .le. 0 ) stop
                       OPEN( UNIT=1, FILE=env_val(1:length)//'/my.dat',... )
                       ......
                       END


   Hints       : Compilation and linking proposition:

                      cc -c -Aa igetenv.c                ; -Aa: ANSI-Standard
                      f77  <my_opts>  MyFortranProg.f  igetenv.o

   Autor/Date: Johannes Reetz, Dirk Husfeld, 12-MAR-1993
   History   :


*/
#include  <stdlib.h>
#include  <string.h>

int igetenv( char * name, char *value, size_t nam_len, size_t val_len)
{
    int    i, env_len;
    char   *cname, *ep;

    for ( i=0; i < val_len; value[i++] = ' ' )
        ;    /* blank out the value area */

    /* convert the env variable's name into a C string (cname) */
    if ( (cname = (char *) malloc( nam_len+1 )) == NULL )
        return 0;    /* failed, no memory */
    strncpy( cname, name, nam_len );
    cname[nam_len] = '\0';
    cname[strcspn(cname," ")] = '\0';   /* drop trailing blanks */

    ep = getenv( cname );    /* enquire the environment for this name */

    free( (void *) cname );    /* free the memory */

    if ( ep == NULL )
        return 0;            /* name not found */

    if ( (env_len = strlen( ep )) > val_len )
        return (-env_len);     /* too long to fit into value */
    else {
       strncpy( value, ep, env_len );
       return env_len;         /* success */
       }
}
