#include  <stdlib.h>
#include  <string.h>

typedef struct {		/* Define string descriptor */
       unsigned short slen;		/* Length of string, 0 for null */
       int stype;			/* type of string, static or dynamic */
       char *s;			/* Addr of string */
      } STRING;

#define STRLEN(__str)    ((long)(__str)->slen)

void call_lineform(argc, argv)
int      argc;
void    *argv[];
   {                           /* corresponding variable */
     void lineform();
     STRING    *linefile,
               *head_atmdate,
               *head_atmuser,
               *head_atmref,
               *head_atmcmt,
               *head_lindate,
               *head_linuser,
               *head_linref,
               *head_lincmt;

     linefile	   = (STRING *) argv[12]; 
     head_atmdate  = (STRING *) argv[19];   
     head_atmuser  = (STRING *) argv[20];
     head_atmref   = (STRING *) argv[21];
     head_atmcmt   = (STRING *) argv[22];
     head_lindate  = (STRING *) argv[23];
     head_linuser  = (STRING *) argv[24];
     head_linref   = (STRING *) argv[25];
     head_lincmt   = (STRING *) argv[26];

     lineform( argv[0], argv[1], argv[2], argv[3], argv[4], argv[5], 
               argv[6],argv[7], argv[8], argv[9], argv[10],argv[11],
               linefile->s,
               argv[13],argv[14],argv[15],
               argv[16],argv[17],argv[18],
               head_atmdate->s,
               head_atmuser->s,
               head_atmref->s,
               head_atmcmt->s,
               head_lindate->s,
               head_linuser->s,
               head_linref->s,
               head_lincmt->s,
               STRLEN(linefile),
               STRLEN(head_atmdate),
               STRLEN(head_atmuser),
               STRLEN(head_atmref),
               STRLEN(head_atmcmt),
               STRLEN(head_lindate),
               STRLEN(head_linuser),
               STRLEN(head_linref),
               STRLEN(head_lincmt)    );
   }
