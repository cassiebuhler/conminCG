/*********************************************************************/
/***    Copyright (c) Robert J. Vanderbei, 1994                    ***/
/***    All Rights Reserved                                        ***/
/*********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

#include "hash.h"
#include "myalloc.h"

static int	hashsize, theta;
static NLIST	**hashtab;  /* pointer table for labels	*/

/***********************************************************
 * Prototype static functions 
 **********************************************************/

static unsigned	hash(
    char *s
);

char *my_strdup(
    char *s
);

/***********************************************************
 * hash():  form hash value for	string s 
 **********************************************************/

static unsigned	hash(
    char *s
)
{
    unsigned hashval;

    for	(hashval = 0; *s != '\0'; s++)
	hashval	= *s + theta * hashval;
    return hashval % hashsize;
}

/***********************************************************
 * hash_init():	 setup hash table 
 **********************************************************/

void hash_init (
    int	len
) 
{
    hashsize = len;
    theta = (int) (0.3874576628*hashsize);
    CALLOC( hashtab, hashsize, NLIST * );
}

/***********************************************************
 * getindex(): look for	s in hashtab and return	associated
 * index
 **********************************************************/

int getindex(
    char *s
)
{
    NLIST *np;

    np = lookup(s);
    if ( np == NULL ) {
	return -1;
    } else {
	return np->index;
    }
}

/***********************************************************
 * lookup(): look for s	in hashtab and return pointer to
 * structure
 **********************************************************/

NLIST *lookup(
    char *s
)
{
    NLIST *np;

    if (hashtab == NULL) hash_init(1);

    for	(np = hashtab[hash(s)];	np != NULL; np = np->next) {
	if (strcmp(s, np->name)	== 0) {
	    return np;		      /* found it  */
	}
    }
    return NULL;		      /* not found */
}

/***********************************************************
 * install(): put (name, index)	in hashtab 
 **********************************************************/

NLIST *install(
    char *s,
    int	 num
)
{
    NLIST *np;
    unsigned hashval;

    np = lookup(s);
    if ( np == NULL ) {			   /* not found	*/
	np = (NLIST *) malloc(sizeof(*np));
	if (np == NULL || (np->name = my_strdup(s)) == NULL) {
	    return NULL;		   /* no room */
	}
	hashval	= hash(s);
	np->next = hashtab[hashval];
	hashtab[hashval] = np;
    }
    np->index =	num;
    return np;
}

#ifdef __cplusplus
}
#endif
