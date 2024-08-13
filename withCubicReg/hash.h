/*********************************************************************/
/***    Copyright (c) Robert J. Vanderbei, 1994                    ***/
/***    All Rights Reserved                                        ***/
/*********************************************************************/

struct nlist {		/* table entry */
    struct nlist *next;	/* next	entry in linked	list */
    char	 *name;	/* defining label */
    int		 index;	/* pointer to numerical	index */
};

typedef	struct nlist NLIST;

void hash_init(
    int	len
);

int getindex(
    char *s
);

NLIST *lookup(
    char *s
);

NLIST *install(
    char *s,
    int	 num
);
