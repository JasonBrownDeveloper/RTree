#ifndef _RTREE_
#define _RTREE_

/*
 * Antonin Guttman: R-Trees: A Dynamic Index Structure for Spatial Searching
 * Proc. 1984 ACM SIGMOD International Conference on Management of Data, pp. 47-57.
 * ISBN 0-89791-128-8
 * http://www-db.deis.unibo.it/courses/SI-LS/papers/Gut84.pdf
 */

#include <stdbool.h> /*bool*/
#include <stddef.h> /*size_t*/
#include <stdint.h> /*uint_fast8_t, int_fast32_t*/

typedef int_fast32_t RTdimension;
typedef uint_fast8_t RTdimensionindex;
typedef uint_fast8_t RTchildindex;

#define RTn 2     /*dimensions; sizeof dimensionindex*/
#define RTPS 4096 /*Pagesize*/

struct RTNodeList {
   struct RTNodeList *Next;
   void *Tuple;
   RTdimension I[RTn*2];
};

struct RTNode;
typedef struct RTNode * RTreePtr;

bool RTNew(RTreePtr *T, struct RTNodeList *list);
bool RTFree(RTreePtr *T);

bool RTSelect(RTreePtr *T, RTdimension S[], struct RTNodeList **list, size_t *count);
bool RTInsert(RTreePtr *T, RTdimension I[], void *Tuple);
bool RTDelete(RTreePtr *T, RTdimension I[], void *Tuple);
bool RTUpdateTuple(RTreePtr *T, RTdimension I[], void *Tuple, void *New);
bool RTUpdateDimension(RTreePtr *T, RTdimension I[], void *Tuple, RTdimension New[]);

#ifdef RTREE_DEBUG
bool RTTrace(struct RTNode *Start, size_t Level, size_t AbsChild, struct RTNode **Out);
bool RTDump(struct RTNode *Start, const char *filename);
#endif

#endif /* _RTREE_ */
