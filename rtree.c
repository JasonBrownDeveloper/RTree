#include <stdio.h>  /*fprintf, fputs*/
#include <stdlib.h> /*malloc, free, NULL, exit*/
#include <string.h> /*memcpy, memmove, memset*/
#include <float.h>  /*LDBL_MAX*/
#include "rtree.h"

#define m 2                        /*minimum children; sizeof childindex*/
#define M (RTPS / sizeof(struct RTNode)) /*maximum children; sizeof childindex*/

#define LEVEL_TOP -1
#define LEVEL_LEAF 1
#define LEVEL_TUPLE 0

#define IS_BRANCH(N) ((N)->Child && (N)->Child[0].Child && (N)->Child[0].Tuple == NULL)
#define IS_LEAF(N)   ((N)->Child && (N)->Child[0].Tuple && (N)->Child[0].Child == NULL)
#define IS_TUPLE(N)  ((N)->Tuple && (N)->Child == NULL)
#define IS_EMPTY(N)  ((N).Parent == NULL) /* Will return true for root node */

static struct RTNode {
   struct RTNode *Parent;
   struct RTNode *Child;
   void *Tuple;
   RTdimension I[RTn*2]; /*{x1,y1,...,x2,y2...}*/
} EMPTY_NODE = {
   .Parent = NULL,
   .Child = NULL,
   .Tuple = NULL,
   .I = {0}
};

static void *mem_alloc(size_t _Size);

static bool Overlap(RTdimension *S1, RTdimension *S2);
static bool Within(RTdimension *S1, RTdimension *S2);
static long double safe_multiply(long double left, long double right);
static long double Volume(RTdimension *S);

static bool InitNodes(struct RTNode *ptr, RTchildindex size);
static bool FreeNodes(struct RTNode *T);

static bool Search(struct RTNode *T, RTdimension S[], struct RTNodeList **list, size_t *count);
static bool Insert(struct RTNode **N, size_t Level, RTdimension I[], void *Tuple, struct RTNode *Branch);
static bool ChooseLeaf(struct RTNode *N, size_t Start, size_t Stop, RTdimension *I, struct RTNode **leaf);
static bool AdjustTree(struct RTNode *N, struct RTNode *NN, struct RTNode **root, struct RTNode **split);
static bool Delete(struct RTNode **T, RTdimension I[], void *Tuple);
static bool FindLeaf(struct RTNode *T, RTdimension I[], void *Tuple, struct RTNode **L, RTchildindex *position);
static bool CondenseTree(struct RTNode *N, struct RTNode **root);
static bool LinearSplit(struct RTNode *L, RTdimension I[], void *Tuple, struct RTNode *Child, struct RTNode **split);
static bool LinearPickSeeds(struct RTNode NL[], long double *width, struct RTNode **hbest, struct RTNode **lbest);

/*Wrapper for malloc checks for out of memory*/
static void *mem_alloc(size_t size) {
   void *mem = malloc(size);
   if (!mem) {
      fputs("fatal: out of memory.\n", stderr);
      exit(EXIT_FAILURE);
   }
   memset(mem, 0, size);
   return mem;
}

/*Returns false if the two shapes don't overlap and TRUE if they do*/
static bool Overlap(RTdimension *S1, RTdimension *S2) {
   RTdimensionindex j = 0, k = 0;

   for (j = 0, k = RTn; j < RTn; ++j, ++k)
      if (S1[j] > S2[k] || S2[j] > S1[k])
         return false;

   return true;
};

/*Returns false if the shape 1 isn't within shape 2 and TRUE if it is*/
static bool Within(RTdimension *S1, RTdimension *S2) {
   RTdimensionindex j = 0, k = 0;

   for (j = 0, k = RTn; j < RTn; ++j, ++k)
      if (S1[j] < S2[j] || S2[k] < S1[k])
         return false;

   return true;
};

static long double safe_multiply(long double left, long double right) {
   int sign = 1;
   if (left == 0 || right == 0) return 0;
   if (left < 0) {left = -left; sign = -sign;}
   if (right < 0) {right = -right; sign = -sign;}
   if (LDBL_MAX / right < left) {
      fputs("fatal: long double overflow\n", stderr);
      exit(EXIT_FAILURE);
   }

   return sign * left * right;
}

/*Returns volume of a shape*/
static long double Volume(RTdimension *S) {
   RTdimensionindex j = 0, k = 0;
   long double volume = 1;

   for (j = 0, k = RTn; j < RTn; ++j, ++k)
      /*1 added to each side to make lines and dots have volume*/
      volume = safe_multiply(volume, (long double)S[k] - S[j] + 1);

   return volume;
}

/*Sets the Tuple of a Node*/
/*In: Parent Node, Size, Tuple, New Tuple */
bool UpdateTuple(RTreePtr *T, RTdimension I[], void *Tuple, void *New) {
   struct RTNode *L = NULL;
   RTchildindex pos = 0;

   FindLeaf(*T, I, Tuple, &L, &pos);

   if (!L)
      return false;

   ((L)->Child+pos)->Tuple = New;
   return true;
}

/*Sets the Dimensions of a node*/
/*In: Parent Node, Size, Tuple, New Dimension */
bool RTUpdateDimensions(RTreePtr *T, RTdimension I[], void *Tuple, RTdimension New[]) {
   struct RTNode *L = NULL;
   RTchildindex pos = 0;

   FindLeaf(*T, I, Tuple, &L, &pos);

   if (!L)
      return false;

   if (Within((RTdimension *)New, (L)->I))
      memcpy(((L)->Child+pos)->I, New, sizeof(((L)->Child+pos)->I));
   else {
      Delete(T, I, Tuple);
      Insert(T, LEVEL_LEAF, New, Tuple, NULL);
   }

   return true;
}

/*Gets the Dimensions of an RTree*/
/*In: Parent Node  Out: Size */
bool RTSelectDimensions(RTreePtr *T, RTdimension I[]) {
   memcpy(I, (*T)->I, sizeof((*T)->I));
   return true;
}

/*Creates a new Tree*/
bool RTNewTree(struct RTNode **T, struct RTNodeList *list) {
   RTchildindex h;
   size_t i, j, k;
   struct RTNodeList *nodelist;
   struct RTNode *node;

   struct RTNode *stack, *stack_next, *branch, *branch_next;

   stack = NULL;
   nodelist = list;
   while (nodelist != NULL) {
      node = (struct RTNode *)mem_alloc(sizeof(struct RTNode));
      InitNodes(node, 1);
      node->Child = (struct RTNode *)mem_alloc(M * sizeof(node->Child[0]));
      InitNodes(node->Child, M);

      memcpy(node->I, nodelist->I, sizeof(node->I));
      for (i = 0; nodelist != NULL && i < M; ++i) {
         memcpy(node->Child[i].I, nodelist->I, sizeof(node->Child[i].I));
         node->Child[i].Tuple = nodelist->Tuple;
         node->Child[i].Parent = node;

         for (j = 0, k = RTn; j < RTn; ++j, ++k) {
            node->I[j] = node->I[j] < nodelist->I[j] ? node->I[j] : nodelist->I[j];
            node->I[k] = node->I[k] > nodelist->I[k] ? node->I[k] : nodelist->I[k];
         }

         nodelist = nodelist->Next;
      }

      node->Parent = stack;
      stack = node;
   }

   while (stack && stack->Parent != NULL) {
      branch = stack;
      stack_next = NULL;
      while(branch != NULL) {
         node = (struct RTNode *)mem_alloc(sizeof(struct RTNode));
         InitNodes(node, 1);
         node->Child = (struct RTNode *)mem_alloc(M * sizeof((*T)->Child[0]));
         InitNodes(node->Child, M);

         memcpy(node->I, branch->I, sizeof(node->I));
         for (i = 0; branch != NULL && i < M; ++i) {
            branch_next = branch->Parent;

            memcpy(node->Child+i, branch, sizeof(node->Child[i]));
            node->Child[i].Parent = node;

            if (node->Child[i].Child != NULL)
               for (h = 0; h < M && !IS_EMPTY(node->Child[i].Child[h]); ++h)
                  node->Child[i].Child[h].Parent = node->Child+i;

            for (j = 0, k = RTn; j < RTn; ++j, ++k) {
               node->I[j] = node->I[j] < branch->I[j] ? node->I[j] : branch->I[j];
               node->I[k] = node->I[k] > branch->I[k] ? node->I[k] : branch->I[k];
            }

            free(branch);
            branch = branch_next;
         }
         node->Parent = stack_next;
         stack_next = node;
      }
      stack = stack_next;
   }

   (*T) = stack;

   if ((*T) == NULL) {
      (*T) = (struct RTNode *)mem_alloc(sizeof(struct RTNode));
      InitNodes(*T, 1);
      (*T)->Child = (struct RTNode *)mem_alloc(M * sizeof((*T)->Child[0]));
      InitNodes((*T)->Child, M);
   }

   return true;
}

/*Initializes an array of struct Nodes*/
static bool InitNodes(struct RTNode *ptr, RTchildindex size) {
   RTchildindex i;

   for (i = 0; i < size; ++i) {
      ptr[i] = EMPTY_NODE;
   }

   return true;
}

/*Frees the children of a Tree*/
static bool FreeNodes(struct RTNode *T) {
   if (IS_BRANCH(T)) {
      RTchildindex i;
      for (i = 0; i < M && !IS_EMPTY(T->Child[i]); ++i)
         FreeNodes(T->Child+i);

      free(T->Child);
      return true;
   /*Property (5) - Root and Leaf*/
   } else if (IS_LEAF(T) || IS_EMPTY(T->Child[0])) {
      free(T->Child);
      return true;
   }

   fputs("rtree on fire!\n", stderr);
   return false;
}

/*Frees a Tree*/
bool RTFreeTree(RTreePtr *T) {
   if (!T || !*T)
      return true;

   if (!FreeNodes(*T))
      return false;

   free(*T);
   *T = NULL;
   return true;
}

/*3.1 Searching*/
/*Algorithm Search*/
/*In: Parent Node, Search Box  Out: Hit List, Hit Count*/
bool RTSelectTuple(RTreePtr *T, RTdimension S[], struct RTNodeList **list, size_t *count) {
   if (!T || !*T) {
      fputs("RTree cannot be NULL.\n", stderr);
      return false;
   }

   if (!list && !count) {
      fputs("Must have List and/or Count.\n", stderr);
      return false;
   }

   return Search(*T, S, list, count);
}

/*In: Parent Node, Search Box  Out: Hit List, Hit Count*/
static bool Search(struct RTNode *T, RTdimension S[], struct RTNodeList **list, size_t *count) {
   struct RTNodeList *curr = NULL;
   size_t cnt = 0;
   RTchildindex i = 0;

   if (count) *count = 0;

   /*S1 [Search subtrees]*/
   if (IS_BRANCH(T)) {
      for (i = 0; i < M && !IS_EMPTY(T->Child[i]); ++i)
         if (Overlap(T->Child[i].I, S)) {
            if (!Search(T->Child+i, S, list, &cnt)) {
               if (list) *list = NULL;
               if (count) *count = 0;
               return false;
            }
            if (count) *count += cnt;
         }

      return true;

   /*S2 [Search leaf node]*/
   /*Property (5) - Root and Leaf*/
   } else if (IS_LEAF(T) || IS_EMPTY(T->Child[0])) {
      for (i = 0; i < M && !IS_EMPTY(T->Child[i]); ++i)
         if (Overlap(T->Child[i].I, S)) {
            if (list) {
               curr = (struct RTNodeList *)mem_alloc(sizeof(struct RTNodeList));
               memcpy(curr->I, T->Child[i].I, sizeof(curr->I));
               curr->Tuple = (T->Child+i)->Tuple;
               curr->Next = *list;
               *list = curr;
            }
            if (count) (*count)++;
         }

      return true;
   }

   fputs("rtree on fire!\n", stderr);
   if (list) *list = NULL;
   if (count) *count = 0;
   return false;
}

/*3.2 Insertion*/
/*Algorithm Insert*/
/*In: Parent Node, Size, Tuple */
bool RTInsertTuple(RTreePtr *N, RTdimension I[], void *Tuple) {
   if (!N || !*N) {
      fputs("RTree cannot be NULL.\n", stderr);
      return false;
   }

   if (!I || !Tuple) {
      fputs("Size and Tuple cannot be NULL.\n", stderr);
      return false;
   }

   /*TODO Guarantee x1 < x2 && y1 < y2 && ...*/

   return Insert(N, LEVEL_LEAF, I, Tuple, NULL);
}

/*In: Parent Node, Level, [Size, Tuple || TupleNode]*/
static bool Insert(struct RTNode **N, size_t Level, RTdimension I[], void *Tuple, struct RTNode *TupleNode) {
   struct RTNode *L = NULL, *LL = NULL, *splitL = NULL, *splitR = NULL, *newRoot = NULL;
   RTdimension *Size;
   size_t Start;
   RTchildindex h, i;

   Size = I;
   if (TupleNode)
      Size = TupleNode->I;
#ifdef RTREE_DEBUG
   else if (!Tuple) {
      fputs("rtree on fire!\n", stderr);
      *N = NULL;
      return false;
   }
#endif

   for (Start = 0, L = *N; L->Child; ++Start, L = L->Child);

   /*Property (5) - Root and Leaf*/
   if (Level == LEVEL_TOP)
      Level = Start;

   /*I1 [Find position for new record]*/
   if (!ChooseLeaf(*N, Start, Level, Size, &L)) {
      *N = NULL;
      return false;
   }

   /*I2 [Add record to leaf node]*/
   for (i = 0; i < M; ++i)
      if (IS_EMPTY(L->Child[i])) {
         if (Tuple) {
            memcpy(L->Child[i].I, I, sizeof(L->Child[i].I));
            L->Child[i].Child = NULL;
            L->Child[i].Tuple = Tuple;
#ifdef RTREE_DEBUG
         } else if (TupleNode) {
#else
         } else {
#endif
            memcpy(L->Child+i, TupleNode, sizeof(L->Child[i]));
         }
#ifdef RTREE_DEBUG
         else {
            fputs("rtree on fire!\n", stderr);
            *N = NULL;
            return false;
         }
#endif

         /* update the newly added child's parent */
         L->Child[i].Parent = L;

         /* update the newly added child's children's parent because of memcpy */
         if (L->Child[i].Child != NULL)
            for (h = 0; h < M && !IS_EMPTY(L->Child[i].Child[h]); ++h)
               L->Child[i].Child[h].Parent = L->Child+i;

         break;
      }

   if (i == M) {
      if (!LinearSplit(L, I, Tuple, TupleNode, &LL)) {
         *N = NULL;
         return false;
      }
   }

   /*I3 [Propagate changes upward]*/
   if (!AdjustTree(L, LL, &splitL, &splitR)) {
      *N = NULL;
      return false;
   }

   /*I4 [Grow tree taller]*/
   if (splitR) {
      newRoot = (struct RTNode *)mem_alloc(sizeof(struct RTNode));
      InitNodes(newRoot, 1);
      newRoot->Child = (struct RTNode *)mem_alloc(M * sizeof(struct RTNode));
      InitNodes(newRoot->Child, M);

      if (!Insert(&newRoot, LEVEL_TOP, NULL, NULL, splitL) || !Insert(&newRoot, LEVEL_TOP, NULL, NULL, splitR)) {
         *N = NULL;
         return false;
      }

      free(splitL);
      free(splitR);
      *N = newRoot;
      return true;
   }

   *N = splitL;
   return true;
}

/*Algorithm ChooseLeaf*/
/*CL1 [Initialize]*/
/*In: Parent Node, Tuple Size  Out: Chosen Leaf*/
static bool ChooseLeaf(struct RTNode *N, size_t Start, size_t Stop, RTdimension *I, struct RTNode **leaf) {
   RTchildindex i = 0;
   RTdimensionindex j = 0, k = 0;
   long double area = 0, increase = 0, min = 0, minarea = 0;
   RTdimension expanded[RTn*2];
   struct RTNode *F = NULL;

   /*Property (5) - Root and Leaf*/
   /*CL2 [Leaf check]*/
   if (Start == Stop) {
      *leaf = N;
      return true;

   /*CL3 [Choose subtree]*/
   } else if (IS_BRANCH(N)) {
      min = N->I[RTn];
      minarea = Volume(N->Child[0].I);
      F = N->Child;

      for (i = 0; i < M && !IS_EMPTY(N->Child[i]); ++i) {
         area = Volume(N->Child[i].I);
         for (j = 0, k = RTn; j < RTn; ++j, ++k) {
            expanded[j] = ( I[j] < N->Child[i].I[j] ) ? I[j] : N->Child[i].I[j];
            expanded[k] = ( I[k] > N->Child[i].I[k] ) ? I[k] : N->Child[i].I[k];
         }
         increase = Volume(expanded) - area;

         if (increase < min || (increase == min && area < minarea)) {
            min = increase;
            minarea = area;
            F = N->Child+i;
         }
      }

      /*CL4 [Descend until a leaf is reached]*/
      return ChooseLeaf(F, Start-1, Stop, I, leaf);
   }

   fputs("rtree on fire!\n", stderr);
   *leaf = NULL;
   return false;
}

/*Algorithm AdjustTree*/
/*AT1 [Initialize]*/
/*In: Start Node, Split Node  Out: Root Node, Split Node*/
static bool AdjustTree(struct RTNode *N, struct RTNode *NN, struct RTNode **root, struct RTNode **split) {
   struct RTNode *P = NULL, *LS = NULL;
   RTchildindex h, i;
   RTdimensionindex j, k;

   /*AT2 [Check if done]*/
   if (N == NULL) {
      *split = NN;
      return true;
   }

   /*AT3 [Adjust covering rectangle in parent entry]*/
   P = N->Parent;
   memcpy(N->I, N->Child[0].I, sizeof(N->I));
   for (i = 0; i < M && !IS_EMPTY(N->Child[i]); ++i)
      for (j = 0, k = RTn; j < RTn; ++j, ++k) {
         if (N->I[j] > N->Child[i].I[j])
            N->I[j] = N->Child[i].I[j];
         if (N->I[k] < N->Child[i].I[k])
            N->I[k] = N->Child[i].I[k];
      }

   /*AT4 [Propagate node split upward]*/
   if (P && NN) {
      for (i = 0; i < M; ++i)
         if (IS_EMPTY(P->Child[i])) {
            memcpy(P->Child+i, NN, sizeof(P->Child[i]));

            /* update the newly added child's parent */
            P->Child[i].Parent = P;

            /* update the newly added child's children's parent because of memcpy */
            if (P->Child[i].Child != NULL)
               for (h = 0; h < M && !IS_EMPTY(P->Child[i].Child[h]); ++h)
                  P->Child[i].Child[h].Parent = P->Child+i;

            free(NN);
            NN = NULL;
            break;
         }

      if (i == M) {
         if (!LinearSplit(P, NULL, NULL, NN, &LS)) {
            *root = NULL;
            *split = NULL;
            return false;
         }
         free(NN);
         NN = LS;
      }
   } else if (P == NULL && NN) {
      *root = N;
      *split = NN;
      return true;
   }

   /*AT5 [Move up to next level]*/
   *root = N;
   return AdjustTree(P, NN, root, split);
}

/*In: Parent Node, Dead Size, Dead Tuple */
bool RTDeleteTuple(RTreePtr *T, RTdimension I[], void *Tuple) {
   if (I == NULL || Tuple == NULL) {
      fputs("Must have Size and Tuple.\n", stderr);
      return false;
   }

   return Delete(T, I, Tuple);
}

/*3.3 Deletion*/
/*Algorithm Delete*/
/*In: Parent Node, Dead Size, Dead Tuple */
static bool Delete(struct RTNode **T, RTdimension I[], void *Tuple) {
   struct RTNode *L = NULL, *newRoot = NULL;
   RTchildindex i, pos = 0;

   /*D1 [Find node containing record]*/
   if (!FindLeaf(*T, I, Tuple, &L, &pos)) {
      *T = NULL;
      return false;
   }

   /*If the tuple wasn't found*/
   if (pos == M || L == NULL)
      return false;

   /*D2 [Delete record]*/
   memmove(L->Child+pos, L->Child+pos+1, (M - 1 - pos) * sizeof(L->Child[pos]));
   InitNodes(L->Child+M-1, 1);

   /*D3 [Propagate changes]*/
   if (!CondenseTree(L, T)) {
      *T = NULL;
      return false;
   }

   /*Property (5) - Root and Leaf*/
   /*D4 [Shorten tree]*/
   if (IS_EMPTY((*T)->Child[1]) && IS_BRANCH(*T)) {
      newRoot = (struct RTNode *)mem_alloc(sizeof(struct RTNode));
      memcpy(newRoot, (*T)->Child+0, sizeof(*newRoot));

      /* update the newly added child's children's parent because of memcpy */
      for (i = 0; i < M && !IS_EMPTY(newRoot->Child[i]); ++i)
         newRoot->Child[i].Parent = newRoot;

      free((*T)->Child);
      free(*T);
      newRoot->Parent = NULL;
      *T = newRoot;
   }

   return true;
}

/*Algorithm FindLeaf*/
/*In: Parent Node, Size, Tuple  Out: Leaf Node, Position*/
static bool FindLeaf(struct RTNode *T, RTdimension I[], void *Tuple, struct RTNode **L, RTchildindex *position) {
   RTchildindex i;

   *L = NULL;
   *position = M;

   /*FL1 [Search subtrees]*/
   if (IS_BRANCH(T)) {
      for (i = 0; i < M && !IS_EMPTY(T->Child[i]); ++i)
         if (Overlap(T->Child[i].I, I)) {
            if (!FindLeaf(T->Child+i, I, Tuple, L, position)) {
               *position = M;
               *L = NULL;
               return false;
            }
            if (*position < M) {
               return true;
            }
         }
      return true;

   /*FL2 [Search leaf node for record]*/
   /*Property (5) - Root and Leaf*/
   } else if (IS_LEAF(T) || IS_EMPTY(T->Child[0])) {
      for (i = 0; i < M && !IS_EMPTY(T->Child[i]); ++i)
         if (T->Child[i].Tuple == Tuple && !memcmp(T->Child[i].I, I, sizeof(T->Child[i].I))) {
            *L = T;
            *position = i;
            return true;
         }
      return true;
   }

   fputs("rtree on fire!\n", stderr);
   *L = NULL;
   *position = M;
   return false;
}

/*Algorithm CondenseTree*/
/*CT1 [Initialize]*/
/*In: Shrunk Node  Out: root*/
static bool CondenseTree(struct RTNode *N, struct RTNode **root) {
   struct CTNodeList
   {
      int level;
      struct RTNode *Node;
      struct CTNodeList *Next;
   };

   struct CTNodeList *Q = NULL, *ptr = NULL, *next = NULL;
   struct RTNode *P = NULL;
   RTchildindex g = 0, h = 0, i = 0;
   RTdimensionindex j = 0, k = 0;
   size_t level = 1;

   /*CT2 [Find parent entry]*/
   while (N->Parent) {
      P = N->Parent;

      /*CT3 [Eliminate under-full node]*/
      for (i = 0; i < m && N->Child[i].Parent; ++i) ;

      if (i < m) {
         for (i = 0; i < m && !IS_EMPTY(N->Child[i]); ++i) {
            ptr = (struct CTNodeList *)mem_alloc(sizeof(struct CTNodeList));
            ptr->Node = (struct RTNode *)mem_alloc(sizeof(struct RTNode));
            ptr->level = level;
            memcpy(ptr->Node, N->Child+i, sizeof(*ptr->Node));
            ptr->Next = Q;
            Q = ptr;
         }

         for (i = 0; i < M && !IS_EMPTY(P->Child[i]); ++i)
            if (P->Child+i == N)
               break;

#ifdef RTREE_DEBUG
         if (i == M) {
            fputs("rtree on fire!\n", stderr);
            if (root) *root = NULL;
            return false;
         }
#endif

         free(P->Child[i].Child);
         memmove(P->Child+i, P->Child+i+1, (M - 1 - i) * sizeof(P->Child[i]));
         InitNodes(P->Child+M-1, 1);

         /* update the children's children's parent because of memmove */
         for (g = 0; g < M && !IS_EMPTY(P->Child[g]); ++g)
            for (h = 0; h < M && !IS_EMPTY(P->Child[g].Child[h]); ++h)
               P->Child[g].Child[h].Parent = P->Child+g;

      /*CT4 [Adjust covering rectangle]*/
      } else {
         memcpy(N->I, N->Child[0].I, sizeof(N->I));
         for (i = 0; i < M && !IS_EMPTY(N->Child[i]); ++i)
            for (j = 0, k = RTn; j < RTn; ++j, ++k) {
               if (N->I[j] > N->Child[i].I[j])
                  N->I[j] = N->Child[i].I[j];
               if (N->I[k] < N->Child[i].I[k])
                  N->I[k] = N->Child[i].I[k];
            }
      }

      /*CT5 [Move up one level in tree]*/
      N = P;
      ++level;
   }

   /*CT6 [Re-insert orphaned entries]*/
   for ( ; Q != NULL; Q = next) {
#ifdef RTREE_DEBUG
      if (IS_TUPLE(Q->Node) || IS_LEAF(Q->Node) || IS_BRANCH(Q->Node)) {
#endif
         if (!Insert(&N, Q->level, NULL, NULL, Q->Node))
            return false;
#ifdef RTREE_DEBUG
      } else {
         fputs("rtree on fire!\n", stderr);
         if (root) *root = NULL;
         return false;
      }
#endif

      next = Q->Next;
      free(Q->Node);
      free(Q);
   }

   if (root) *root = N;
   return true;
}

/*3.5.3 A Linear-Cost Algorithm*/
/*Algorithm LinearSplit*/
/*In: Full Node, [Extra Size, Extra Tuple || Extra Branch]  Out: Split Node*/
static bool LinearSplit(struct RTNode *L, RTdimension I[], void *Tuple, struct RTNode *Child, struct RTNode **split) {
   struct RTNode NL[M+1], *seedA = NULL, *seedB = NULL, *LL = NULL, *Parent = NULL;
   RTchildindex i = 0, A = 1, B = 1;
   RTdimensionindex j = 0, k = 0;
   long double Larea = 0, Lincrease = 0, LLarea = 0, LLincrease = 0, width[RTn];
   RTdimension Lexpanded[RTn*2], LLexpanded[RTn*2];

   /*Copy children into bigger array*/
   memcpy(NL, L->Child, M * sizeof(NL[0]));

   /*Copy extra to the end of the array*/
   if (Tuple) {
      NL[M].Parent = NULL;
      memcpy(NL[M].I, I, sizeof(NL[M].I));
      NL[M].Child = NULL;
      NL[M].Tuple = Tuple;
#ifdef RTREE_DEBUG
   } else if (Child) {
#else
   } else {
#endif
      memcpy(NL+M, Child, sizeof(NL[M]));
   }
#ifdef RTREE_DEBUG
   else {
      fputs("rtree on fire!\n", stderr);
      *split = NULL;
      return false;
   }
#endif

   Parent = L->Parent;

   /*Find the width of all sides of L if E was a child*/
   for (j = 0, k = RTn; j < RTn; ++j, ++k) {
      width[j] = (long double)((NL[M].I[k] > L->I[k]) ? NL[M].I[k] : L->I[k]) - ((NL[M].I[j] < L->I[j]) ? NL[M].I[j] : L->I[j]);
   }

   /*
    * From here all the children are in an array and we know how much area is covered
    */

   LL = (struct RTNode *)mem_alloc(sizeof(struct RTNode));
   InitNodes(LL, 1);
   LL->Child = (struct RTNode *)mem_alloc(M * sizeof(struct RTNode));
   InitNodes(LL->Child, M);

   /*LS1 [Pick first entry for each group]*/
   if (!LinearPickSeeds(NL, width, &seedA, &seedB)) {
      *split = NULL;
      return false;
   }

   /*Last check passed clear L*/
   InitNodes(L->Child, M);
   memset(L->I, 0, sizeof(L->I));
   L->Parent = NULL;
   L->Tuple = NULL;

   if(!Insert(&L, LEVEL_TOP, NULL, NULL, seedA) || !Insert(&LL, LEVEL_TOP, NULL, NULL, seedB)) {
      *split = NULL;
      return false;
   }

   /*LS2 [Check if done]*/
   /*for loop checks 'if all entries have been assigned' for LS2 and acts as the PickNext algorithm*/
   for (i = 0; i < M+1; ++i) {
      if (NL+i != seedA && NL+i != seedB) {
         if ((M + 1) - (A + B) == m - A) {
            if(!Insert(&L, LEVEL_TOP, NULL, NULL, NL+i)) {
               *split = NULL;
               return false;
            }
         } else if ((M + 1) - (A + B) == m - B) {
            if(!Insert(&LL, LEVEL_TOP, NULL, NULL, NL+i)) {
               *split = NULL;
               return false;
            }
         } else {
            /*LS3 [Select entry to assign]*/
            Larea = Volume(L->I);
            LLarea = Volume(LL->I);

            for (j = 0, k = RTn; j < RTn; ++j, ++k) {
               Lexpanded[j] = ( NL[i].I[j] < L->I[j] )  ? NL[i].I[j] : L->I[j];
               Lexpanded[k] = ( NL[i].I[k] > L->I[k] )  ? NL[i].I[k] : L->I[k];
               LLexpanded[j] = ( NL[i].I[j] < LL->I[j] ) ? NL[i].I[j] : LL->I[j];
               LLexpanded[k] = ( NL[i].I[k] > LL->I[k] ) ? NL[i].I[k] : LL->I[k];
            }

            Lincrease = Volume(Lexpanded) - Larea;
            LLincrease = Volume(LLexpanded) - LLarea;

            if (Lincrease < LLincrease || (Lincrease == LLincrease && Larea < LLarea)) {
               if(!Insert(&L, LEVEL_TOP, NULL, NULL, NL+i)) {
                  *split = NULL;
                  return false;
               }
               ++A;
            } else if (Lincrease > LLincrease || (Lincrease == LLincrease && Larea > LLarea)) {
               if(!Insert(&LL, LEVEL_TOP, NULL, NULL, NL+i)) {
                  *split = NULL;
                  return false;
               }
               ++B;
            } else if (A < B) {
               if(!Insert(&L, LEVEL_TOP, NULL, NULL, NL+i)) {
                  *split = NULL;
                  return false;
               }
               ++A;
            } else {
               if(!Insert(&LL, LEVEL_TOP, NULL, NULL, NL+i)) {
                  *split = NULL;
                  return false;
               }
               ++B;
            }
         }
      }
   }

   L->Parent = Parent;
   *split = LL;
   return true;
}

/*Algorithm LinearPickSeeds*/
/*In: Node List, Node Width  Out: High Best, Low Best*/
static bool LinearPickSeeds(struct RTNode NL[], long double *width, struct RTNode **hbest, struct RTNode **lbest) {
   long double separation = 0.0, sbest = -1.0;
   RTchildindex i = 0, low = 0, high = 1;
   RTdimensionindex j = 0, k = 0;

   /*LPS1 [Find extreme rectangles along all dimensions]*/
   for (j = 0, k = RTn; j < RTn; ++j, ++k) {
      for (i = 0; i < M; ++i) {
         if (NL[i].I[j] > NL[low].I[j])
            if (i != high)
               low = i;
         if (NL[i].I[k] < NL[high].I[k])
            if (i != low)
               high = i;
      }

      /*LPS2 [Adjust for shape of the rectangle]*/
      separation = ((long double)NL[low].I[j] - NL[high].I[k]) / width[j];

      /*LPS3 [Select the most extreme pair]*/
      if (separation > sbest) {
         *hbest = NL+high;
         *lbest = NL+low;
         sbest = separation;
      }
   }

   if (*hbest != *lbest)
      return true;

   fputs("rtree on fire!\n", stderr);
   *hbest = NULL;
   *lbest = NULL;
   return false;
}

#ifdef RTREE_DEBUG
#include <math.h>   /*pow*/

bool RTTrace(struct RTNode *Start, size_t Level, size_t AbsChild, struct RTNode **Out) {
   size_t count, group, place, last, factor;

   last = 0;
   for ( ; Level > 0; Level--) {
      count = (size_t)pow(M, Level);
      group = count / M;
      place = group + last;
      for (factor = 0; factor < M; ++factor)
         if (AbsChild + 1 <= place + group * factor)
            break;

      if (factor == M) {
         fputs("Bad Trace.\n", stderr);
         *Out = NULL;
         return false;
      }

      Start = Start->Child+factor;
      if (Start->Parent == NULL) {
         *Out = NULL;
         return true;
      }
      last += group * factor;
   }

   *Out = Start;
   return true;
}

bool RTDump(struct RTNode *Start, const char *filename) {
   size_t i, j, level, max, size, cnt;
   RTdimensionindex k;
   char *pad = NULL;
   struct RTNode *curr;
   FILE *log = NULL;

   if ((log = fopen(filename, "w")) == NULL)
      return false;

   for (level = 0, curr = Start; curr; ++level, curr = curr->Child);

   max = (size_t)pow(M, level - 1) * 16;
   pad = (char *)mem_alloc((max / 2 - 8) * sizeof(char) + 1);
   memset(pad, ' ', (max / 2 - 8) * sizeof(pad[0]));
   pad[(max / 2 - 8)] = '\0';

   fprintf(log, "%s %14p %s\n", pad, Start, pad);
   for(k = 0; k < RTn*2; ++k) {
      fprintf(log, "%s %d-%12d %s\n", pad, k, Start->I[k], pad);
   }
   fprintf(log, "%s %14p %s\n", pad, Start->Parent, pad);
   fprintf(log, "%s %14p %s\n", pad, Start->Child, pad);
   fprintf(log, "%s %14p %s\n\n", pad, Start->Tuple, pad);

   for (i = 1; i < level; ++i) {
      cnt = (size_t)pow(M, i);
      size = (max / cnt - 16);
      pad[size] = '\0';

      for (j = 0; j < cnt; ++j) {
         if (!trace(Start, i, j, &curr)) {
            fclose(log);
            return false;
         }
         if (j == 0 && size > 0)
            pad[size/2] = '\0';

         if(curr)
            fprintf(log, "%s %14p ", pad, curr);
         else
            fprintf(log, "%s XXXXXXXXXXXXXX ", pad);

         if (j == 0 && size > 0)
            pad[size/2] = ' ';

         if (j == cnt - 1 && size > 0) {
            pad[size/2] = '\0';
            fputs(pad, log);
            pad[size/2] = ' ';
         }
      }
      fputs("\n", log);

      for (k = 0; k < RTn*2; ++k) {
         for (j = 0; j < cnt; ++j) {
            if (!trace(Start, i, j, &curr)) {
               fclose(log);
               return false;
            }
            if (j == 0 && size > 0)
               pad[size/2] = '\0';

            if(curr)
               fprintf(log, "%s %d-%12d ", pad, k, curr->I[k]);
            else
               fprintf(log, "%s XXXXXXXXXXXXXX ", pad);

            if (j == 0 && size > 0)
               pad[size/2] = ' ';

            if (j == cnt - 1 && size > 0) {
               pad[size/2] = '\0';
               fputs(pad, log);
               pad[size/2] = ' ';
            }
         }
      fputs("\n", log);
      }

      for (j = 0; j < cnt; ++j) {
         if (!trace(Start, i, j, &curr)) {
            fclose(log);
            return false;
         }
         if (j == 0 && size > 0)
            pad[size/2] = '\0';

         if(curr)
            fprintf(log, "%s %14p ", pad, curr->Parent);
         else
            fprintf(log, "%s XXXXXXXXXXXXXX ", pad);

         if (j == 0 && size > 0)
            pad[size/2] = ' ';

         if (j == cnt - 1 && size > 0) {
            pad[size/2] = '\0';
            fputs(pad, log);
            pad[size/2] = ' ';
         }
      }
      fputs("\n", log);

      for (j = 0; j < cnt; ++j) {
         if (!trace(Start, i, j, &curr)) {
            fclose(log);
            return false;
         }
         if (j == 0 && size > 0)
            pad[size/2] = '\0';

         if(curr)
            fprintf(log, "%s %14p ", pad, curr->Child);
         else
            fprintf(log, "%s XXXXXXXXXXXXXX ", pad);

         if (j == 0 && size > 0)
            pad[size/2] = ' ';

         if (j == cnt - 1 && size > 0) {
            pad[size/2] = '\0';
            fputs(pad, log);
            pad[size/2] = ' ';
         }
      }
      fputs("\n", log);

      for (j = 0; j < cnt; ++j) {
         if (!trace(Start, i, j, &curr)) {
            fclose(log);
            return false;
         }
         if (j == 0 && size > 0)
            pad[size/2] = '\0';

         if(curr)
            fprintf(log, "%s %14p ", pad, curr->Tuple);
         else
            fprintf(log, "%s XXXXXXXXXXXXXX ", pad);

         if (j == 0 && size > 0)
            pad[size/2] = ' ';

         if (j == cnt - 1 && size > 0) {
            pad[size/2] = '\0';
            fputs(pad, log);
            pad[size/2] = ' ';
         }
      }
      fputs("\n\n", log);
   }

   free(pad);
   fclose(log);
   return true;
}
#endif /* RTREE_DEBUG */
