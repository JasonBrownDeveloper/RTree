I needed an efficient way to index spatial data in a small application I was writing.  R-Trees described by Antonin Guttman (http://www-db.deis.unibo.it/courses/SI-LS/papers/Gut84.pdf) fit the requirements.  Unfortunately, I could not find a library or implementation not tied to a database.

Here is my implementation in C.  I used the same names and notations as described in the paper.  The comments describe what rule in the paper is currently being implemented.  I deviated from the paper slightly in the CondenseTree method by merging under-full nodes into siblings instead of re-inserting leaf nodes.  I've found that this performs faster and doesn't cause too much fragmentation.

I also added a bulk loading feature when creating a new tree, but it is not optimized and may create a very fractured tree. 
