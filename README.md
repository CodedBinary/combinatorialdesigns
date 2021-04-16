# combinatorialdesigns

A script to generate combinatorial designs. Currently working on those based on affine and projective geometries and difference systems. Can also check for various simple properties of combinatorial designs with some small work.

Uses Sagemath. At the moment I'm running it by running `sage` and then `load("design.py")`

# What is a design?

There is one particular kind of design called a balanced incomplete block design. A balanced incomplete block design is formed from a set of points V and a collection of k-subsets of V, called B, with the property that every two points in V occur in the same number of elements of B. 

For instance, if V={0,1,2,3,4,5,6} and B is given by
```
{{1,2,4},
{2,3,5},
{3,4,6},
{4,5,0},
{5,6,1},
{6,0,2},
{0,1,3}}
```
Then we have any two elements of V (for instance, 2 and 5) occurring in exactly one block ({2,3,5}).
