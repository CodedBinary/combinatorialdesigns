# combinatorialdesigns

A script to generate combinatorial designs. Currently working on those based on affine and projective geometries and difference systems. Can also check for various simple properties of combinatorial designs with some small work.

Uses Sagemath. At the moment I'm running it by executing `sage design.sage`

# Features
 - Calculation of parameters for designs
 - Computation of designs of projective and affine geometries
 - Basic existence checking, including obvious necessary conditions, attempts at using the Bruck Ryser Chowler theorem
 - Basic checking for ability to compute a design from Hadamard matrix, affine or projective geometry
 - Computation of derived, residual, complement designs
 - Computation of cyclic designs given a set of sets to act on
 - Generation of difference set via quadratic residues

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

There are lots of ways to calculate designs, whether it is from a finite geometry, a different design, or a difference set.

# Examples

Suppose we want to look at the d flats of the affine geometry AG(n,q) where n=3, q=2, d=2.
```
sage: x = AG(3,2,2)
sage: x.parameters
(8, 4, 3)
```

At this point, the blocks of the design haven't been calculated yet. This lets us calculate parameters of other designs without having to wait for a useless calculation. Suppose we want to see the blocks.

```
sage: x.generate()
sage: x.list_blocks()
[[(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)],
 [(0, 0, 1), (1, 0, 1), (0, 1, 1), (1, 1, 1)],
 [(0, 0, 0), (1, 0, 0), (0, 1, 1), (1, 1, 1)],
 [(0, 1, 0), (1, 1, 0), (0, 0, 1), (1, 0, 1)],
 [(0, 0, 0), (1, 0, 1), (0, 1, 0), (1, 1, 1)],
 [(1, 0, 0), (0, 0, 1), (1, 1, 0), (0, 1, 1)],
 [(0, 0, 0), (1, 0, 1), (0, 1, 1), (1, 1, 0)],
 [(1, 0, 0), (0, 0, 1), (1, 1, 1), (0, 1, 0)],
 [(0, 0, 0), (1, 0, 0), (0, 0, 1), (1, 0, 1)],
 [(0, 1, 0), (1, 1, 0), (0, 1, 1), (1, 1, 1)],
 [(0, 0, 0), (1, 1, 0), (0, 0, 1), (1, 1, 1)],
 [(1, 0, 0), (0, 1, 0), (1, 0, 1), (0, 1, 1)],
 [(0, 0, 0), (0, 1, 0), (0, 0, 1), (0, 1, 1)],
 [(1, 0, 0), (1, 1, 0), (1, 0, 1), (1, 1, 1)]]
```

But at this point, the blocks and points are given as points and subspaces of a vector space. That isn't nice for design theory. Affine geometries have an easy way to turn these into integers, though! If we had a different design, we can call biject_obvious(), which always works.
```
sage: x.biject()
sage: x.list_blocks()
[[1, 2, 3, 4],
 [5, 6, 7, 8],
 [1, 2, 7, 8],
 [3, 4, 5, 6],
 [1, 6, 3, 8],
 [2, 5, 4, 7],
 [1, 6, 7, 4],
 [2, 5, 8, 3],
 [1, 2, 5, 6],
 [3, 4, 7, 8],
 [1, 4, 5, 8],
 [2, 3, 6, 7],
 [1, 3, 5, 7],
 [2, 4, 6, 8]]
sage: x.list_points()
[1, 2, 3, 4, 5, 6, 7, 8]
sage: x.list_resolution_classes()
[[[1, 2, 3, 4], [5, 6, 7, 8]],
 [[1, 2, 7, 8], [3, 4, 5, 6]],
 [[1, 6, 3, 8], [2, 5, 4, 7]],
 [[1, 6, 7, 4], [2, 5, 8, 3]],
 [[1, 2, 5, 6], [3, 4, 7, 8]],
 [[1, 4, 5, 8], [2, 3, 6, 7]],
 [[1, 3, 5, 7], [2, 4, 6, 8]]]
```
Notice that the point set and resolution classes updated themselves. Let's explore that a bit.
```
sage: x  = AG(3,2,2)
sage: x.generate()
sage: x.list_points()
[(0, 0, 0),
 (1, 0, 0),
 (0, 1, 0),
 (1, 1, 0),
 (0, 0, 1),
 (1, 0, 1),
 (0, 1, 1),
 (1, 1, 1)]
sage: y = x.complement()
sage: y.generate()
sage: y.list_points()
[(0, 0, 0),
 (1, 0, 0),
 (0, 1, 0),
 (1, 1, 0),
 (0, 0, 1),
 (1, 0, 1),
 (0, 1, 1),
 (1, 1, 1)]
sage: x.biject()
sage: y.list_points()
[1, 2, 3, 4, 5, 6, 7, 8]
sage: x.list_points()
[1, 2, 3, 4, 5, 6, 7, 8]
```

So they automatically update. If this behaviour is undesired, we have a method to fix that.

```
sage: x = AG(3,2,2)
sage: x.generate()
sage: x.list_points()
[(0, 0, 0),
 (1, 0, 0),
 (0, 1, 0),
 (1, 1, 0),
 (0, 0, 1),
 (1, 0, 1),
 (0, 1, 1),
 (1, 1, 1)]
sage: y = x.complement()
sage: y.generate()
sage: y.decouple_points()
sage: x.biject()
sage: y.list_points()
[(0, 0, 0),
 (1, 0, 0),
 (0, 1, 0),
 (1, 1, 0),
 (0, 0, 1),
 (1, 0, 1),
 (0, 1, 1),
 (1, 1, 1)]
sage: y.list_blocks()
[[(0, 0, 1), (1, 0, 1), (0, 1, 1), (1, 1, 1)],
 [(0, 0, 0), (1, 0, 0), (0, 1, 0), (1, 1, 0)],
 [(0, 1, 0), (1, 1, 0), (0, 0, 1), (1, 0, 1)],
 [(0, 0, 0), (1, 0, 0), (0, 1, 1), (1, 1, 1)],
 [(1, 0, 0), (1, 1, 0), (0, 0, 1), (0, 1, 1)],
 [(0, 0, 0), (0, 1, 0), (1, 0, 1), (1, 1, 1)],
 [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1)],
 [(0, 0, 0), (1, 1, 0), (1, 0, 1), (0, 1, 1)],
 [(0, 1, 0), (1, 1, 0), (0, 1, 1), (1, 1, 1)],
 [(0, 0, 0), (1, 0, 0), (0, 0, 1), (1, 0, 1)],
 [(1, 0, 0), (0, 1, 0), (1, 0, 1), (0, 1, 1)],
 [(0, 0, 0), (1, 1, 0), (0, 0, 1), (1, 1, 1)],
 [(1, 0, 0), (1, 1, 0), (1, 0, 1), (1, 1, 1)],
 [(0, 0, 0), (0, 1, 0), (0, 0, 1), (0, 1, 1)]]
```

As you can see, we bijected x, but we decoupled y first. As a result, the points of y are what they were at the time of decoupling.

We can create designs from designs, though. Let's look at the complement.

```
sage: y = x.complement()
sage: x.parameters
(8, 4, 3)
sage: y.parameters
(8, 4, 3)
sage: y.generate()
sage: y.list_blocks()
[[5, 6, 7, 8],
 [1, 2, 3, 4],
 [3, 4, 5, 6],
 [1, 2, 7, 8],
 [2, 4, 5, 7],
 [1, 3, 6, 8],
 [2, 3, 5, 8],
 [1, 4, 6, 7],
 [3, 4, 7, 8],
 [1, 2, 5, 6],
 [2, 3, 6, 7],
 [1, 4, 5, 8],
 [2, 4, 6, 8],
 [1, 3, 5, 7]]
```

Perhaps we want to know if a design exists at all.

```
sage: x = BIBD([53,13,3])
sage: x.existence()
(False, 'common divisibility argument mod 5')
sage: x = BIBD([7,3,1])
sage: x.existence()
(True, 'Hadamard design of order 2')
```

At the moment, we are working on a way to check if the complement exists smoothly.
