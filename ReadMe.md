# Generalized HSP (GHSP) Search
The GHSP Search algorithm uses ball-partioning-based index to accelerate the search for the HSP neighbors of a query 
point. The main idea is to allow pivots, the keys of the index, to account for groups of points within their pivot 
domain (the ball centered at the pivot). The two main processes of the HSP algorithm are adapted to use a list of pivots
rather than a list of points:

1. Active Nearest Neighbor Search: finding the next HSP neighbor within the list of active points
2. Active Point Validation: testing all active points against the HSP inequalities between the query and the new HSP 
neighbor. 


#### Two-Layer GHSP Search
The two-layer GHSP search algorithm uses a simple, 2-layer index of pivots. This can be thought of as a list of pivots, 
with each pivot containing a pointer to its pivot domain. The algorithm proceeds with two lists: a list of active 
pivots, and a list of active points. The list of active pivots are able to rule out groups of points by their pivot 
through the GHSP inequalities, while the list of active points are reserved for those points that can no longer be ruled 
out by their pivot, and must be tested against the HSP inequalities. 

#### Recursive GHSP Search
The effectiveness of the 2-Layer GHSP is limited by this trade-off on radius size: a radius too small and there are too
many pivots (less savings), a radius too large and those domains closeby that must be examined by brute force will have
a lot of points in them. The natural solution is to extend this to multiple layers: The top layer has coarse radius 
pivots, and those pivot domains contain fine-radius pivots of the second layer, and those pivot domains contain groups 
of points. This GHSP algorithm could contain three separate lists of active pivots/points, one for each layer. 

However, consider the reason for the name Generalized HSP: when $r=0$ in the GHSP inequalities, we have the HSP 
inequalities! We are performing the same exact operations on each pivot/point in every list: thus, it is possible to use
a single list $A$ that contains all active pivots/points. Each member of the list will store its radius. The active 
pivots will contain pointers to their pivot domain which can be added to the list when necessary, and the active points 
will not contain any pointers.

For this, a new struct called `Pivot` was created in order to store all information about a single pivot. This includes
the index (within the dataset), the radius, the set of children, and the max distance to a child (smaller radius). The
pivot index thus creates a recursive list of pivots, and each pivot contains a pointer to a list of pivots below it. 