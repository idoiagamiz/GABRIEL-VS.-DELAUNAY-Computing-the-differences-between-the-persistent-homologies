# GABRIEL VS. DELAUNAY: Computing the differences between the persistent homologies
Persistent homology is a very useful tool in Topological  Data Analysis. In the last years, the Delaunay triangulation has been used to get the persistence diagram, but what if we could find a more efficient manner of getting the same persistent homology?

In this project we will explain the existence of two sub-complexes of the Delaunay triangulation that can help us to answer this problem: the alpha complexes and the Gabriel complexes.

The alpha complex is the most widely used method to compute persistent homology for large low-dimensional data sets. It is more common than the Gabriel simplices and it has been studied for more time. In fact, it has been used a lot because it works.  

Gabriel complex has some problems, but we want to discover in which way can give us information that we can use, even if it does not totally work. That is why we are specially interested in this type of complex, a field to investigate that could give us great benefits. 
On the one hand, an isomorphism between the 0-persistent homology of both complexes will be proved on the plane. On the other hand, we will see that in three dimensions how even despite of not getting the same 1-persistent homology  there exist some information for this homologies.

Suppose that S is a set of points in the plane and Del(S) the corresponding Delaunay triangulation. Let K be the simplicial complex consisting on those points and the edges and triangles of Del(S).
Let K1 be another simplicial complex, which consists just on those simplices with Gabriel edges. By an isomorphism we will see how both of this simplicial complexes end up giving us the same 0-persistent homology.

