# GABRIEL VS. DELAUNAY: Computing the differences between the persistent homologies
Persistent homology is a very useful tool in Topological  Data Analysis. In the last years, the Delaunay triangulation has been used to get the persistence diagram, but what if we could find a more efficient manner of getting the same persistent homology?

In this project we will explain the existence of two sub-complexes of the Delaunay triangulation that can help us to answer this problem: the alpha complexes and the Gabriel complexes.

The alpha complex is the most widely used method to compute persistent homology for large low-dimensional data sets. It is more common than the Gabriel simplices and it has been studied for more time. In fact, it has been used a lot because it works.  

Gabriel complex has some problems, but we want to discover in which way can give us information that we can use, even if it does not totally work. That is why we are specially interested in this type of complex, a field to investigate that could give us great benefits. 

In dimension two, suppose that S is a set of points in the plane and Del(S) the corresponding Delaunay triangulation. Let K be the simplicial complex consisting on those points and the edges and triangles of Del(S).
Let K1 be another simplicial complex, which consists just on those simplices with Gabriel edges. By an isomorphism we will see how both of this simplicial complexes end up giving us the same 0-persistent homology. The entire code is available in ApppendixA.

In dimension three, there is not an isomorphism between this two 1-homologies, so first we will compute an experiment that is going to fail. After that, we will try to define an extended version of the Gabriel 1- and 2-simplices such that the isomorphism is fulfilled.  This process is based on  try and error until we reach the desired isomorphism. When deciding which of the Delaunay 1-simplices are also Gabriel we usually check that the number of points of S that are closer to the midpoint of the edge pq than p and q is zero. Now we will allow this number of points to be less or equal to 1,2,3... until we get what we want. We do the same with the number of points that are in the cicumsphere of the triangles. It takes sense because in each step we are getting more and more closer to the real Delaunay simplices. The entire code is available in AppendixB2. 

