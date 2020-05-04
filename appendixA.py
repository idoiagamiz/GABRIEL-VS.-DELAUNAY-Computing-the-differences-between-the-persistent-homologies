'''------------------R^2---------------------'''
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
import itertools
from scipy.spatial.distance import cdist, pdist
from itertools import combinations


#Coordinates of the points
S= np.random.rand(100, 2)
#S=np.array([[0.03251924, 0.26965579], [0.52825821, 0.99552744], [0.83619923, 0.41257578], [0.94820483, 0.38014814]])
print('S:', S)

'''--------------------Del(P)-----------------------'''
#The 2-simplices(the triangles)
Del = Delaunay(S)
simplices=Del.simplices
print('Triangles Delaunay:',simplices)

#Plot the triangulation
plt.triplot(S[:, 0], S[:, 1], simplices,c='black')
plt.plot(S[:, 0], S[:, 1], 'o')
plt.gca().set_aspect('equal', adjustable='box')
plt.show()

#The 1-simplices(the edges)
lines = set()
for simplex in simplices: #a simplex is something like this [9 0 8], each number represent a point of the triangulation
    for i in range(len(simplex)): #for each point
        for j in range(i+1, len(simplex)):
            if (simplex[i], simplex[j]) and (simplex[j], simplex[i]) not in lines:
                lines.add((simplex[i], simplex[j]))


'''--------------GABRIEL 1-simplices, Del1(P)--------------------------------'''
#We know that they are some of those lines
def is_gabriel(line, P, epsilon=1e-6):
    midpoint = 0.5 * (P[line[0]] + P[line[1]]) #midpoint of the edge
    radius = 0.5 * np.sqrt(np.sum((P[line[0]] - P[line[1]])**2))
    return np.sum(cdist([midpoint], P) < radius - epsilon) == 0
#A point is gabriel if the midpoint is not closer to any other element of P than P[line[0]] and P[line[1]]
#that means that the distance can not be smaller than the radius

#Gabriel 1-simplices
gabriel_lines = set()
for line in lines:
    if is_gabriel(line, S):
        gabriel_lines.add(line)


#Plot of the gabriel 1-simplices
plt.scatter(S[:, 0], S[:, 1])
for line in gabriel_lines:
    plt.plot([S[line[0], 0], S[line[1], 0]], [S[line[0], 1], S[line[1], 1]],c='black')
    plt.gca().set_aspect('equal', adjustable='box')
plt.show()



'''--------------PERSISTENT HOMOLOGY of Del and Del1--------------------'''
#Simplicial complex for Delaunay, K: triangles, lines and points
K =[]
for face in range(len(S)): #the vertices (0-simplex)
    K.append(frozenset((face,)))
for face in lines: #the lines (1-simplex)
    K.append(frozenset(face))
for face in simplices: #the triangles (2-simplex)
    K.append(frozenset(face))

#Simplicial complex for Gabriel, K1:triangles that have Gabriel edges, Gabriel edges and points
def check_gabriel_edges(face):
    i = 0
    for facet in itertools.combinations(face, len(face)-1):
        if frozenset(facet) in K1:
            i = i + 1
    return i==len(face)

K1 = []
for face in range(len(S)):  # the vertices (0-simplex)
    K1.append(frozenset((face,)))
for face in gabriel_lines:  # the gabriel lines (1-simplex)
    K1.append(frozenset(face))
for face in simplices:
    if check_gabriel_edges(face):
        K1.append(frozenset(face))


"-----FILTRATION-----"
def filt_value_fct(face):
    if len(face) == 1:
        return 0
    return np.max(pdist(S[list(face)]))

#K
filtration_values = []
for face in K:
    filtration_values.append(filt_value_fct(face))
order = np.argsort(filtration_values,kind='mergesort')

K = np.array(K)
filtration_values = np.array(filtration_values)
filtration_values = filtration_values[order]
K = K[order] #We reorder K
face_index_dictionary = {face: index for index, face in enumerate(K)}


#K1
filtration_values1 = []
for face in K1:
    filtration_values1.append(filt_value_fct(face))
order1 = np.argsort(filtration_values1,kind='mergesort')

K1 = np.array(K1)
filtration_values1 = np.array(filtration_values1)
filtration_values1 = filtration_values1[order1]
K1 = K1[order1] #We reorder K1
face_index_dictionary1 = {face: index for index, face in enumerate(K1)}

'''-----BOUNDARY MATRIX------'''
def boundary(K,face_index_dictionary):
    n=len(K)
    partial = np.zeros((n, n))
    #Rows: elements of K in order, columns:elements of K in order
    for j in range(0, n):
        if len(K[j]) == 1:  #This columns are always 0, a point doesn't have boundary
            continue
        else:  # If it is not a point we calculate the boundary
            boundary = list(face_index_dictionary[frozenset(facet)] for facet in combinations(K[j], len(K[j]) - 1))
            for i in range(0, n):
                if face_index_dictionary[K[i]] in boundary:
                    partial[i][j] = 1
                else:
                    continue
    return partial

'''---PERSISTENT HOMOLOGY MATRIX---'''
def persistent_homology(D): #D is the partial, the boundary matrix
    R=D
    low=[]
    n=len(D)
    for j in range(0, n):
        #First we define the low
        count=0
        for i in range(n-1,-1,-1):
            if D[i][j]==1:
                low.append(i)
                count=1
                break
            else:
                continue
        if count==0: #All the rows are 0
            low.append('undefined')

        for j0 in range(len(low)-1):
            if low[j0]==low[j] and low[j]!='undefined':
                for i in range(n):
                    R[i][j]=(R[i][j0]+R[i][j])% 2
                persistent_homology(R)

    #Recalculate the low at the end:
    low=[]
    for j in range(0, n):
        count=0
        for i in range(n-1,-1,-1):
            if D[i][j]==1:
                low.append(i)
                count=1
                break
            else:
                continue
        if count==0:
            low.append('undefined')

    return R,low

'''------PERSISTENT k-HOMOLOGY-----'''
def pers0_homology(K,partial):
    low=persistent_homology(partial)[1]
    pers0_homology=dict()
    for i in range(len(S)):
        if i in low:
            x = low.index(i)
            pers0_homology[i] = K[x]
        else:
            pers0_homology[i] = None
    return pers0_homology

def pers1_homology(K,partial):
    low = persistent_homology(partial)[1]
    pers1_homology = dict()
    for element in K:
        if len(element) == 2:  # the births are the edges
            i = np.where(K == element)[0][0]
            if i in low:  # not all the edges are a birth
                x = low.index(i)
                pers1_homology[element] = K[x]
            else:
                pers1_homology[element] = None
    return pers1_homology

'''------------RESULTS-------------------'''
print('K ordered:',K)
partial=boundary(K,face_index_dictionary)
print('Boundary of Del: \n',partial)
persistentDel=persistent_homology(partial)[0]
print('persistent_homology(partial) of Del: \n',persistentDel)

print('K1 ordered:',K1)
partial1=boundary(K1,face_index_dictionary1)
print('boundary of Del1: \n',partial1)
persistentDel1=persistent_homology(partial1)[0]
print('persistent_homology(partial1)of DEl1: \n',persistentDel1)

#THE KEYS OF THE DICTIONARY ARE THE BIRTHS AND THE VALUES THE DIES.
p0_homologyDel=pers0_homology(K,partial)
print('p0_homologyDel:',p0_homologyDel)
p0_homologyDel1=pers0_homology(K1,partial1)
print('p0_homologyDel1:',p0_homologyDel1)

p1_homologyDel=pers1_homology(K,partial)
print('p1_homologyDel:',p1_homologyDel)
p1_homologyDel1=pers1_homology(K1,partial1)
print('p1_homologyDel1:',p1_homologyDel1)




