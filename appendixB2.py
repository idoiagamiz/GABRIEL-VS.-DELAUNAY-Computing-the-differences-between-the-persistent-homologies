#R^3 Here I do NOT look if the triangles are Gabriel! I just take from the DELAUNAY triangles the ones that has all  edges gabriel.
'''------------------R^3---------------------'''
import numpy as np
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import itertools
from scipy.spatial.distance import cdist, pdist
from itertools import combinations
'''For the Gabriel extension we change line 68'''

#Coordinates of the points
S = np.random.rand(10, 3)
'''S=np.array([[0.33291177, 0.55892937, 0.47344252],
 [0.46523839, 0.949858,   0.66214016],
 [0.51876311, 0.92521649, 0.09629331],
 [0.9883815,  0.14287922, 0.15805861]])'''
print('S:', S)

'''--------------------Del(P)-----------------------'''
#The 3-simplices(the tetrahedrons)
Del = Delaunay(S)
tetrahedrons=Del.simplices
print('Tetrahedrons Delaunay:',tetrahedrons)

#The 1-simplices(the edges)
lines=set()
for simplex in tetrahedrons:
    for subset in itertools.combinations(simplex,2):
        if subset not in lines:
            lines.add(subset)
print('lines:',lines)

#Plot the Delaunay tetrahedralization
ax = plt.axes(projection='3d')
ax.scatter(S[:, 0], S[:, 1], S[:, 2], c='black')
for line in lines:
    plt.plot([S[line[0], 0], S[line[1], 0]], [S[line[0], 1], S[line[1], 1]], [S[line[0], 2], S[line[1], 2]], c='black')
plt.show()


#the triangles
triangles=set()
for simplex in tetrahedrons:
    for subset in itertools.combinations(simplex, 3):
        if subset not in triangles:
            triangles.add(subset)
print('triangles:',triangles)


'''--------------GABRIEL 2-simplices, Del2(P)--------------------------------'''
#We know that they are some of those lines and some of those triangles
def is_gabriel(line, P):
    midpoint = 0.5 * (P[line[0]] + P[line[1]])  #we calculate te midpoint of the line
    radius = 0.5 * np.sqrt(np.sum((P[line[0]] - P[line[1]])**2))
    num_inner_points=0
    for point in P:
        if (point==P[line[0]]).all() or (point==P[line[1]]).all():
            continue
        else:
            if cdist([midpoint],[point])<=radius:
                num_inner_points=num_inner_points+1
    return num_inner_points

#Gabriel 1-simplices
gabrielext_lines = set()
for line in lines:
    if is_gabriel(line, S)<=0:  #THIS IS WHAT WE CHANGE: <=1,2,3,4,5
        gabrielext_lines.add(line)
print('gabriel_lines:', gabrielext_lines)


#We use this function to decide if the Delaunay triangle is inside K2 or not
def is_gabrieltriext(triangle, P):
    sum=0
    for edge in itertools.combinations(triangle, 2):
        if edge in gabrielext_lines:
            sum=sum+1
    if sum==3:
        return 1
    else:
        return 0

gabrielext_triangles = set()
for triangle in triangles:
    if is_gabrieltriext(triangle, S):
        gabrielext_triangles.add(triangle)
print('Gabrielext triangles:', gabrielext_triangles)

plt.figure()
#Plot of the gabriel 1 and 2-simplices
ax = plt.axes(projection='3d')
ax.scatter(S[:, 0], S[:, 1], S[:, 2], c='black') #plot the points
for line in gabrielext_lines:
    plt.plot([S[line[0], 0], S[line[1], 0]], [S[line[0], 1], S[line[1], 1]], [S[line[0], 2], S[line[1], 2]], c='black')
for triangle in gabrielext_triangles:
    plt.plot([S[triangle[0], 0], S[triangle[1], 0],S[triangle[2], 0]], [S[triangle[0], 1], S[triangle[1], 1],S[triangle[2], 1]], [S[triangle[0], 2], S[triangle[1], 2],S[triangle[2], 2]], c='blue')
plt.show()

'''--------------PERSISTENT HOMOLOGY of Del and Del2--------------------'''
#K is a simplicial complex: tetrahaedrons, triangles, lines and points
K =[]
for face in range(len(S)): #the vertices (0-simplex)
    K.append(frozenset((face,)))
for face in lines: #the lines (1-simplex)
    K.append(frozenset(face))
for face in triangles: #the three points(the triangles) (2-simplex)
    K.append(frozenset(face))
for face in tetrahedrons: #the tetrahedrons (3-simplex)
    K.append(frozenset(face))

#Simplicial complex for Gabriel, K2: tetrahedrons and triangles that have Gabriel edges, Gabriel edges and points
def check_gabriel_edges(face):
    i = 0
    for facet in itertools.combinations(face, len(face)-1):
        if frozenset(facet) in K2:
            i = i + 1
    return i==len(face)

K2 =[]
for face in range(len(S)): #the vertices (0-simplex)
    K2.append(frozenset((face,)))
for face in gabrielext_lines: #the gabriel lines (1-simplex)
    K2.append(frozenset(face))
for face in gabrielext_triangles:
    K2.append(frozenset(face))
for face in tetrahedrons:
    if check_gabriel_edges(face):
        K2.append(frozenset(face))

"-----FILTRATION-----"
def filt_value_fct(face):
    if len(face) == 1:
        return 0
    return np.max(pdist(S[list(face)]))
#K
filtration_values = []
for face in K:
    filtration_values.append(filt_value_fct(face))
order = np.argsort(filtration_values,kind='mergesort') #Returns the indices that would sort an array.

K = np.array(K)
filtration_values = np.array(filtration_values)
filtration_values = filtration_values[order]
K = K[order]  #we reorder K
face_index_dictionary = {face: index for index, face in enumerate(K)}

#K2
filtration_values2 = []
for face in K2:
    filtration_values2.append(filt_value_fct(face))
order2 = np.argsort(filtration_values2, kind='mergesort')

K2 = np.array(K2)
filtration_values2 = np.array(filtration_values2)
filtration_values2 = filtration_values2[order2]
K2 = K2[order2]  #we reorder K
face_index_dictionary2 = {face: index for index, face in enumerate(K2)}

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
def persistent_homology(D): #D is the partial and it is a matrix, with numpy
    R=D
    low=[]
    n=len(D)
    for j in range(0, n):
        count=0
        for i in range(n-1,-1,-1):
            if D[i][j]==1:
                low.append(i)
                count=1
                break
            else:
                continue
        if count==0:#all the rows are 0
            low.append('undefined')

        for j0 in range(len(low)-1): #len(low) is j+1 now
            if low[j0]==low[j] and low[j]!='undefined':
                #add column j0 to column j
                for i in range(n):
                    R[i][j]=(R[i][j0]+R[i][j])% 2
                persistent_homology(R)
    #calculate the low at the end:
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
        if count==0:#all the rows are 0
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

print('K2 ordered:',K2)
partial2=boundary(K2,face_index_dictionary2)
print('boundary of Del2: \n',partial2)
persistentDel2=persistent_homology(partial2)[0]
print('persistent_homology(partial1)of DEl2: \n',persistentDel2)

#THE KEYS OF THE DICTIONARY ARE THE BIRTHS AND THE VALUES THE DIES.
p0_homologyDel=pers0_homology(K,partial)
print('p0_homologyDel:',p0_homologyDel)
p0_homologyDel2=pers0_homology(K2,partial2)
print('p0_homologyDel1:',p0_homologyDel2)

p1_homologyDel=pers1_homology(K,partial)
print('p1_homologyDel:',p1_homologyDel)
p1_homologyDel2=pers1_homology(K2,partial2)
print('p1_homologyDel1:',p1_homologyDel2)

