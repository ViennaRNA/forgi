import forgi.threedee.utilities.vector as ftuv
import collections as col
import numpy as np
import itertools as it
import networkx as nx

class Projection2D(object):
    """A 2D Projection of a CoarseGrainRNA unto a 2D-plain"""
    def __init__(self, cg, proj_direction):
        """
        @param cg:  an CoarseGrainRNA object with 3D coordinates for every element
                    .. note:: The projection is generated from this cg, but it is not associated with it after construction.
                              Thus future changes of the cg are not reflected in the projection.
        @param proj_direction: a vector (in 3D space) in the direction of projection. 
                               The length of this vector is not used.
        
        """        
        #: The projected coordinates of all stems
        self._coords=dict()

        self._cross_points=None
        self._proj_graph=None

        #Calculate orthonormal basis of projection plane.
        proj_direction=np.array(proj_direction)
        _, unit_vec1, unit_vec2=ftuv.create_orthonormal_basis(proj_direction)
        self._unit_vec1=unit_vec1
        self._unit_vec2=unit_vec2
        self._proj_direction=proj_direction
        self._project(cg)


  
    @property
    def proj_direction(self):
        """A vector describing the direction of the projection"""
        return self._proj_direction
    """
    @proj_direction.setter
    def proj_direction(self, val):
        '''Set a new projection direction and recalculate the projected points'''
        self._cross_points=None
        self._proj_graph=None
        self._proj_direction=proj_direction
        self._project()
    """
    @property
    def crossingPoints(self):
        """
        All points, where 2 segments intersect
        
        A list of triples (key1, key2, coordinate), where coordinate is a 2D vector.
        """
        if self._cross_points is None:
            self._cross_points=col.defaultdict(list)
            for i, key1 in enumerate(self._coords):
                for j, key2 in enumerate(self._coords):
                    if j>i: #_coords is not modified during the iteration, thus the order is the same for i and j.
#                        import matplotlib.pyplot as plt
#                        plt.plot(self._coords[key1][0][0], self._coords[key1][0][1], 'bo')
#                        plt.plot(self._coords[key1][1][0], self._coords[key1][1][1], 'bo')
#                        plt.plot(self._coords[key2][0][0], self._coords[key2][0][1], 'ro')
#                        plt.plot(self._coords[key2][1][0], self._coords[key2][1][1], 'ro')
#                        plt.axis(self.get_bounding_square(0))
                        for cr in ftuv.seg_intersect(self._coords[key1], self._coords[key2]):
#                            plt.plot(cr[0], cr[1], 'go')
                            self._cross_points[key1].append((cr, key2))
                            self._cross_points[key2].append((cr, key1))
#                        plt.text(self._coords[key1][0][0], self._coords[key1][0][1], key1)
#                       plt.text(self._coords[key2][0][0], self._coords[key2][0][1], key2)        
#                        plt.show()  
        return self._cross_points
    @property
    def proj_graph(self):
        """A graph describing the projected RNA.

        This graph is stored as a networkx graph object.
        """
        if self._proj_graph is None:
            self._build_proj_graph()
        return self._proj_graph
    @property
    def convex_hull(self):
          points=[ x for x in a for a in self._coords.values()]
          return scipy.spatial.ConvexHull(points)        
    '''def _build_proj_graph(self):
        """
        Generate a graph from the 2D projection. 
  
        This is implemented as a dictionary {coordinate1:coordinate2 for all edges coordinate1-coordinate2}
        As edges are undirected, each edge is present twice (in both directions) in this dictionary.
        """
        proj_graph=col.defaultdict(set)            
        for key, element in self._coords.items():
            crs=self.crossingPoints
            sortedCrs=ftuv.sortAlongLine(element[0], element[1], [x[0] for x in crs[key]])     
            oldpoint=None
            for point in sortedCrs:
                point=(point[0], point[1])
                if oldpoint is not None:                    
                    proj_graph[oldpoint].add(point)
                    proj_graph[point].add(oldpoint)
                oldpoint=point
        self._proj_graph=proj_graph'''
    def _build_proj_graph(self):
        """
        Generate a graph from the 2D projection. 
  
        This is implemented as a networkx.Graph with the coordinates as nodes.
        """
        proj_graph=nx.Graph()
        for key, element in self._coords.items():
            crs=self.crossingPoints
            sortedCrs=ftuv.sortAlongLine(element[0], element[1], [x[0] for x in crs[key]])     
            oldpoint=None
            for point in sortedCrs:
                point=(point[0], point[1]) #Tuple, to be hashable
                if oldpoint is not None:
                    proj_graph.add_edge(oldpoint, point)                 
                oldpoint=point
        self._proj_graph=proj_graph
    def _project(self, cg):
        """
        Calculates the 2D coordinates of all coarse grained elements by vector rejection.
        Stores them inside self._coords
        """
        self._coords=dict()
        #Project all coordinates to this plane
        for key, val in cg.coords.items():
            start=np.array([np.dot(self._unit_vec1,val[0]), np.dot(self._unit_vec2,val[0])])
            end=np.array([np.dot(self._unit_vec1,val[1]), np.dot(self._unit_vec2,val[1])])
            self._coords[key]=(start, end)
    def get_bounding_box(self, margin=0.):
        """
        Returns the coordinates for a box that contains all points of the 2D projection.

        :param margin: increase the bounding box in every direction by this margin.
        :returns: left, right, bottom, top
        """
        points=list(p for x in self._coords.values() for p in x)
        #print("P", points)
        left=min(x[0] for x in points)-margin
        right=max(x[0] for x in points)+margin
        bottom=min(x[1] for x in points)-margin
        top=max(x[1] for x in points)+margin
        #print "BB",  left, right, bottom, top
        return left, right, bottom, top
    def get_bounding_square(self, margin=0.):
        """
        Returns the coordinates for a square that contains all points of the 2D projection.

        :param margin: increase the bounding box in every direction by this margin.
        :returns: left, right, bottom, top
        """
        bb=self.get_bounding_box()
        length=max([bb[1]-bb[0],bb[3]-bb[2]])/2+margin
        x=(bb[0]+bb[1])/2
        y=(bb[2]+bb[3])/2
        return x-length, x+length, y-length, y+length
    def plot(self):
        """Plots the 2D projection"""
        import matplotlib.pyplot as plt
        import matplotlib.lines as lines
        fig, ax = plt.subplots()
        for s,e in self.proj_graph.edges():
            line=lines.Line2D([s[0], e[0]],[s[1],e[1]])                
            ax.add_line(line)
        for key,p in self._coords.items():
            plt.plot(p[0][0], p[0][1], 'ro')
            plt.plot(p[1][0], p[1][1], 'ro')
            ax.text((p[0][0]+p[1][0])/2, (p[0][1]+p[1][1])/2, key)
        for key,l in self.crossingPoints.items():
            for p in l:
                #print("P", p)
                plt.plot(p[0][0], p[0][1], 'go')
        plt.axis(self.get_bounding_square(5))
        plt.show()


