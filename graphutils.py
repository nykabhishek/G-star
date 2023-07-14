import sys
sys.path.append('DubinsLineSegToLineSeg/')
from shapely.geometry import Point, Polygon, MultiPolygon, LineString, MultiLineString, mapping
import numpy as np
from gstar_utils import is_route_feasible, find_slope, retrieve_gate, retrieve_minEdge, HeadingsToSectors
from DubinsLineSegToLineSeg import DubinsL2L, utils, dubutils
from shapely.ops import split, nearest_points
from copy import copy, deepcopy
import matplotlib.pyplot as plt
from types import SimpleNamespace
import networkx as nx
import dubins as du
import math
# import sys
# sys.path.append('DubinsLineSegToLineSeg/')
# from DubinsLineSegToLineSeg import DubinsLineToLine_V2, utils

# sectors = 4
step_size = 0.1
# continuity_constraint = 0.1
# angle_continuity_constraint = np.pi/4
# Tolerance = 0.2 #


class Graph_Gstar:
    """
    Class implementing networkX graph to find G* paths.

    Attributes
    ----------
    StartConf, EndConf : 3-tuple
        List of (x, y, rho) coordinates in the frame of the environment representing the start and end congfiguration of the Dubins vehicle.

    graph

    self.ObstacleList

    bounding_box : 4-tuple
        Coordinates of the lower left and upper right corners of the bounding box containing the obstacle.

    Methods
    -------
    None

    """

    def __init__(self, StartPoint, EndPoint, rho, Sectors, tolerances):

        self.GateCounter = 0
        self.graph = nx.Graph()
        self.Sectors = Sectors

        self.Start = Point(StartPoint[0], StartPoint[1])
        self.End = Point(EndPoint[0], EndPoint[1])
        self.rho = rho
        self.tolerances = tolerances

        self.gateGeomList = []

        self.LineStart = [(StartPoint[0], StartPoint[1]-0.01),
                          (StartPoint[0], StartPoint[1]+0.01)]
        self.gateGeomList.append(LineString(self.LineStart))
        # self.LineStart = DubinsL2L.LineSegment((StartPoint[0], StartPoint[1]-0.01), (StartPoint[0], StartPoint[1]+0.01))
        # self.gateGeomList.append(self.LineStart)
        self.LineEnd = [(EndPoint[0], EndPoint[1]-0.01),
                        (EndPoint[0], EndPoint[1]+0.01)]
        self.gateGeomList.append(LineString(self.LineEnd))
        # self.LineEnd = DubinsL2L.LineSegment((EndPoint[0], EndPoint[1]-0.01), (EndPoint[0], EndPoint[1]+0.01))
        # self.gateGeomList.append(self.LineEnd)

        self.graph.add_node(self.GateCounter, id='s', geom=LineString(
            self.LineStart), point=self.Start, sector=None, x_dist=0)
        StartId = self.GateCounter
        self.GateCounter += 1

        self.graph.add_node(self.GateCounter, id='e', geom=LineString(
            self.LineEnd), point=self.End, sector=None, x_dist=self.End.distance(self.Start))
        EndId = self.GateCounter
        self.GateCounter += 1

        eucPath = LineString([self.Start, self.End])
        self.graph.add_edge(
            StartId, EndId, weight=eucPath.length, geom=eucPath)

        self.eucLB_free, self.dubLB_free = None, None
        self.eucLowerPath, self.dubLowerPath = None, None
        self.eucLowerBound, self.dubLowerBound, self.dubUpperBound = None, None, None
        self.eucLB_time, self.dubLB_time, self.dubUB_time = None, None, None

        self.eucGraph = None
        self.dubUB_graph = None
        self.dubUB_nodeCounter = 0
        self.dubUB_headingList = None


def node_id_list(Id, G):
    id_list = []
    for n in list(G.graph.nodes):
        if G.graph.nodes[n]['id'] == Id:
            id_list.append(n)
    return id_list


def get_shortest_path(start_id, end_id, G):
    ShortestPath, ShortestLength = [], np.inf
    start_id_list, end_id_list = node_id_list(
        start_id, G), node_id_list(end_id, G)
    for start in start_id_list:
        for end in end_id_list:
            path = nx.shortest_path(
                G.graph, source=start, target=end, weight='weight', method='dijkstra')
            length = nx.shortest_path_length(
                G.graph, source=start, target=end, weight='weight', method='dijkstra')
            if length < ShortestLength:
                ShortestPath, ShortestLength = path, length
    return ShortestPath, ShortestLength


def is_path_exists(StartId, EndId, G):
    start_id_list, end_id_list = node_id_list(
        StartId, G), node_id_list(EndId, G)
    for s in start_id_list:
        for e in end_id_list:
            if not nx.has_path(G.graph, source=s, target=e):
                return False
    return True


def ConstructGates(G, Map, Path, graphType='euclidean'):
    ''' 1. Draw gates for obstacles intesecting the path '''
    for i in range(len(Path)-1):
        u, v = Path[i], Path[i+1]
        Edge = G.graph.get_edge_data(u, v)
        # if is_feasible_connection(Edge, ObstacleMap) == False:
        #
        for obstacle in Map.ObstacleList:
            if mapping(Edge['geom'])['type'] == 'LineString':
                # if Edge['geom'].crosses(obstacle):
                if not is_feasible_connection(Edge, obstacle, G.tolerances):
                    # nx.G.graph.remove_edge()
                    G = add_gate(G, Edge['geom'], obstacle, Map, graphType)

            elif mapping(Edge['geom'])['type'] == 'MultiLineString':
                for segment in list(Edge['geom']):
                    # if segment.crosses(obstacle):
                    if not is_feasible_connection(segment, obstacle, G.tolerances):
                        G = add_gate(G, segment, obstacle, Map, graphType)
    return G


def add_gate(G, Edge, Obstacle, Map, Type):
    ''' List of points intersecting the obstacle'''
    points = list(Edge.intersection(Obstacle).boundary)
    ''' Draw gates at midway of points'''
    avg_X, avg_Y = 0, 0
    for p in points:
        avg_X = avg_X + (p.x)/2
        avg_Y = avg_Y + (p.y)/2

    gate_geom = LineString([(avg_X, Map.map_maxy), (avg_X, Map.map_miny)])

    for obs in Map.ObstacleList:
        gate_geom = gate_geom - obs

    # ExistingGateList = []
    # try:
    #     ExistingGateList = [G.graph.nodes[n]['geom'] for n in list(G.graph.nodes)]
    # except AttributeError:
    #     print("Error adding new gates")

    if Type == 'Euclidean':
        if mapping(gate_geom)['type'] == 'LineString':
            # if gate_geom not in ExistingGateList:
            if gate_geom not in G.gateGeomList:
                G.graph.add_node(G.GateCounter, id=G.GateCounter,
                                 geom=gate_geom, sector=None, x_dist=avg_X)
                G.gateGeomList.append(gate_geom)
                G.GateCounter += 1
        elif mapping(gate_geom)['type'] == 'MultiLineString':
            for gate_segment_geom in list(gate_geom):
                # if gate_segment_geom not in ExistingGateList:
                if gate_segment_geom not in G.gateGeomList:
                    G.graph.add_node(G.GateCounter, id=G.GateCounter,
                                     geom=gate_segment_geom, sector=None, x_dist=avg_X)
                    G.gateGeomList.append(gate_segment_geom)
                    G.GateCounter += 1

    elif Type == 'Dubins':
        sector_list = [(2*np.pi*i/G.Sectors, 2*np.pi*(i+1)/G.Sectors)
                       for i in range(G.Sectors)]
        if mapping(gate_geom)['type'] == 'LineString':
            # if gate_geom not in ExistingGateList:
            if gate_geom not in G.gateGeomList:
                for sector in sector_list:
                    G.graph.add_node(G.GateCounter, id=G.GateCounter,
                                     geom=gate_geom, sector=sector, x_dist=avg_X)
                    G.gateGeomList.append(gate_geom)
                    G.GateCounter += 1
        elif mapping(gate_geom)['type'] == 'MultiLineString':
            for gate_segment_geom in list(gate_geom):
                # if gate_segment_geom not in ExistingGateList:
                if gate_segment_geom not in G.gateGeomList:
                    for sector in sector_list:
                        G.graph.add_node(
                            G.GateCounter, id=G.GateCounter, geom=gate_segment_geom, sector=sector, x_dist=avg_X)
                        G.gateGeomList.append(gate_segment_geom)
                        G.GateCounter += 1

    else:
        print('TypeError: Check selected path type')
    return G


def is_feasible_path(Path, G, obstacleList):
    feasible = True
    for i in range(len(Path)-1):
        u = Path[i]
        v = Path[i+1]
        Edge = G.graph.get_edge_data(u, v)

        for obstacle in obstacleList:
            feasible = is_feasible_connection(Edge, obstacle, G.tolerances)
            if not feasible:
                return False

    return feasible


def is_feasible_connection(Edge, obstacle, Tolerances):
    feasible = True
    ''' if the path is a single line segment '''
    if mapping(Edge['geom'])['type'] == 'LineString':
        if Edge['geom'].within(obstacle) or Edge['geom'].overlaps(obstacle):
            return False
        if Edge['geom'].crosses(obstacle):
            if mapping(obstacle)['type'] == 'Polygon':
                s = Edge['geom'].intersection(obstacle).length
                if s > Tolerances['polygon_intersection']:
                    return False
            else:
                s = Edge['geom'].intersection(obstacle).length
                d = math.sqrt(4*obstacle.area/math.pi)
                if (s/d) > Tolerances['circle_intersection_ratio']:
                    return False
    elif mapping(Edge['geom'])['type'] == 'MultiLineString':
        for segment in list(Edge['geom']):
            if segment.within(obstacle) or segment.overlaps(obstacle):
                return False
            if segment.crosses(obstacle):
                if mapping(obstacle)['type'] == 'Polygon':
                    s = segment.intersection(obstacle).length
                    if s > Tolerances['polygon_intersection']:
                        return False
                else:
                    s = segment.intersection(obstacle).length
                    d = math.sqrt(4*obstacle.area/math.pi)
                    if (s/d) > (s/d) > Tolerances['circle_intersection_ratio']:
                        return False
    return feasible


def is_route_continuous(Path, G):
    for i in range(1, len(Path)-1):
        _, b = G.graph.get_edge_data(Path[i-1], Path[i])['geom'].boundary
        c, _ = G.graph.get_edge_data(Path[i], Path[i+1])['geom'].boundary
        if b.distance(c) > G.tolerances['continuity']:
            # print('points: ', b, c)
            # print('distance', b.distance(c))
            return False
    return True


def is_angle_matching(Path, G):
    for i in range(1, len(Path)-1):
        # TO DO: Add boundary point sorting based on x-distace from start point
        AngIn = G.graph.get_edge_data(Path[i-1], Path[i])['pathConf'][1][2]
        AngOut = G.graph.get_edge_data(Path[i], Path[i+1])['pathConf'][0][2]
        pathConfList = sorted([AngIn, AngOut])
        if utils.Angdiff(pathConfList[0], pathConfList[1]) > G.tolerances['angular']*np.pi/180:
            print('Angles', (pathConfList[0], pathConfList[1]))
            print('Angle diff',  utils.Angdiff(
                pathConfList[0], pathConfList[1])*180/np.pi)
            return False
    return True


def EuclideanToDubinsGraph(G, heading_restricted, startAngle, goalAngle):
    sector_list = [(2*np.pi*i/G.Sectors, 2*np.pi*(i+1)/G.Sectors)
                   for i in range(G.Sectors)]
    for n in list(G.graph.nodes):
        if G.graph.nodes[n]['id'] == 's':
            if heading_restricted:
                G.graph.add_node(G.GateCounter, id='s', point=G.graph.nodes[n]['point'], geom=G.graph.nodes[n]['geom'], sector=(
                    startAngle, startAngle+0.01), x_dist=G.graph.nodes[n]['x_dist'])
                G.GateCounter += 1
            else:
                for sector in sector_list:
                    G.graph.add_node(
                        G.GateCounter, id='s', point=G.graph.nodes[n]['point'], geom=G.graph.nodes[n]['geom'], sector=sector, x_dist=G.graph.nodes[n]['x_dist'])
                    G.GateCounter += 1
        elif G.graph.nodes[n]['id'] == 'e':
            if heading_restricted:
                G.graph.add_node(G.GateCounter, id='e', point=G.graph.nodes[n]['point'], geom=G.graph.nodes[n]['geom'], sector=(
                    goalAngle, goalAngle+0.01), x_dist=G.graph.nodes[n]['x_dist'])
                G.GateCounter += 1
            else:
                for sector in sector_list:
                    G.graph.add_node(
                        G.GateCounter, id='e', point=G.graph.nodes[n]['point'], geom=G.graph.nodes[n]['geom'], sector=sector, x_dist=G.graph.nodes[n]['x_dist'])
                    G.GateCounter += 1
        else:
            for sector in sector_list:
                G.graph.add_node(G.GateCounter, id=G.GateCounter,
                                 geom=G.graph.nodes[n]['geom'], sector=sector, x_dist=G.graph.nodes[n]['x_dist'])
                G.GateCounter += 1
        G.graph.remove_node(n)
    G.graph = G.graph.to_directed()
    return G


def ConnectGates_Dubins(G):
    """
    Function to connect Dubins edges between gates.

    Attributes
    ----------
    Graph : nx.graph()
        A graph containing Gates as nodes.

    Output
    -------
    Graph : nx.graph()
        A graph containing Gates connected using Dubins paths.

    """
    G.graph.clear_edges()
    smallest_gate_x_dist, largest_gate_x_dist = np.inf, 0
    for n in list(G.graph.nodes):
        if G.graph.nodes[n]['x_dist'] < smallest_gate_x_dist:
            smallest_gate_x_dist = G.graph.nodes[n]['x_dist']
        if G.graph.nodes[n]['x_dist'] > largest_gate_x_dist:
            largest_gate_x_dist = G.graph.nodes[n]['x_dist']

    ''' 4/1. Start path'''
    nearest_gates, nearest_gate_dist = [], np.inf
    for n in list(G.graph.nodes):
        if G.graph.nodes[n]['x_dist'] == smallest_gate_x_dist:
            nearest_gate_dist = G.graph.nodes[n]['x_dist']
            nearest_gates.append(n)

    prev_x_dist = nearest_gate_dist
    prev_gates = nearest_gates

    ''' 4/2. Gate connections '''
    while prev_x_dist < largest_gate_x_dist:
        nearest_gates, nearest_gate_dist = [], np.inf
        for n in list(G.graph.nodes):
            if G.graph.nodes[n]['x_dist'] < nearest_gate_dist and G.graph.nodes[n]['x_dist'] > prev_x_dist:
                nearest_gate_dist = G.graph.nodes[n]['x_dist']
                nearest_gates = [n]
            elif G.graph.nodes[n]['x_dist'] == nearest_gate_dist:
                nearest_gates.append(n)

        for a in prev_gates:
            for b in nearest_gates:
                line1 = list(zip(*G.graph.nodes[a]['geom'].coords.xy))
                line2 = list(zip(*G.graph.nodes[b]['geom'].coords.xy))
                sector1 = G.graph.nodes[a]['sector']
                sector2 = G.graph.nodes[b]['sector']
                # try:
                # minLength, minConfStart_1, minConfGoal_1, minPathType, minPathSegLengths = DubinsLineToLine_V2.DubinsLineToLineV2(line1, sector1, line2, sector2, G.rho)
                L2LDub_1 = DubinsL2L.Line2LineDubins(
                    line1, sector1, line2, sector2, G.rho)
                minLength_1, minPath_1 = L2LDub_1.MinDub_L2L()
                minConfStart_1, minConfGoal_1 = list(
                    minPath_1.iniPos)+[minPath_1.iniHead], list(minPath_1.finalPos)+[minPath_1.finalHead]
                path_1 = du.shortest_path(minConfStart_1, minConfGoal_1, G.rho)
                configurations_1, _ = path_1.sample_many(step_size)
                temp_path_1 = list(
                    map(lambda c: (c[0], c[1]), configurations_1))
                temp_path_1.append((minConfGoal_1[0], minConfGoal_1[1]))
                temp_path_1 = LineString(temp_path_1)
                # G.graph.add_edge(a, b, weight=temp_path_1.length, geom=temp_path_1, pathConf=[minConfStart_1, minConfGoal_1], parent_gates=[a, b])
                G.graph.add_edge(a, b, weight=minLength_1, geom=temp_path_1, pathConf=[
                                 minConfStart_1, minConfGoal_1], parent_gates=[a, b])

                # NewLine
                # minLength, minConfStart_2, minConfGoal_2, minPathType, minPathSegLengths = DubinsLineToLine_V2.DubinsLineToLineV2(line2, sector2, line1, sector1, G.rho)
                L2LDub_2 = DubinsL2L.Line2LineDubins(
                    line2, sector2, line1, sector1, G.rho)
                minLength_2, minPath_2 = L2LDub_2.MinDub_L2L()
                minConfStart_2, minConfGoal_2 = list(
                    minPath_2.iniPos)+[minPath_2.iniHead], list(minPath_2.finalPos)+[minPath_2.finalHead]
                path_2 = du.shortest_path(minConfStart_2, minConfGoal_2, G.rho)
                configurations_2, _ = path_2.sample_many(step_size)
                temp_path_2 = list(
                    map(lambda c: (c[0], c[1]), configurations_2))
                temp_path_2.append((minConfGoal_2[0], minConfGoal_2[1]))
                temp_path_2 = LineString(temp_path_2)
                # G.graph.add_edge(b, a, weight=temp_path_2.length, geom=temp_path_2, pathConf=[minConfStart_2, minConfGoal_2], parent_gates=[b, a])
                G.graph.add_edge(b, a, weight=minLength_2, geom=temp_path_2, pathConf=[
                                 minConfStart_2, minConfGoal_2], parent_gates=[b, a])

                # except:
                #     None
        prev_x_dist = nearest_gate_dist
        prev_gates = nearest_gates

    return G


def ConnectGates_Euclidean(G):
    """
    Function to connect Euclidean edges between gates.

    Attributes
    ----------
    Graph : nx.graph()
        A graph containing Gates as nodes.

    Output
    -------
    Graph : nx.graph()
        A graph containing Gates connected using Euclidean paths.

    """

    G.graph.clear_edges()
    smallest_gate_x_dist, largest_gate_x_dist = np.inf, 0
    for n in list(G.graph.nodes):
        if G.graph.nodes[n]['x_dist'] < smallest_gate_x_dist:
            smallest_gate_x_dist = G.graph.nodes[n]['x_dist']
        if G.graph.nodes[n]['x_dist'] > largest_gate_x_dist:
            largest_gate_x_dist = G.graph.nodes[n]['x_dist']

    ''' 4/1. Start path'''
    nearest_gates, nearest_gate_dist = [], np.inf
    for n in list(G.graph.nodes):
        if G.graph.nodes[n]['x_dist'] == smallest_gate_x_dist:
            nearest_gate_dist = G.graph.nodes[n]['x_dist']
            nearest_gates.append(n)

    prev_x_dist = nearest_gate_dist
    prev_gates = nearest_gates

    ''' 4/2. Gate connections '''

    while prev_x_dist < largest_gate_x_dist:
        nearest_gates, nearest_gate_dist = [], np.inf
        for n in list(G.graph.nodes):
            if G.graph.nodes[n]['x_dist'] < nearest_gate_dist and G.graph.nodes[n]['x_dist'] > prev_x_dist:
                nearest_gate_dist = G.graph.nodes[n]['x_dist']
                nearest_gates = [n]
            elif G.graph.nodes[n]['x_dist'] == nearest_gate_dist:
                nearest_gates.append(n)

        for a in prev_gates:
            for b in nearest_gates:
                near_points = nearest_points(
                    G.graph.nodes[a]['geom'], G.graph.nodes[b]['geom'])
                temp_path = LineString(near_points)
                minConfStart = [near_points[0].x,
                                near_points[0].y, find_slope(temp_path)]
                minConfGoal = [near_points[1].x,
                               near_points[1].y, find_slope(temp_path)]
                G.graph.add_edge(a, b, weight=temp_path.length, geom=temp_path, pathConf=[
                                 minConfStart, minConfGoal], parent_gates=[a, b])
        prev_x_dist = nearest_gate_dist
        prev_gates = nearest_gates

    return G


def UpdatePath_Dubins(G):
    # Graph.clear_edges()
    smallest_gate_x_dist, largest_gate_x_dist = np.inf, 0
    for n in list(G.graph.nodes):
        if G.graph.nodes[n]['x_dist'] < smallest_gate_x_dist:
            smallest_gate_x_dist = G.graph.nodes[n]['x_dist']
        if G.graph.nodes[n]['x_dist'] > largest_gate_x_dist:
            largest_gate_x_dist = G.graph.nodes[n]['x_dist']

    nearest_gates, nearest_gate_dist = [], np.inf
    for n in list(G.graph.nodes):
        if G.graph.nodes[n]['x_dist'] == smallest_gate_x_dist:
            nearest_gate_dist = G.graph.nodes[n]['x_dist']
            nearest_gates.append(n)

    prev_x_dist = nearest_gate_dist
    prev_gates = nearest_gates

    while prev_x_dist < largest_gate_x_dist:
        nearest_gates, nearest_gate_dist = [], np.inf
        for n in list(G.graph.nodes):
            if G.graph.nodes[n]['x_dist'] < nearest_gate_dist and G.graph.nodes[n]['x_dist'] > prev_x_dist:
                nearest_gate_dist = G.graph.nodes[n]['x_dist']
                nearest_gates = [n]
            elif G.graph.nodes[n]['x_dist'] == nearest_gate_dist:
                nearest_gates.append(n)

        for a in prev_gates:
            for b in nearest_gates:
                if not G.graph.has_edge(a, b):
                    line1 = list(zip(*G.graph.nodes[a]['geom'].coords.xy))
                    line2 = list(zip(*G.graph.nodes[b]['geom'].coords.xy))
                    sector1 = G.graph.nodes[a]['sector']
                    sector2 = G.graph.nodes[b]['sector']
                    try:
                        # minLength, minConfStart, minConfGoal, minPathType, minPathSegLengths = DubinsLineToLine_V2.DubinsLineToLineV2(line1, sector1, line2, sector2, G.rho)
                        # minLength, minConfStart, minConfGoal, minPathType, minPathSegLengths = DubinsL2L.Line2LineDubins(line1, sector1, line2, sector2, G.rho)
                        L2LDub = DubinsL2L.Line2LineDubins(
                            line1, sector1, line2, sector2, G.rho)
                        minLength, minPath = L2LDub.MinDub_L2L()
                        minConfStart, minConfGoal = list(
                            minPath.iniPos)+[minPath.iniHead], list(minPath.finalPos)+[minPath.finalHead]
                        path = du.shortest_path(
                            minConfStart, minConfGoal, G.rho)
                        configurations, _ = path.sample_many(step_size)
                        temp_path = list(
                            map(lambda c: (c[0], c[1]), configurations))
                        temp_path.append((minConfGoal[0], minConfGoal[1]))
                        temp_path = LineString(temp_path)
                        # G.graph.add_edge(a, b, weight=temp_path.length, geom=temp_path, pathConf=[minConfStart, minConfGoal], parent_gates=[G.graph.nodes[a]['id'], G.graph.nodes[b]['id']])
                        G.graph.add_edge(a, b, weight=minLength, geom=temp_path, pathConf=[
                                         minConfStart, minConfGoal], parent_gates=[G.graph.nodes[a]['id'], G.graph.nodes[b]['id']])

                    except:
                        None
        prev_x_dist = nearest_gate_dist
        prev_gates = nearest_gates

    return G


def UpdatePath_Euclidean(G):
    # Graph.clear_edges()
    smallest_gate_x_dist, largest_gate_x_dist = np.inf, 0
    for n in list(G.graph.nodes):
        if G.graph.nodes[n]['x_dist'] < smallest_gate_x_dist:
            smallest_gate_x_dist = G.graph.nodes[n]['x_dist']
        if G.graph.nodes[n]['x_dist'] > largest_gate_x_dist:
            largest_gate_x_dist = G.graph.nodes[n]['x_dist']

    nearest_gates, nearest_gate_dist = [], np.inf
    for n in list(G.graph.nodes):
        if G.graph.nodes[n]['x_dist'] == smallest_gate_x_dist:
            nearest_gate_dist = G.graph.nodes[n]['x_dist']
            nearest_gates.append(n)

    prev_x_dist = nearest_gate_dist
    prev_gates = nearest_gates

    while prev_x_dist < largest_gate_x_dist:
        nearest_gates, nearest_gate_dist = [], np.inf
        for n in list(G.graph.nodes):
            if G.graph.nodes[n]['x_dist'] < nearest_gate_dist and G.graph.nodes[n]['x_dist'] > prev_x_dist:
                nearest_gate_dist = G.graph.nodes[n]['x_dist']
                nearest_gates = [n]
            elif G.graph.nodes[n]['x_dist'] == nearest_gate_dist:
                nearest_gates.append(n)

        for a in prev_gates:
            for b in nearest_gates:
                if not G.graph.has_edge(a, b):
                    try:
                        near_points = nearest_points(
                            G.graph.nodes[a]['geom'], G.graph.nodes[b]['geom'])
                        temp_path = LineString(near_points)
                        minConfStart = [near_points[0].x,
                                        near_points[0].y, find_slope(temp_path)]
                        minConfGoal = [near_points[1].x,
                                       near_points[1].y, find_slope(temp_path)]
                        G.graph.add_edge(a, b, weight=temp_path.length, geom=temp_path, pathConf=[
                                         minConfStart, minConfGoal], parent_gates=[G.graph.nodes[a]['id'], G.graph.nodes[b]['id']])
                    except:
                        None

        prev_x_dist = nearest_gate_dist
        prev_gates = nearest_gates

    return G


def GraphUpdate_BreakGate(G: Graph_Gstar, Path) -> Graph_Gstar:
    """
    Function to split the gates not satifying the continuty constraint.

    Attributes
    ----------
    Graph : nx.graph()
        A graph containing Gates as nodes and paths between gates as Edges
    Path : List
        Coordinates of the lower left and upper right corners of the bounding box containing the obstacle.
    PathConstraint : variable
        Coordinates of the center of the bounding box.
    Counter : int
        The polygon representing the obstacle.

    Output
    -------
    Graph : nx.graph()
        A graph containing Gates as nodes and paths between gates as Edges

    """

    NodesToDelete = []

    for i in range(1, len(Path)-1):
        # TO DO: Add boundary point sorting based on x-distace from start point
        # print('path', Path[i-1], Path[i])
        # print(G.graph.get_edge_data(Path[i-1], Path[i]))
        _, b = G.graph.get_edge_data(Path[i-1], Path[i])['geom'].boundary
        c, _ = G.graph.get_edge_data(Path[i], Path[i+1])['geom'].boundary

        if b.distance(c) > G.tolerances['continuity']:
            gate_bound_i, gate_bound_j = G.graph.nodes[Path[i]
                                                       ]['geom'].boundary
            avg_X = gate_bound_i.x
            avg_Y = (gate_bound_i.y + gate_bound_j.y)/2

            # if G.graph.nodes[i]['id'] == 's':
            #     G.graph.add_node(G.GateCounter, id='s', geom = LineString([gate_bound_i, (gate_bound_i.x,avg_Y)]), heading=G.graph.nodes[i]['headings'], x_dist=avg_X)
            #     G.GateCounter+=1
            #     G.graph.add_node(G.GateCounter, id='s', geom = LineString([(gate_bound_i.x,avg_Y), gate_bound_j]), heading=G.graph.nodes[i]['headings'], x_dist=avg_X)
            #     G.GateCounter+=1

            # elif G.graph.nodes[i]['id'] == 'e':
            #     G.graph.add_node(G.GateCounter, id='e', geom = LineString([gate_bound_i, (gate_bound_i.x,avg_Y)]), heading=G.graph.nodes[i]['headings'], x_dist=avg_X)
            #     G.GateCounter+=1
            #     G.graph.add_node(G.GateCounter, id='e', geom = LineString([(gate_bound_i.x,avg_Y), gate_bound_j]), heading=G.graph.nodes[i]['headings'], x_dist=avg_X)
            #     G.GateCounter+=1

            # else:
            if not G.graph.nodes[Path[i]]['id'] == 's' or not G.graph.nodes[Path[i]]['id'] == 'e':
                G.graph.add_node(G.GateCounter, id=G.GateCounter, geom=LineString(
                    [gate_bound_i, (avg_X, avg_Y)]), sector=G.graph.nodes[Path[i]]['sector'], x_dist=avg_X)
                G.gateGeomList.append(LineString(
                    [gate_bound_i, (avg_X, avg_Y)]))
                G.GateCounter += 1
                G.graph.add_node(G.GateCounter, id=G.GateCounter, geom=LineString(
                    [(avg_X, avg_Y), gate_bound_j]), sector=G.graph.nodes[Path[i]]['sector'], x_dist=avg_X)
                G.gateGeomList.append(LineString(
                    [(avg_X, avg_Y), gate_bound_j]))
                G.GateCounter += 1

                if not Path[i] in NodesToDelete:
                    NodesToDelete.append(Path[i])

            else:
                print("Graph Error")

    # print('Delete Nodes: ', NodesToDelete)
    for n in NodesToDelete:
        try:
            G.gateGeomList.remove(G.graph.nodes[n]['geom'])
            # print('removed gate')
        except ValueError:
            None
        G.graph.remove_node(n)  # check correctness

    return G


def GraphUpdate_BreakAngle(G: Graph_Gstar, Path: list) -> Graph_Gstar:
    """
    Updates the graph by breaking angles in the given path based on a tolerance threshold.

    Args:
        G (Graph_Gstar): The input graph.
        Path (list): The path in the graph.

    Returns:
        Graph_Gstar: The updated graph.

    """
    NodesToDelete = []

    for i in range(1, len(Path) - 1):
        # Get the angle information from the graph edges
        AngIn = G.graph.get_edge_data(Path[i - 1], Path[i])['pathConf'][1][2]
        AngOut = G.graph.get_edge_data(Path[i], Path[i + 1])['pathConf'][0][2]
        
        pathConfList = [AngIn, AngOut]
        
        # Check if the angle difference exceeds the tolerance threshold
        if utils.Angdiff(pathConfList[0], pathConfList[1]) > G.tolerances['angular'] * np.pi/180:
            MidAng = utils.MidAng(pathConfList[0], pathConfList[1])
            
            if not (G.graph.nodes[Path[i]]['id'] == 's' or G.graph.nodes[Path[i]]['id'] == 'e'):
                # Add two new nodes with updated sectors
                G.graph.add_node(G.GateCounter, id=G.GateCounter, geom=G.graph.nodes[Path[i]]['geom'], sector=(
                    pathConfList[0], MidAng), x_dist=G.graph.nodes[Path[i]]['x_dist'])
                G.GateCounter += 1
                G.graph.add_node(G.GateCounter, id=G.GateCounter, geom=G.graph.nodes[Path[i]]['geom'], sector=(
                    MidAng, pathConfList[1]), x_dist=G.graph.nodes[Path[i]]['x_dist'])
                G.GateCounter += 1
                
                if Path[i] not in NodesToDelete:
                    NodesToDelete.append(Path[i])
            else:
                print("Graph Error")

    # Remove nodes that were marked for deletion
    for n in NodesToDelete:
        G.graph.remove_node(n)

    return G


def PlotGstar(shortest_path, G, Map, label, save_path='./images', action='save'):
    """Plots the shortest path on a map, including obstacles and edges.

    Args:
        shortest_path (list): A list of node IDs representing the shortest path.
        G (GstarGraph): The graph to plot.
        Map (Map): The map containing obstacles and start/end points.
        label (str): The title of the plot.
        save_path (str, optional): The path to save the image to. Defaults to './images'.
        action (str, optional): Either 'save' to save the image or 'display' to show it. Defaults to 'save'.
    """
    plt.ioff()
    fig = plt.figure()

    # Plot obstacles
    for obstacle in Map.ObstacleList:
        x, y = obstacle.exterior.xy
        plt.fill(x, y, c="blue")

    # Plot start/end points
    plt.plot(Map.StartPoint.x, Map.StartPoint.y, 'x', c='black')
    plt.plot(Map.EndPoint.x, Map.EndPoint.y, 'x', c='black')

    # Plot edges
    for n in list(G.graph.nodes):
        node = G.graph.nodes[n]
        if mapping(node['geom'])['type'] == 'LineString':
            x, y = node['geom'].xy
            plt.plot(x, y, c="green")
        elif mapping(node['geom'])['type'] == 'MultiLineString':
            for line in list(node['geom']):
                x, y = line.xy
                plt.plot(x, y, c="green")

    # Plot shortest path
    for i in range(len(shortest_path)-1):
        u = shortest_path[i]
        v = shortest_path[i+1]
        edge = G.graph.get_edge_data(u, v)
        if mapping(edge['geom'])['type'] == 'LineString':
            x, y = edge['geom'].xy
            plt.plot(x, y, c="black", linewidth=3)
        elif mapping(edge['geom'])['type'] == 'MultiLineString':
            for line in list(edge['geom']):
                x, y = line.xy
                plt.plot(x, y, c="black", linewidth=3)

    # Set equal axis scaling and title
    plt.axis('equal')
    plt.title(label)

    # Either display or save the plot
    if action == 'display':
        plt.show(block=True)
    elif action == 'save':
        plt.savefig(save_path, format='svg')
    plt.close('all')


def PlotMap(Map, label, save_path='./images', action='save'):

    for obstacle in Map.ObstacleList:
        x, y = obstacle.exterior.xy
        plt.fill(x, y, c="blue")

    plt.plot(Map.StartPoint.x, Map.StartPoint.y, 'x', c='black')
    plt.plot(Map.EndPoint.x, Map.EndPoint.y, 'x', c='black')
    plt.axis('equal')
    plt.title(label)
    if action == 'display':
        plt.show(block=True)
    elif action == 'save':
        plt.savefig(save_path, format='svg')
    plt.close('all')