# Author: Abhishek Nayak
# from re import S
import sys
sys.path.append('DubinsLineSegToLineSeg/')

from graphutils import is_path_exists, is_route_continuous, is_feasible_path, is_angle_matching, get_shortest_path
from graphutils import ConstructGates, ConnectGates_Dubins, GraphUpdate_BreakGate, UpdatePath_Dubins, GraphUpdate_BreakAngle
from graphutils import ConnectGates_Euclidean, UpdatePath_Euclidean, EuclideanToDubinsGraph
from graphutils import PlotGstar
from graphutils import Graph_Gstar
from DubinsLineSegToLineSeg import DubinsL2L
from map_utils import Map
from copy import deepcopy
import pickle as pkl
import dubins as du
import time
import networkx as nx
import numpy as np
import csv
import os
import yaml
import logging

from networkx.algorithms.shortest_paths.generic import shortest_path
# from networkx.classes.graph import Graph
sys.path.append('DubinsLineSegToLineSeg/')

step_size = 0.1
bound_resolution = 1


def GstarPaths(Map, G, timeLimit, imgPath, heading_restricted, headingAngles):

    '''
    This function computes the Euclidean and Dubins lower and upper bounds of the shortest path
    between the start and goal points in the input Map. If the shortest path is feasible, the
    function terminates. Otherwise, it constructs and connects gates to update the graph and
    continues the computation until the shortest path is feasible or the time limit is reached.

    Args:
    - Map: An object of the Map class, containing the obstacles and start and goal points
    - G: An object of the Graph class, containing the nodes and edges of the graph
    - timeLimit: The maximum computation time allowed for finding the shortest path
    - imgPath: The path where the images of the computed paths will be saved
    - heading_restricted: A boolean indicating whether the heading angles at the start and goal
                          points are restricted or not
    - headingAngles: A dictionary containing the restricted heading angles at the start and goal
                     points (if any)

    Returns:
        A tuple containing the lower bounds of the G* algorithm, the map object, and the GStarGraph object.
    '''

    # Compute the Euclidean upper bound of the shortest path
    eucBoundStartTime = time.monotonic()
    short_path_list, G.eucLB_free = get_shortest_path('s', 'e', G)

    if heading_restricted:
        outAngle = headingAngles['startAngle']*np.pi/180
        inAngle = headingAngles['goalAngle']*np.pi/180
    else:
        outAngle, inAngle = 0, 0

    # Compute the Dubins upper bound of the shortest path
    delta = 0.001
    dubLB_free = DubinsL2L.Line2LineDubins(G.LineStart, (outAngle-delta, outAngle+delta), G.LineEnd, (inAngle-delta, inAngle+delta), G.rho)
    G.dubLB_free, minPath = dubLB_free.MinDub_L2L()

    # Print the Euclidean and Dubins upper bounds
    print("Euclidean Upper Bound (No Obstacles): ", G.eucLB_free)
    print("Dubins Upper Bound (No Obstacles): ", G.dubLB_free)

    ''' Step 1: Construct Euclidean Lower bounds '''
    if not is_feasible_path(short_path_list, G, Map.ObstacleList):
        improve_bounds = True
    else:
        improve_bounds = False

    print('### Solving Euclidean Lower Bound ###')
    loop = 0
    G.eucLB_time, eucBoundStartTime = 0, time.monotonic()
    improve_bounds = True
    
    while improve_bounds and (G.eucLB_time < 0.001*timeLimit):
        path_exists = is_path_exists('s', 'e', G)
        path_feasible = is_feasible_path(short_path_list, G, Map.ObstacleList)
        path_continuous = is_route_continuous(short_path_list, G)

        if path_exists and path_continuous and path_feasible:
            improve_bounds = False
        else:
            if not path_feasible:
                G = ConstructGates(G, Map, short_path_list,graphType='Euclidean')
                G = ConnectGates_Euclidean(G)
                short_path_list, short_path_length = get_shortest_path('s', 'e', G)

            if not path_continuous:
                G = GraphUpdate_BreakGate(G, short_path_list)
                G = UpdatePath_Euclidean(G)
                short_path_list, short_path_length = get_shortest_path('s', 'e', G)
        
        print('Path: ', short_path_list)
        G.eucLB_time = time.monotonic()-eucBoundStartTime
        loop += 1
        plotTitle = "Euclidean Lower Bound (length = " + str(round(short_path_length, 3)) + ")"
        PlotGstar(short_path_list, G, Map, plotTitle, save_path=imgPath+'/step_'+str(loop)+'_euc', action='save')

    G.eucLowerPath, G.eucLowerBound = short_path_list, short_path_length
    G.eucGraph = deepcopy(G.graph)

    print(nx.info(G.graph))
    print('Euclidean Lower Bound: ', G.eucLowerBound)
    print("Time (Euclindean Lower Bound): ", G.eucLB_time)
    plotTitle = "Euclidean Lower Bound (length = " + str(round(G.eucLowerBound, 3)) + ")"
    PlotGstar(G.eucLowerPath, G, Map, plotTitle, save_path=imgPath+'/eucLB', action='save')

    """ Step 2:  Compute Dubins Lower Bound """
    print('### Solving Dubins Lower Bound ###')
    times, Gstar_path_lengths = [], []
    dubBoundStartTime = time.monotonic()
    improve_bounds = True

    G = EuclideanToDubinsGraph(G, heading_restricted, startAngle=outAngle, goalAngle=inAngle)
    G = ConnectGates_Dubins(G)

    short_path_list, short_path_length = get_shortest_path('s', 'e', G)
    times.append(time.monotonic() - dubBoundStartTime)
    Gstar_path_lengths.append(short_path_length)
    loop += 1
    
    plotTitle = "G* (length = " + str(round(short_path_length, 3)) + ")"
    PlotGstar(short_path_list, G, Map, plotTitle, save_path=imgPath +'/step_'+str(loop)+'_dub', action='save')

    print(short_path_list)
    print('G* lower bound distance: ', short_path_length)
    # PlotGstar(short_path_list, G.graph, self.Map, loop, title)
    G.dubLB_time = time.monotonic()-dubBoundStartTime

    '''### G* lower bounds ###'''
    '''### REPLACE EUCLIDEAN PATHS WITH DubLineToLine ###'''
    while improve_bounds and (G.dubLB_time < timeLimit):
        path_exists = is_path_exists('s', 'e', G)
        path_feasible = is_feasible_path(short_path_list, G, Map.ObstacleList)
        path_continuous = is_route_continuous(short_path_list, G)
        angle_matching = is_angle_matching(short_path_list, G)

        if path_exists and path_continuous and path_feasible and angle_matching:
            improve_bounds = False
            # break
        else:
            if not path_feasible:
                logging.info('loop {}', loop)
                logging.info('Constructing New Gates...')
                G = ConstructGates(G, Map, short_path_list, graphType='Dubins')
                G = ConnectGates_Dubins(G)
                short_path_list, short_path_length = get_shortest_path('s', 'e', G)
                # PlotGstar(short_path_list, G, Map, title, save_path='./images', action='save')

            if not path_continuous:
                logging.info('Breaking Gates...')
                G = GraphUpdate_BreakGate(G, short_path_list)
                # G = ConnectGates_Dubins(G)
                G = UpdatePath_Dubins(G)
                short_path_list, short_path_length = get_shortest_path('s', 'e', G)
                # PlotGstar(short_path_list, G, Map, title, save_path='./images', action='display')

            if not angle_matching:
                logging.info('Optimizing Angle...')
                G = GraphUpdate_BreakAngle(G, short_path_list)
                ## G = ConnectGates_Dubins(G)
                G = UpdatePath_Dubins(G)
                short_path_list, short_path_length = get_shortest_path('s', 'e', G)
                # PlotGstar(short_path_list, G, Map, title, save_path='./images', action='display')

            G.dubLB_time = time.monotonic()-dubBoundStartTime
            times.append(G.dubLB_time)
            Gstar_path_lengths.append(short_path_length)
            loop += 1
            plotTitle = "G* (length = " + str(round(short_path_length, 3)) + ")"
            PlotGstar(short_path_list, G, Map, plotTitle, save_path=imgPath+'/step_'+str(loop)+'_dub', action='save')
            
    
    # short_path_list, short_path_length = get_shortest_path('s', 'e', G)
    G.dubLowerPath, G.dubLowerBound = short_path_list, short_path_length
    print("G* Bound Time: ", G.dubLB_time)
    print("G* Path: ", G.dubLowerPath)
    print('G* lower bound distance: ', G.dubLowerBound)
    plotTitle = "G* (length = " + str(round(G.dubLowerBound, 3)) + ")"

    return G.dubLowerPath, Map, G


if __name__ == "__main__":

    timeLimit = 1800
    radius_list = [1, 2, 3]
    initial_Sectors = 3
    start_conf = (0, 4.5, 0)
    end_conf = (16, 4.5, 0)
    heading_restricted = True

    instance_paths = ['./test_folder/final']

    for path in instance_paths:

        with open(path+'/tolerances.yaml') as toleances_yaml:
            tolerances = yaml.load(toleances_yaml, Loader=yaml.FullLoader)

        with open(path+'/heading.yaml') as heading_yaml:
            headingAngles = yaml.load(heading_yaml, Loader=yaml.FullLoader)

        result_fields = ['name', 'path', 'Obstacles', 'turning_radius', 'continuity_tolerance', 'angle_tolerance', 'node_count',
                         'edge_count', 'eucLB_NoObstacles', 'dubLB_NoObstacles', 'eucLB', 'eucLB_time', 'dubLB', 'dubLB_time',  'dubUB', 'dubUB_time']
        result_filename = path+'/results_' + \
            time.strftime("%Y%m%d-%H%M%S")+'.csv'

        # writing to csv file
        with open(result_filename, 'w') as csvfile:
            # creating a csv writer object
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(result_fields)

            # writing the data rows
            instance_results = []
            for (instance_path, dirs, files) in os.walk(path):
                for f in files:
                    if f.endswith(".pkl"):
                        obstacle_count = instance_path.split("/")[-2]
                        obstacle_count = obstacle_count.split('_')[1]

                        # M = Map(start_conf, end_conf, (15,8), obstacle_count, shape=5)
                        filepath = instance_path+'/' + f
                        if f.split('_')[0] == 'Map':
                            m = open(filepath, 'rb')
                            Map = pkl.load(m)

                            for rho in radius_list:
                                print(
                                    "\n--------------------------------------------------------------------------")
                                print("Instance: ", instance_path +
                                      '/'+f, " Radius: ", rho)
                                print("Tolerances: ", tolerances)
                                G = Graph_Gstar(
                                    start_conf, end_conf, rho, initial_Sectors, tolerances)
                                imgPath = instance_path+'/img/r_'+str(rho)
                                if not os.path.exists(imgPath):
                                    os.mkdir(imgPath)

                                try:

                                    start_time = time.monotonic()
                                    # if not os.path.exists(instance_path+'/graphs/G_r'+str(rho)+'.pkl'):
                                    try:
                                        dubLB_path, Map, G = GstarPaths(
                                            Map, G, timeLimit, imgPath, heading_restricted, headingAngles)
                                        PlotGstar(dubLB_path, G, Map, 'G* Lower Bound Path',
                                                    save_path=instance_path+'/dub_LB_r'+str(rho), action='save')
                                        print('test')

                                        total_time = time.monotonic() - start_time
                                        print('Time of execution: ',
                                                total_time)
                                    except:
                                        None

                                    if not os.path.exists(instance_path+'/graphs'):
                                        os.mkdir(instance_path+'/graphs')
                                    G_filename = instance_path + \
                                        '/graphs/G_r'+str(rho)+'.pkl'
                                    G_file = open(G_filename, 'wb')
                                    pkl.dump(G, G_file)
                                    G_file.close()

                                    instance_data = [
                                        f, 
                                        filepath, 
                                        obstacle_count, 
                                        rho, 
                                        tolerances['continuity'], 
                                        tolerances['angular'], 
                                        G.graph.number_of_nodes(), 
                                        G.graph.size(), 
                                        G.eucLB_free, 
                                        G.dubLB_free, 
                                        G.eucLowerBound, 
                                        G.eucLB_time, 
                                        G.dubLowerBound, 
                                        G.dubLB_time, 
                                        G.dubUpperBound, 
                                        G.dubUB_time
                                        ]
                                    csvwriter.writerow(instance_data)

                                except:
                                    None