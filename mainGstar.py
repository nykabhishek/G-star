import yaml
import time
import csv
import os
import pickle as pkl
import numpy as np

from graphutils import Graph_Gstar
from Gstar import GstarPaths, PlotGstar

if __name__ == "__main__":

    # Set the time limit for the algorithm to execute
    timeLimit = 600

    # Set the initial number of sectors and starting and ending configurations
    initial_Sectors = 3
    start_conf = (0, 4.5, 0)
    end_conf = (16, 4.5, 0)

    # Specify whether the heading is restricted or not
    heading_restricted = True

    # Specify the list of radii to use for the algorithm
    radius_list = [1, 2]

    # Specify the paths to the instance files
    instance_paths = ['./instances/map1']

    # Loop through each instance path
    for path in instance_paths:
        # Load the tolerances for this instance from the YAML file
        with open(path+'/tolerances.yaml') as toleances_yaml:
            tolerances = yaml.load(toleances_yaml, Loader=yaml.FullLoader)

        # Load the heading angles for this instance from the YAML file
        with open(path+'/heading.yaml') as heading_yaml:
            headingAngles = yaml.load(heading_yaml, Loader=yaml.FullLoader)


        # Set up the result fields for the CSV file
        result_fields = [
            'name', 
            'path', 
            'Obstacles', 
            'turning_radius', 
            'continuity_tolerance', 
            'angle_tolerance', 
            'node_count', 
            'edge_count', 
            'eucLB_NoObstacles', 
            'dubLB_NoObstacles', 
            'eucLB', 
            'eucLB_time', 
            'dubLB', 
            'dubLB_time',  
            'dubUB', 
            'dubUB_time']
        result_filename = path+'/results_' + time.strftime("%Y%m%d-%H%M%S")+'.csv'

        # Open the CSV file and write the result fields to the first row
        with open(result_filename, 'w') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow(result_fields)

            # Loop through each instance path and file
            instance_results = []
            for (instance_path, dirs, files) in os.walk(path):
                for f in files:
                    # Check if the file is a pickled Map object
                    if f.endswith(".pkl") and f.split('_')[0] == 'Map':
                        obstacle_count = instance_path.split("/")[-2]
                        obstacle_count = obstacle_count.split('_')[1]

                        # Load the Map object from the pickled file
                        # filepath = instance_path + '/'+f
                        filepath = os.path.join(instance_path, f)
                        m = open(filepath, 'rb')
                        Map = pkl.load(m)

                        for rho in radius_list:
                            print("\n--------------------------------------------------------------------------")
                            # print("Instance: ", instance_path + '/'+ f, " Radius: ", rho)
                            print("Tolerances: ", tolerances)

                            # Instantiate the Graph_Gstar object
                            G = Graph_Gstar(start_conf, end_conf, rho, initial_Sectors, tolerances)

                            # Set the image save path
                            imgPath = instance_path+'/img'+'/r_'+str(rho)

                            # Record the start time
                            start_time = time.time()

                            # if not os.path.exists(instance_path+'/graphs/G_r'+str(rho)+'.pkl'):
                            
                            try:
                                # Calculate the GstarPaths
                                dubLB_path, Map, G = GstarPaths(Map, G, timeLimit, imgPath, heading_restricted, headingAngles)

                                # Plot the GstarPaths
                                PlotGstar(dubLB_path, G, Map, 'G* Lower Bound Path', save_path=instance_path+'/dub_LB_r'+str(rho), action='save')

                                total_time = time.time() - start_time
                                print('Time of execution: ', total_time)
                            except:
                                None


                            # Ensure the 'graphs' directory exists
                            graphs_path = os.path.join(instance_path, 'graphs')
                            if not os.path.exists(graphs_path):
                                os.makedirs(graphs_path)

                            # Save the graph to a file
                            graph_filename = os.path.join(graphs_path, f'G_r{rho}.pkl')
                            with open(graph_filename, 'wb') as f:
                                pkl.dump(G, f)

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
