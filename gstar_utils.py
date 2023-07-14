import numpy as np
from typing import List
import networkx as nx
from shapely.geometry import Point, LineString

def HeadingsToSectors(headings):
    '''This function takes a list of headings and returns a list of sectors
    
    Parameters
    ----------
    headings
        a list of headings, in degrees, that you want to convert to sectors.
    
    '''
    SectorList = [(headings[k], headings[k+1]) for k in range(len(headings)-1)]
    # SectorList.append((headings[-1], 2*np.pi))
    SectorList.append((headings[-1], 0))
    return SectorList

def is_route_feasible(route_dict) -> bool:
    """Check if the given routes are feasible.

    Parameters
    ----------
    route_dict : dict
        A dictionary of the form {'route_id': [list of stops]}.

    Returns
    -------
    bool
        True if all routes are feasible, False otherwise.

    """
    for path in route_dict.values():
        if not path['feasible']:
            return False
    return True


def retrieve_gate(index, gate_dict):
    '''
    Retrieves a gate object with a matching index from a list of gate dictionaries.

    Parameters:
        index (int): The index of the gate to be retrieved.
        gate_dict (list[dict]): A list of gate dictionaries, where each dictionary represents a gate.

    Returns:
        dict: The gate dictionary that matches the given index.
    '''
    for item in gate_dict:
        if item['id'] == index:
            return item


def is_path_exists(StartId, EndId, G) -> bool:
    exists = True
    start_id_list, end_id_list = node_id_list(StartId, G.graph), node_id_list(EndId, G.graph)
    for s in start_id_list:
        for e in end_id_list:
            exists = nx.has_path(G.graph, source=s, target=e)
            if not exists:
                break
    return exists

    # """
    # Determines if a path exists between the start and end nodes in the given graph.

    # Args:
    #     StartId (str): The ID of the start node.
    #     EndId (str): The ID of the end node.
    #     G (GstarGraph): The graph to search for a path.

    # Returns:
    #     bool: True if a path exists, False otherwise.
    # """
    # start_id_list, end_id_list = node_id_list(StartId, G.graph), node_id_list(EndId, G.graph)

    # return any(nx.has_path(G.graph, source=s, target=e) for s in start_id_list for e in end_id_list)

def node_id_list(id: int, graph: nx.Graph) -> List[int]:
    """Returns a list of node IDs that match a given ID in the given graph.

    Args:
        id (int): The ID to match.
        graph (nx.Graph): The graph to search.

    Returns:
        List[int]: A list of node IDs in the graph that match the given ID.
    """
    id_list = []
    for node in graph.nodes:
        # Check if the node's 'id' attribute matches the given ID
        if graph.nodes[node]['id'] == id:
            id_list.append(node)
    return id_list

def valid_connection(connection, pathway):
    valid = True
    for path in pathway:
        if connection==path or connection.within(path) or path.contains(connection) or path.overlaps(connection):
            valid = False
    return valid

def find_projection(point, line):
    """
    Finds the projection of a point onto a line using NumPy.

    Args:
        point: A shapely.geometry.Point object representing the point to project.
        line: A shapely.geometry.LineString object representing the line to project onto.

    Returns:
        A tuple containing the (x, y) coordinates of the projection.
    """
    if not isinstance(point, Point):
        raise ValueError("The 'point' argument must be a shapely.geometry.Point object.")
    if not isinstance(line, LineString):
        raise ValueError("The 'line' argument must be a shapely.geometry.LineString object.")

    x = np.array(point.coords[0])
    start_point = np.array(line.coords[0])
    end_point = np.array(line.coords[len(line.coords) - 1])
    n = end_point - start_point
    n /= np.linalg.norm(n, 2)
    projection_coords = start_point + n * np.dot(x - start_point, n)
    return projection_coords

def find_slope(line):
    """
    Compute the slope of a LineString object.

    Parameters:
    line (LineString): a LineString object representing a line segment

    Returns:
    float: the slope of the line segment, in radians, measured from the positive x-axis in counterclockwise direction

    The slope of a line is the angle it makes with the positive x-axis, measured counterclockwise from the x-axis.

    """
    a, b = line.boundary
    
    # linear equation: y = k*x + m
    # k = (b.y - a.y) / (b.x - a.x)
    # m = a.y - k * a.x
    # slope = np.arctan(k)

    slope = np.arctan((b.y - a.y)/ (b.x - a.x))
    if slope < 0:
        slope = 2*np.pi + slope

    return slope

def retrieve_minEdge(graph, node_id_1, node_id_2):
    """
    Retrieve the minimum weight edge between two nodes in a graph.

    Parameters:
    graph (networkx.Graph): The graph to search for the edge.
    node_id_1 (hashable): The ID of the first node.
    node_id_2 (hashable): The ID of the second node.

    Returns:
    dict: The edge dictionary representing the minimum weight edge between the two nodes.
    """
    edge_dict = graph.get_edge_data(node_id_1, node_id_2)
    min_edge = min(edge_dict.values(), key=lambda e: e['weight'])
    return min_edge


# def extend_gate(self, temp_gate):
#     #NEEDS REFINEMENT

#     a, b = temp_gate.boundary
#     if a.x == b.x:  # vertical line
#         extended_gate = LineString([(a.x, self.map_miny), (a.x, self.map_maxy)])
#     elif a.y == b.y:  # horizonthal line
#         extended_gate = LineString([(self.map_minx, a.y), (self.map_maxx, a.y)])
#     else:
#         # linear equation: y = k*x + m
#         k = (b.y - a.y) / (b.x - a.x)
#         m = a.y - k * a.x
#         y0 = k * self.map_minx + m
#         y1 = k * self.map_maxx + m
#         x0 = (self.map_miny - m) / k
#         x1 = (self.map_maxy - m) / k
#         points_on_boundary_lines = [Point(self.map_minx, y0), Point(self.map_maxx, y1), 
#                                     Point(x0, self.map_miny), Point(x1, self.map_maxy)]
#         points_sorted_by_distance = sorted(points_on_boundary_lines, key=self.bounding_box.distance)
#         extended_gate = LineString(points_sorted_by_distance[:2])

#     print('gate:', extended_gate)

#     '''
#     print(MultiPolygon(self.map_obstacle).wkt)
#     gate = split(gate, MultiPolygon(self.map_obstacle)) 
#     '''

#     refined_gate = extended_gate-MultiPolygon(self.map_obstacle)
#     print('difference:', refined_gate)

#     self.gate_list.append(refined_gate)
#     self.gates.append(refined_gate)

#     return

def find_dubins_path_points(path_dict):
    return


if __name__ == '__main__':

    n_disc = 3
    headings = [k*2*np.pi/(n_disc) for k in range(n_disc)]
    sectors = HeadingsToSectors(headings)
    print('Sectors ;', sectors)
