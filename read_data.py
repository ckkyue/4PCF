import numpy as np
import redshift_distance as rd

def read(data_path):
    data = []
    with open(data_path, "r") as file:
        lines = file.readlines()
        for line in lines[1:]:  # Exclude the first line
            data.append(line.strip().split(","))
    return data

def get_weights(data):
    weights = []
    for i in range(len(data)):
        row = data[i]
        w_fkp = np.float64(row[4])
        w_sys = np.float64(row[9])
        w_noz = np.float64(row[6])
        w_cp = np.float64(row[5])
        w = w_fkp * w_sys * (w_noz + w_cp - 1)
        weights.append(w)
    return weights

def get_carts(data):
    carts = []
    for i in range(len(data)):
        row = data[i]
        RA = np.float64(row[1])
        DEC = np.float64(row[2])
        redshift = np.float64(row[3])
        r = rd.redshift_to_dist(redshift)
        DEC = DEC * np.pi / 180
        RA = RA * np.pi / 180
        x = r * np.cos(DEC) * np.cos(RA)
        y = r * np.cos(DEC) * np.sin(RA)
        z = r * np.sin(DEC)
        cart = [x, y, z]
        carts.append(cart)
    return carts