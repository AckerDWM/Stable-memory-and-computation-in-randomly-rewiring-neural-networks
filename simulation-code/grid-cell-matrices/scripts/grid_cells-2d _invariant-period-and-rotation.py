from __future__ import division
import numpy as np

def grids(pixels):
    pixels = pixels / 100
    c = np.random.uniform(low=0.0, high=99.0, size=2)
    _lambda = .6 # inter-vertex spacing
    theta = 0.0 # rotation
    rotation_matrix = np.array([
        [np.cos(theta), -np.sin(theta)], 
        [np.sin(theta), np.cos(theta)]
        ])
    w_k = np.array(
        [[np.cos(np.pi/2),np.sin(np.pi/2)], 
         [np.cos(np.pi/6),np.sin(np.pi/6)], 
         [np.cos(np.pi/6),np.sin(-np.pi/6)]])
    w_k = np.dot(w_k, rotation_matrix)
    rate = np.zeros(100)
    for i in range(100):
        r = np.array([pixels[i], 1])
        vals = np.zeros(3)
        for j in range(3):
            spacing = (4*np.pi)/(3**0.5*_lambda)
            vals[j] = np.cos(spacing*np.dot(w_k[j,:], (r-c)))
        rate[i] = np.sum(vals)
    return rate

g = lambda x: np.exp(0.3 * (x+3/2))-1
G = lambda pixels: g(grids(pixels))
pixels = np.arange(100)
grid_cells = np.array([G(pixels) for _ in range(10000)])

# save the result
fname = "/Users/danielacker/Desktop/grid_cells-2d-invariant-period-and-rotation.npy"
np.save(file=fname, arr=grid_cells)


