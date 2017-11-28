import numpy as np


def eccentric_annomaly(period, epsilon, time, tol=0.0001):
    ji = 2*np.pi*time/period

    curr_u = np.random.uniform(0, np.pi, 1)
    num = (-curr_u*np.cos(curr_u) + np.sin(curr_u))*epsilon + ji
    next_u = num/(1-epsilon*np.cos(curr_u))

    while next_u < 0 or next_u > np.pi:
        curr_u = np.random.uniform(0, np.pi, 1)
        num = (-curr_u*np.cos(curr_u) + np.sin(curr_u))*epsilon + ji
        next_u = num/(1-epsilon*np.cos(curr_u))

    dist = abs(curr_u - next_u)

    while dist > tol:
        curr_u = next_u
        numer = (-curr_u*np.cos(curr_u) + np.sin(curr_u))*epsilon + ji
        next_u = numer/(1-epsilon*np.cos(curr_u))
        dist = abs(curr_u - next_u)
    
    return next_u


def runge_kutta(func, t_0, x_0, t_final, steps):
    h = (t_0 + t_final)/steps
    curr_x = x_0

    for i in range(steps):
        k_1 = func(t_0 + i*h, curr_x)
        k_2 = func(t_0 + i*h + h/2, curr_x + k_1*h/2)
        k_3 = func(t_0 + i*h + h/2, curr_x + k_2*h/2)
        k_4 = func(t_0 + (i+1)*h + h/2, curr_x + k_3*h)
        curr_x = curr_x + h*(k_1 + 2*k_2 + 2*k_3 + 1*k_4)/6

    return curr_x
