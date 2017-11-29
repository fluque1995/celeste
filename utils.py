import numpy as np


def get_x1(period, epsilon, time):
    ji = 2*np.pi*time/period
    u_1 = -1
    if time < period/2:
        while u_1 < 0 or u_1 > np.pi:
            u_0 = np.random.uniform(0, np.pi, 1)
            num = (-u_0*np.cos(u_0) + np.sin(u_0))*epsilon + ji
            u_1 = num/(1-epsilon*np.cos(u_0))
    else:
        while u_1 < np.pi or u_1 > 2*np.pi:
            u_0 = np.random.uniform(np.pi, 2*np.pi, 1)
            num = (-u_0*np.cos(u_0) + np.sin(u_0))*epsilon + ji
            u_1 = num/(1-epsilon*np.cos(u_0))

    return u_1


def eccentric_annomaly(period, epsilon, time, tol=0.001):
    ji = 2*np.pi*time/period

    curr_u = next_u = get_x1(period, epsilon, time)
    dist = np.infty
    
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
