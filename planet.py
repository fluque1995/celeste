import numpy as np
import utils


class Planet:
    def __init__(self, name, epsilon, a, period):
        self.name = name
        self.epsilon = epsilon
        self.a = a
        self.period = period
        
        self.mu = 4*np.pi**2*self.a**3/self.period**2
        self.c = np.sqrt(self.mu*self.a*(1-epsilon**2))

    def __str__(self):
        return self.name

    def position(self, time):
        eccentric_annomaly = utils.eccentric_annomaly(
            self.period,
            self.epsilon,
            time
        )
        
        position = [
            self.a*(np.cos(eccentric_annomaly) - self.epsilon)[0],
            self.a*np.sqrt(1-self.epsilon**2)*np.sin(eccentric_annomaly)[0]
        ]
        
        return position

    def distance_to_sun(self, time):
        return np.linalg.norm(self.position(time))

    def speed(self, time):
        eccentric_annomaly = utils.eccentric_annomaly(
            self.period,
            self.epsilon,
            time
        )

        denom = self.a**2
        denom *= np.sqrt(1 - self.epsilon**2)
        denom *= 1 - self.epsilon*np.cos(eccentric_annomaly)
        quotient = self.c/denom

        speed = [
            (-self.a*quotient*np.sin(eccentric_annomaly))[0],
            (self.a*quotient *
             np.sqrt(1 - self.epsilon**2)*np.cos(eccentric_annomaly))[0]
        ]

        return speed

    def speed_module(self, time):
        return np.linalg.norm(self.speed(time))

    def real_annomaly_deriv(self, t, theta):
        num = self.c*(1 + self.epsilon*np.cos(theta))**2
        denom = (self.a**2*(1 - self.epsilon**2)**2)
        return num/denom

    def real_annomaly(self, time):
        return utils.runge_kutta(self.real_annomaly_deriv, 0, 0, time, 5000)

    def energy(self):
        return -self.c**2/(2*self.a**2*(1-self.epsilon**2))

    def energy_from_time(self, time):
        return (self.speed_module(time)**2/2 -
                self.mu/self.distance_to_sun(time))

    def angular_moment_from_time(self, time):
        real_annomaly = self.real_annomaly(time)
        return (self.distance_to_sun(time)**2 *
                self.real_annomaly_deriv(time, real_annomaly))


planets_dict = {
    "Mercurio": Planet("Mercurio", 0.206, 0.387, 87.97),
    "Venus": Planet("Venus", 0.007, 0.723, 224.7),
    "La Tierra": Planet("La tierra", 0.017, 1, 365.26),
    "Marte": Planet("Marte", 0.093, 1.524, 686.98),
    "Júpiter": Planet("Júpiter", 0.048, 5.203, 4332.6),
    "Saturno": Planet("Saturno", 0.056, 9.546, 10759),
    "Urano": Planet("Urano", 0.047, 19.2, 30687),
    "Neptuno": Planet("Neptuno", 0.009, 30.09, 60784),
}
