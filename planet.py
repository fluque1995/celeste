import numpy as np
import utils
import scipy.special as sp


class Planet:
    def __init__(self, name, epsilon, a, period):
        self.name = name
        self.epsilon = epsilon
        self.a = a
        self.period = period
        self.orbit = None
        
        self.mu = 4*np.pi**2*self.a**3/self.period**2
        self.c = np.sqrt(self.mu*self.a*(1-epsilon**2))

    def __str__(self):
        return self.name

    def get_time_in_period(self, time):
        return time % self.period

    def position(self, time):

        time = self.get_time_in_period(time)
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

        time = self.get_time_in_period(time)
        eccentric_annomaly = utils.eccentric_annomaly(
            self.period,
            self.epsilon,
            time
        )

        denom = self.period*(1-self.epsilon*np.cos(eccentric_annomaly))
        quotient = 2*np.pi*self.a/denom

        speed = [
            (-quotient*np.sin(eccentric_annomaly))[0],
            (quotient *
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
        time = self.get_time_in_period(time)
        return utils.runge_kutta(self.real_annomaly_deriv, 0, 0, time, 5000)

    def energy(self):
        return -self.c**2/(2*self.a**2*(1-self.epsilon**2))

    def energy_from_time(self, time):
        return (self.speed_module(time)**2/2 -
                self.mu/self.distance_to_sun(time))

    def angular_moment_from_time(self, time):
        time = self.get_time_in_period(time)
        real_annomaly = self.real_annomaly(time)
        return (self.distance_to_sun(time)**2 *
                self.real_annomaly_deriv(time, real_annomaly))

    def real_annomaly_from_eccentric(self, time):
        ecc_an = utils.eccentric_annomaly(self.period,
                                          self.epsilon,
                                          time)
        real_an = np.arccos((np.cos(ecc_an)[0] - self.epsilon) /
                            (1 - self.epsilon*np.cos(ecc_an)[0]))
        real_an = real_an if 2*time < self.period else 2*np.pi - real_an
        return real_an

    def get_orbit(self, npoints):
        if self.orbit is None:
            ts = list(np.arange(0, self.period, self.period/npoints))
            self.orbit = [self.position(t) for t in ts]
            self.orbit.append(self.orbit[0])
            self.orbit = np.asarray(self.orbit)
        return self.orbit

    def eccentric_annomaly(self, time):
        time = self.get_time_in_period(time)
        return utils.eccentric_annomaly(self.period,
                                        self.epsilon,
                                        time)[0]
        
    def eccentric_annomaly_bessel(self, time, nfuncs):
        time = self.get_time_in_period(time)
        ji = 2*np.pi*time/self.period
        ecc_an = ji

        for i in range(nfuncs):
            bessel = sp.jv((i+1), (i+1)*self.epsilon)
            ecc_an += ((2/(i+1))*bessel*np.sin(2*np.pi*time/self.period))

        return ecc_an


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
