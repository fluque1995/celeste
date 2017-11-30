import numpy as np
import utils
import plotly.offline as py
import plotly.graph_objs as go


class Planet:
    def __init__(self, name, epsilon, a, period):
        self.name = name
        self.epsilon = epsilon
        self.a = a
        self.period = period
        
        self.mu = 4*np.pi**2*self.a**3/self.period**2
        self.c = np.sqrt(self.mu*self.a*(1-epsilon**2))

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

    def print_information(self, time):
        print("Posición de {} en {}: {}".format(self.name,
                                                time,
                                                self.position(time)))
        print("Distancia al sol de {} en {}: {}".format(self.name,
                                                        time,
                                                        self.
                                                        distance_to_sun(time)))
        print("Velocidad de {} en {}: {}".format(self.name,
                                                 time,
                                                 self.speed(time)))
        print("Módulo de la velocidad de {} en {}: {}"
              .format(self.name,
                      time,
                      self.
                      speed_module(time)))
        print("Anomalía real de {} en {}: {}"
              .format(self.name,
                      time,
                      self.real_annomaly(time)))
        print("Energía de {} en {}: {}"
              .format(self.name,
                      time,
                      self.energy_from_time(time)))
        print("Energía (constante) de {}: {}"
              .format(self.name,
                      self.energy()))

    def display_orbit(self, time):
        ts = list(np.arange(0, self.period, self.period/500))
        xs = [self.position(t) for t in ts]
        xs.append(xs[0])
        xs = np.asarray(xs)

        orbit = go.Scattergl(x=xs[:, 0], y=xs[:, 1],
                             name='órbita')

        position = self.position(time)
        planet = go.Scattergl(x=[position[0]], y=[position[1]],
                              mode='markers',
                              marker=dict(
                                  size=15,
                                  color='rgba(152, 0, 0, .8)',
                                  line=dict(
                                      width=2,
                                      color='rgb(0, 0, 0)'
                                  )
                              ),
                              name='{}: día {}'
                              .format(self.name, time))

        sun = go.Scattergl(x=[0], y=[0],
                           mode='markers',
                           marker=dict(
                               size=20,
                               color='rgba(230, 230, 0, .9)',
                               line=dict(width=1, color='rgb(100,100,0)')
                           ),
                           name='Sol')
        rng = int(self.a) + 1
        layout = go.Layout(
            width=700, height=600,
            xaxis=dict(
                anchor='y',
                range=[-rng, rng]
            ),
            yaxis=dict(
                anchor='x',
                autorange=False,
                range=[-rng, rng],
            )
        )
        data = [orbit, planet, sun]
        fig = go.Figure(data=data, layout=layout)
        py.plot(fig, filename='planet-orbit.html')


planets_list = [
    Planet("Mercurio", 0.206, 0.387, 87.97),
    Planet("Venus", 0.007, 0.723, 224.7),
    Planet("La tierra", 0.017, 1, 365.26),
    Planet("Marte", 0.093, 1.524, 686.98),
    Planet("Júpiter", 0.048, 5.203, 4332.6),
    Planet("Saturno", 0.056, 9.546, 10759),
    Planet("Urano", 0.047, 19.2, 30687),
    Planet("Neptuno", 0.009, 30.09, 60784),
]

planets_list[7].print_information(25)

planets_list[7].display_orbit(25)

