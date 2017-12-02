import plotly.offline as py
import plotly.graph_objs as go


class Displayer:

    def compare_eccentric_anomalies(self, planet, time):
        print("Anomalía excéntrica de {} el día {},"
              .format(planet.name, time) +
              "calculada usando funciones de Bessel: {}"
              .format(planet.eccentric_annomaly_bessel(time, 20)))

        print("Anomalía excéntrica de {} el día {},"
              .format(planet.name, time) +
              "calculada por el método de Newton: {}"
              .format(planet.eccentric_annomaly(time)))

        print("Diferencia entre ambos valores: {}"
              .format(abs(planet.eccentric_annomaly_bessel(time, 20) -
                          planet.eccentric_annomaly(time))))

    def print_information(self, planet, time):

        print("Posición de {} en el día {}: {}"
              .format(planet.name, time, planet.position(time)))

        print("Distancia al sol de {} en el día {}: {}\n"
              .format(planet.name, time, planet.distance_to_sun(time)))

        print("Velocidad de {} en el día {}: {}"
              .format(planet.name, time, planet.speed(time)))

        print("Módulo de la velocidad de {} en el día {}: {}\n"
              .format(planet.name, time, planet.speed_module(time)))

        print("Anomalía real de {} en el día {}: {}"
              .format(planet.name, time, planet.real_annomaly(time)))

        print("Anomalía real de {} en el día {}"
              .format(planet.name, time) +
              " (cálculo a partir de la anomalía excéntrica): {}\n"
              .format(planet.real_annomaly_from_eccentric(time)))

        print("Energía de {} en el día {}: {}"
              .format(planet.name, time, planet.energy_from_time(time)))

        print("Energía (constante) de {}: {}\n"
              .format(planet.name, planet.energy()))

        print("Módulo del momento angular de {} en el día {}: {}"
              .format(planet.name, time,
                      planet.angular_moment_from_time(time)))

        print("Módulo del momento angular de {} (constante): {}\n"
              .format(planet.name, planet.c))

        self.compare_eccentric_anomalies(planet, time)

    def display_orbit(self, planet, time):

        xs = planet.get_orbit(200)

        orbit = go.Scattergl(x=xs[:, 0], y=xs[:, 1],
                             name='órbita')

        position = planet.position(time)
        planet_pos = go.Scattergl(x=[position[0]], y=[position[1]],
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
                                  .format(planet.name, time))

        sun = go.Scattergl(x=[0], y=[0],
                           mode='markers',
                           marker=dict(
                               size=20,
                               color='rgba(230, 230, 0, .9)',
                               line=dict(width=1, color='rgb(100,100,0)')
                           ),
                           name='Sol')

        rng = int(planet.a) + 1

        layout = go.Layout(
            width=600, height=400,
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
        data = [orbit, planet_pos, sun]
        fig = go.Figure(data=data, layout=layout)
        py.iplot(fig)

