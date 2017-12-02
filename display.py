import plotly.offline as py
import plotly.graph_objs as go


class Displayer:

    def print_information(self, planet, time):

        print("Posición de {} en el día {}: {}"
              .format(planet.name, time, planet.position(time)))

        print("Distancia al sol de {} en el día {}: {}"
              .format(planet.name, time, planet.distance_to_sun(time)))

        print("Velocidad de {} en el día {}: {}"
              .format(planet.name, time, planet.speed(time)))

        print("Módulo de la velocidad de {} en el día {}: {}"
              .format(planet.name, time, planet.speed_module(time)))

        print("Anomalía real de {} en el día {}: {}"
              .format(planet.name, time, planet.real_annomaly(time)))

        print("Energía de {} en el día {}: {}"
              .format(planet.name, time, planet.energy_from_time(time)))

        print("Energía (constante) de {}: {}"
              .format(planet.name, planet.energy()))

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
            width=700, height=500,
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
