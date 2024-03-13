import plotly.graph_objects as go
# import my_constants as constants
import Analysis.constants as constants

class GraphObject:

    def __init__(self, fname, xAxis, yAxis, title, dataPairs):
        self.path = constants.PATH
        self.fname = fname
        self.xAxis = xAxis
        self.yAxis = yAxis
        self.title = title
        self.colors = constants.COLORS
        self.colorCounters = {}
        for color in self.colors.keys():
            self.colorCounters[color] = 0
        self.dataPairs = dataPairs

    def _update_fig(self, figure):
        figure.update_layout(
            xaxis=dict(
                showline=True,
                showgrid=False,
                showticklabels=True,
                linecolor='rgb(51, 51, 51)',
                linewidth=2,
                ticks='outside',
                tickfont=dict(
                    family='Helvetica',
                    size=12,
                    color='rgb(51, 51, 51)',
                ),
            ),
            yaxis=dict(
                showgrid=False,
                zeroline=False,
                showline=True,
                showticklabels=True,
                linecolor='rgb(51, 51, 51)',
                linewidth=2,
                ticks='outside',
                tickfont=dict(
                    family='Helvetica',
                    size=12,
                    color='rgb(51, 51, 51)',
                ),
            ),
            autosize=False,
            margin=dict(
                # autoexpand=True,
                l=100,
                r=20,
                t=110,
            ),
            # showlegend=False,
            plot_bgcolor='white'
        )

    def _get_next_color(self, palette):
        result = self.colors[palette][self.colorCounters[palette] % len(self.colors[palette])]
        self.colorCounters[palette] += 1
        return result

    def _create_trace(self, xVals, yVals, label, mode, palette):
        return go.Scatter(x=xVals, y=yVals, mode=mode, 
            name=label, 
            line=dict(color=self._get_next_color(palette), width=4), 
            marker=dict(color=self._get_next_color(palette),
                        size=10,
                        line=dict(
                            color='rgba(51, 51, 51, 1)',
                            width=1
                        )),
            connectgaps=False,
            # text=texts
        )
    
    def _plot_data(self, traces):
        fig = go.Figure()
        for trace in traces:
            fig.add_trace(trace)
        
        self._update_fig(fig)
        fig.update_layout(
            title={
                'text': self.title,
                'y':0.9,
                'x':0.5,
                'xanchor': 'center',
                'yanchor': 'top',
                'font': dict(color='#333333'), },
            xaxis_title=dict(text=self.xAxis, font=dict(color='#333333')),
            yaxis_title=dict(text=self.yAxis, font=dict(color='#333333')))

        with open(f'{self.path}scripts/regression/exports/evaluation/html/{self.fname}.html', 'w') as f:
            f.write(fig.to_html(include_plotlyjs='cdn'))
        fig.write_image(f'{self.path}scripts/regression/exports/evaluation/svg/{self.fname}.svg')
        fig.write_image(f'{self.path}scripts/regression/exports/evaluation/jpg/{self.fname}.jpg')


    def plot(self):
        traces = []
        for xVals, yVals, label, mode, palette in self.dataPairs:
            traces.append(self._create_trace(xVals, yVals, label, mode, palette))
        self._plot_data(traces)