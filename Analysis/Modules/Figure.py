import os
import plotly.graph_objects as go

class Styler:

    def __init__(self):
        self.BLACK = "rgb(51, 51, 51)"
        self.GREY = "rgb(236, 236, 236)"
        self.TITLE_SIZE = 12
        self.FONT = "Helvetica"
        self.BG_COLOR = "white"
        self.TICK_SIZE = 14
        self.AXIS_TITLE_SIZE = 18
        self.ANNOTATION_SIZE = 14

    def _update_fig(self, figure, title=""):
        """A wrapper function with common settings for Figure object from plotly

        Args:
            figure (obj): a plotly Figure graph object
        """
        figure.update_layout(
            title=title,
            font=dict(
                family=self.FONT,
                size=self.TITLE_SIZE,
                color=self.BLACK,
            ),
            # showlegend=False,
            plot_bgcolor=self.BG_COLOR,
        )

    def _update_axes(self, figure, xdtick=1, ydtick=1, xtitle="", ytitle=""):

        figure.update_xaxes(
            title=dict(
                text=xtitle,
                font=dict(
                    family=self.FONT,
                    size=self.AXIS_TITLE_SIZE,
                    color=self.BLACK,
                ),
            ),
            showline=True,
            showgrid=True,
            mirror=True,
            dtick=xdtick,
            showticklabels=True,
            linecolor=self.BLACK,
            linewidth=2,
            ticks="outside",
            tickfont=dict(
                family=self.FONT,
                size=self.TICK_SIZE,
                color=self.BLACK,
            ),
        )
        figure.update_yaxes(
            title=dict(
                text=ytitle,
                font=dict(
                    family=self.FONT,
                    size=self.AXIS_TITLE_SIZE,
                    color=self.BLACK,
                ),
            ),
            showgrid=True,
            gridcolor=self.GREY,
            zeroline=True,
            mirror=True,
            dtick=ydtick,
            zerolinecolor=self.BLACK,
            showline=True,
            showticklabels=True,
            linecolor=self.BLACK,
            linewidth=2,
            ticks="outside",
            tickfont=dict(
                family=self.FONT,
                size=self.TICK_SIZE,
                color=self.BLACK,
            ),
        )

    def _save_fig(
        self, figure:go.Figure, path:str, fname:str, html:bool=False, jpg:bool=False, svg:bool=False, pdf: bool = True
    ):
        if html:
            figure.write_html(
                os.path.join(path, "html", f"{fname}.html"),
                include_plotlyjs="cdn",
            )
        if jpg:
            figure.write_image(os.path.join(path, "jpg", f"{fname}.jpg"), scale=4.0)
        if svg:
            figure.write_image(os.path.join(path, "svg", f"{fname}.svg"))
        if pdf:
            figure.write_image(os.path.join(path, "pdf", f"{fname}.pdf"))
