from dash import html, register_page
import dash_bootstrap_components as dbc

from definitions.general_funcs import bold_it
from definitions.clnm_funcs import display_column2

register_page(__name__, path='/clnm')

layout = dbc.Row([
    dbc.Col(width={'size': 10, 'offset': 1},
            children=[
                html.Br(),
                html.Div(['Results of the', bold_it('cross-lag network analysis'), 'performed using the variables listed below.']),
                html.Hr(),

                dbc.Row([
                    # Temporal network
                    display_column2('t', '6vw'),
                    # Contemporaneous network
                    display_column2('c', '36vw'),
                    # Between person network
                    display_column2('b', '65vw')
                ])
            ])
])
