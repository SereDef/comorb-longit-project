from dash import Dash, html, page_container
import dash_bootstrap_components as dbc

import definitions.layout_styles_local as styles

from callbacks import clpm_callbacks, clnm_callbacks, csnm_callbacks

FONT_AWESOME = "https://use.fontawesome.com/releases/v5.13.0/css/all.css"
external_stylesheets = [dbc.themes.BOOTSTRAP, FONT_AWESOME]

app = Dash(__name__,
           use_pages=True,
           external_stylesheets=external_stylesheets,
           meta_tags=[{'name': 'viewport', 'content': 'width=device-width, initial-scale=1'}],
           suppress_callback_exceptions=True)
server = app.server

navbar_content = [
    dbc.NavItem(dbc.NavLink('Data overview', href='/', active='exact')),
    dbc.NavItem(dbc.NavLink('Cross-lag panel model (1)', href='/clpm', active='exact')),
    dbc.NavItem(dbc.NavLink('Cross-lag panel model (2)', href='/riclpm', active=False)),
    dbc.NavItem(dbc.NavLink('Cross-lag network model', href='/clnm', active='exact')),
    dbc.NavItem(dbc.NavLink('Cross-sectional network models', href='/csnm', active='exact'))
]

app.layout = html.Div([
    # TITLE
    html.H1(id='title',
            style=styles.TITLE,
            children='Longitudinal modelling of the co-development of depression and cardio-metabolic risk \
            from childhood to young adulthood'),
    html.Br(),

    # CONTENT
    dbc.Row(
        dbc.Col(width={'size': 10, 'offset': 1},
                children=[dbc.Nav(navbar_content, pills=True, fill=True,
                                  style=styles.NAVBAR)])
    ),
    html.Br(),

    page_container
])


if __name__ == '__main__':
    app.run_server(debug=True)
