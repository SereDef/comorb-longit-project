from dash import Dash, html, page_container
import dash_bootstrap_components as dbc

from callbacks import clpm_callbacks, clnm_callbacks, csnm_callbacks

FONT_AWESOME = "https://use.fontawesome.com/releases/v5.13.0/css/all.css"
external_stylesheets = [dbc.themes.BOOTSTRAP, FONT_AWESOME]

app = Dash(__name__,
           use_pages=True,
           external_stylesheets=external_stylesheets,
           suppress_callback_exceptions=True)

navbar_content = [
    dbc.NavItem(dbc.NavLink('Data overview', href='/', active='exact')),
    dbc.NavItem(dbc.NavLink('Cross-lag panel model', href='/clpm', active='exact')),
    dbc.NavItem(dbc.NavLink('Cross-lag network model', href='/clnm', active='exact')),
    dbc.NavItem(dbc.NavLink('Cross-sectional network model', href='/csnm', active='exact'))
]

app.layout = html.Div([
    # TITLE
    html.H1(id='title',
            style={'textAlign': 'center', 'font-weight': 'bold'},
            children='Longitudinal modelling of the co-development of depression and cardio-metabolic risk \
            from childhood to young adulthood'),
    html.Br(),

    # CONTENT
    dbc.Row(
        dbc.Col(width={'size': 10, 'offset': 1},
                children=[dbc.Nav(navbar_content, pills=True, fill=True)])
    ),
    html.Br(),

    page_container
])


if __name__ == '__main__':
    app.run_server(debug=True)
