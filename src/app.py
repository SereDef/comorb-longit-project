from dash import Dash, html, page_container
import dash_bootstrap_components as dbc

from callbacks import clpm_callbacks, clnm_callbacks, csnm_callbacks

FONT_AWESOME = "https://use.fontawesome.com/releases/v5.13.0/css/all.css"
external_stylesheets = [dbc.themes.BOOTSTRAP, FONT_AWESOME]

app = Dash(__name__,
           use_pages=True,
           external_stylesheets=external_stylesheets,
           suppress_callback_exceptions=True)

app.layout = html.Div([
    # TITLE
    html.H1(id='title',
            style={'textAlign': 'center', 'font-weight': 'bold'},
            children='Longitudinal modelling of the co-development of depression and cardio-metabolic risk \
            from childhood to young adulthood'),

    # CONTENT
    dbc.Row(
        style={'textAlign': 'center'},
        children=[dbc.Col(dbc.NavLink('Data overview', href='/')),
                  dbc.Col(dbc.NavLink('Cross-lag panel model', href='/clpm')),
                  dbc.Col(dbc.NavLink('Cross-lag network model', href='/clnm')),
                  dbc.Col(dbc.NavLink('Cross-sectional network model', href='/csnm'))]
    ),
    page_container
])


if __name__ == '__main__':
    app.run_server(debug=True)
