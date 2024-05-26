from plotly.subplots import make_subplots
import plotly.graph_objects as go
import numpy as np


r = np.linspace(0, 10, 100)

def func1(r):
    return r**1

def func2(r):
    return r**2
    
def func3(r):
    return r**3

def func4(r):
    return r**4

def func5(r):
    return r**5

def func6(r):
    return r**6


fig = make_subplots(rows=1, cols=3, subplot_titles=("(a)energias simples", "(b)energias complexas", "(c) Segundo Coeficiente Virial(B) em função da Temperatura(T)"))






fig.add_trace(
go.Scatter(x=r, y=func1(r)),
row=1, col=1
)
fig.add_trace(
go.Scatter(x=r, y=func2(r)),
row=1, col=1
)
fig.add_trace(
go.Scatter(x=r, y=func3(r)),
row=1, col=1
)
fig.add_trace(
go.Scatter(x=r, y=func4(r)),
row=1, col=1
)
fig.add_trace(
go.Scatter(x=r, y=func5(r)),
row=1, col=1
)
fig.add_trace(
go.Scatter(x=r, y=func6(r)),
row=1, col=1
)



fig.add_trace(
go.Scatter(x=r, y=func1(r)),
row=1, col=2
)
fig.add_trace(
go.Scatter(x=r, y=func2(r)),
row=1, col=2
)
fig.add_trace(
go.Scatter(x=r, y=func3(r)),
row=1, col=2
)
fig.add_trace(
go.Scatter(x=r, y=func4(r)),
row=1, col=2
)
fig.add_trace(
go.Scatter(x=r, y=func5(r)),
row=1, col=2
)
fig.add_trace(
go.Scatter(x=r, y=func6(r)),
row=1, col=2
)




fig.add_trace(
go.Scatter(x=r, y=func1(r), name= 'B sem c2'),
row=1, col=3
)
fig.add_trace(
go.Scatter(x=r, y=func2(r)),
row=1, col=3
)
fig.add_trace(
go.Scatter(x=r, y=func3(r)),
row=1, col=3
)

fig.add_trace(
go.Scatter(x=r, y=func4(r)),
row=1, col=3
)


fig.update_layout(height=600, width=800, title_text="Side By Side Subplots")



fig.update_yaxes(title_text="yaxis 1 title", row=1, col=1)
fig.update_yaxes(title_text="yaxis 2 title", range=[0, 100], row=1, col=2)
fig.update_yaxes(title_text="yaxis 3 title", showgrid=False, row=2, col=1)
fig.update_yaxes(title_text="yaxis 4 title", row=2, col=2)

fig.update_xaxes(title_text="yaxis 3 title", range=[0, 5], showgrid=False, row=1, col=3)



fig.update_layout(height=800, width=600, template='plotly_white')
fig.show()





