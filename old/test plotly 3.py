import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px




x = np.linspace(0,10,25)
y = np.sin(x)
y[-1] = -200
cond = y > -25


fig = px.scatter(x[cond],y[cond])

xlim = fig.update_xaxes()
ylim = fig.update_yaxes()



fig = px.scatter(x[cond == False], y[cond == False])
fig.update_xaxes(xlim)
fig.update_yaxes(ylim)

fig.show()


