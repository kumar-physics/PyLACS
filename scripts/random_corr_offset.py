import random
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import statsmodels.api as sm
import numpy as np

# Generate N random float numbers between -5 and +5
N = 50  # change if needed
numbers = [random.uniform(0, 5) for _ in range(N)]

# Compute sum and diff only for consecutive pairs (i, i+1)
sums = []
diffs = []
for i in range(N - 1):  # up to N-2
    s = numbers[i] + numbers[i + 1]
    d = numbers[i] - numbers[i + 1]
    sums.append(s)
    diffs.append(d)

# Create DataFrame for plotting
df = pd.DataFrame({
    "Difference": diffs,
    "Sum": sums
})

# Robust linear regression: y ~ x
X = sm.add_constant(df["Difference"])  # add intercept
y = df["Sum"]

rlm_model = sm.RLM(y, X, M=sm.robust.norms.HuberT())
rlm_results = rlm_model.fit()

slope = rlm_results.params["Difference"]
intercept = rlm_results.params["const"]

print("Robust Linear Fit:")
print(f"y = {slope:.3f} * x + {intercept:.3f}")

# Generate fitted line for plotting
x_fit = np.linspace(df["Difference"].min(), df["Difference"].max(), 200)
y_fit = slope * x_fit + intercept

# Scatter plot
fig = px.scatter(df, x="Difference", y="Sum",
                 title="Sum vs Difference for Consecutive Pairs (Robust Fit)",
                 opacity=0.7)

# Add regression line
fig.add_trace(go.Scatter(x=x_fit, y=y_fit, mode="lines",
                         name="Robust Linear Fit", line=dict(color="red")))

fig.show()
