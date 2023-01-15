import plotly.graph_objects as go

# create a nested pie chart to visualise the results from the meta-assembly
fig =go.Figure(go.Sunburst(
    labels=["Meta-assembly", "Transcripts<br>220'907", "Genes<br>62'314", "Novel Genes<br>5.89%", "Single Exon<br>13.29%", 
            "Novel transcripts<br>5.70%"],
    parents=["", "Meta-assembly", "Meta-assembly", "Genes<br>62'314", "Transcripts<br>220'907", "Transcripts<br>220'907" ],
    values=[100, 220907, 62314*2, 3669*2, 29354, 12597],
))
fig.update_layout(margin = dict(t=0, l=0, r=0, b=0))

fig.show()
