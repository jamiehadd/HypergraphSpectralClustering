figs: fig/4-heatmaps.png fig/algorithm-demo.png fig/contact-primary-school-classes-2.png fig/clustering-math.png

fig/4-heatmaps.png: throughput/binary-detection/exp-6.csv scripts/binary-detection-viz.R
	Rscript scripts/binary-detection-viz.R

fig/algorithm-demo.png: scripts/pipeline-viz.jl
	julia --project=. scripts/pipeline-viz.jl

throughput/binary-detection/exp-6.csv: scripts/binary-detection-experiments.jl
	julia --project=. scripts/binary-detection-experiments.jl

fig/contact-primary-school-classes-2.png: scripts/data-clustering-viz.jl throughput/data-throughput
	julia --project=. scripts/data-clustering-viz.jl

fig/clustering-math.png: scripts/mat-sx-viz.jl
	julia --project=. scripts/mat-sx-viz.jl
