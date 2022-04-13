.PHONY : figs
.PHONY : clean

JULIA = julia --project=.

figs: fig/4-heatmaps.png fig/algorithm-demo.png fig/contact-primary-school-classes-1.png fig/clustering-math.png

clean: 
	rm -f fig/SN-congress-bills-1.png
	rm -f fig/SN-congress-bills-2.png
	rm -f fig/clustering-math.png
	rm -f fig/algorithm-demo.png 
	rm -f fig/4-heatmaps.png
	rm -f fig/contact-primary-school-classes-1.png
	rm -f fig/contact-primary-school-classes-2.png
	rm -f fig/contact-high-school-classes-1.png
	rm -f fig/contact-high-school-classes-2.png

fig/4-heatmaps.png: scripts/binary-detection-viz.R throughput/binary-detection/exp-6.csv throughput/binary-detection/exp-vanilla.csv
	Rscript $<

fig/algorithm-demo.png: scripts/pipeline-viz.jl
	$(JULIA) $<

throughput/bulk-throughput/exp-6.csv: scripts/binary-detection-experiments.jl
	$(JULIA) $<

throughput/bulk-throughput/exp-vanilla.csv: scripts/binary-detection-vanilla.jl
	$(JULIA) $<

fig/contact-primary-school-classes-1.png: scripts/data-clustering-viz.jl throughput/data-throughput
	$(JULIA) $<

fig/clustering-math.png: scripts/math-sx-viz.jl
	$(JULIA) $<
