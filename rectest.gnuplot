set data style xyerrorbars
set key right bottom
set xlabel "Empirical frequency that reconstruction is correct"
set ylabel "Posterior probability of ML reconstruction (to one decimal place)"
set terminal postscript landscape color

set xrange [0:1]
set yrange [0:1]

set title "Simulated data"

set output "figs/hammer_paired.eps"
plot x t "", "Hammerhead_1.gpdata" index 0 t "Covariant base-pair model", "Hammerhead_1.gpdata" index 1 t "Non-covariant point substitution model"

set output "figs/hammer_unpaired.eps"
plot x t "", "Hammerhead_1.gpdata" index 2 t "PFOLD single-stranded model", "Hammerhead_1.gpdata" index 3 t "General reversible model"


set title "Cross-validation on RFAM data"

set output "figs/xval_paired.eps"
plot x t "", "Published.xval-gpdata" index 0 t "Covariant base-pair model", "Published.xval-gpdata" index 1 t "Non-covariant point substitution model"

set output "figs/xval_unpaired.eps"
plot x t "", "Published.xval-gpdata" index 2 t "PFOLD single-stranded model", "Published.xval-gpdata" index 3 t "General reversible model"


set title "RFAM cross-validation (rates estimated from RFAM)"

set output "figs/pubxval_paired.eps"
plot x t "", "Published.pubxval-gpdata" index 0 t "Covariant base-pair model", "Published.pubxval-gpdata" index 1 t "Non-covariant point substitution model"

set output "figs/pubxval_unpaired.eps"
plot x t "", "Published.pubxval-gpdata" index 2 t "PFOLD single-stranded model", "Published.pubxval-gpdata" index 3 t "General reversible model"




set data style histogram
set style fill solid

set key right top
set xlabel "Posterior probability of ML reconstruction"
set ylabel "Number of basepair reconstructions"

set xrange [-.5:5.75]

set title "Simulated data"
set yrange [0:7000]
set output "figs/hammer_hist.eps"

plot "Hammerhead_1.histdata" using 2 t "Covariant model, incorrect" ls 1, "Hammerhead_1.histdata" using 3 t "Non-covariant model, incorrect" ls 3, "Hammerhead_1.histdata" using 4 t "Covariant model, incorrect AND non-canonical" ls 4, "Hammerhead_1.histdata" using 5:xticlabels(1) t "Non-covariant model, incorrect AND non-canonical" ls 5

set title "Cross-validation on RFAM data"
set yrange [0:2500]
set output "figs/xval_hist.eps"

plot "Published.xval-histdata" using 2 t "Covariant model, incorrect" ls 1, "Published.xval-histdata" using 3 t "Non-covariant model, incorrect" ls 3, "Published.xval-histdata" using 4 t "Covariant model, incorrect AND non-canonical" ls 4, "Published.xval-histdata" using 5:xticlabels(1) t "Non-covariant model, incorrect AND non-canonical" ls 5

set title "RFAM cross-validation (rates estimated from RFAM)"
set yrange [0:2500]
set output "figs/pubxval_hist.eps"

plot "Published.pubxval-histdata" using 2 t "Covariant model, incorrect" ls 1, "Published.pubxval-histdata" using 3 t "Non-covariant model, incorrect" ls 3, "Published.pubxval-histdata" using 4 t "Covariant model, incorrect AND non-canonical" ls 4, "Published.pubxval-histdata" using 5:xticlabels(1) t "Non-covariant model, incorrect AND non-canonical" ls 5
