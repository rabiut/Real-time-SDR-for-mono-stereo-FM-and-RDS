# example.gnuplut : configuration for plotting (change as needed)

reset                                   # reset
set size ratio 0.2                      # set relative size of plots
set grid xtics ytics                    # grid: enable both x and y lines
set grid lt 1 lc rgb '#cccccc' lw 1     # grid: thin gray lines
set multiplot layout 3,1 scale 1.0,1.0  # set two plots for this figure

# time domain
set ylabel 'Q'               # set y-axis label
set xlabel 'I'                   # set x-axis label
                     # set x plot range
plot '../data/IQ.dat' with points pt 7 ps 0.2

# freq domain (Fourier)
set ylabel 'Amplitude'              # set y-axis label
set xlabel 'Time (s)'               # set x-axis label
set yrange [-1.0:1.0]                      # set y plot range
set xrange [0.006:0.012]                       # set x plot range
plot '../data/signal3.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#0000ff' title 'In-Phase','../data/signal1.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#ff0000' title 'Quadrature','../data/signal2.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#00ff00' title 'Pre RRC filter'

# freq domain (PSD)
set ylabel 'Bit'            # set y-axis label
set xlabel 'bit #'             # set x-axis label                    # set y plot range
set yrange [0.0:1.0]
set xrange [0.000:50]                       # set x plot range
# add your own .dat file for PSD as part of the take-home
plot '../data/symbPairing.dat' using 1:2 with lines lt 1 lw 2 lc rgb '#0000ff' notitle



unset multiplot
