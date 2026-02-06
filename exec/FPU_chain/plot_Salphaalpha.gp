# Set output format and file
set terminal pdf font 'Helvetica,12'
set output 'plots.pdf'

# Set multiplot configuration for 2x3 layout
set multiplot layout 2, 3 title "Plots of S_alphaalpha_avg1000000" font 'Helvetica,14'

# Use larger font for titles and labels
set label font 'Helvetica,10'
set key font 'Helvetica,10'
set xlabel font 'Helvetica,10'
set ylabel font 'Helvetica,10'
set title font 'Helvetica,12'

# Looping over rows 5 to 10 for y-axis
do for [i=5:10] {
    # Set individual plot titles
    set title sprintf("Plot for Column %d", i)

    # Set x-axis tics to show only first and last points
    stats 'S_alphaalpha_avg1000000' using 4 nooutput
    set xtics (STATS_min, STATS_max)

    # Plot with row 4 as x and current row i as y using lines only
    plot 'S_alphaalpha_avg1000000' using 4:i title sprintf("Column %d", i) with lines lc rgb 'purple'
}

# End multiplot
unset multiplot
