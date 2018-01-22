using Plots
using StatPlots

function nice_plots(df)
    c = cor(df[:sttr], df[:mtld])
    @df df scatter(:sttr, :mtld;
                   xlab="STTR", ylab="MTLD",
                   title=string("Correlation = ", c))
end
