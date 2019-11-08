

# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800


# --------------------------------------------------------------------------
# If you don't have the ggplot2 package installed, uncomment the following line:
# install.packages( "ggplot2" );
library( ggplot2 );
# --------------------------------------------------------------------------
# If you don't have the scales package installed, uncomment the following line:
# install.packages( "scales" );
library( scales );


# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0031929","TOR signaling", 0.065,-5.033,-3.146, 3.922,-2.5850,0.657,0.000),
c("GO:0043094","cellular metabolic compound salvage", 0.314, 4.477,-5.405, 4.605,-1.3635,0.770,0.025),
c("GO:0008213","protein alkylation", 0.343, 0.426, 6.368, 4.643,-1.1062,0.759,0.052),
c("GO:0015748","organophosphate ester transport", 0.144,-5.034, 3.704, 4.268,-1.4976,0.771,0.059),
c("GO:0046834","lipid phosphorylation", 0.177, 0.013,-6.921, 4.357,-1.0691,0.764,0.063),
c("GO:0043101","purine-containing compound salvage", 0.132, 6.218,-1.374, 4.230,-1.3536,0.747,0.136),
c("GO:0006265","DNA topological change", 0.268, 4.278, 2.615, 4.536,-1.2790,0.667,0.155),
c("GO:0044271","cellular nitrogen compound biosynthetic process",22.502, 5.123,-2.746, 6.460,-1.0438,0.741,0.232),
c("GO:0048585","negative regulation of response to stimulus", 0.344,-5.289,-2.018, 4.645,-1.0742,0.650,0.284),
c("GO:0038179","neurotrophin signaling pathway", 0.007,-4.380,-4.020, 2.963,-1.7721,0.689,0.300),
c("GO:0034404","nucleobase-containing small molecule biosynthetic process", 1.408, 3.601,-2.839, 5.257,-1.3468,0.672,0.303),
c("GO:0016570","histone modification", 0.373, 2.608, 5.237, 4.680,-1.0768,0.684,0.423),
c("GO:0006869","lipid transport", 0.270,-4.618, 4.283, 4.539,-1.0255,0.768,0.458),
c("GO:0007033","vacuole organization", 0.102, 4.526, 4.616, 4.119,-1.0039,0.737,0.501));

one.data <- data.frame(revigo.data);
names(one.data) <- revigo.names;
one.data <- one.data [(one.data$plot_X != "null" & one.data$plot_Y != "null"), ];
one.data$plot_X <- as.numeric( as.character(one.data$plot_X) );
one.data$plot_Y <- as.numeric( as.character(one.data$plot_Y) );
one.data$plot_size <- as.numeric( as.character(one.data$plot_size) );
one.data$log10_p_value <- as.numeric( as.character(one.data$log10_p_value) );
one.data$frequency <- as.numeric( as.character(one.data$frequency) );
one.data$uniqueness <- as.numeric( as.character(one.data$uniqueness) );
one.data$dispensability <- as.numeric( as.character(one.data$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

p1 <- ggplot( data = one.data );
p1 <- p1 + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) + scale_size_area();
p1 <- p1 + scale_colour_gradientn( colours = c(goodCols[1], goodCols[2], goodCols[3], goodCols[5]), limits = c( min(one.data$log10_p_value), 0) );
p1 <- p1 + geom_point( aes(plot_X, plot_Y, size = plot_size), shape = 21, fill = "transparent", colour = I (alpha ("black", 0.6) )) + scale_size_area();
p1 <- p1 + scale_size( range=c(5, 30)) + theme_bw(); # + scale_fill_gradientn(colours = heat_hcl(7), limits = c(-300, 0) );
ex <- one.data [ one.data$dispensability < 0.15, ]; 
p1 <- p1 + geom_label( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 4 );
p1 <- p1 + labs (y = "semantic space x", x = "semantic space y");
p1 <- p1 + theme(legend.key = element_blank()) ;
one.x_range = max(one.data$plot_X) - min(one.data$plot_X);
one.y_range = max(one.data$plot_Y) - min(one.data$plot_Y);
p1 <- p1 + xlim(min(one.data$plot_X)-one.x_range/10,max(one.data$plot_X)+one.x_range/10);
p1 <- p1 + ylim(min(one.data$plot_Y)-one.y_range/10,max(one.data$plot_Y)+one.y_range/10);



# --------------------------------------------------------------------------
# Output the plot to screen

p1;

# Uncomment the line below to also save the plot to a file.
# The file type depends on the extension (default=pdf).

# ggsave("C:/Users/path_to_your_file/revigo-plot.pdf");
dev.print(pdf,"/Users/lukas/sciebo/inquilineGenomics18/GOenrichment/revigo/RELAX.BP.intensified.revigo.pdf",width=10,height=8);