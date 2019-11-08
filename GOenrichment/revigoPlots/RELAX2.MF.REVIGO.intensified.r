

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
revigo.data <- rbind(c("GO:0000287","magnesium ion binding", 1.785, 4.238,-3.495, 5.400,-1.6364,0.796,0.000),
c("GO:0004001","adenosine kinase activity", 0.008,-5.014, 4.789, 3.060,-2.3279,0.570,0.000),
c("GO:0004427","inorganic diphosphatase activity", 0.064, 3.605, 6.126, 3.954,-1.7077,0.677,0.017),
c("GO:0061783","peptidoglycan muralytic activity", 0.218, 0.191,-4.794, 4.486,-1.3206,0.778,0.021),
c("GO:0003676","nucleic acid binding",21.226,-4.102,-4.279, 6.475,-1.7011,0.832,0.093),
c("GO:0003796","lysozyme activity", 0.045, 5.185, 1.859, 3.806,-1.6635,0.706,0.155),
c("GO:0017171","serine hydrolase activity", 1.238, 4.786, 4.062, 5.241,-1.0443,0.689,0.199),
c("GO:0003998","acylphosphatase activity", 0.024, 2.705, 7.121, 3.528,-1.5884,0.683,0.311),
c("GO:0016779","nucleotidyltransferase activity", 1.954,-5.103, 2.596, 5.439,-1.2441,0.589,0.338),
c("GO:0019205","nucleobase-containing compound kinase activity", 0.345,-5.804, 3.027, 4.686,-1.2366,0.587,0.470),
c("GO:0019206","nucleoside kinase activity", 0.065,-5.126, 4.153, 3.962,-1.1141,0.550,0.675));

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
dev.print(pdf,"/Users/lukas/sciebo/inquilineGenomics18/GOenrichment/revigo/RELAX.MF.intensified.revigo.pdf",width=10,height=8);