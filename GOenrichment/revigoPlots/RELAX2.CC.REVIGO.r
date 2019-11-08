

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
source("~/sciebo/librarySchrader.R")

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0005623","cell",65.001, 3.798,-3.825, 3.856,-1.0560,0.947,0.000),
c("GO:0016459","myosin complex", 0.253,-3.270, 0.941, 1.462,-1.9626,0.584,0.000),
c("GO:0043226","organelle",45.392, 0.494,-6.545, 3.700,-2.1871,0.918,0.000),
c("GO:0044424","intracellular part",54.273, 4.570, 2.924, 3.778,-1.4437,0.785,0.078),
c("GO:0044422","organelle part",27.186,-1.520, 6.985, 3.478,-1.2958,0.752,0.118),
c("GO:0005942","phosphatidylinositol 3-kinase complex", 0.072,-4.309,-4.725, 0.954,-1.3778,0.750,0.158),
c("GO:0005774","vacuolar membrane", 0.616,-5.241, 4.466, 1.839,-1.1759,0.561,0.175),
c("GO:0071203","WASH complex", 0.063,-6.216,-2.338, 0.903,-1.1141,0.665,0.234),
c("GO:0005643","nuclear pore", 0.425,-5.855, 0.980, 1.681,-1.3161,0.547,0.273),
c("GO:0044464","cell part",64.711, 5.156, 1.133, 3.854,-1.0560,0.839,0.288),
c("GO:0016592","mediator complex", 0.299,-6.254, 0.136, 1.531,-1.2299,0.578,0.358),
c("GO:0005856","cytoskeleton", 6.391, 0.622, 5.159, 2.849,-1.7696,0.668,0.361),
c("GO:0005635","nuclear envelope", 1.086,-6.255, 3.146, 2.083,-1.0022,0.585,0.611));

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
dev.print(pdf,"/Users/lukas/sciebo/inquilineGenomics18/GOenrichment/revigo/RELAX.CC.revigo.pdf",width=10,height=8);