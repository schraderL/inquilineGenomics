

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
revigo.data <- rbind(c("GO:0005085","guanyl-nucleotide exchange factor activity", 0.693,-0.120,-1.304, 1.892,-2.7986,0.908,0.000),
c("GO:0005158","insulin receptor binding", 0.126,-5.993, 5.128, 1.176,-1.7540,0.863,0.000),
c("GO:0005217","intracellular ligand-gated ion channel activity", 0.099,-0.477,-5.327, 1.079,-1.3617,0.843,0.000),
c("GO:0019787","ubiquitin-like protein transferase activity", 2.268, 6.830, 0.697, 2.403,-1.4531,0.735,0.000),
c("GO:0044877","macromolecular complex binding", 2.960,-3.974,-4.754, 2.519,-1.5462,0.848,0.071),
c("GO:0005539","glycosaminoglycan binding", 0.315,-2.523, 5.836, 1.556,-1.0372,0.800,0.078),
c("GO:0008135","translation factor activity, RNA binding", 0.675,-5.660, 0.987, 1.881,-1.2850,0.825,0.085),
c("GO:0004828","serine-tRNA ligase activity", 0.027, 1.135,-6.864, 0.602,-1.1903,0.825,0.101),
c("GO:0016854","racemase and epimerase activity", 0.081, 3.954,-5.046, 1.000,-1.1663,0.818,0.111),
c("GO:0016638","oxidoreductase activity, acting on the CH-NH2 group of donors", 0.108, 6.258,-3.627, 1.114,-1.2269,0.816,0.114),
c("GO:0016298","lipase activity", 0.576, 4.079, 4.423, 1.813,-1.3833,0.762,0.136),
c("GO:0005515","protein binding",20.274,-5.049,-2.532, 3.353,-2.5287,0.847,0.141),
c("GO:0016740","transferase activity",12.580, 4.866, 0.679, 3.146,-1.2287,0.775,0.207),
c("GO:0016801","hydrolase activity, acting on ether bonds", 0.072, 3.281, 6.428, 0.954,-1.1644,0.783,0.238),
c("GO:0101005","ubiquitinyl hydrolase activity", 0.360, 4.542, 5.471, 1.613,-1.0993,0.767,0.274),
c("GO:1901265","nucleoside phosphate binding", 9.574,-5.151,-0.250, 3.027,-2.0550,0.805,0.282),
c("GO:0003774","motor activity", 0.648, 5.256, 4.445, 1.863,-1.0796,0.761,0.290),
c("GO:0008201","heparin binding", 0.144,-2.130, 6.487, 1.230,-1.0602,0.806,0.338),
c("GO:0004672","protein kinase activity", 2.286, 6.544, 0.169, 2.407,-1.2741,0.709,0.479),
c("GO:0005524","ATP binding", 6.713,-3.758, 2.744, 2.873,-1.0777,0.742,0.492));

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

dev.print(pdf,"/Users/lukas/sciebo/inquilineGenomics18/GOenrichment/revigo/RELAX.BP.revigo.pdf",width=10,height=8);
