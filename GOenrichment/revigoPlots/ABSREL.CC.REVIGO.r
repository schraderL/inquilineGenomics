

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
revigo.data <- rbind(c("GO:0005576","extracellular region",10.049, 3.174, 5.427, 3.046,-1.3059,0.903,0.000),
c("GO:0030496","midbody", 0.235,-0.512, 4.930, 1.431,-1.7373,0.848,0.000),
c("GO:0032991","macromolecular complex",24.860,-4.991, 0.314, 3.439,-1.0562,0.918,0.000),
c("GO:0000346","transcription export complex", 0.072, 5.740,-3.195, 0.954,-1.5574,0.679,0.023),
c("GO:0042995","cell projection", 4.508,-5.680,-2.971, 2.698,-1.4401,0.837,0.034),
c("GO:0038201","TOR complex", 0.063,-3.323, 3.062, 0.903,-1.4602,0.836,0.038),
c("GO:0044463","cell projection part", 2.345, 7.152, 2.706, 2.415,-1.4401,0.665,0.045),
c("GO:0043228","non-membrane-bounded organelle",15.526,-2.692,-6.851, 3.235,-1.8636,0.765,0.095),
c("GO:0001518","voltage-gated sodium channel complex", 0.018, 2.021, 0.058, 0.477,-1.0969,0.785,0.132),
c("GO:0031932","TORC2 complex", 0.036, 6.508,-5.235, 0.699,-1.2565,0.736,0.206),
c("GO:0000347","THO complex", 0.054, 6.939,-3.448, 0.845,-1.5574,0.683,0.274),
c("GO:0005875","microtubule associated complex", 3.377, 2.836,-5.105, 2.573,-1.4598,0.485,0.294),
c("GO:0034464","BBSome", 0.054, 6.437,-0.739, 0.845,-1.5519,0.561,0.525),
c("GO:0016459","myosin complex", 0.253, 3.324,-5.557, 1.462,-1.2248,0.564,0.575),
c("GO:0043232","intracellular non-membrane-bounded organelle",15.526, 0.153,-7.432, 3.235,-2.1918,0.582,0.599),
c("GO:0005856","cytoskeleton", 6.391,-0.002,-7.122, 2.849,-1.1034,0.605,0.682));

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
dev.print(pdf,"/Users/lukas/sciebo/inquilineGenomics18/GOenrichment/revigo/RELAX.CC.intensified.revigo.pdf",width=10,height=8);