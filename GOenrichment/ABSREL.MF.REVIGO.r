

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
revigo.data <- rbind(c("GO:0003774","motor activity", 0.648, 6.224, 0.715, 1.863,-1.6383,0.801,0.000),
c("GO:0005085","guanyl-nucleotide exchange factor activity", 0.693, 5.515,-3.674, 1.892,-1.1367,0.898,0.000),
c("GO:0005215","transporter activity", 7.730,-6.273,-3.065, 2.934,-1.3768,0.905,0.000),
c("GO:0005488","binding",49.132, 2.742,-5.608, 3.737,-1.1739,0.947,0.000),
c("GO:0038024","cargo receptor activity", 0.144, 1.209,-1.979, 1.230,-1.0223,0.898,0.000),
c("GO:0043395","heparan sulfate proteoglycan binding", 0.054,-5.519, 3.773, 0.845,-1.3665,0.694,0.000),
c("GO:0020037","heme binding", 1.197,-1.434,-5.369, 2.127,-1.2441,0.821,0.061),
c("GO:0008276","protein methyltransferase activity", 0.441, 3.574, 5.981, 1.699,-1.5528,0.719,0.116),
c("GO:0005515","protein binding",20.274,-3.570,-0.858, 3.353,-1.4559,0.833,0.120),
c("GO:0016787","hydrolase activity",16.917, 5.301, 3.343, 3.274,-1.1675,0.807,0.178),
c("GO:0019899","enzyme binding", 3.158,-3.667, 4.870, 2.547,-1.1367,0.761,0.232),
c("GO:0046906","tetrapyrrole binding", 1.206,-2.588,-4.964, 2.130,-1.0655,0.821,0.234),
c("GO:0042393","histone binding", 0.351,-3.115, 6.728, 1.602,-1.0915,0.781,0.279),
c("GO:0004707","MAP kinase activity", 0.063, 2.687, 7.031, 0.903,-1.0410,0.780,0.283),
c("GO:0017137","Rab GTPase binding", 0.585,-4.205, 6.049, 1.820,-1.1024,0.776,0.296),
c("GO:0008061","chitin binding", 1.044,-6.806, 1.410, 2.068,-1.2676,0.783,0.344),
c("GO:0016741","transferase activity, transferring one-carbon groups", 1.197, 3.927, 5.241, 2.127,-1.1938,0.752,0.368),
c("GO:0016817","hydrolase activity, acting on acid anhydrides", 5.102, 5.692, 1.921, 2.754,-1.5376,0.787,0.368));

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
dev.print(pdf,"/Users/lukas/sciebo/inquilineGenomics18/GOenrichment/revigo/ABSREL.MF.revigo.pdf",width=10,height=8);