

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
revigo.data <- rbind(c("GO:0006022","aminoglycan metabolic process", 1.257,-2.460,-6.355, 2.158,-1.7799,0.784,0.000),
c("GO:0007017","microtubule-based process", 4.404, 3.005, 6.458, 2.701,-1.5935,0.877,0.000),
c("GO:0023052","signaling",18.423, 3.067, 5.458, 3.322,-1.6968,0.962,0.000),
c("GO:0030705","cytoskeleton-dependent intracellular transport", 0.580,-1.318, 5.422, 1.826,-1.4134,0.913,0.000),
c("GO:0050896","response to stimulus",27.608,-4.502, 0.583, 3.497,-1.1624,0.967,0.000),
c("GO:0007265","Ras protein signal transduction", 1.600, 6.657,-0.783, 2.262,-1.5784,0.616,0.058),
c("GO:0007154","cell communication",19.021,-5.850, 1.549, 3.335,-1.6882,0.891,0.064),
c("GO:0017144","drug metabolic process", 5.810,-0.952,-3.988, 2.821,-1.2733,0.804,0.087),
c("GO:0006928","movement of cell or subcellular component", 5.968,-2.512, 2.199, 2.833,-2.2076,0.874,0.163),
c("GO:0006383","transcription from RNA polymerase III promoter", 0.246, 1.670,-7.623, 1.462,-1.7570,0.716,0.184),
c("GO:0016226","iron-sulfur cluster assembly", 0.132,-3.815,-2.428, 1.204,-1.1537,0.783,0.199),
c("GO:0019538","protein metabolic process",22.756,-0.314,-6.667, 3.413,-1.8894,0.790,0.208),
c("GO:0099111","microtubule-based transport", 0.554,-1.723, 4.618, 1.806,-1.3242,0.850,0.281),
c("GO:0048583","regulation of response to stimulus", 8.183, 7.279,-1.112, 2.969,-2.1427,0.679,0.281),
c("GO:1901564","organonitrogen compound metabolic process",10.864, 0.061,-7.642, 3.092,-2.0223,0.777,0.290),
c("GO:0080135","regulation of cellular response to stress", 0.958, 6.268, 0.407, 2.041,-1.4473,0.654,0.296),
c("GO:0006359","regulation of transcription from RNA polymerase III promoter", 0.079, 3.555,-6.077, 1.000,-1.9101,0.635,0.313),
c("GO:0031163","metallo-sulfur cluster assembly", 0.132,-5.678,-1.760, 1.204,-1.1831,0.904,0.317),
c("GO:0043170","macromolecule metabolic process",37.690,-1.329,-6.571, 3.632,-1.0737,0.822,0.321),
c("GO:0015748","organophosphate ester transport", 0.202,-0.643, 5.422, 1.380,-1.2916,0.892,0.360),
c("GO:0055085","transmembrane transport", 7.032,-1.790, 5.160, 2.904,-1.0975,0.909,0.372),
c("GO:0023014","signal transduction by protein phosphorylation", 1.503, 4.578,-2.001, 2.236,-1.2503,0.578,0.381),
c("GO:0010646","regulation of cell communication", 7.251, 6.311,-2.691, 2.917,-1.8894,0.639,0.398),
c("GO:0023051","regulation of signaling", 7.278, 7.010,-2.654, 2.919,-1.8996,0.661,0.409),
c("GO:0045892","negative regulation of transcription, DNA-templated", 2.918, 3.014,-4.519, 2.522,-1.4318,0.435,0.431),
c("GO:0051716","cellular response to stimulus",19.768, 6.631, 2.008, 3.352,-1.2668,0.743,0.452),
c("GO:0035556","intracellular signal transduction", 5.502, 6.412,-1.442, 2.797,-1.1925,0.597,0.500),
c("GO:0051174","regulation of phosphorus metabolic process", 2.646, 5.031,-4.019, 2.480,-1.0830,0.594,0.526),
c("GO:0080134","regulation of response to stress", 2.312, 6.172, 0.737, 2.422,-1.0386,0.669,0.546),
c("GO:0031098","stress-activated protein kinase signaling cascade", 0.668, 6.438,-0.183, 1.886,-1.4609,0.626,0.599),
c("GO:0006040","amino sugar metabolic process", 1.064,-4.486,-5.376, 2.086,-1.5969,0.846,0.608),
c("GO:0007165","signal transduction",15.338, 6.101,-1.317, 3.242,-1.7721,0.559,0.614));

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
dev.print(pdf,"/Users/lukas/sciebo/inquilineGenomics18/GOenrichment/revigo/ABSREL.BP.revigo.pdf",width=10,height=8);