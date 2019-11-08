

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
revigo.data <- rbind(c("GO:0009987","cellular process",64.762,-1.238,-3.584, 3.867,-1.3270,0.990,0.000),
c("GO:0023052","signaling",18.423,-1.395, 1.649, 3.322,-1.1029,0.976,0.000),
c("GO:0040011","locomotion", 5.704,-0.546, 3.290, 2.813,-1.0590,0.972,0.000),
c("GO:0046903","secretion", 2.259,-5.933,-0.053, 2.412,-2.2007,0.726,0.000),
c("GO:0065007","biological regulation",38.323,-0.613, 5.968, 3.640,-2.1487,0.982,0.000),
c("GO:0070646","protein modification by small protein removal", 0.554, 5.914,-3.623, 1.806,-1.7352,0.828,0.000),
c("GO:0097485","neuron projection guidance", 2.312,-4.374,-5.764, 2.422,-1.4145,0.834,0.056),
c("GO:0048167","regulation of synaptic plasticity", 0.220,-1.257,-7.059, 1.415,-1.1904,0.781,0.098),
c("GO:0006629","lipid metabolic process", 4.140, 3.576, 3.536, 2.674,-1.4389,0.799,0.108),
c("GO:0006793","phosphorus metabolic process",10.117, 7.161, 0.489, 3.061,-2.3468,0.816,0.135),
c("GO:0050789","regulation of biological process",34.895, 0.238,-6.388, 3.599,-2.4202,0.821,0.165),
c("GO:0043412","macromolecule modification",13.501, 6.484,-3.376, 3.187,-1.8508,0.827,0.169),
c("GO:0016310","phosphorylation", 6.082, 6.341, 2.461, 2.841,-1.1308,0.818,0.198),
c("GO:0071704","organic substance metabolic process",47.288, 7.000,-1.782, 3.731,-1.5800,0.844,0.228),
c("GO:0008614","pyridoxine metabolic process", 0.026, 3.833, 4.607, 0.602,-1.3904,0.748,0.244),
c("GO:0031503","protein complex localization", 0.281,-5.456, 2.720, 1.519,-1.0306,0.779,0.251),
c("GO:0006434","seryl-tRNA aminoacylation", 0.026, 3.751, 1.326, 0.602,-1.1759,0.757,0.271),
c("GO:0006766","vitamin metabolic process", 0.105, 2.524, 6.199, 1.114,-1.0872,0.843,0.295),
c("GO:0035556","intracellular signal transduction", 5.502,-1.136,-6.475, 2.797,-2.6198,0.742,0.306),
c("GO:0048869","cellular developmental process",21.350,-1.428, 6.568, 3.386,-1.5560,0.851,0.308),
c("GO:0019216","regulation of lipid metabolic process", 0.343, 2.338,-4.591, 1.602,-1.2941,0.743,0.314),
c("GO:0051168","nuclear export", 0.615,-5.542, 1.373, 1.851,-1.1302,0.755,0.327),
c("GO:0006730","one-carbon metabolic process", 0.211, 3.047, 5.370, 1.398,-1.0110,0.801,0.343),
c("GO:0032879","regulation of localization", 3.780,-3.645,-2.270, 2.634,-1.2993,0.708,0.343),
c("GO:0009605","response to external stimulus", 8.596,-4.324,-7.011, 2.991,-1.0097,0.945,0.348),
c("GO:0019219","regulation of nucleobase-containing compound metabolic process",10.961, 3.157,-2.745, 3.096,-2.5528,0.603,0.379),
c("GO:0044260","cellular macromolecule metabolic process",31.942, 5.837,-1.309, 3.561,-1.4283,0.750,0.387),
c("GO:0051049","regulation of transport", 2.558,-3.953,-1.844, 2.465,-1.7447,0.639,0.389),
c("GO:0051641","cellular localization", 7.709,-5.759, 0.484, 2.943,-1.3556,0.769,0.415),
c("GO:0009914","hormone transport", 0.237,-4.940,-1.412, 1.447,-1.6108,0.714,0.417),
c("GO:0016192","vesicle-mediated transport", 5.661,-5.917, 0.272, 2.810,-1.2321,0.749,0.446),
c("GO:0071166","ribonucleoprotein complex localization", 0.457,-5.192, 2.240, 1.724,-1.0092,0.796,0.451),
c("GO:0051236","establishment of RNA localization", 0.888,-5.972, 1.347, 2.009,-1.0255,0.745,0.452),
c("GO:0036211","protein modification process",12.490, 6.078,-2.701, 3.153,-2.0000,0.800,0.476),
c("GO:0010646","regulation of cell communication", 7.251, 0.432,-6.972, 2.917,-1.1024,0.774,0.487),
c("GO:0032989","cellular component morphogenesis", 7.568,-1.883, 6.683, 2.936,-1.0921,0.859,0.497),
c("GO:0023051","regulation of signaling", 7.278,-0.202,-6.865, 2.919,-1.0640,0.808,0.500),
c("GO:0044271","cellular nitrogen compound biosynthetic process",16.208, 5.542, 0.364, 3.266,-1.4401,0.699,0.522),
c("GO:0032940","secretion by cell", 2.057,-4.982, 1.528, 2.371,-1.7399,0.645,0.535),
c("GO:0050794","regulation of cellular process",32.557, 0.571,-6.197, 3.569,-1.6383,0.757,0.544),
c("GO:0008654","phospholipid biosynthetic process", 0.650, 4.420, 3.043, 1.875,-1.1772,0.731,0.570),
c("GO:0051649","establishment of localization in cell", 6.206,-5.579, 0.811, 2.849,-1.0414,0.726,0.601),
c("GO:0016070","RNA metabolic process",15.575, 5.595,-0.491, 3.249,-1.1261,0.693,0.614),
c("GO:0090087","regulation of peptide transport", 0.202,-4.708,-1.879, 1.380,-1.4283,0.669,0.639),
c("GO:0018130","heterocycle biosynthetic process",12.516, 5.636, 0.850, 3.154,-1.1457,0.703,0.645),
c("GO:0019438","aromatic compound biosynthetic process",12.578, 5.642, 0.678, 3.156,-1.2396,0.708,0.646),
c("GO:0010468","regulation of gene expression",12.314, 2.901,-4.172, 3.147,-1.5072,0.669,0.661),
c("GO:0031326","regulation of cellular biosynthetic process",11.620, 2.972,-3.134, 3.122,-1.1226,0.606,0.687),
c("GO:0032774","RNA biosynthetic process",10.205, 5.381,-0.207, 3.065,-1.0809,0.672,0.694));

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
p1 <- p1 + geom_label( data = ex, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)),fill="white", size = 4 );
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

dev.print(pdf,"/Users/lukas/sciebo/inquilineGenomics18/GOenrichment/revigo/RELAX.MF.revigo.pdf",width=10,height=8);
