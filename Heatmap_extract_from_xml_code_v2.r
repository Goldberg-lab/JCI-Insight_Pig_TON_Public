#################Comparison of spatial analysis methods to construct topographic cell densities maps. Xml version.
################# Garza-Gisholt, E., Hemmi, J.M., Hart, N.S. and Collin, S.P.


###This version extracts the information from an xml file exported from StereoInvestigator. If you have your information in illustrator or other graphical program, then run the script to extract information from svg file
###We recommend using R Studio especially if you are novice in the use of R
### To download R Studio for any platform, go to http://www.rstudio.com.
### You have to press Ctrl + Enter to run a command. 

###The symbol "#" before any line in R corresponds to notes. Some commands that do not need to run every time have the # symbol; if you need to run them just delete the # and run the command 

### If you have an error message don't panic, instead read what the error says and try to figure out what is wrong with the file or with the commands

### You can always go back and run the commands until you find where the problem is. Also you can reset the values with the button "Clear All" in the Workspace.

### To look for help with specific commands you can run the help from R with ? and the function you want to know for example ?library

### The first part of the script extracts the information from the xml file
### It is good idea to open the xml file to familiarize yourself with the structure of the file and the data that have to be extracted.
### The extraction of the information from the xml file was done with the collaboration of Duncan Temple Lang, Dec 14, 2011. He is the creator of the XML package for R

### It is necessary that you save a file with the name retina.xml exporting the tracing points from Stereo Investigator. This is done by File->Export tracing-> then name the file and select the xml extension.

###We suggest you to use the same contour to count rods, cones and subsampling. We also recommend you to have a directory for each retina and copy this script into it. When you run the script all the images will be saved in the directory where the script is located.

###First set your directory to where you opened the file. A copy of the script should be copied to the directory of your retina, open that copy and then in the menu "Session" select "Set Working Directory" the option "To source file location"

### The first part recalls the x and y coordinates for the contour, optic nerve and each of the markers. It extracts the "Site" where the marker was placed.
###The packages needed to run the script can be installed for the first time using the next command deleting the "#" sign and pressing ctrl+Enter .

#install.packages(c("XML", "plyr", "stringr", "spatstat", "ggplot2", "fields", "Akima", "sp", "RColorBrewer", "raster", "maptools"))


library(XML)  
doc = xmlParse("File S1 retina.xml")
nsURI = c(n = "http://www.mbfbioscience.com/2007/neurolucida")
#contours = getNodeSet(doc, "//n:contour", nsURI)
contour = getNodeSet(doc, "//n:contour[1]/n:point", nsURI)
opticnerve = getNodeSet(doc, "//n:contour[2]/n:point", nsURI) 
markers = getNodeSet(doc, "//n:marker/n:point", nsURI)
getXY =
	function(node)
	{
		as(xmlAttrs(node)[c("x", "y")], "numeric")
	}

contour.xy=as.data.frame(t(sapply(contour, getXY)))
opticnerve.xy=as.data.frame(t(sapply(opticnerve, getXY)))
markers.xy = as.data.frame(t(sapply(markers, getXY)))

names(markers.xy) = c("x", "y")
names(contour.xy)=c("x", "y")
names(opticnerve.xy)=c("x", "y")

###If the retina needs to be rotated, it is better to do it mathematically from the coordinates with the following commands. In this example, the retina is rotated 180 degrees but some cases you will not need rotation. To select a different angle just change "angled=180" to the degrees that you need.  
###In that case it is possible to select all the rotate data commands and comment on them using Ctr+Shift+C

###Rotate the data

#names(markers.xy) = c("xn", "yn")
#names(contour.xy)=c("xn", "yn")
#names(opticnerve.xy)=c("xn", "yn")

#angled=180
#angle=angled/180*pi
#rotcentre=c((mean(opticnerve.xy$xn)),(mean(opticnerve.xy$yn)))

#contour.xy[, "xn"] <- contour.xy$xn+rotcentre[1]
#contour.xy[, "yn"]<- contour.xy$yn+rotcentre[2]
#opticnerve.xy[, "xn"] <- opticnerve.xy$xn+rotcentre[1]
#opticnerve.xy[, "yn"]<- opticnerve.xy$yn+rotcentre[2]
#markers.xy[, "xn"] <- markers.xy$xn+rotcentre[1]
#markers.xy[, "yn"]<- markers.xy$yn+rotcentre[2]

#contour.xy[,"x"] <- contour.xy$xn*cos(angle)-contour.xy$yn*sin(angle)
#contour.xy[,"y"]<- contour.xy$xn*sin(angle)+contour.xy$yn*cos(angle)
#opticnerve.xy[,"x"]<- opticnerve.xy$xn*cos(angle)-opticnerve.xy$yn*sin(angle)
#opticnerve.xy[,"y"]<- opticnerve.xy$xn*sin(angle)+opticnerve.xy$yn*cos(angle)
#markers.xy[,"x"] <- markers.xy$xn*cos(angle)-markers.xy$yn*sin(angle)
#markers.xy[,"y"]<- markers.xy$xn*sin(angle)+markers.xy$yn*cos(angle)

########End of rotate 

markers.xy$Site = factor(xpathSApply(doc, "//n:marker/n:property", xmlValue, namespaces = nsURI))
#detach("package:raster")
library(plyr)
library(ggplot2)
library(stringr)

###The ggplot helps to graphically observe if the contours and the markers are correct. If you do not get any plot, check the xml structure again to see that you are extracting the right information and you do not have any errors.

ggplot(data = contour.xy, aes(x, y)) + 
	geom_path(data= contour.xy, aes(x, y)) +
	geom_path(data= opticnerve.xy, aes(x, y)) +
	geom_point(data=markers.xy, col = "blue") +
	coord_equal()

###The next commands configure the data to be analysed. The markers are counted and the x and y coordinates for each Site is obtained by averaging all the markers.

counts<-data.frame(count(markers.xy, "Site")) 
counts$x<-with(markers.xy, tapply(x, Site, mean))
counts$y<-with(markers.xy, tapply(y, Site, mean))

counts$Site= substring (counts$Site, 2) #KH remove the first part of the site field

counts$freq=as.numeric(counts$freq)
counts$Site=as.numeric(counts$Site)
sapply(counts, class) #KH spits out a tibble of the data classes
counts[, "Site"] <- counts$Site+1 #KH adds +1 to each site so that they are non-zero
counts<-counts[with(counts, order(Site)), ] #KH orders the sites numerically

### Now that the data frame is complete, the next step is to convert the number of cells counted to a standard value, in this case the number of cells per square millimetre. First, you need to change the counting frame value to the counting frame you used. Remember that this is a variable number and you have to change it for each retina analysed. In some cases, for example, when you count photoreceptors it is necessary to express the result in thousands of cells per square millimetre; in this case you need to divide the number of cells by 1000 and express the values in the right units.
counting.frame<- 400*400
counts = transform(counts, cells = (counts$freq * (1000000/counting.frame)))
counts$cells<- round(counts$cells, digits = 0)
head(counts)

###If you want to manually remove any outliers that you identified before, it is possible to do this with the command:

#counts<- subset(counts, !(Site %in% c(1,2,3,4)))
 
###This way it is possible to delete more than one row.


library(spatstat)
library(fields)
library(akima)
library(RColorBrewer)
library(sp)

###The next part of the script will set up the graphical parameters to construct the maps. The package spatstat creates a "window" with the function owin that is the area that will be analysed. 
### The contour of the retina should be drawn in anticlockwise direction and the optic nerve in a clockwise direction. Otherwise, if the owin command marks an error then the nod direction should be reversed as "list(x=rev(xp), y=rev(yp))"


xp<- as.vector(contour.xy$x)
yp<- as.vector(contour.xy$y)
xd<-as.vector(opticnerve.xy$x)
yd<-as.vector(opticnerve.xy$y)

retina <- owin(poly=list(list(x=rev(xp), y=rev(yp)), list(x=xd, y=yd)))
par(mar=c(0.6, 0.6, 0.6, 0.6)) #sets graphical parameters

###Other possible sources of error are if the contour self intersects. In this case, it is recommend identifying the nods where it intersects and deleting them.

plot(retina, hatch=TRUE)

retinamask<-as.mask(retina)

### The next commands set up the graphical parameters for the maps like the colors, the mask and the scale bar. 
### The three color gradients that we use are grey gradient, rainbow gradient but cutting the darker blues from the spectrum or the heat gradient (maps published).


#bw<-rev(grey.level(256))
color<-designer.colors( 256, tim.colors(5), x= c(-0.2, 0.2, 0.4, 0.7, 1.2))
heat<-rev(heat.colors(256))

xs<-as.vector(counts$x)
ys<-as.vector(counts$y)
cells<-as.vector(counts$cells)
samcells<-ppp(xs, ys, window=retina, marks=cells) #KH makes a mysterious object
plot(unmark(samcells), main='', pch=".") #KH plots the mysterious object, retina contour without optic nerve
text(samcells, labels=marks(samcells), cex=0.7) #KH writes in cell density at each coordinate

xrange <- range(xp, na.rm=TRUE)
yrange <- range(yp, na.rm=TRUE)
zrange <- c(30, 1.04*max(cells))

###The xbox and ybox ranges are used for the mask of the maps. Sometimes it does not cover enough area and in that case you can increase the value that extends the range.

xbox<-xrange + c((if(xrange[1]<0) (0.02*xrange[1]) else (-0.02*xrange[1])),
								 (if(xrange[2]>0) (0.02*xrange[2]) else (-0.02*xrange[2]))) 

ybox<-yrange + c((if(yrange[1]<0) (0.02*yrange[1]) else (-0.02*yrange[1])),
								 (if(yrange[2]>0) (0.02*yrange[2]) else (-0.02*yrange[2])))


###The scale bar can be modified to millimetres changing unit="mm" and then reducing the scale to 0.001. The size of the scale bar is specified with the size at the end of the function.

scalebar<-function (size, unit="cm", scale=.0001, t.cex= 0.8)
{
	x=0.98*xrange[2]-size
	y=yrange[1]+(0.06*(yrange[2]-yrange[1]))
	xvals=size * c(0, 0.5, 1) + x
	yvals=c(0, 0.01*(yrange[2]-yrange[1]), 0.03*(yrange[2]-yrange[1]), 0.04*(yrange[2]-yrange[1]))+ y
	for (i in 1:2) rect(xvals[i], yvals[3], xvals[i + 2], yvals[4], 
											col = "black")
	labels <- c(paste(size*scale, unit))
	text(xvals[c(2)], yvals[1], labels = labels, adj = 0.5, 
			 cex = t.cex)
}

size<-10000

mask<-function()
{
	polypath(c(xp, NA, c(xbox, rev(xbox))),
					 c(yp, NA, rep(ybox, each=2)),
					 col="white", rule="evenodd", lty=0)
	polypath(xd, yd, col="black")
	plot(retina, main='', add=TRUE, lwd=2, scalebar(size))
}

###The first map is the Gaussian Kernel Smoother from the spatstat package. The sigma value can be adjusted to the distance between points. If it is omitted, the smoothing kernel bandwidth is chosen by least squares cross-validation.
###It is possible to change graphic parameters in the plot and contour functions. The col= can be changed to bw for black and white. nlevels is the number of contours but levels=c() and allows you to specify what contours will be plotted. For more options look ?contour and ?plot

dens<-Smooth(samcells, sigma = 2000)
plot(dens, main='', col=heat, win=retina, zlim=zrange)
#contour(dens, add=TRUE, nlevels=5, asp=1, drawlabels=TRUE, levels=c(50, 100, 150, 200, 250), labcex=0.7, lwd=2)
mask()

### The second map is the akima linear interpolation. The sequence of values can be modified in the "by=". In the example, a value is calculated every 200 microns. In this case, the command to plot the map is surface.

akimalin<-interp(xs, ys, cells, 
								 xo=seq(xrange[1], xrange[2], by=200),
								 yo=seq(yrange[1], yrange[2], by=200), 
								 linear=TRUE)
surface(akimalin, asp=1, col=heat, axes=FALSE, levels=c(50, 100, 150, 200, 250), ylim=yrange, xlim=xrange, zlim=zrange) 
mask()

### The third and fourth maps both work with the function Tps from the package fields. It gives a krig object that allows predicting values with the function. To calculate values every 200 microns, a grid is created with the following command.

grid<- make.surface.grid( list( seq((xrange[1]), (xrange[2]), by=200), seq((yrange[1]), (yrange[2]), by=200)))
coord<-cbind(xs, ys)

### The third map is a spline cubic interpolation. It uses the Tps function with a lambda value of 0.

kinterp<-Tps(coord, cells, lambda=0)
look<- predict(kinterp, grid) 
out.p<-as.surface( grid, look)
surface(out.p, asp=1, col=heat, axes=FALSE, levels=c(50, 100, 150, 200, 250), zlim=zrange)
mask()

### The fourth map is the Tps with the generalized cross validation (GCV) smoothing value. It is possible to change the smoothing with the degrees of freedom (df=) in the Tps function. 

k<-Tps(coord, cells)
look2<- predict(k, grid) 
out.p2<-as.surface( grid, look2)
surface(out.p2, asp=1, col=heat, lwd.poly=0.8, levels=c(50, 100, 150, 200, 250), axes=FALSE, zlim=zrange)
mask()

################RESIDUALS###################
###The next list of commands will analyse the residuals of the two smoothing models Gks and Tps comparing the observed values to the modelled values. The maps show the position of the variation and the plots show the variation in the x and y axes. 


library(raster)
library(maptools)

denssp<-as.SpatialGridDataFrame.im(dens)
densspras <- raster(denssp)

coords<-as.data.frame(coord)
coords$observed<-counts$cells
coords$tpsinterp<-predict(kinterp)
coords$tps<-predict(k)
coords$gks<-extract(densspras, as.data.frame(coord))

coords$tps.res<-(abs(coords$observed-coords$tps)*100/(coords$observed))
coords$gks.res<-(abs(coords$observed-coords$gks)*100/(coords$observed))

par(mar=c(2.5, 2.5, 2.5, 2.5))

res.tps.diff<-as.vector(coords$tps.res)
kinterp.tps<-Tps(coord, res.tps.diff, lambda=0)                    
surface(kinterp.tps,  asp=1, col=color, axes=TRUE, labcex=0.8, ylim=yrange, levels=c(10,50,100), zlim=c(0,200))
mask()

res.gks.diff<-as.vector(coords$gks.res)
kinterp.gks<-Tps(coord, res.gks.diff, lambda=0)                    
surface(kinterp.gks,  asp=1, col=color, axes=TRUE, labcex=0.8, ylim=yrange, levels=c(10,50,100))
mask()


#########TRANSECTS###########
### It is very useful to draw transects in the maps and extract the data from those transects.  

y= rep(akimalin$y, each = length(akimalin$x))
x= rep(akimalin$x, length(akimalin$y))
z = as.numeric(akimalin$z)
akimalinsp = data.frame(x, y, z)

coords.long<-as.data.frame(grid)
coords.long$akimalin<-akimalinsp$z
coords.long$tpsinterp<-predict(kinterp, grid)
coords.long$tps<-predict(k, grid)
coords.long$gks<-extract(densspras, grid)
names(coords.long)<-c("xs", "ys", "observed", "tpslinear", "tps", "gks")
coords.long$xs<- round(coords.long$xs, digits = 0)
coords.long$ys<- round(coords.long$ys, digits = 0)
coords.long<-coords.long[!duplicated(coords.long$gks),]

par(mar=c(1.5, 1.5, 1.5, 1.5))
plot(dens, main='', col=bw, win=retina, bty="n", axes=TRUE)
plot(retina, main='', add=TRUE, lwd=2)

### The next commands help to decide where to place transects. The first command will provide the coordinates with the highest x and y values. The tables show the coordinates and the numbers of sites per coordinate.

coords.long[which.max(coords.long$tpslinear), ]

table(coords.long$xs)
table(coords.long$ys)

transectx1.val<--32078
transecty1.val<--19550

transectx1<-subset(coords.long, ys==transectx1.val)
transecty1<-subset(coords.long, xs==transecty1.val)

transectx1.observations<-subset(coords, ys>=(transectx1.val-500)& ys<=(transectx1.val+500)) 
transecty1.observations<-subset(coords, xs>=(transecty1.val-500)& xs<=(transecty1.val+500)) 


plot(dens, main='', col=bw, win=retina, bty="n", axes=TRUE)
plot(retina, main='', add=TRUE, lwd=2)
lines(transectx1$ys~transectx1$xs, lwd=2, col="red")
lines(transecty1$ys~transecty1$xs, lwd=2, col="blue")

ggplot(transectx1, aes(xs)) +
	geom_line(aes(y = observed))+
	geom_point(data=transectx1.observations, aes(y=observed), col = "blue") +
	scale_x_continuous(limits=c(min(coords$xs), max(coords$xs)), "")+
	scale_y_continuous(limits=c(0, 300), "")+
	theme_bw()+
	theme(legend.position="top")

ggplot(transectx1, aes(xs)) +
	geom_line(aes(y = tpslinear))+
	geom_point(data=transectx1.observations, aes(y=observed), col = "blue") +
	scale_x_continuous(limits=c(min(coords$xs), max(coords$xs)), "")+
	scale_y_continuous(limits=c(0, 300), "")+
	theme_bw()+
	theme(legend.position="top")

ggplot(transectx1, aes(xs)) +
	geom_line(aes(y = tps))+
	geom_point(data=transectx1.observations, aes(y=observed), col = "blue") +
	scale_x_continuous(limits=c(min(coords$xs), max(coords$xs)), "")+
	scale_y_continuous(limits=c(0, 300), "")+
	theme_bw()+
	theme(legend.position="top")

ggplot(transectx1, aes(xs)) +
	geom_line(aes(y = gks))+
	geom_point(data=transectx1.observations, aes(y=observed), col = "blue") +
	scale_x_continuous(limits=c(min(coords$xs), max(coords$xs)), "")+
	scale_y_continuous(limits=c(0, 300), "")+
	theme_bw()+
	theme(legend.position="top")

#############DISTRIBUTION FUNCTIONS############

###It is possible to create the denstity distribution curves to compare the functions. Also it is possible to create the empirical cumulative distribution function (ecdf). 


ggplot(coords.long)+
	geom_density(aes(x=gks, colour="Gks"))+
	geom_density(aes(x=observed, colour="Akimainterp"))+
	geom_density(aes(x=tpslinear, colour="Tpsinterp"))+
	geom_density(aes(x=tps, colour = "Tps"))+
 	scale_colour_manual("",values=c("Akimainterp"="black", "Tpsinterp"="orange","Tps"="blue", "Gks"= "red"), breaks = c("Akimainterp", "Tpsinterp", "Tps", "Gks"))+
	scale_x_continuous(expand = c(0, 0), limits=c(0, max(coords$tps)) )+
	xlab(expression("cells per "*mm^2))+	
	scale_y_continuous(expand = c(0, 0), "")+
	theme_bw()+
	theme(legend.position="top")+
	theme(axis.text.y = element_text(angle = 90, hjust = 1))

coords.long<-na.omit(coords.long)

tpsinterpolation.ecdf<-ecdf(coords.long$tpslinear)
akima.ecdf<-ecdf(coords.long$observed)
tps.ecdf<-ecdf(coords.long$tps)
gks.ecdf<-ecdf(coords.long$gks)

tps.95<-quantile(tps.ecdf, c(.95)) 
gks.95<-quantile(gks.ecdf, c(.95))
tpsinterp.95<-quantile(tpsinterpolation.ecdf, c(.95)) 
akima.95<-quantile(akima.ecdf, c(.95))


table.ecdf.obs<-as.data.frame(coords.long$observed)
names(table.ecdf.obs)<-"cells"
table.ecdf.obs$model<-"akima"

table.ecdf.tpsint<-as.data.frame(coords.long$tpslinear)
names(table.ecdf.tpsint)<-"cells"
table.ecdf.tpsint$model<-"tpsint"

table.ecdf.tps<-as.data.frame(coords.long$tps)
names(table.ecdf.tps)<-"cells"
table.ecdf.tps$model<-"tps"

table.ecdf.gks<-as.data.frame(coords.long$gks)
names(table.ecdf.gks)<-"cells"
table.ecdf.gks$model<-"gks"

table.ecdf<-rbind(table.ecdf.obs, table.ecdf.tpsint, table.ecdf.tps, table.ecdf.gks)

ecdf <- ddply(table.ecdf, .(model), summarize,
							cells = unique(cells),
							ecdf = ecdf(cells)(unique(cells)))

ggplot(ecdf, aes(cells, ecdf, color = model)) + 
	geom_hline(yintercept=0.95, linetype = "longdash")+
	geom_segment(aes(x=tpsinterp.95, y=0, xend=tpsinterp.95, yend=0.95), 
							 colour="orange", linetype = "longdash")+
	geom_segment(aes(x=tps.95, y=0, xend=tps.95, yend=0.95), 
							 colour="blue", linetype = "longdash")+
	geom_segment(aes(x=gks.95, y=0, xend=gks.95, yend=0.95), 
							 colour="red", linetype = "longdash")+
	geom_segment(aes(x=akima.95, y=0, xend=akima.95, yend=0.95), 
							 colour="black", linetype = "longdash")+
	scale_colour_manual("", values=c("akima"="black", "tpsint"="orange","tps"="blue", "gks"="red"), breaks = c("akima", "tpsint", "tps", "gks"))+
	scale_x_continuous(expand = c(0, 0), limits=c(0, max(coords$tps)) )+
	xlab(expression("cells per "*mm^2))+
	scale_y_continuous(expand = c(0, 0), "")+
  theme_bw()+
	theme(legend.position="top")+
	geom_step()+
	theme(axis.text.y = element_text(angle = 90, hjust = 1))

###Finally, maps can be exported in pdf, jpeg, tiff and other formats. For pdf:

#pdf('name.pdf')
#All the lines of the plots that want to be added (can be more than one plot)
#dev.off()

###For publication, using Arial font is good to follow instructions from the blog: http://r.789695.n4.nabble.com/How-to-enable-Arial-font-for-postcript-pdf-figure-on-Windows-td3017809.html

