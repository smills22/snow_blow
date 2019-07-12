##    snow_blow_v1.py (corresponds to working version 0.48)
##
##    Snow_blow model: A simple model for calculating snow erosion and
##    distribution from a DEM and wind direction.  The program is designed for
##    use in palaeo-environmental settings where meteorological observations
##    are not available.
##
##    If you use this code, please cite Mills et al. (in review) "Testing and
##    application of a model for snow redistribution 1 (Snow_Blow) in the
##    Ellsworth Mountains, Antarctica. Journal of Glaciology.
##
##    The default parameters are set to the values used in Mills et al. paper.
##
##    Contact Stephanie Mills: smills@uow.edu.au for further information.
##
##    Copyright (C) 2019 Stephanie Mills
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <https://www.gnu.org/licenses/>.


# ---START ROUTING FUNCTION---------------------------------------------------------------------------------------
def routing(ns,ew):
    global rows,columns,sd,se,wd,ns_orig,ew_orig,se_index,se_check

    dirpp = ([3,4,5,2,-1,6,1,0,7])
    pp = numpy.zeros(9)
    rcount = 0

    #check that it is possible to do the window on this cell...
    #note this has to be rows-1 because of 0 array dimension
    if (rows-1 > ns > 0) and (cols-1 > ew > 0):

        #Check we haven't been here already... 
        #If we have, then can just return...
        if se_check[ns,ew] == 0:
            se_check[ns,ew] = 1

            for i in range(ns-1, ns+2):
                for j in range(ew-1, ew+2):

                    pp = numpy.zeros(9)
                    if rcount != 4:

                        if wd[i,j] < 0:
                            wd[i,j] = windDir

                        #convert wind direction to flow direction
                        if wd[i,j] < 180:
                            fd = wd[i,j] + 180
                        else:
                            fd = wd[i,j] - 180
                        
                        #calculate sector
                        sector = int(fd/45.0)
                        if sector == 8:
                            sector = 0

                        alpha2 = ((sector + 1) * 45.0) - fd

                        if sector != 7:
                            pp[sector] = alpha2/45.0
                            pp[sector + 1] = 1- pp[sector]
                        else:
                            pp[sector] = alpha2/45.0
                            pp[0] = 1 - pp[sector]

                        if pp[dirpp[rcount]] > 0:
                            routing(i,j)
                            # set up a new list for [ns,ew] using the upwind cell but shifted by 1.
                            # (don't want to overwrite if this is the second call)
                            for k in range(1, vh_dist):
                                se_index[ns,ew,k] = se_index[ns,ew,k] + pp[dirpp[rcount]] * se_index[i,j,k-1]
                            # add upstream neighbour snow eroded to cell 0.
                            se_index[ns,ew,0]= se_index[ns,ew,0] + pp[dirpp[rcount]] * se[i,j]
                      
                    rcount = rcount + 1

# ---END ROUTING FUNCTION-----------------------------------------------------------------------------------------

# ---START MAIN PROGRAM-------------------------------------------------------------------------------------------

# Import system modules
print "Importing System Modules..."
import arcpy, os, math, sys, string, numpy
from arcpy.sa import *
from arcpy import env

# Set up workspace
print "Setting up workspace..."
arcpy.env.overwriteOutput = True
arcpy.CheckOutExtension("Spatial")
Workspacefolder = "C:\\Example\\folder\\"
env.workspace = Workspacefolder

# Input values
print "Defining parameters..."
#Output text file option:
paramsfile = 1 # 1 is write to output text file, 0 is no output
#DEM file:
DEM = "C:\\Example\\dem.tif"  
#Wind options:
windDir = 122.5 # wind direction - don't use a value of 360, use 0 instead.
windSpeed = 15.0   # wind speed m/s 
tWindSpeed = 5  # threshold wind speed m/s 
#Routing options:
snow_routing = 1 # 1 is routing on, 0 is off
winddir_option = 1 # 1 is blowing down amended wind direction, 0 is blowing down background wind direction
iterations = 8 # number of iterations for calculation of amount of snow accumated in each cell. 
#Deposition options: 
dep_dist = 690.0 # distance over which deposition occurs (m). There is NO bug testing for an odd size kernel
mean_dist = 150.0 # Exponential: uses mean distance, not max distance
#Erosion options:
erosion_scaling = 1 # 1 is erosion scaling on, 0 is off
erosion_scale = 1 # maximum erosion value
#Curvature options:
curvature = 1 # 1 is curvature on, 0 is curvature off
curvatureval = 0.5
#Slope options:
minslopeval = 5.0 # slope value at which slope index is 0(degrees)
mslope = 1 # 1 is maxslopeval on, 0 is maxslopeval off (uses highest value in DEM as maxslopeval)
maxslopeval = 20.0 # slope value at which slope index saturates at 1(degrees)
#Other options:
boundary_snow = 1 # 1 is boundary snow on, 0 is off

# Write parameters to an output file
if paramsfile == 1:
    #open params.txt
    f = open(Workspacefolder + "params.txt", 'w')
    f.write("DEM: " + DEM + "\n")
    f.write("windDir: " + str(windDir) + "\n")
    f.write("windSpeed: " + str(windSpeed) + "\n")
    f.write("tWindSpeed: " + str(tWindSpeed) + "\n")
    f.write("snow_routing: " + str(snow_routing) + "\n")
    f.write("winddir_option: " + str(winddir_option) + "\n")
    f.write("iterations: " + str(iterations) + "\n")
    f.write("dep_dist: " + str(dep_dist) + "\n")
    f.write("mean_dist: " + str(mean_dist) + "\n")
    f.write("erosion_scaling: " + str(erosion_scaling) + "\n")
    f.write("erosion_scale: " + str(erosion_scale) + "\n")
    f.write("curvature: " + str(curvature)  + "\n")
    f.write("curvatureval: " + str(curvatureval) + "\n")
    f.write("minslopeval: "  + str(minslopeval)  + "\n")
    f.write("mslope: " + str(mslope) + "\n")
    f.write("maxslopeval: "  + str(maxslopeval)  + "\n")
    f.write("boundary_snow: " + str(boundary_snow) + "\n")
    #close params.txt
    f.close()

# Calculable parameters 
rows = int(str(arcpy.GetRasterProperties_management(DEM, "ROWCOUNT")))
cols = int(str(arcpy.GetRasterProperties_management(DEM, "COLUMNCOUNT")))
cellsize =  str(arcpy.GetRasterProperties_management(DEM, "CELLSIZEX"))
ymax = float(str(arcpy.GetRasterProperties_management(DEM, "TOP")))
xmin = float(str(arcpy.GetRasterProperties_management(DEM, "LEFT")))
xmax = float(str(arcpy.GetRasterProperties_management(DEM, "RIGHT")))
ymin = float(str(arcpy.GetRasterProperties_management(DEM, "BOTTOM")))
vh_dist = int(round(dep_dist/float(cellsize),0)) #Round up when depositing snow  across cells
d_dist = int(round(dep_dist/math.sqrt(float(cellsize)**2+float(cellsize)**2),0))
L = 1.0/(mean_dist) #Exponential - Lamba
A = 1.0 # Units of snow to shift
max_dist = math.log(1-0.99)/-L

# Print metadata
print "  Columns:", cols
print "  Rows:", rows
print "  Cellsize:", cellsize

# Set Python recursive limit to number of cells in grid (default is 1000)
sys.setrecursionlimit(rows*cols)
print "Python recursion limit set to:",rows*cols

# Create raster layers with variables set used in the model (constant across the grid)
windDirRas = CreateConstantRaster(windDir, "FLOAT", cellsize, Extent(xmin, ymin, xmax, ymax))
windDirRas.save("windDir")
windSpeedRas = CreateConstantRaster(windSpeed, "FLOAT", cellsize, Extent(xmin, ymin, xmax, ymax))
windSpeedRas.save("windSpeed")
tWindSpeedRas = CreateConstantRaster(tWindSpeed, "FLOAT", cellsize, Extent(xmin, ymin, xmax, ymax))
tWindSpeedRas.save("tWindSpeed")

# Calculate deposition weights-----------------------------------------------------------------------------------
print "Calculating weights..."
#Initialise weights
#vh_weights. Use for loop to apply weights function. Currently set up for hardcoded 3x3 kernel
#Applies exponential using B = A [exp(-L*x1) - exp(-L*x2)] but A is multiplied in the iterative loop below
vh_weights=[] #initialise list
B_cum1 = 0 #initialise variable to store cumulative weight
for i in range(0, int(vh_dist)):
        x1 = 0 + (i* int(round(float(cellsize))) ) # Exponential - closest distance to source for a target cell
        x2 = int(round(float(cellsize)))  + (i* int(round(float(cellsize))) ) # Exponential - farthest distance to source for a target cell
        B = math.exp(-L*x2) - math.exp(-L*x1)
        vh_weights.append(-B)
        B_cum1 = B+B_cum1

#scale vh_weights calculated above to sum to 1.
vh_sum = sum(vh_weights)
for i in range(0, int(vh_dist)):
    vh_weights[i] = vh_weights[i] / float(vh_sum)

#d_weights
d_weights=[] #initialise list
B_cum2 = 0 #initialise variable to store cumulative weight
for i in range(0, int(d_dist)):
        x1 = 0 + (i* int(round(math.sqrt(float(cellsize)**2+float(cellsize)**2))) ) # Exponential - closest distance to source for a target cell
        x2 = int(round(math.sqrt(float(cellsize)**2+float(cellsize)**2))  + (i* int(round(math.sqrt(float(cellsize)**2+float(cellsize)**2))))) # Exponential - farthest distance to source for a target cell
        B = math.exp(-L*x2) - math.exp(-L*x1)
        d_weights.append(-B)
        B_cum2 = B+B_cum2

# Calculate slope in % and degrees, and aspect ------------------------------------------------------------------
print "Calculating slope & aspect..."
slope_deg = Slope(DEM, "DEGREE")
slope_deg.save("slope_deg")
slope_perc = Slope(DEM, "PERCENT_RISE")
slope_perc.save("slope_perc")
aspect = Aspect(DEM)
aspect.save("aspect")
maxSlope = float(str(arcpy.GetRasterProperties_management("slope_deg", "MAXIMUM")))

# Calculate amended wind direction for all cells depending on slope and aspect-----------------------------------
print "Calculating wind direction..."
#AS PER PURVES p315 and Ryan 1977
#first calculate the deflection
WindDefl = -0.225 * slope_perc * Sin(2 * ((aspect * 0.01745) - (windDirRas * 0.01745)))
#then we need to know whether to add the deflection (between -90 and 90) or subtract the deflection (outside -90 and 90)
#first calculate the absolute difference:
WindAspectDiffAbs = aspect - windDirRas    
#then scale this to between -180 and 180:
WindAspectDiffSca = Con(WindAspectDiffAbs <= -180,WindAspectDiffAbs + 360,      #if diff <= -180 then add 360
                        Con(WindAspectDiffAbs > 180,WindAspectDiffAbs - 360,    #if diff > 180 then substract 360
                            WindAspectDiffAbs))                                 #otherwise no change

#add the deflection if between -90 and 90, otherwise subtract the deflection
amWindDir=Con(WindAspectDiffSca == 180, windDirRas,
              Con(WindAspectDiffSca == 0, windDirRas,
                  Con(WindAspectDiffSca > 90,(windDirRas - WindDefl),           #between 90 and 180 then subtract difference
                      Con(WindAspectDiffSca >- 90,(windDirRas + WindDefl),      #between -90 and 90 then add difference 
                          Con(WindAspectDiffSca >- 180,(windDirRas - WindDefl), #between -180 and -90 then subtract difference
                              windDirRas)))))                                   #otherwise no change
amWindDir.save("amWindDir")

#calculate amended wind direction taking into account circular nature
amWindDirFin = Con(amWindDir < 0,amWindDir+360,                             #if negative value then add 360
                   Con(amWindDir >= 360, amWindDir-360,                     #if greater than or equal to 360 (don't want a winddir of 360)
                       amWindDir))                                          #otherwise no change
amWindDirFin.save("amWindDirFin")

# Calculate components of the shelter index ----------------------------------------------------------------------
print "Calculating shelter index..."

# Calculate Aspect index - continous index values
#calculate lee aspect - opposite direction to aspect - taking 0-360 degrees limits into account, i.e. add 180 deg to wind direction value
#if it is between 0 and 180 deg, and subtract 180 degrees if it is between 180 and 360 deg.
leeAspect = Con(amWindDirFin < 180,amWindDirFin + 180,   #if wind between 0 and 180 then add 180
                    amWindDirFin - 180)                     #otherwise subtract 180
leeAspect.save("leeAspect")
#calculate difference between aspect and lee aspect of the prevailing wind
aspectdiff = Abs(leeAspect - aspect)
aspectdiff.save("aspectdiff")
#convert this to a 0-1 scale
aspectIndex = Con(aspectdiff > 315.0, 1.0-((360.0- aspectdiff)/45.0), Con(aspectdiff < 45, 1.0-(aspectdiff/45.0),0))
aspectIndex.save("aspectIndex")

# Calculate slope index
if mslope == 0:
    maxslopeval = maxSlope
slopeIndex = Con(slope_deg >= minslopeval,Con(slope_deg <= maxslopeval, (slope_deg-minslopeval)/(maxslopeval-minslopeval),1),0)
slopeIndex.save("slopeIndex")

#calculate curvature index
if curvature == 1:
    Curvature(DEM, "#", "#","plancurve")
    #need to use the "plancurve" output from above:
    curveIndex = Con(Raster("plancurve") >= 0,Con(Raster("plancurve") <= curvatureval, 1.0-(Raster("plancurve")/curvatureval),0),1)
    curveIndex.save("curveIndex")

# Calculate shelter index
if curvature == 1:
    shelterIndex = aspectIndex * slopeIndex * curveIndex
else:
    shelterIndex = aspectIndex * slopeIndex
shelterIndex.save("shelterIndex")

# Amend wind speed accordingly - Wind speed reduced based upon sheter index --------------------------------------
amWindSpeed = windSpeed - (windSpeed * shelterIndex)
amWindSpeed.save("amWindSpeed")

# Calculate the amount of snow to be eroded from each cell--------------------------------------------------------
# ...if the amended wind speed is greater than the threshold wind speed. (it will stay constant with each iteration below)
print "Calculating erosion..."
snowEroded = Con(amWindSpeed > tWindSpeed,pow(amWindSpeed, 3) - pow(tWindSpeed, 3),0)
snowEroded.save("snoweroded")

# Switch to NumPy to do the snow routing--------------------------------------------------------------------------
if snow_routing == 1:

    #convert rasters with amended wind direction and amount of snow eroded from each cell to arrays for
    #processing across the grid further down
    if winddir_option == 1:
        wd = arcpy.RasterToNumPyArray("amWindDirFin")
    elif winddir_option == 0:
        wd = arcpy.RasterToNumPyArray("windDir")
    se = arcpy.RasterToNumPyArray("snowEroded")

    #display max erosion
    max_erosion = numpy.max(se)
    print "  Max erosion: ", max_erosion

    #Snow movement
    print "Snow routing..."
    print "  No. Iterations:", iterations
    if boundary_snow == 1:
        print "  Boundary snow: On"
    else:
        print "  Boundary snow: Off"

    # Create empty arrays with the same number of rows and columns as other rasters - to be filled in with
    # values for:
    # Snow depth (sd):
    sd = numpy.zeros((rows, cols))
    # Snow accumulation grid (sa):
    sa = numpy.ones((rows, cols))
    # set up array to store the upwind contributing totals
    se_index = numpy.zeros((rows, cols,vh_dist))
    # set up array to record which cells we have already visited
    se_check = numpy.zeros((rows, cols))

    #-0.0 values were messing up the routing - so here I set them all to 0.
    for ns_orig in range(0, rows):
        for ew_orig in range(0, cols):
            if (se[ns_orig,ew_orig] < 0.0):
                se[ns_orig,ew_orig] = 0.0

    # normalising snow eroded grid to be between 0 and 1
    # (recording initial erosion in se_orig)
    if erosion_scaling == 1:
        se_orig = (se.copy()/max_erosion)*erosion_scale
        # set se to these initial values
        se = se_orig.copy()
        # scale sa to be equal to max erosion
        sa = sa * erosion_scale
    else:
        se_orig = se.copy()
        # scale sa to be equal to max erosion
        sa = sa * max_erosion
        # need to reset erosion_scale for the purpose of setting the boundary conditions
        erosion_scale = max_erosion

    print "  Max erosion scaled to: ", erosion_scale

    #set up the boundary values here if boundary_snow is on 
    if boundary_snow == 1: 

        for ns_orig in range(0, rows): 
            for ew_orig in range(0, cols):    
                if (ns_orig == rows-1) or (ns_orig == 0) or (ew_orig == cols-1) or(ew_orig == 0):  
                    se_index[ns_orig,ew_orig] = erosion_scale 
                    se[ns_orig,ew_orig] = erosion_scale 

    # variable to count iterations:
    count = 1
    # ITERATIONS LOOP
    while count <= iterations:

        print "\nIteration:", count

        # **RESET ARRAYS IF NOT FIRST ITERATIONS
        if count > 1:
            print "  Resetting arrays..."
             
            for ns_orig in range(0, rows):
                for ew_orig in range(0, cols):
                    
                    # Calculate snow that can be redistributed in iterations > 1

                    # case 1: if all the snow erodes from the cell
                    if sa[ns_orig,ew_orig] < se_orig[ns_orig,ew_orig]:
                        # set se to be the snow that was deposited in the previous iteration
                        se[ns_orig,ew_orig] = sa[ns_orig,ew_orig]
                    # case 2: need to reset se as it might have been changed in a previous iteration to be less
                    else:
                        se[ns_orig,ew_orig] = se_orig[ns_orig,ew_orig]
                    
                    # reset the index grid (the one that stores the upwind snow input)
                    se_index[ns_orig,ew_orig] = 0
                    # reset the check grid
                    se_check[ns_orig,ew_orig] = 0
                    # reset grid to store snow 
                    sd[ns_orig,ew_orig] = 0

                    # Now need to reset the boundary cells if boundary conditions are on.
                    if boundary_snow == 1:
                        if (ns_orig == rows-1) or (ns_orig == 0) or (ew_orig == cols-1) or (ew_orig == 0):
                            se_index[ns_orig,ew_orig] = erosion_scale
                            se[ns_orig,ew_orig] = erosion_scale

        # SNOW ROUTING
        print "  Routing snow..."
        # rows loop
        for ns_orig in range(0, rows):
            # columns loop
            for ew_orig in range(0, cols):

                # call routing algorithm for cell ns_orig,ew_orig
                routing(ns_orig,ew_orig)

                # calculate sd for cell ns_orig,ew_orig
                for k in range(0, vh_dist):
                    sd[ns_orig,ew_orig] = sd[ns_orig,ew_orig] + se_index[ns_orig,ew_orig,k]*vh_weights[k]
                    
                # Calculate snow accumulation (sa)
                sa[ns_orig,ew_orig] = sa[ns_orig,ew_orig] + sd[ns_orig,ew_orig] - se[ns_orig,ew_orig]

                # check for negatives
                if sa[ns_orig,ew_orig] < 0.0:
                    sa[ns_orig,ew_orig] = 0.0

        print "  End of iteration", count

        #Convert sa array to raster with the same extent and cell size as DEM and save it
        snowAcc = arcpy.NumPyArrayToRaster(sa, arcpy.Point(xmin, ymin), DEM, DEM)
        #snowAccExp.save(Workspacefolder + "/snow_acc.tif")
        snowAcc.save("snow_acc"+str(count))
        
        count = count + 1

    # END ITERATIONS LOOP
    # End snow routing------------------------------------------------------------------------------------------------

    # Save final files to rasters-------------------------------------------------------------------------------------
    # Convert sd array to raster with the same extent and cell size as DEM and save it
    snowDep = arcpy.NumPyArrayToRaster(sd, arcpy.Point(xmin, ymin), DEM, DEM)
    snowDep.save("snowDep")
    #Convert se_orig array to raster with the same extent and cell size as DEM and save it
    snowErodedInit = arcpy.NumPyArrayToRaster(se_orig, arcpy.Point(xmin, ymin), DEM, DEM)
    snowErodedInit.save("snowerodedIni")
    #Convert last se array instance to raster with the same extent and cell size as DEM and save it
    snowErodedFin = arcpy.NumPyArrayToRaster(se, arcpy.Point(xmin, ymin), DEM, DEM)
    snowErodedFin.save("snowerodedfin")
    #Convert sa array to raster with the same extent and cell size as DEM and save it
    snowAcc = arcpy.NumPyArrayToRaster(sa, arcpy.Point(xmin, ymin), DEM, DEM)
    snowAcc.save("snow_acc")
    #Calculate final index field and save it
    IndexFinal = arcpy.NumPyArrayToRaster(((sa-erosion_scale)/erosion_scale), arcpy.Point(xmin, ymin), DEM, DEM)
    IndexFinal.save("index_final")

print "\nEnd"

# ---END MAIN PROGRAM---
                                   



   




    

