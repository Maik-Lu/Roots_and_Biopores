#@ Integer(label="Length filter in voxel", value=60,min=0, max=200, style="slider") t_s
#@ Double(value=0.6, min=0.0, max=1.0, label="Vesselness filter") t_v
#@ Double(value=0.7, min=0.0, max=1.0, label="Mask for unsharp") mask
#@ Integer(value=1, min=0, max=10, label="Radius for unsharp") radius
#@ String(label="Auto threshold?",choices={"Otsu", "Own value"},style="radioButtonHorizontal") Threshold_value
#@ Integer(label="Upper_Threshold if own value checked", value=120) Upper_Threshold
#@ Double(label="Ratio for lower threshold", value=3) Threshold_ratio
#@OUTPUT Dataset output
#@ OpService ops
#@ DatasetService ds
#@ ConvertService convertService
#@ UIService ui
#@ImageJ ij
#@OUTPUT ImagePlus output
# This script segments 3D roots (based on their tubeness like shape) in 3D CT 
#MLucas 01.2022

import ij.IJ
from net.imagej import Dataset
from ij import ImagePlus
from ij.measure import ResultsTable
from ij.plugin import ImageCalculator
import cmath
from ij import IJ, WindowManager
from ij.plugin.frame import ThresholdAdjuster
def Tubeness(data,sigma1, low, high):
	data1 = convertService.convert(data, Dataset)
	slices = data1.getDepth()
	tube = ops.filter().tubeness(data1,sigma1)
	output = ds.create(tube)
	IJ.run("Conversions...", " ");
	tube8bit = ops.convert().uint8(output)
	tube8bit = ds.create(tube8bit)	
	ui.show(tube8bit)
	result2 = IJ.getImage()
	result2.setTitle("result2")
	tube8bit3=IJ.run(result2 ,"Properties...", "channels=1 slices="+str(slices)+" frames=1 pixel_width=1.0000 pixel_height=1.0000 voxel_depth=1.0000");	
	result=IJ.run(tube8bit3, "3D Hysteresis Thresholding", "high="+str(high)+" low="+str(low)+" connectivity");
	result2.close();
	return result
	

def segmentsmallBiopores(scale05):
	results = []
	for i in range(1,5):
		Tubeness(scale05,i,55,45)
		a = IJ.getImage()
		results.append(a)
		
	imp12 = ImageCalculator.run(results[0],results[1], "Max create stack");
	imp13 = ImageCalculator.run(imp12, results[2], "Max create stack");
	imp14 = ImageCalculator.run(imp13, results[3], "Max create stack");
	imp14.show()
	smallBiopores = IJ.getImage() 
	smallBiopores.setTitle("smallBiopores")
	for i in range(0,4):
		results[i].close()
	imp13.close()
	imp12.close()		
	return smallBiopores

def segmentlargeBiopores(scale02):
	width, height, ch , slices, fr = scale02.getDimensions()
	scale02.setSlice(slices);
	IJ.run(scale02,"Add Slice", "");
	IJ.run(scale02, "Flip Z", "");
	scale02.setSlice(slices+1);
	IJ.run(scale02,"Add Slice", "");
	IJ.run(scale02, "Flip Z", "");
	scale02.show()
	mask = IJ.getImage()
	mask.setTitle("mask");
	results = []
	maxStack = []
	for i in range(2,5):
		Tubeness(mask,i,55,30)
		a = IJ.getImage()
		results.append(a)
	for i in range(4,31):
		Tubeness(mask,i, 50, 30)
		a = IJ.getImage()
		results.append(a)
	a = ImageCalculator.run(results[0], results[1], "Max create stack");	
	maxStack.append(a)
	for i in range(1,28):
		a = ImageCalculator.run(results[i+1],maxStack[i-1], "Max create stack");
		maxStack.append(a)
	maxStack[27].show()
	largeBiopores2 = IJ.getImage() 
	largeBiopores2.setTitle("largeBiopores")
	for i in range(0,30):
		results[i].close()	
	for i in range(0,26):
		maxStack[i].close()
	largeBiopores2.setSlice(1);
	IJ.run(largeBiopores2,"Delete Slice", "");
	largeBiopores2.setSlice(slices+1);
	IJ.run(largeBiopores2,"Delete Slice", "");
	largeBiopores2.show();
	largeBioporesRemoved4 =largeBiopores2.resize(total_width/2,total_height/2,  total_slices/2, "none");
	largeBioporesRemoved4.show();
	largeBioporesfinal1= IJ.getImage();
	largeBioporesfinal1.setTitle("largeBiioporesfinal1");
	largeBiopores2.changes= False;
	largeBiopores2.close();
	mask.changes= False;
	mask.close();
	return largeBioporesfinal1

def RemoveBloobs(image,ves, length, dyn):
	image = IJ.getImage()	
	image.setTitle("bio");
	smallBioporesRemoved=IJ.run(image,"Distance Transform Watershed 3D", "distances=[Borgefors (3,4,5)] output=[32 bits] normalize dynamic="+str(dyn)+" connectivity=6");
	smallBioporeswatershed = IJ.getImage()	
	IJ.run(smallBioporesRemoved,"Analyze Regions 3D", "equivalent_ellipsoid ellipsoid_elongations surface_area_method=[Crofton (13 dirs.)] euler_connectivity=26");	
	IJ.renameResults("biodist-watershed-morpho", "Results");
	rt = ResultsTable.getResultsTable();
	ElliR3 = rt.getColumn(5)
	ElliR2 = rt.getColumn(4)
	ElliR1 = rt.getColumn(3)
	Rb = []
	Ra = []
	Vesselness = []	
	for i in range(len(ElliR3)):
		Rb.append(ElliR3[i]/cmath.sqrt(ElliR2[i]*ElliR1[i])); 
		Ra.append(ElliR2[i]/ElliR1[i]); 
		Vesselness.append(abs(cmath.exp(-Rb[i]*Rb[i])*cmath.exp(-Ra[i]*Ra[i]))); 
		if (Vesselness[i] >ves)and (ElliR1[i] > length):
			rt.setValue("flag",i ,255) 
		else:
			rt.setValue("flag",i ,0) 	
	rt.updateResults()
	rt.show("flaged")	
	IJ.run("Assign Measure to Label", "results="+"flaged"+" column="+"flag"+" min=0 max=255"); 	
	IJ.run("8-bit", "");
	smallBioporesRemoved = IJ.getImage(smallBioporesRemoved)
	smallBioporesRemoved.setTitle("smallBioporesRemoved")	
	smallBioporeswatershed.close()
	image.close()
	return smallBioporesRemoved

gray= IJ.getImage()
IJ.run(gray, "Unsharp Mask...", "radius="+str(radius)+" mask="+str(mask)+" stack");
if Threshold_value == "Otsu":
	IJ.setAutoThreshold(gray,"Otsu dark stack")
	gray.show()
	ImProc = gray.getProcessor()
	upper_thr = ImProc.getMinThreshold()
else:
	upper_thr = Threshold
IJ.run("Threshold...");
ThresholdAdjuster.setMode("Over/Under");	
lower_thr = round(upper_thr/Threshold_ratio);
IJ.setThreshold(gray,lower_thr,upper_thr);
ij.Prefs.blackBackground = False;
IJ.run("Convert to Mask", "method=Default background=Dark black");
thresh = IJ.getImage()
thresh.setTitle("thresh");
total_width, total_height, ch , total_slices, fr = thresh.getDimensions()
scale05 = thresh.resize(total_width/2,total_height/2,  total_slices/2, "none");
dyn =2;
smallBiopores = segmentsmallBiopores(scale05)
scale02 = thresh.resize(total_width/5,total_height/5,  total_slices/5, "none");

largeBioporesfinal1 = segmentlargeBiopores(scale02)
Bioporesfinal3 = ImageCalculator.run(largeBioporesfinal1,smallBiopores, "ADD create stack");
Bioporesfinal3.show();
Bioporesfinal3 = IJ.getImage();
Bioporesfinal3.setTitle("Bioporesfinal3");
t_s_small=t_s*0.5
smallBioporesRemoved= RemoveBloobs(Bioporesfinal3, t_v, t_s_small, dyn)
smallBioporesRemoved.show()
Roots = IJ.getImage();
Roots.setTitle("Roots");

smallBiopores.close()
largeBioporesfinal1.close()

output = (Roots)
#IJ.run("Close All", "")




#print Tube


