import os, sys, shutil, arcpy
import traceback, time
from arcpy.sa import *
from arcpy import env


arcpy.CheckOutExtension("Spatial")

DEBUGGING = False

def log(message):
    arcpy.AddMessage(message)
    with file(sys.argv[0]+".log", 'a') as logFile:
        logFile.write("%s:\t%s\n" % (time.asctime(), message))

class Toolbox(object):
    def __init__(self):
        """Define the toolbox (the name of the toolbox is the name of the
        .pyt file)."""
        self.label = "Toolbox"
        self.alias = ""

        # List of tool classes associated with this toolbox
        self.tools = [Test, TopoHydro,ImpCov,Runoff,GetNEXRAD,ScenarioAnalysis]

class Scenarioforeachfield(object):
    def __init__(self):
        self.label = "Scenarioforeachfield"
        self.description = "Scenario analysis for each field"
        self.canRunInBackground = False
        
        arcpy.env.Workspace = self.Workspace = os.path.split(__file__)[0]
        log("Workspace = " + arcpy.env.Workspace)
        arcpy.env.overwriteOutput = True       

    def getParameterInfo(self):
        """Define parameter definitions"""
        
        param0 = arcpy.Parameter(
            displayName="BMP Points",
            name="bmppts",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input",
            multiValue=False)  
        
        param1 = arcpy.Parameter(
            displayName="Status Field",
            name="statusField",
            datatype="Field",
            parameterType="Required",
            direction="Input",
            multiValue=False)  

        
        params = [ param0, param1 ]
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return
            
    def execute(self, parameters, messages):
        try:
            BMPpath = parameters[0].valueAsText
            BMPs = "BMPs"
            scen = parameters[1].valueAsText
            arcpy.MakeFeatureLayer_management(BMPpath, BMPs)


            log("Parameter is %s" % (parameters[0].valueAsText))
            Lakes = "Lakes"
            LanduseExisting = "LanduseExisting"
            Impervious = "Impervious"
            BMPs = parameters[0].valueAsText
            inFeatures = [Lakes, LanduseExisting, Impervious]
            unionE = "unionE"
            arcpy.Union_analysis (inFeatures,unionE,"ALL", "", "GAPS")

            arcpy.AddField_management("unionE", "TotalNitrogen", "FLOAT")
            arcpy.AddField_management("unionE", "TotalPhosphorus", "FLOAT")
            arcpy.AddField_management("unionE", "Sediment", "FLOAT")
            arcpy.AddField_management("unionE", "Copper", "FLOAT")
            arcpy.AddField_management("unionE", "FecalColiform", "FLOAT")
            arcpy.AddField_management("unionE", "Zinc", "FLOAT")    

            arcpy.SelectLayerByAttribute_management(unionE, "NEW_SELECTION", "Reclassify = 'Roadways'")
            arcpy.CalculateField_management(unionE, "TotalNitrogen", "12.2", "VB", "")
            arcpy.CalculateField_management(unionE, "TotalPhosphorus", "1.8", "VB", "")
            arcpy.CalculateField_management(unionE, "Sediment", "405", "VB", "")
            arcpy.CalculateField_management(unionE, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionE, "Zinc", "0.1", "VB", "")
            arcpy.CalculateField_management(unionE, "FecalColiform", "2.8", "VB", "")


            arcpy.SelectLayerByAttribute_management(unionE, "NEW_SELECTION", "Reclassify in ( 'Commercial' )")
            arcpy.CalculateField_management(unionE, "TotalNitrogen", "14", "VB", "")
            arcpy.CalculateField_management(unionE, "TotalPhosphorus", "2.7", "VB", "")
            arcpy.CalculateField_management(unionE, "Sediment", "400", "VB", "")
            arcpy.CalculateField_management(unionE, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionE, "Zinc", "0", "VB", "")
            arcpy.CalculateField_management(unionE, "FecalColiform", "8.4", "VB", "")


            arcpy.SelectLayerByAttribute_management(unionE, "NEW_SELECTION", "Reclassify IN ( 'Industrial' )")
            arcpy.CalculateField_management(unionE, "TotalNitrogen", "10.4", "VB", "")
            arcpy.CalculateField_management(unionE, "TotalPhosphorus", "1.9", "VB", "")
            arcpy.CalculateField_management(unionE, "Sediment", "372.5", "VB", "")
            arcpy.CalculateField_management(unionE, "Copper", "0.2", "VB", "")
            arcpy.CalculateField_management(unionE, "Zinc", "0.1", "VB", "")
            arcpy.CalculateField_management(unionE, "FecalColiform", "8.4", "VB", "")

            arcpy.SelectLayerByAttribute_management(unionE, "NEW_SELECTION", "Reclassify in ( 'Institutional' , 'Research Triangle Park' )")
            arcpy.CalculateField_management(unionE, "TotalNitrogen", "9.4", "VB", "")
            arcpy.CalculateField_management(unionE, "TotalPhosphorus", "0.7", "VB", "")
            arcpy.CalculateField_management(unionE, "Sediment", "50", "VB", "")
            arcpy.CalculateField_management(unionE, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionE, "Zinc", "0", "VB", "")
            arcpy.CalculateField_management(unionE, "FecalColiform", "14.8", "VB", "")

            arcpy.SelectLayerByAttribute_management(unionE, "NEW_SELECTION", "Reclassify IN ( 'Parks and Open Space' , 'Agricultural' )")
            arcpy.CalculateField_management(unionE, "TotalNitrogen", "2.3", "VB", "")
            arcpy.CalculateField_management(unionE, "TotalPhosphorus", "0.1", "VB", "")
            arcpy.CalculateField_management(unionE, "Sediment", "10", "VB", "")
            arcpy.CalculateField_management(unionE, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionE, "Zinc", "0", "VB", "")
            arcpy.CalculateField_management(unionE, "FecalColiform", "12", "VB", "")

            arcpy.SelectLayerByAttribute_management(unionE, "NEW_SELECTION", "Reclassify IN ( 'High Density Residential' , 'Medium Density Residential' )")
            arcpy.CalculateField_management(unionE, "TotalNitrogen", "11.2", "VB", "")
            arcpy.CalculateField_management(unionE, "TotalPhosphorus", "1.6", "VB", "")
            arcpy.CalculateField_management(unionE, "Sediment", "242.5", "VB", "")
            arcpy.CalculateField_management(unionE, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionE, "Zinc", "0", "VB", "")
            arcpy.CalculateField_management(unionE, "FecalColiform", "30.2", "VB", "")

            arcpy.SelectLayerByAttribute_management(unionE, "NEW_SELECTION", "Reclassify IN( 'Very Low Density Residential' , 'Low Density Residential' )")
            arcpy.CalculateField_management(unionE, "TotalNitrogen", "6.4", "VB", "")
            arcpy.CalculateField_management(unionE, "TotalPhosphorus", "0.7", "VB", "")
            arcpy.CalculateField_management(unionE, "Sediment", "150", "VB", "")
            arcpy.CalculateField_management(unionE, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionE, "Zinc", "0", "VB", "")
            arcpy.CalculateField_management(unionE, "FecalColiform", "16.2", "VB", "")

            arcpy.SelectLayerByAttribute_management(unionE, "CLEAR_SELECTION", "")

            LanduseFuture = "LanduseFuture"

            inFeatures1 = [Lakes, LanduseFuture, Impervious]
            unionF = "unionF"
            arcpy.Union_analysis (inFeatures1,unionF,"ALL", "", "GAPS")
            arcpy.AddField_management("unionF", "TotalNitrogen", "FLOAT")
            arcpy.AddField_management("unionF", "TotalPhosphorus", "FLOAT")
            arcpy.AddField_management("unionF", "Sediment", "FLOAT")
            arcpy.AddField_management("unionF", "Copper", "FLOAT")
            arcpy.AddField_management("unionF", "FecalColiform", "FLOAT")
            arcpy.AddField_management("unionF", "Zinc", "FLOAT")

            arcpy.SelectLayerByAttribute_management(unionF, "NEW_SELECTION", "Reclassify = 'Roadways'")
            arcpy.CalculateField_management(unionF, "TotalNitrogen", "12.2", "VB", "")
            arcpy.CalculateField_management(unionF, "TotalPhosphorus", "1.8", "VB", "")
            arcpy.CalculateField_management(unionF, "Sediment", "405", "VB", "")
            arcpy.CalculateField_management(unionF, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionF, "Zinc", "0.1", "VB", "")
            arcpy.CalculateField_management(unionF, "FecalColiform", "2.8", "VB", "")


            arcpy.SelectLayerByAttribute_management(unionF, "NEW_SELECTION", "Reclassify in ( 'Commercial' )")
            arcpy.CalculateField_management(unionF, "TotalNitrogen", "14", "VB", "")
            arcpy.CalculateField_management(unionF, "TotalPhosphorus", "2.7", "VB", "")
            arcpy.CalculateField_management(unionF, "Sediment", "400", "VB", "")
            arcpy.CalculateField_management(unionF, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionF, "Zinc", "0", "VB", "")
            arcpy.CalculateField_management(unionF, "FecalColiform", "8.4", "VB", "")


            arcpy.SelectLayerByAttribute_management(unionF, "NEW_SELECTION", "Reclassify IN ( 'Industrial' )")
            arcpy.CalculateField_management(unionF, "TotalNitrogen", "10.4", "VB", "")
            arcpy.CalculateField_management(unionF, "TotalPhosphorus", "1.9", "VB", "")
            arcpy.CalculateField_management(unionF, "Sediment", "372.5", "VB", "")
            arcpy.CalculateField_management(unionF, "Copper", "0.2", "VB", "")
            arcpy.CalculateField_management(unionF, "Zinc", "0.1", "VB", "")
            arcpy.CalculateField_management(unionF, "FecalColiform", "8.4", "VB", "")

            arcpy.SelectLayerByAttribute_management(unionF, "NEW_SELECTION", "Reclassify in ( 'Institutional' , 'Research Triangle Park' )")
            arcpy.CalculateField_management(unionF, "TotalNitrogen", "9.4", "VB", "")
            arcpy.CalculateField_management(unionF, "TotalPhosphorus", "0.7", "VB", "")
            arcpy.CalculateField_management(unionF, "Sediment", "50", "VB", "")
            arcpy.CalculateField_management(unionF, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionF, "Zinc", "0", "VB", "")
            arcpy.CalculateField_management(unionF, "FecalColiform", "14.8", "VB", "")

            arcpy.SelectLayerByAttribute_management(unionF, "NEW_SELECTION", "Reclassify IN ( 'Parks and Open Space' , 'Agricultural' )")
            arcpy.CalculateField_management(unionF, "TotalNitrogen", "2.3", "VB", "")
            arcpy.CalculateField_management(unionF, "TotalPhosphorus", "0.1", "VB", "")
            arcpy.CalculateField_management(unionF, "Sediment", "10", "VB", "")
            arcpy.CalculateField_management(unionF, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionF, "Zinc", "0", "VB", "")
            arcpy.CalculateField_management(unionF, "FecalColiform", "12", "VB", "")

            arcpy.SelectLayerByAttribute_management(unionF, "NEW_SELECTION", "Reclassify IN ( 'High Density Residential' , 'Medium Density Residential' )")
            arcpy.CalculateField_management(unionF, "TotalNitrogen", "11.2", "VB", "")
            arcpy.CalculateField_management(unionF, "TotalPhosphorus", "1.6", "VB", "")
            arcpy.CalculateField_management(unionF, "Sediment", "242.5", "VB", "")
            arcpy.CalculateField_management(unionF, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionF, "Zinc", "0", "VB", "")
            arcpy.CalculateField_management(unionF, "FecalColiform", "30.2", "VB", "")

            arcpy.SelectLayerByAttribute_management(unionF, "NEW_SELECTION", "Reclassify IN( 'Very Low Density Residential' , 'Low Density Residential' )")
            arcpy.CalculateField_management(unionF, "TotalNitrogen", "6.4", "VB", "")
            arcpy.CalculateField_management(unionF, "TotalPhosphorus", "0.7", "VB", "")
            arcpy.CalculateField_management(unionF, "Sediment", "150", "VB", "")
            arcpy.CalculateField_management(unionF, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionF, "Zinc", "0", "VB", "")
            arcpy.CalculateField_management(unionF, "FecalColiform", "16.2", "VB", "")

            arcpy.SelectLayerByAttribute_management(unionF, "CLEAR_SELECTION", "")
           
            Raster_Nitrogen = "Raster_Nitrogen"
            Raster_Phosphorus = "Raster_Phosphorus"
            Raster_Sediment = "Raster_Sediment"
            Raster_Copper = "Raster_Copper"
            Raster_Zinc = "Raster_Zinc"
            Raster_Fecal = "Raster_Fecal"
            arcpy.PolygonToRaster_conversion(unionE, "TotalNitrogen", Raster_Nitrogen, "CELL_CENTER", "NONE", "40")
            arcpy.PolygonToRaster_conversion(unionE, "TotalPhosphorus", Raster_Phosphorus, "CELL_CENTER", "NONE", "40")
            arcpy.PolygonToRaster_conversion(unionE, "Sediment", Raster_Sediment, "CELL_CENTER", "NONE", "40")
            arcpy.PolygonToRaster_conversion(unionE, "Copper", Raster_Copper, "CELL_CENTER", "NONE", "40")
            arcpy.PolygonToRaster_conversion(unionE, "Zinc", Raster_Zinc, "CELL_CENTER", "NONE", "40")
            arcpy.PolygonToRaster_conversion(unionE, "FecalColiform", Raster_Fecal, "CELL_CENTER", "NONE", "40")

            RasterF_Nitrogen = "RasterF_Nitrogen"
            RasterF_Phosphorus = "RasterF_Phosphorus"
            RasterF_Sediment = "RasterF_Sediment"
            RasterF_Copper = "RasterF_Copper"
            RasterF_Zinc = "RasterF_Zinc"
            RasterF_Fecal = "RasterF_Fecal"

            arcpy.PolygonToRaster_conversion(unionF, "TotalNitrogen", RasterF_Nitrogen, "CELL_CENTER", "NONE", "40")
            arcpy.PolygonToRaster_conversion(unionF, "TotalPhosphorus", RasterF_Phosphorus, "CELL_CENTER", "NONE", "40")
            arcpy.PolygonToRaster_conversion(unionF, "Sediment", RasterF_Sediment, "CELL_CENTER", "NONE", "40")
            arcpy.PolygonToRaster_conversion(unionF, "Copper", RasterF_Copper, "CELL_CENTER", "NONE", "40")
            arcpy.PolygonToRaster_conversion(unionF, "Zinc", RasterF_Zinc, "CELL_CENTER", "NONE", "40")
            arcpy.PolygonToRaster_conversion(unionF, "FecalColiform", RasterF_Fecal, "CELL_CENTER", "NONE", "40")

            demPath = "DEM"
            dem = Raster(demPath)
            fill = Fill(dem)
            flowDirection = FlowDirection(fill)
            flowAccEN = FlowAccumulation(flowDirection,Raster_Nitrogen, "FLOAT")
            flowAccEP = FlowAccumulation(flowDirection,Raster_Phosphorus, "FLOAT")
            flowAccEN.save("C:/ss/output/flowAccEN1.tif")
            flowAccEP.save("C:/ss/output/flowAccEP1.tif")
            flowAccES = FlowAccumulation(flowDirection,Raster_Sediment, "FLOAT")
            flowAccEC = FlowAccumulation(flowDirection,Raster_Copper, "FLOAT")
            flowAccES.save("C:/ss/output/flowAccES1.tif")
            flowAccEC.save("C:/ss/output/flowAccEC1.tif")
            flowAccEZ = FlowAccumulation(flowDirection,Raster_Zinc, "FLOAT")
            flowAccEF = FlowAccumulation(flowDirection,Raster_Fecal, "FLOAT")
            flowAccEZ.save("C:/ss/output/flowAccEZ1.tif")
            flowAccEF.save("C:/ss/output/flowAccEF1.tif")

            FlowAcc3_Nitrogen = FlowAccumulation(flowDirection,RasterF_Nitrogen, "FLOAT")
            FlowAcc3_Phosphorus = FlowAccumulation(flowDirection,RasterF_Phosphorus, "FLOAT")
            FlowAcc3_Sediment = FlowAccumulation(flowDirection,RasterF_Sediment, "FLOAT")
            FlowAcc3_Copper = FlowAccumulation(flowDirection,RasterF_Copper, "FLOAT")
            FlowAcc3_Zinc = FlowAccumulation(flowDirection,RasterF_Zinc, "FLOAT")
            FlowAcc3_Fecal = FlowAccumulation(flowDirection,RasterF_Fecal, "FLOAT")



            BMP_TN = "BMP_TN"
            BMP_TP = "BMP_TP"
            BMP_FC = "BMP_FC"
            BMP_CU = "BMP_CU"
            BMP_Zn = "BMP_Zn"
            BMP_Sed = "BMP_Sed"

            BMPpath = parameters[0].valueAsText
            BMPs = "BMPs"
            scen = parameters[1].valueAsText
            arcpy.MakeFeatureLayer_management(BMPpath, BMPs)
            arcpy.SelectLayerByAttribute_management(BMPs, "NEW_SELECTION", " %s = 'TRUE'" %scen)
            arcpy.PointToRaster_conversion(BMPs, "TN_Eff_Ex", BMP_TN, "MOST_FREQUENT", "NONE", "40")
            arcpy.PointToRaster_conversion(BMPs, "TP_Eff_Ex", BMP_TP, "MOST_FREQUENT", "NONE", "40")
            arcpy.PointToRaster_conversion(BMPs, "FC_Eff_Ex", BMP_FC, "MOST_FREQUENT", "NONE", "40")
            arcpy.PointToRaster_conversion(BMPs, "CU_Eff_Ex", BMP_CU, "MOST_FREQUENT", "NONE", "40")
            arcpy.PointToRaster_conversion(BMPs, "Zn_Eff_Ex", BMP_Zn, "MOST_FREQUENT", "NONE", "40")
            arcpy.PointToRaster_conversion(BMPs, "Sed_Eff_Ex", BMP_Sed, "MOST_FREQUENT", "NONE", "40")

            FlowAcc3_Nitrogen_BMP = FlowAccumulation(flowDirection, BMP_TN, "FLOAT")
            FlowAcc3_Phosphorus_BMP = FlowAccumulation(flowDirection, BMP_TP, "FLOAT")
            FlowAcc3_Fecal_BMP = FlowAccumulation(flowDirection, BMP_FC, "FLOAT")
            FlowAcc3_Sediment_BMP = FlowAccumulation(flowDirection, BMP_Sed, "FLOAT")
            FlowAcc3_Copper_BMP = FlowAccumulation(flowDirection, BMP_CU, "FLOAT")
            FlowAcc3_Zinc_BMP = FlowAccumulation(flowDirection, BMP_Zn, "FLOAT")



            maxN = Minus(692.2000122070313,FlowAcc3_Nitrogen_BMP)
            maxP = Minus(1082.5,FlowAcc3_Phosphorus_BMP)
            maxS = Minus(2205.800048828125,FlowAcc3_Sediment_BMP)
            maxF = Minus(1503,FlowAcc3_Fecal_BMP)
            maxC = Minus(1074,FlowAcc3_Copper_BMP)
            maxZ = Minus(1382,FlowAcc3_Zinc_BMP)

            BMPred_N = Divide(maxN, 692.2000122070313)
            BMPred_P = Divide(maxP, 1082.5)
            BMPred_se = Divide(maxS, 2205.800048828125)
            BMPred_fe = Divide(maxF, 1503)
            BMPred_Copper = Divide(maxC, 1074)
            BMPred_Z = Divide(maxZ,1382)

            futureSe = Times(FlowAcc3_Sediment,BMPred_se)
            futureFe = Times(FlowAcc3_Fecal, BMPred_fe)
            futureZ = Times(FlowAcc3_Zinc, BMPred_Z)
            futureC = Times(FlowAcc3_Copper, BMPred_Copper)
            futureP = Times(BMPred_P, FlowAcc3_Phosphorus)
            futureN = Times(BMPred_N, FlowAcc3_Nitrogen)
            futureSe.save("C:/ss/output/futureSe1.tif")
            futureFe.save("C:/ss/output/futureFe1.tif")
            futureZ.save("C:/ss/output/futureZ1.tif")
            futureC.save("C:/ss/output/futureC1.tif")
            futureP.save("C:/ss/output/futureP1.tif")
            futureN.save("C:/ss/output/futureN1.tif")

            log("final")
            Nfinal = Con(IsNull(flowAccEN),0,Divide(futureN, flowAccEN))
            Pfinal = Con(IsNull(flowAccEP),0,Divide(futureP, flowAccEP))
            Sfinal = Con(IsNull(flowAccES),0,Divide(futureSe, flowAccES))
            Zfinal = Con(IsNull(flowAccEZ),0,Divide(futureZ, flowAccEZ))
            Cfinal = Con(IsNull(flowAccEC),0,Divide(futureC, flowAccEC))
            Ffinal = Con(IsNull(flowAccEF),0,Divide(futureFe, flowAccEF))
            Nfinal.save("C:/ss/output/Nfinal.tif")
            Pfinal.save("C:/ss/output/Pfinal.tif")
            Sfinal.save("C:/ss/output/Sfinal.tif")
            Zfinal.save("C:/ss/output/Zfinal.tif")
            Cfinal.save("C:/ss/output/Cfinal.tif")
            Ffinal.save("C:/ss/output/Ffinal.tif")
            log("reclassify")
            maxN = 1146.75
            maxP = 3889
            maxS = 8749
            maxF = 299.334
            maxZ = 3685
            maxC = 494
            Reclass_N = Reclassify(Nfinal, "Value", "0 %f NODATA;%f %f 1;%f %f 2;%f %f 3;%f %f 4" %(0.05*maxN, 0.05*maxN, 0.1 * maxN, 0.1 * maxN, 0.25 * maxN, 0.25 * maxN, 0.5 * maxN, 0.5 * maxN, maxN), "DATA")
            Reclass_P = Reclassify(Pfinal, "Value", "0 %f NODATA;%f %f 1;%f %f 2;%f %f 3;%f %f 4" %(0.05*maxP, 0.05*maxP, 0.1 * maxP, 0.1 * maxP, 0.25 * maxP, 0.25 * maxP, 0.5 * maxP, 0.5 * maxP, maxP), "DATA")
            Reclass_S = Reclassify(Sfinal, "Value", "0 %f NODATA;%f %f 1;%f %f 2;%f %f 3;%f %f 4" %(0.05*maxS, 0.05*maxS, 0.1 * maxS, 0.1 * maxS, 0.25 * maxS, 0.25 * maxS, 0.5 * maxS, 0.5 * maxS, maxS), "DATA")
            Reclass_F = Reclassify(Ffinal, "Value", "0 %f NODATA;%f %f 1;%f %f 2;%f %f 3;%f %f 4" %(0.05*maxF, 0.05*maxF, 0.1 * maxF, 0.1 * maxF, 0.25 * maxF, 0.25 * maxF, 0.5 * maxF, 0.5 * maxF, maxF), "DATA")
            Reclass_Z = Reclassify(Zfinal, "Value", "0 %f NODATA;%f %f 1;%f %f 2;%f %f 3;%f %f 4" %(0.05*maxZ, 0.05*maxZ, 0.1 * maxZ, 0.1 * maxZ, 0.25 * maxZ, 0.25 * maxZ, 0.5 * maxZ, 0.5 * maxZ, maxZ), "DATA")
            Reclass_C = Reclassify(Cfinal, "Value", "0 %f NODATA;%f %f 1;%f %f 2;%f %f 3;%f %f 4" %(0.05*maxC, 0.05*maxC, 0.1 * maxC, 0.1 * maxC, 0.25 * maxC, 0.25 * maxC, 0.5 * maxC, 0.5 * maxC, maxC), "DATA")
            StreamT_N = "StreamT_N"
            StreamT_P = "StreamT_P"
            StreamT_S = "StreamT_S"
            StreamT_F = "StreamT_F"
            StreamT_Z = "StreamT_Z"
            StreamT_C = "StreamT_C"
            log("stream to feature")
            arcpy.gp.StreamToFeature_sa(Reclass_N, flowDirection, StreamT_N, "SIMPLIFY")
            arcpy.gp.StreamToFeature_sa(Reclass_P, flowDirection, StreamT_P, "SIMPLIFY")
            arcpy.gp.StreamToFeature_sa(Reclass_S, flowDirection, StreamT_S, "SIMPLIFY")
            arcpy.gp.StreamToFeature_sa(Reclass_F, flowDirection, StreamT_F, "SIMPLIFY")
            arcpy.gp.StreamToFeature_sa(Reclass_Z, flowDirection, StreamT_Z, "SIMPLIFY")
            arcpy.gp.StreamToFeature_sa(Reclass_C, flowDirection, StreamT_C, "SIMPLIFY")


			
        except Exception as err:
            log(traceback.format_exc())
            log(err)
            raise err
        return

class TopoHydro(object):
    def __init__(self):
        self.label = "Topography and Hydrology Analysis"
        self.description = "Establishes the watershed and stream network"
        self.canRunInBackground = False
        
        arcpy.env.Workspace = self.Workspace = os.path.split(__file__)[0]
        log("Workspace = " + arcpy.env.Workspace)
        arcpy.env.overwriteOutput = True       

    def getParameterInfo(self):
        """Define parameter definitions"""
        
        param0 = arcpy.Parameter(
            displayName="Input Digital Elevation Model",
            name="DEM",
            datatype="DERasterDataset",
            parameterType="Required",
            direction="Input",
            multiValue=False)  
            
        param1 = arcpy.Parameter(
            displayName="Analysis Mask",
            name="Mask",
            datatype="DEFeatureClass",
            parameterType="Optional",
            direction="Input",
            multiValue=False)  
        
        param2 = arcpy.Parameter(
            displayName="Threshold accumulation for Stream formation (acres)",
            name="StreamFormation",
            datatype="GPDouble",
            parameterType="Required",
            direction="Input",
            multiValue=False)  
            
        param3 = arcpy.Parameter(
            displayName="Existing vector stream to use to modify drainage",
            name="ExistingStreams",
            datatype="DEFeatureClass",
            parameterType="Optional",
            direction="Input",
            multiValue=False)
        
        params = [ param0, param1, param2, param3 ]
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return
            
    def execute(self, parameters, messages):
        try:
            demPath = parameters[0].valueAsText
            
            arcpy.env.extent = demPath
            arcpy.env.snapRaster = parameters[1].valueAsText
            arcpy.env.cellSize = demPath

            dem = Raster(demPath)
            fill = Fill(dem)
            flowd = FlowDirection(fill)
            flowacc = FlowAccumulation(flowd,"", "FLOAT")
            reclassify = Reclassify(flowacc, "Value", "0 51689 1;51689 169257 2;169257 411422 3;411422 770795 4","DATA")
            ## should set the workspace first, as geoprocessing env variable
            reclassify.save("C:/ss/output/reclassify.tif")
            if DEBUGGING: flowAccumulation.save("flowaccumulation.tif") 

			
        except Exception as err:
            log(traceback.format_exc())
            log(err)
            raise err
        return

class ImpCov(object):
    def __init__(self):
        self.label = "Imperviousness Analysis"
        self.description = "Impervious area contributions"
        self.canRunInBackground = False
        
        arcpy.env.Workspace = self.Workspace = os.path.split(__file__)[0]
        log("Workspace = " + arcpy.env.Workspace)
        arcpy.env.overwriteOutput = True       

    def getParameterInfo(self):
        """Define parameter definitions"""
        
        param0 = arcpy.Parameter(
            displayName="Impervious Areas",
            name="ImperviousAreas",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input",
            multiValue=False)  
            
        param1 = arcpy.Parameter(
            displayName="Lakes",
            name="Lakes",
            datatype="DEFeatureClass",
            parameterType="Optional",
            direction="Input",
            multiValue=False)  
       
        params = [ param0, param1 ]
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return
            
    def execute(self, parameters, messages):
        try:
            log("Parameters are %s, %s" % (parameters[0].valueAsText, parameters[1].valueAsText))

            dem = "DEM"
            fill = Fill(dem)
            flowd = FlowDirection(fill)
            flowaccum = FlowAccumulation(flowd,"", "FLOAT")
            rastercalc = Times(flowaccum, 0.036730458)
            reclass = Reclassify(rastercalc, "Value", "0 6997 NODATA;6997 28312.03125 1", "DATA")

            # Local variables: define these by finding them in your workspace or getting them as tool parameters
            Impervious = parameters[0].valueAsText
	    arcpy.CalculateField_management(Impervious, "LENGTH", "1", "VB", "")
            Imper = "Imper"
            arcpy.FeatureToRaster_conversion(Impervious, "LENGTH", Imper, "4")
            rasterblock = BlockStatistics(Imper, "Rectangle 10 10 CELL", "SUM", "DATA")
            rasteragg = Aggregate(rasterblock, "10", "MEAN", "EXPAND", "DATA")

            iaccum = FlowAccumulation(flowd, rasteragg, "FLOAT")
            Divide1 = Divide(iaccum, flowaccum)
            impervreclass = Reclassify(Divide1, "Value", "0 10 1;10 20 2;20 30 3;30 40 4;40 50 5;50 60 6;60 70 7;70 80 8;80 90 9;90 100 10", "DATA")
            imper_mult1 = Times(impervreclass, reclass)		
            StreamToFeature(imper_mult1, flowd, "Stream_Task3Imp", "SIMPLIFY")
            imper_mult1.save("C:/ss/output/imper_mult1.tif") 

        except Exception as err:
            log(traceback.format_exc())
            log(err)
            raise err
        return

class Runoff(object):
    def __init__(self):
        self.label = "Runoff Calculations"
        self.description = "Calculation of standard storm flows via USGS regression equations"
        self.canRunInBackground = False
        
        arcpy.env.Workspace = self.Workspace = os.path.split(__file__)[0]
        log("Workspace = " + arcpy.env.Workspace)
        arcpy.env.overwriteOutput = True       

    def getParameterInfo(self):
        """Define parameter definitions"""
        
        param0 = arcpy.Parameter(
            displayName="Impervious Areas",
            name="ImperviousAreas",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input",
            multiValue=False)  
            
        param1 = arcpy.Parameter(
            displayName="Lakes",
            name="Lakes",
            datatype="DEFeatureClass",
            parameterType="Optional",
            direction="Input",
            multiValue=False)  

        param2 = arcpy.Parameter(
            displayName="Curve Number",
            name="Landuse",
            datatype="DEFeatureClass",
            #parameterType="Required",
            parameterType="Optional",
            direction="Input",
            multiValue=False)  
               
        params = [ param0, param1, param2 ]
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return
            
    def execute(self, parameters, messages):
        try:
            log("Parameters are %s, %s" % (parameters[0].valueAsText, parameters[1].valueAsText))

            dem = "DEM"
            fill = Fill(dem)
            flowd = FlowDirection(fill)
            flowaccum = FlowAccumulation(flowd,"", "FLOAT")
            rastercalc = Times(flowaccum, 0.036730458)
            reclass = Reclassify(rastercalc, "Value", "0 6997 NODATA;6997 28312.03125 1", "DATA")

            Impervious = parameters[0].valueAsText
	        arcpy.CalculateField_management(Impervious, "LENGTH", "1", "VB", "")
            Imper_rast = "Imper_rast"
            arcpy.FeatureToRaster_conversion(Impervious, "LENGTH", Imper_rast, "4")
            block_rast = BlockStatistics(Imper_rast, "Rectangle 10 10 CELL", "SUM", "DATA")
            agg_rast = Aggregate(block_rast, "10", "MEAN", "EXPAND", "DATA")

            imper_accum = FlowAccumulation(flowd, agg_rast, "FLOAT")
            Divide1 = Divide(imper_accum, flowaccum)
            reclass_imperv = Reclassify(Divide1, "Value", "0 10 1;10 20 2;20 30 3;30 40 4;40 50 5;50 60 6;60 70 7;70 80 8;80 90 9;90 100 10", "DATA")
            imper_mult1 = Times(reclass_imperv, reclass)		
            StreamToFeature(imper_mult1, flowd, "Stream_Task3Imp", "SIMPLIFY")

            flowaccum_sqmi=Times(rastercalc, 0.0015625)
            P_2=Power(flowaccum_sqmi,0.691)
            P_5=Power(flowaccum_sqmi,0.670)
            P_10=Power(flowaccum_sqmi,0.665)
            P_25=Power(flowaccum_sqmi,0.655)
            P_50=Power(flowaccum_sqmi,0.650)
            P_100=Power(flowaccum_sqmi,0.643)
            T_2=Times(P_2,144)
            T_5=Times(P_5,248)
            T_10=Times(P_10,334)
            T_25=Times(P_25,467)
            T_50=Times(P_50,581)
            T_100=Times(P_100,719)
            PO_2=Power(T_2,0.338)
            PO_5=Power(T_5,0.338)
            PO_10=Power(T_10,0.338)
            PO_25=Power(T_25,0.338)
            PO_50=Power(T_50,0.338)
            PO_100=Power(T_100,0.338)
            p1=Power(Divide1,0.436)
            p2=Power(flowaccum_sqmi,0.390)
            t1=Times(p2,28.5)
            t2=Times(t1,p1)
            recur_2I=Times(t2,PO_2)
            recur_5I=Times(t2,PO_5)
            recur_10I=Times(t2,PO_10)
            recur_25I=Times(t2,PO_25)
            recur_50I=Times(t2,PO_50)
            recur_100I=Times(t2,PO_100)


        except Exception as err:
            log(traceback.format_exc())
            log(err)
            raise err
        return

class GetNEXRAD(object):
    def __init__(self):
        self.label = "Get NEXRAD rainfall"
        self.description = "Get a raster of rainfall for a specific rain event from NEXRAD weather radar"
        self.canRunInBackground = False
        
        arcpy.env.Workspace = self.Workspace = os.path.split(__file__)[0]
        log("Workspace = " + arcpy.env.Workspace)
        arcpy.env.overwriteOutput = True       

    def getParameterInfo(self):
        """Define parameter definitions"""
        
        param0 = arcpy.Parameter(
            displayName="Start Date",
            name="startDate",
            datatype="GPDate",
            parameterType="Required",
            direction="Input",
            multiValue=False)  
        
        param1 = arcpy.Parameter(
            displayName="End Date",
            name="endDate",
            datatype="GPDate",
            parameterType="Required",
            direction="Input",
            multiValue=False)  
        
        param2 = arcpy.Parameter(
            displayName="Radar Station ID",
            name="radarID",
            datatype="GPString",
            parameterType="Required",
            direction="Input",
            multiValue=False)  
            
        params = [ param0, param1, param2 ]
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return
            
    def execute(self, parameters, messages):
        try:
            log("Parameter is %s" % (parameters[0].valueAsText))
            
            # code for vector -> raster from Tyler Pitts
            lie=arcpy.CheckOutExtension ('Spatial')
            log ( lie )
            shapefiles=[]
            stormfolder=1
            while stormfolder<=4:
                #log ( 'checking folder',stormfolder )
                shapefiles=[]
                for root, dirs, files in os.walk('storm'+str(stormfolder)):
                    for file in files:
                        if file.endswith('.shp'):
                                shapefiles.append('storm'+str(stormfolder)+'/'+file)
                log ( 'done creating an array of the shapefiles' )
                log ( 'converting to rasters' )
                rasters=[]
                for x in range(len(shapefiles)):
                    #Slog ( 'converting',shapefiles[x] )
                    raster=arcpy.PolygonToRaster_conversion(shapefiles[x], 'value', 'storm'+str(stormfolder)+'/raster'+str(x), 'CELL_CENTER', 'NONE',0.00012196015)
                    rasters.append(raster)
                log ( 'completed raster conversion' )
                log ( 'calculating cell statistics' )
                maxreflect=CellStatistics (rasters, 'MAXIMUM', 'DATA')
                maxreflect.save('storm'+str(stormfolder)+'/reflect'+str(stormfolder)+'.tif')
                lowerLeft = arcpy.Point(maxreflect.extent.XMin,maxreflect.extent.YMin)
                cellSize = maxreflect.meanCellWidth
                reflectence=arcpy.RasterToNumPyArray(maxreflect)
                rows=len(reflectence)
                cols=len(reflectence[0])
                rainfallraster=numpy.zeros((rows,cols))
                for row in range(rows):
                    for col in range(cols):
                        if reflectence[row][col]<0:
                            rainfallraster[row][col]=0
                        rainfallraster[row][col]=(reflectence[row][col]/300)**(1/1.4)
                where_are_NaNs = numpy.isnan(rainfallraster)
                rainfallraster[where_are_NaNs]=0
                newraster=arcpy.NumPyArrayToRaster(rainfallraster,lowerLeft,cellSize)
                newraster.save('storm'+str(stormfolder)+'/rainfall'+str(stormfolder)+'.tif')
                stormfolder=stormfolder+1
                log ( 'completed rainfall calc' )
                #log ( 'complete with folder',stormfolder )
            log ( 'finished making max reflectance rasters' )

        except Exception as err:
            log(traceback.format_exc())
            log(err)
            raise err
        return

class ScenarioAnalysis(object):
    def __init__(self):
        self.label = "Scenario Analysis"
        self.description = "Compute a quantification of Watershed-wide Improvement based on BMP Buildout Scenario"
        self.canRunInBackground = False
        
        arcpy.env.Workspace = self.Workspace = os.path.split(__file__)[0]
        log("Workspace = " + arcpy.env.Workspace)
        arcpy.env.overwriteOutput = True       

    def getParameterInfo(self):
        """Define parameter definitions"""
        
        param0 = arcpy.Parameter(
            displayName="BMP Points",
            name="bmppts",
            datatype="DEFeatureClass",
            parameterType="Required",
            direction="Input",
            multiValue=False)  
        
        param1 = arcpy.Parameter(
            displayName="Status Field",
            name="statusField",
            datatype="Field",
            parameterType="Required",
            direction="Input",
            multiValue=False)  

        
        params = [ param0, param1 ]
        return params

    def isLicensed(self):
        return True

    def updateParameters(self, parameters):
        return

    def updateMessages(self, parameters):
        return
            
    def execute(self, parameters, messages):
        try:
            log("Parameter is %s" % (parameters[0].valueAsText))


            Lakes = "Lakes"
            LanduseExisting = "LanduseExisting"
            Impervious = "Impervious"
            BMPs = parameters[0].valueAsText
            inFeatures = [Lakes, LanduseExisting, Impervious]
            unionE = "unionE"





            arcpy.Union_analysis (inFeatures,unionE,"ALL", "", "GAPS")

            arcpy.AddField_management("unionE", "TotalNitrogen", "FLOAT")
            arcpy.AddField_management("unionE", "TotalPhosphorus", "FLOAT")
            arcpy.AddField_management("unionE", "Sediment", "FLOAT")
            arcpy.AddField_management("unionE", "Copper", "FLOAT")
            arcpy.AddField_management("unionE", "FecalColiform", "FLOAT")
            arcpy.AddField_management("unionE", "Zinc", "FLOAT")    

            arcpy.SelectLayerByAttribute_management(unionE, "NEW_SELECTION", "Reclassify = 'Roadways'")
            arcpy.CalculateField_management(unionE, "TotalNitrogen", "12.2", "VB", "")
            arcpy.CalculateField_management(unionE, "TotalPhosphorus", "1.8", "VB", "")
            arcpy.CalculateField_management(unionE, "Sediment", "405", "VB", "")
            arcpy.CalculateField_management(unionE, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionE, "Zinc", "0.1", "VB", "")
            arcpy.CalculateField_management(unionE, "FecalColiform", "2.8", "VB", "")


            arcpy.SelectLayerByAttribute_management(unionE, "NEW_SELECTION", "Reclassify in ( 'Commercial' )")
            arcpy.CalculateField_management(unionE, "TotalNitrogen", "14", "VB", "")
            arcpy.CalculateField_management(unionE, "TotalPhosphorus", "2.7", "VB", "")
            arcpy.CalculateField_management(unionE, "Sediment", "400", "VB", "")
            arcpy.CalculateField_management(unionE, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionE, "Zinc", "0", "VB", "")
            arcpy.CalculateField_management(unionE, "FecalColiform", "8.4", "VB", "")


            arcpy.SelectLayerByAttribute_management(unionE, "NEW_SELECTION", "Reclassify IN ( 'Industrial' )")
            arcpy.CalculateField_management(unionE, "TotalNitrogen", "10.4", "VB", "")
            arcpy.CalculateField_management(unionE, "TotalPhosphorus", "1.9", "VB", "")
            arcpy.CalculateField_management(unionE, "Sediment", "372.5", "VB", "")
            arcpy.CalculateField_management(unionE, "Copper", "0.2", "VB", "")
            arcpy.CalculateField_management(unionE, "Zinc", "0.1", "VB", "")
            arcpy.CalculateField_management(unionE, "FecalColiform", "8.4", "VB", "")

            arcpy.SelectLayerByAttribute_management(unionE, "NEW_SELECTION", "Reclassify in ( 'Institutional' , 'Research Triangle Park' )")
            arcpy.CalculateField_management(unionE, "TotalNitrogen", "9.4", "VB", "")
            arcpy.CalculateField_management(unionE, "TotalPhosphorus", "0.7", "VB", "")
            arcpy.CalculateField_management(unionE, "Sediment", "50", "VB", "")
            arcpy.CalculateField_management(unionE, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionE, "Zinc", "0", "VB", "")
            arcpy.CalculateField_management(unionE, "FecalColiform", "14.8", "VB", "")

            arcpy.SelectLayerByAttribute_management(unionE, "NEW_SELECTION", "Reclassify IN ( 'Parks and Open Space' , 'Agricultural' )")
            arcpy.CalculateField_management(unionE, "TotalNitrogen", "2.3", "VB", "")
            arcpy.CalculateField_management(unionE, "TotalPhosphorus", "0.1", "VB", "")
            arcpy.CalculateField_management(unionE, "Sediment", "10", "VB", "")
            arcpy.CalculateField_management(unionE, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionE, "Zinc", "0", "VB", "")
            arcpy.CalculateField_management(unionE, "FecalColiform", "12", "VB", "")

            arcpy.SelectLayerByAttribute_management(unionE, "NEW_SELECTION", "Reclassify IN ( 'High Density Residential' , 'Medium Density Residential' )")
            arcpy.CalculateField_management(unionE, "TotalNitrogen", "11.2", "VB", "")
            arcpy.CalculateField_management(unionE, "TotalPhosphorus", "1.6", "VB", "")
            arcpy.CalculateField_management(unionE, "Sediment", "242.5", "VB", "")
            arcpy.CalculateField_management(unionE, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionE, "Zinc", "0", "VB", "")
            arcpy.CalculateField_management(unionE, "FecalColiform", "30.2", "VB", "")

            arcpy.SelectLayerByAttribute_management(unionE, "NEW_SELECTION", "Reclassify IN( 'Very Low Density Residential' , 'Low Density Residential' )")
            arcpy.CalculateField_management(unionE, "TotalNitrogen", "6.4", "VB", "")
            arcpy.CalculateField_management(unionE, "TotalPhosphorus", "0.7", "VB", "")
            arcpy.CalculateField_management(unionE, "Sediment", "150", "VB", "")
            arcpy.CalculateField_management(unionE, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionE, "Zinc", "0", "VB", "")
            arcpy.CalculateField_management(unionE, "FecalColiform", "16.2", "VB", "")

            arcpy.SelectLayerByAttribute_management(unionE, "CLEAR_SELECTION", "")

            LanduseFuture = "LanduseFuture"

            inFeatures1 = [Lakes, LanduseFuture, Impervious]
            unionF = "unionF"
            arcpy.Union_analysis (inFeatures1,unionF,"ALL", "", "GAPS")
            arcpy.AddField_management("unionF", "TotalNitrogen", "FLOAT")
            arcpy.AddField_management("unionF", "TotalPhosphorus", "FLOAT")
            arcpy.AddField_management("unionF", "Sediment", "FLOAT")
            arcpy.AddField_management("unionF", "Copper", "FLOAT")
            arcpy.AddField_management("unionF", "FecalColiform", "FLOAT")
            arcpy.AddField_management("unionF", "Zinc", "FLOAT")

            arcpy.SelectLayerByAttribute_management(unionF, "NEW_SELECTION", "Reclassify = 'Roadways'")
            arcpy.CalculateField_management(unionF, "TotalNitrogen", "12.2", "VB", "")
            arcpy.CalculateField_management(unionF, "TotalPhosphorus", "1.8", "VB", "")
            arcpy.CalculateField_management(unionF, "Sediment", "405", "VB", "")
            arcpy.CalculateField_management(unionF, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionF, "Zinc", "0.1", "VB", "")
            arcpy.CalculateField_management(unionF, "FecalColiform", "2.8", "VB", "")


            arcpy.SelectLayerByAttribute_management(unionF, "NEW_SELECTION", "Reclassify in ( 'Commercial' )")
            arcpy.CalculateField_management(unionF, "TotalNitrogen", "14", "VB", "")
            arcpy.CalculateField_management(unionF, "TotalPhosphorus", "2.7", "VB", "")
            arcpy.CalculateField_management(unionF, "Sediment", "400", "VB", "")
            arcpy.CalculateField_management(unionF, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionF, "Zinc", "0", "VB", "")
            arcpy.CalculateField_management(unionF, "FecalColiform", "8.4", "VB", "")


            arcpy.SelectLayerByAttribute_management(unionF, "NEW_SELECTION", "Reclassify IN ( 'Industrial' )")
            arcpy.CalculateField_management(unionF, "TotalNitrogen", "10.4", "VB", "")
            arcpy.CalculateField_management(unionF, "TotalPhosphorus", "1.9", "VB", "")
            arcpy.CalculateField_management(unionF, "Sediment", "372.5", "VB", "")
            arcpy.CalculateField_management(unionF, "Copper", "0.2", "VB", "")
            arcpy.CalculateField_management(unionF, "Zinc", "0.1", "VB", "")
            arcpy.CalculateField_management(unionF, "FecalColiform", "8.4", "VB", "")

            arcpy.SelectLayerByAttribute_management(unionF, "NEW_SELECTION", "Reclassify in ( 'Institutional' , 'Research Triangle Park' )")
            arcpy.CalculateField_management(unionF, "TotalNitrogen", "9.4", "VB", "")
            arcpy.CalculateField_management(unionF, "TotalPhosphorus", "0.7", "VB", "")
            arcpy.CalculateField_management(unionF, "Sediment", "50", "VB", "")
            arcpy.CalculateField_management(unionF, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionF, "Zinc", "0", "VB", "")
            arcpy.CalculateField_management(unionF, "FecalColiform", "14.8", "VB", "")

            arcpy.SelectLayerByAttribute_management(unionF, "NEW_SELECTION", "Reclassify IN ( 'Parks and Open Space' , 'Agricultural' )")
            arcpy.CalculateField_management(unionF, "TotalNitrogen", "2.3", "VB", "")
            arcpy.CalculateField_management(unionF, "TotalPhosphorus", "0.1", "VB", "")
            arcpy.CalculateField_management(unionF, "Sediment", "10", "VB", "")
            arcpy.CalculateField_management(unionF, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionF, "Zinc", "0", "VB", "")
            arcpy.CalculateField_management(unionF, "FecalColiform", "12", "VB", "")

            arcpy.SelectLayerByAttribute_management(unionF, "NEW_SELECTION", "Reclassify IN ( 'High Density Residential' , 'Medium Density Residential' )")
            arcpy.CalculateField_management(unionF, "TotalNitrogen", "11.2", "VB", "")
            arcpy.CalculateField_management(unionF, "TotalPhosphorus", "1.6", "VB", "")
            arcpy.CalculateField_management(unionF, "Sediment", "242.5", "VB", "")
            arcpy.CalculateField_management(unionF, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionF, "Zinc", "0", "VB", "")
            arcpy.CalculateField_management(unionF, "FecalColiform", "30.2", "VB", "")

            arcpy.SelectLayerByAttribute_management(unionF, "NEW_SELECTION", "Reclassify IN( 'Very Low Density Residential' , 'Low Density Residential' )")
            arcpy.CalculateField_management(unionF, "TotalNitrogen", "6.4", "VB", "")
            arcpy.CalculateField_management(unionF, "TotalPhosphorus", "0.7", "VB", "")
            arcpy.CalculateField_management(unionF, "Sediment", "150", "VB", "")
            arcpy.CalculateField_management(unionF, "Copper", "0", "VB", "")
            arcpy.CalculateField_management(unionF, "Zinc", "0", "VB", "")
            arcpy.CalculateField_management(unionF, "FecalColiform", "16.2", "VB", "")

            arcpy.SelectLayerByAttribute_management(unionF, "CLEAR_SELECTION", "")
           
            Raster_Nitrogen = "Raster_Nitrogen"
            Raster_Phosphorus = "Raster_Phosphorus"
            Raster_Sediment = "Raster_Sediment"
            Raster_Copper = "Raster_Copper"
            Raster_Zinc = "Raster_Zinc"
            Raster_Fecal = "Raster_Fecal"
            arcpy.PolygonToRaster_conversion(unionE, "TotalNitrogen", Raster_Nitrogen, "CELL_CENTER", "NONE", "40")
            arcpy.PolygonToRaster_conversion(unionE, "TotalPhosphorus", Raster_Phosphorus, "CELL_CENTER", "NONE", "40")
            arcpy.PolygonToRaster_conversion(unionE, "Sediment", Raster_Sediment, "CELL_CENTER", "NONE", "40")
            arcpy.PolygonToRaster_conversion(unionE, "Copper", Raster_Copper, "CELL_CENTER", "NONE", "40")
            arcpy.PolygonToRaster_conversion(unionE, "Zinc", Raster_Zinc, "CELL_CENTER", "NONE", "40")
            arcpy.PolygonToRaster_conversion(unionE, "FecalColiform", Raster_Fecal, "CELL_CENTER", "NONE", "40")

            RasterF_Nitrogen = "RasterF_Nitrogen"
            RasterF_Phosphorus = "RasterF_Phosphorus"
            RasterF_Sediment = "RasterF_Sediment"
            RasterF_Copper = "RasterF_Copper"
            RasterF_Zinc = "RasterF_Zinc"
            RasterF_Fecal = "RasterF_Fecal"

            arcpy.PolygonToRaster_conversion(unionF, "TotalNitrogen", RasterF_Nitrogen, "CELL_CENTER", "NONE", "40")
            arcpy.PolygonToRaster_conversion(unionF, "TotalPhosphorus", RasterF_Phosphorus, "CELL_CENTER", "NONE", "40")
            arcpy.PolygonToRaster_conversion(unionF, "Sediment", RasterF_Sediment, "CELL_CENTER", "NONE", "40")
            arcpy.PolygonToRaster_conversion(unionF, "Copper", RasterF_Copper, "CELL_CENTER", "NONE", "40")
            arcpy.PolygonToRaster_conversion(unionF, "Zinc", RasterF_Zinc, "CELL_CENTER", "NONE", "40")
            arcpy.PolygonToRaster_conversion(unionF, "FecalColiform", RasterF_Fecal, "CELL_CENTER", "NONE", "40")

            BMP_TN = "BMP_TN"
            BMP_TP = "BMP_TP"
            BMP_FC = "BMP_FC"
            BMP_CU = "BMP_CU"
            BMP_Zn = "BMP_Zn"
            BMP_Sed = "BMP_Sed"


            arcpy.PointToRaster_conversion(BMPs, "TN_Eff_Ex", BMP_TN, "MOST_FREQUENT", "NONE", "40")
            arcpy.PointToRaster_conversion(BMPs, "TP_Eff_Ex", BMP_TP, "MOST_FREQUENT", "NONE", "40")
            arcpy.PointToRaster_conversion(BMPs, "FC_Eff_Ex", BMP_FC, "MOST_FREQUENT", "NONE", "40")
            arcpy.PointToRaster_conversion(BMPs, "CU_Eff_Ex", BMP_CU, "MOST_FREQUENT", "NONE", "40")
            arcpy.PointToRaster_conversion(BMPs, "Zn_Eff_Ex", BMP_Zn, "MOST_FREQUENT", "NONE", "40")
            arcpy.PointToRaster_conversion(BMPs, "Sed_Eff_Ex", BMP_Sed, "MOST_FREQUENT", "NONE", "40")



            demPath = "DEM"
            arcpy.env.extent = demPath
            arcpy.env.cellSize = demPath
            dem = Raster(demPath)
            fill = Fill(dem)
            flowDirection = FlowDirection(fill)
            flowAccEN = FlowAccumulation(flowDirection,Raster_Nitrogen, "FLOAT")
            flowAccEP = FlowAccumulation(flowDirection,Raster_Phosphorus, "FLOAT")
            flowAccEN.save("C:/ss/output/flowAccEN1.tif")
            flowAccEP.save("C:/ss/output/flowAccEP1.tif")
            flowAccES = FlowAccumulation(flowDirection,Raster_Sediment, "FLOAT")
            flowAccEC = FlowAccumulation(flowDirection,Raster_Copper, "FLOAT")
            flowAccES.save("C:/ss/output/flowAccES1.tif")
            flowAccEC.save("C:/ss/output/flowAccEC1.tif")
            flowAccEZ = FlowAccumulation(flowDirection,Raster_Zinc, "FLOAT")
            flowAccEF = FlowAccumulation(flowDirection,Raster_Fecal, "FLOAT")
            flowAccEZ.save("C:/ss/output/flowAccEZ1.tif")
            flowAccEF.save("C:/ss/output/flowAccEF1.tif")

            FlowAcc3_Nitrogen = FlowAccumulation(flowDirection,RasterF_Nitrogen, "FLOAT")
            FlowAcc3_Phosphorus = FlowAccumulation(flowDirection,RasterF_Phosphorus, "FLOAT")
            FlowAcc3_Sediment = FlowAccumulation(flowDirection,RasterF_Sediment, "FLOAT")
            FlowAcc3_Copper = FlowAccumulation(flowDirection,RasterF_Copper, "FLOAT")
            FlowAcc3_Zinc = FlowAccumulation(flowDirection,RasterF_Zinc, "FLOAT")
            FlowAcc3_Fecal = FlowAccumulation(flowDirection,RasterF_Fecal, "FLOAT")

            FlowAcc3_Nitrogen_BMP = FlowAccumulation(flowDirection, BMP_TN, "FLOAT")
            FlowAcc3_Phosphorus_BMP = FlowAccumulation(flowDirection, BMP_TP, "FLOAT")
            FlowAcc3_Fecal_BMP = FlowAccumulation(flowDirection, BMP_FC, "FLOAT")
            FlowAcc3_Sediment_BMP = FlowAccumulation(flowDirection, BMP_Sed, "FLOAT")
            FlowAcc3_Copper_BMP = FlowAccumulation(flowDirection, BMP_CU, "FLOAT")
            FlowAcc3_Zinc_BMP = FlowAccumulation(flowDirection, BMP_Zn, "FLOAT")

            maxN = Minus(692.2000122070313,FlowAcc3_Nitrogen_BMP)
            maxP = Minus(1082.5,FlowAcc3_Phosphorus_BMP)
            maxS = Minus(2205.800048828125,FlowAcc3_Sediment_BMP)
            maxF = Minus(1503,FlowAcc3_Fecal_BMP)
            maxC = Minus(1074,FlowAcc3_Copper_BMP)
            maxZ = Minus(1382,FlowAcc3_Zinc_BMP)

            BMPred_N = Divide(maxN, 692.2000122070313)
            BMPred_P = Divide(maxP, 1082.5)
            BMPred_se = Divide(maxS, 2205.800048828125)
            BMPred_fe = Divide(maxF, 1503)
            BMPred_Copper = Divide(maxC, 1074)
            BMPred_Z = Divide(maxZ,1382)

            futureSe = Times(FlowAcc3_Sediment,BMPred_se)
            futureFe = Times(FlowAcc3_Fecal, BMPred_fe)
            futureZ = Times(FlowAcc3_Zinc, BMPred_Z)
            futureC = Times(FlowAcc3_Copper, BMPred_Copper)
            futureP = Times(BMPred_P, FlowAcc3_Phosphorus)
            futureN = Times(BMPred_N, FlowAcc3_Nitrogen)
            futureSe.save("C:/ss/output/futureSe1.tif")
            futureFe.save("C:/ss/output/futureFe1.tif")
            futureZ.save("C:/ss/output/futureZ1.tif")
            futureC.save("C:/ss/output/futureC1.tif")
            futureP.save("C:/ss/output/futureP1.tif")
            futureN.save("C:/ss/output/futureN1.tif")
            log("final")
            Nfinal = Con(IsNull(flowAccEN),0,Divide(futureN, flowAccEN))
            Pfinal = Con(IsNull(flowAccEP),0,Divide(futureP, flowAccEP))
            Sfinal = Con(IsNull(flowAccES),0,Divide(futureSe, flowAccES))
            Zfinal = Con(IsNull(flowAccEZ),0,Divide(futureZ, flowAccEZ))
            Cfinal = Con(IsNull(flowAccEC),0,Divide(futureC, flowAccEC))
            Ffinal = Con(IsNull(flowAccEF),0,Divide(futureFe, flowAccEF))
            Nfinal.save("C:/ss/output/Nfinal.tif")
            Pfinal.save("C:/ss/output/Pfinal.tif")
            Sfinal.save("C:/ss/output/Sfinal.tif")
            Zfinal.save("C:/ss/output/Zfinal.tif")
            Cfinal.save("C:/ss/output/Cfinal.tif")
            Ffinal.save("C:/ss/output/Ffinal.tif")
            log("reclassify")
            maxN = 1146.41
            maxP = 3889
            maxS = 8749
            maxF = 299.333
            maxZ = 3685
            maxC = 494
            Reclass_N = Reclassify(Nfinal, "Value", "0 %f NODATA;%f %f 1;%f %f 2;%f %f 3;%f %f 4" %(0.05*maxN, 0.05*maxN, 0.1 * maxN, 0.1 * maxN, 0.25 * maxN, 0.25 * maxN, 0.5 * maxN, 0.5 * maxN, maxN), "DATA")
            Reclass_P = Reclassify(Pfinal, "Value", "0 %f NODATA;%f %f 1;%f %f 2;%f %f 3;%f %f 4" %(0.05*maxP, 0.05*maxP, 0.1 * maxP, 0.1 * maxP, 0.25 * maxP, 0.25 * maxP, 0.5 * maxP, 0.5 * maxP, maxP), "DATA")
            Reclass_S = Reclassify(Sfinal, "Value", "0 %f NODATA;%f %f 1;%f %f 2;%f %f 3;%f %f 4" %(0.05*maxS, 0.05*maxS, 0.1 * maxS, 0.1 * maxS, 0.25 * maxS, 0.25 * maxS, 0.5 * maxS, 0.5 * maxS, maxS), "DATA")
            Reclass_F = Reclassify(Ffinal, "Value", "0 %f NODATA;%f %f 1;%f %f 2;%f %f 3;%f %f 4" %(0.05*maxF, 0.05*maxF, 0.1 * maxF, 0.1 * maxF, 0.25 * maxF, 0.25 * maxF, 0.5 * maxF, 0.5 * maxF, maxF), "DATA")
            Reclass_Z = Reclassify(Zfinal, "Value", "0 %f NODATA;%f %f 1;%f %f 2;%f %f 3;%f %f 4" %(0.05*maxZ, 0.05*maxZ, 0.1 * maxZ, 0.1 * maxZ, 0.25 * maxZ, 0.25 * maxZ, 0.5 * maxZ, 0.5 * maxZ, maxZ), "DATA")
            Reclass_C = Reclassify(Cfinal, "Value", "0 %f NODATA;%f %f 1;%f %f 2;%f %f 3;%f %f 4" %(0.05*maxC, 0.05*maxC, 0.1 * maxC, 0.1 * maxC, 0.25 * maxC, 0.25 * maxC, 0.5 * maxC, 0.5 * maxC, maxC), "DATA")
            StreamT_N = "StreamT_N"
            StreamT_P = "StreamT_P"
            StreamT_S = "StreamT_S"
            StreamT_F = "StreamT_F"
            StreamT_Z = "StreamT_Z"
            StreamT_C = "StreamT_C"
            log("stream to feature")
            arcpy.gp.StreamToFeature_sa(Reclass_N, flowDirection, StreamT_N, "SIMPLIFY")
            arcpy.gp.StreamToFeature_sa(Reclass_P, flowDirection, StreamT_P, "SIMPLIFY")
            arcpy.gp.StreamToFeature_sa(Reclass_S, flowDirection, StreamT_S, "SIMPLIFY")
            arcpy.gp.StreamToFeature_sa(Reclass_F, flowDirection, StreamT_F, "SIMPLIFY")
            arcpy.gp.StreamToFeature_sa(Reclass_Z, flowDirection, StreamT_Z, "SIMPLIFY")
            arcpy.gp.StreamToFeature_sa(Reclass_C, flowDirection, StreamT_C, "SIMPLIFY")

            StreamInvPts = "StreamInvPts"

            inRasterList = " 'C:/ss/output/flowAccEN1.tif' flowAccEN; 'C:/ss/output/flowAccEP1.tif' flowAccEP;'C:/ss/output/flowAccES1.tif' flowAccES;'C:/ss/output/flowAccEC1.tif' flowAccEC;'C:/ss/output/flowAccEZ1.tif' flowAccEZ;'C:/ss/output/flowAccEF1.tif' flowAccEF;'C:/ss/output/futureSe1.tif' futureSe;'C:/ss/output/futureFe1.tif' futureFe;'C:/ss/output/futureN1.tif' futureN;'C:/ss/output/futureC1.tif' futureC;'C:/ss/output/futureP1.tif' futureP;'C:/ss/output/futureZ1.tif' futureZ;"


            arcpy.gp.ExtractMultiValuesToPoints_sa(StreamInvPts, inRasterList, "NONE")

            arcpy.AddField_management("StreamInvPts", "N", "FLOAT")
            arcpy.AddField_management("StreamInvPts", "P", "FLOAT")
            arcpy.AddField_management("StreamInvPts", "S", "FLOAT")
            arcpy.AddField_management("StreamInvPts", "C", "FLOAT")
            arcpy.AddField_management("StreamInvPts", "F", "FLOAT")
            arcpy.AddField_management("StreamInvPts", "Z", "FLOAT")
            arcpy.SelectLayerByAttribute_management(StreamInvPts, "NEW_SELECTION", "flowAccEN > 0")
            arcpy.CalculateField_management(StreamInvPts, "N", "[futureN] / [flowAccEN]", "VB", "")

            arcpy.SelectLayerByAttribute_management(StreamInvPts, "NEW_SELECTION", "flowAccEP > 0")
            arcpy.CalculateField_management(StreamInvPts, "P", "[futureP] / [flowAccEP]", "VB", "")
            arcpy.SelectLayerByAttribute_management(StreamInvPts, "NEW_SELECTION", "flowAccEC > 0")
            arcpy.CalculateField_management(StreamInvPts, "C", "[futureC] / [flowAccEC]", "VB", "")
            arcpy.SelectLayerByAttribute_management(StreamInvPts, "NEW_SELECTION", "flowAccEZ > 0")
            arcpy.CalculateField_management(StreamInvPts, "Z", "[futureZ] / [flowAccEZ]", "VB", "")

            arcpy.SelectLayerByAttribute_management(StreamInvPts, "NEW_SELECTION", "flowAccEF = 0")
            arcpy.SelectLayerByAttribute_management(StreamInvPts, "SWITCH_SELECTION", "")
            arcpy.CalculateField_management(StreamInvPts, "F", "[futureFe] / [flowAccEF]", "VB", "")
            arcpy.SelectLayerByAttribute_management(StreamInvPts, "NEW_SELECTION", "flowAccES > 0")
            arcpy.CalculateField_management(StreamInvPts, "S", "[futureSe] / [flowAccES]", "VB", "")

            arcpy.SelectLayerByAttribute_management(StreamInvPts, "CLEAR_SELECTION", "")



        except Exception as err:
            log(traceback.format_exc())
            log(err)
            raise err
        return



