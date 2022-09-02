/* import (5 entries)
var L8sr:ImageCollection"USGS Landsat 8 Surface Reflectance Tier 1"
var L7sr:ImageCollection"USGS Landsat 7 Surface Reflectance Tier 1"
var L5sr:ImageCollection"USGS Landsat 5 Surface Reflectance Tier 1"
var studyarea:Table users/zqqcnu/HC_34Y_studyarea/S1314
var hhkBoundry:Table users/zqqcnu/huanghekou/huanghekou_boundry
*/

var hhkBoundry = ee.FeatureCollection("users/zqqcnu/huanghekou/huanghekou_boundry"),
    studyarea = ee.FeatureCollection("users/zqqcnu/HC_34Y_studyarea/S1314")

var timeField = 'system:time_start';
//Screen image
//----------------------------------------Seagrass green period-----------------------------------------------------
var start = ee.Date('2013-6-1');
var end = ee.Date('2013-9-10');

var start1 = ee.Date('2014-6-1');
var end1 = ee.Date('2014-9-10');

//----------------------------------------- seagrass senescence period---------------------------------------------------------
var start3 = ee.Date('2013-9-21');
var end3 = ee.Date('2013-12-1');

var start4 = ee.Date('2014-9-21');
var end4 = ee.Date('2014-12-1');


//-----------------------------surface reflectance QA band to mask clouds.

var cloudMaskL457 = function(image) {
  var qa = image.select('pixel_qa');
  // If the cloud bit (5) is set and the cloud confidence (7) is high
  // or the cloud shadow bit is set (3), then it's a bad pixel.
  var cloud = qa.bitwiseAnd(1 << 5)
                .and(qa.bitwiseAnd(1 << 7))
                .or(qa.bitwiseAnd(1 << 4))
                .or(qa.bitwiseAnd(1 << 3));
  // Remove edge pixels that don't occur in all bands
  var mask = image.mask().reduce(ee.Reducer.min());
  
  return image.updateMask(cloud.not()).updateMask(mask).divide(10000)
      .select('B1','B2', 'B3', 'B4', 'B5', 'B7')
      .copyProperties(image, ['system:time_start']);
};


var inBands8 = ee.List(['B2','B3','B4','B5','B6','B7'])
var outBands8 = ee.List(['B1','B2','B3','B4','B5','B7']); 
function maskL8sr(image) {
  // Bits 3 and 5 are cloud shadow and cloud, respectively.
  var cloudShadowBitMask = 1 << 3;
  var cloudsBitMask = 1 << 5;
  var snowBitMask = 1 << 4;

  // Get the pixel QA band.
  var qa = image.select('pixel_qa');

  // Both flags should be set to zero, indicating clear conditions.
  var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
      .and(qa.bitwiseAnd(cloudsBitMask).eq(0))
      .and(qa.bitwiseAnd(snowBitMask).eq(0));

  // Return the masked image, scaled to reflectance, without the QA bands.
  return image.updateMask(mask).divide(10000)
      .select(inBands8,outBands8)
      .copyProperties(image, ["system:time_start"]);
}



//-----------------------------L578 Band combinations and spectral indices
function index(image) {
  var date = ee.Date(image.get(timeField));
  var years = date.difference(ee.Date('1970-01-01'), 'year');
  var nir = image.select("B4");
  var red = image.select("B3");
  var green = image.select("B2");
  var blue = image.select("B1");
  var swir = image.select("B5");

    return image
    .addBands(image.normalizedDifference(['B4', 'B3']).rename('NDVI'))
    // Add an LSWI band.
    .addBands(image.normalizedDifference(['B2', 'B5']).rename('mNDWI'))
    // Add a time band.
    //.addBands(ee.Image(years).rename('t'))
    .float()
    // Add a constant band.
    //.addBands(ee.Image.constant(1));
}

//--------------------green period date1--------------------

var filteredLandsat1 = L8sr
  .filterBounds(studyarea)
  .filterDate(start,end)
  .map(maskL8sr)
  .map(index)
  

var filteredLandsat2 = L5sr
  .filterBounds(studyarea)
  .filterDate(start,end)
  .map(cloudMaskL457)
  .map(index)
  
  
var filteredLandsat3 = L7sr
  .filterBounds(studyarea)
  .filterDate(start,end)
  .map(cloudMaskL457)
  .map(index)
  

//-------------------------green period date 2----------------------
var filteredLandsat11 = L8sr
  .filterBounds(studyarea)
  .filterDate(start1,end1)
  .map(maskL8sr)
  .map(index)
  
print(filteredLandsat1)
var filteredLandsat21 = L5sr
  .filterBounds(studyarea)
  .filterDate(start1,end1)
  .map(cloudMaskL457)
  .map(index)
  
var filteredLandsat31 =  L7sr
  .filterBounds(studyarea)
  .filterDate(start1,end1)
  .map(cloudMaskL457)
  .map(index)
  
  
var collection =filteredLandsat1.merge(filteredLandsat2).merge(filteredLandsat3)
                .merge(filteredLandsat11).merge(filteredLandsat21).merge(filteredLandsat31);
//print(collection,"green period collection")

//-----------------------------senescence date1---------------------------

var filteredLandsatKW1 = L8sr
  .filterBounds(studyarea)
  .filterDate(start3,end3)
  .map(maskL8sr)
  .map(index)

print(filteredLandsatKW1)
var filteredLandsatKW2 = L5sr
  .filterBounds(studyarea)
  .filterDate(start3,end3)
  .map(cloudMaskL457)
  .map(index)
  
var filteredLandsatKW3 =  L7sr
  .filterBounds(studyarea)
  .filterDate(start3,end3)
  .map(cloudMaskL457)
  .map(index)
  

//-------------------------senescence date 2----------------------
var filteredLandsatKW11 = L8sr
  .filterBounds(studyarea)
  .filterDate(start4,end4)
  .map(maskL8sr)
  .map(index)
  
var filteredLandsatKW21 = L5sr
  .filterBounds(studyarea)
  .filterDate(start4,end4)
  .map(cloudMaskL457)
  .map(index)
  
  
var filteredLandsatKW31 =  L7sr
  .filterBounds(studyarea)
  .filterDate(start4,end4)
  .map(cloudMaskL457)
  .map(index)
  
  
var collectionKW =filteredLandsatKW1.merge(filteredLandsatKW2).merge(filteredLandsatKW3)
                .merge(filteredLandsatKW11).merge(filteredLandsatKW21).merge(filteredLandsatKW31);
//print(collectionKW,"senescence collection")


var visParams = {
  bands: ['B4', 'B3', 'B2'],
  min: 0,
  max: 0.3,
  gamma: 1.4,
};
Map.centerObject(studyarea,11)

//-----------------------------------Tasseled Cap Transformation-----------------------------

var L5_collection =filteredLandsat2.merge(filteredLandsat21);
//print(L5_collection,"green_L5")

var L7_collection =filteredLandsat3.merge(filteredLandsat31);
//print(L7_collection,"green_L7")  

var L8_collection =filteredLandsat1.merge(filteredLandsat11);
//print(L8_collection,"green_L8")


//Maximum NDVI synthesis during the green perid
var greenestL5 = L5_collection.qualityMosaic('NDVI').clip(studyarea);

var greenestL7 = L7_collection.qualityMosaic('NDVI').clip(studyarea);

var greenestL8 = L8_collection.qualityMosaic('NDVI').clip(studyarea);

 
 //landsat 8----------------------------------------------------------------- 
var landsat8 = greenestL8.select('B1','B2','B3','B4','B5','B7');
//print(landsat8);

  // Define an Array of Tasseled Cap coefficients
  var coefficientsL8 = ee.Array([
    [0.3029, 0.2786, 0.4733, 0.5599, 0.508, 0.1872],
    [-0.2941, -0.243, -0.5424, 0.7276, 0.0713, -0.1608],
    [0.1511, 0.1973, 0.3283, 0.3407, -0.7117, -0.4559],
    [-0.8239, 0.0849, 0.4396, -0.0580, 0.2013, -0.2773],
    [-0.3294, 0.0557, 0.1056, 0.1855, -0.4349, 0.8085],
    [0.1079, -0.9023, 0.4119, 0.0575, -0.0259, 0.0252]
  ]);

  var arrayImage1D = landsat8.toArray();
  // make an Array Image with a 2-D Array per pixel, 6x1
  var arrayImage2D = arrayImage1D.toArray(1);

  // do a matrix multiplication: 6x6 times 6x1
  var image8 = ee.Image(coefficientsL8)
    .matrixMultiply(arrayImage2D)
    .arrayProject([0]) // get rid of the extra dimensions
    .arrayFlatten([
      ['brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth']
  ]).select(['brightness', 'greenness', 'wetness']);


//landsat 7 -------------------------------------------------------------------------
var landsat7 = greenestL7.select('B1','B2','B3','B4','B5','B7');
  print(landsat7);
  // Define an Array of Tasseled Cap coefficients
  var coefficientsL7 = ee.Array([
    [ 0.3561, 0.3972, 0.3904, 0.6966, 0.2286, 0.1596],
    [-0.3344,-0.3544,-0.4556, 0.6966,-0.0242,-0.2630],
    [0.2626, 0.2141, 0.0926, 0.0656,-0.7629,-0.5388],
    [0.0805,-0.0498, 0.1950,-0.1327, 0.5752,-0.7775],
    [-0.7252,-0.0202, 0.6683, 0.0631,-0.1494,-0.0274],
    [0.4000,-0.8172, 0.3832, 0.0602,-0.1095, 0.0985]
  ]);

  var arrayImage1D = landsat7.toArray();
  // make an Array Image with a 2-D Array per pixel, 6x1
  var arrayImage2D = arrayImage1D.toArray(1);

  // do a matrix multiplication: 6x6 times 6x1
  var image7 = ee.Image(coefficientsL7)
    .matrixMultiply(arrayImage2D)
    .arrayProject([0]) // get rid of the extra dimensions
    .arrayFlatten([
      ['brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth']
  ]).select(['brightness', 'greenness', 'wetness']);


//landsat 5------------------------------------------------------------------------------------------
var landsat5 = greenestL5.select('B1','B2','B3','B4','B5','B7');
  print(landsat5);
  // Define an Array of Tasseled Cap coefficients
  var coefficientsL5 = ee.Array([
    [0.3037, 0.2793, 0.4743, 0.5585, 0.5082, 0.1863],
    [-0.2848, -0.2435, -0.5436, 0.7243, 0.0840, -0.1800],
    [0.1509, 0.1973, 0.3279, 0.3406, -0.7112, -0.4572],
    [-0.8242, 0.0849, 0.4392, -0.0580, 0.2012, -0.2768],
    [-0.3280, 0.0549, 0.1075, 0.1855, -0.4357, 0.8085],
    [0.1084, -0.9022, 0.4120, 0.0573, -0.0251, 0.0238]
  ]);

  var arrayImage1D = landsat5.toArray();
  // make an Array Image with a 2-D Array per pixel, 6x1
  var arrayImage2D = arrayImage1D.toArray(1);

  // do a matrix multiplication: 6x6 times 6x1
  var image5 = ee.Image(coefficientsL5)
    .matrixMultiply(arrayImage2D)
    .arrayProject([0]) // get rid of the extra dimensions
    .arrayFlatten([
      ['brightness', 'greenness', 'wetness', 'fourth', 'fifth', 'sixth']
  ]).select(['brightness', 'greenness', 'wetness']);


// //---------------------select brightness--------------------------------------------
   
  var image8_b= image8.select("brightness")
  var image7_b= image7.select("brightness")
// var image5_b= image5.select("brightness");
  var img_brightness = image7_b.min(image8_b).clip(studyarea);
  //print(img_min_brightness);
//Map.addLayer(img_brightness,{'palette':['black','white'],'min':0.03,'max':0.4},'img_brightness')


//---------------------------------------Maximum NDVI synthesis-------------------------


var greenest = collection.qualityMosaic('NDVI').clip(hhkBoundry);
 Map.addLayer(greenest,visParams,'NDVI_max');
 
var NDVImax = greenest.select('NDVI');


 
//------------------------Mean NDVI in senescence perid--------------------
var greenest_KWmean = collectionKW.select('NDVI').mean().clip(studyarea);
//Map.addLayer(greenest_KWmean,visParam,'Se_NDVI_mean');


//-----------------------green NDVImax-senescence NDVI_mean-----------------
var NDVImax_KWmean = greenest.select('NDVI').subtract(greenest_KWmean);
//Map.addLayer(NDVImax_KWmean,{'palette':['white','black'],'min':0,'max':1.3},'生长季NDVImax-枯萎季NDVImean')

//-----------------------When NDVImax, mNDWI at the same time----------------
var greenest_mNDWI = greenest.select('mNDWI').clip(studyarea);
var mNDWI_0 = greenest_mNDWI.gte(0)
Map.addLayer(mNDWI_0,{},'mNDWI_0')

var palettes = require('users/gena/packages:palettes');
var palette = palettes.misc.tol_rainbow[7];

//---------------------------------SGI seagrass index--------------------------
var SGI = (NDVImax.add(NDVImax_KWmean)).divide(img_brightness);
Map.addLayer(SGI, {min: -1, max: 12, palette: palette}, 'SGI');



Export.image.toDrive({
  image: SGI,
  description: "1314SGI",
  fileNamePrefix: "1314SGI",
  scale: 30,
  region: studyarea,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: mNDWI_0,
  description: "1314mNDWI",
  fileNamePrefix: "1314mNDWI",
  scale: 30,
  region: studyarea,
  maxPixels: 1e13
});