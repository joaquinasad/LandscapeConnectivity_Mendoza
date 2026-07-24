//---------------------------------------------------------------------
// 1-INTEGRACIÓN Y RECORTE DE CLASES ORIGINALES DE MAPBIOMAS (1985-2024)
//---------------------------------------------------------------------

// Área de estudio
var roi = ee.FeatureCollection("projects/earthengine-legacy/assets/users/joaquinasad/areadeestudio3");

// Importación de activos de MapBiomas Cuyo y Urbano
var mapbiomasGeneral = ee.Image("projects/mapbiomas-argentina/assets/LAND-COVER/COLLECTION-2/GENERAL/CLASSIFICATION/FINAL_CLASSIFICATION/CUYO/CUYO-FINAL-v1");
var mapbiomasUrbano = ee.Image("projects/mapbiomas-argentina/assets/LAND-COVER/COLLECTION-2/GENERAL/URBAN/CLASSIFICATION/FINAL_CLASSIFICATION/urban_argentina_1985-2024_v1");

// Módulo de paletas oficiales de MapBiomas Argentina
var mapbiomasPalette = require('users/mapbiomas/modules:Palettes.js').get('classification9')
  .concat(["#000000","#000000","#000000","#000000","#000000","#000000","#000000","#a2c830"]);

// Parámetros de visualización utilizando la paleta completa original
var visClass = {
    'min': 0, 
    'max': 77,
    'palette': mapbiomasPalette
};

// Lista de años del periodo de estudio
var years = ee.List.sequence(1985, 2024);

// Integrarción de la capa urbana
var procesarAñoOriginal = function(year) {
  var yearStr = ee.Number(year).format('%04d');
  var bandName = ee.String('classification_').cat(yearStr);
  
  // Seleccionar la banda del año correspondiente para ambas imágenes
  var imgGeneral = mapbiomasGeneral.select([bandName], ['classification']);
  var imgUrbano = mapbiomasUrbano.select([bandName], ['classification']);
  
  // Superposición del área urbana (clase 24) sobre la clasificación general
  var imgIntegrada = imgGeneral.where(imgUrbano.eq(24), 24);
  
  // Recortarte al límite del polígono de estudio
  var imgRecortada = imgIntegrada.clip(roi);
  
  // Retornar la imagen con metadatos temporales asignados
  return imgRecortada
    .set('year', year)
    .set('system:time_start', ee.Date.fromYMD(year, 1, 1).millis());
};

// Creación de la ImageCollection con todas las clases originales mapping de la función
var coleccionOriginales = ee.ImageCollection(years.map(procesarAñoOriginal));

print('Colección MapBiomas Integrada con clases originales (1985-2024):', coleccionOriginales);

// --------------------------------------------------------------------
// 2 Y 3 - MSPA AGRUPADO Y ANÁLISIS ESPACIAL DE SUPERFICIES (1985-2024)
// --------------------------------------------------------------------

//Selección de clases naturales (MapBiomas)
var clasesNaturales = [3, 4, 66, 77, 45, 12, 11, 33, 34, 25]; 

// Distancia de borde (2 píxeles = 60 metros)
var edgeWidthPixels = 2; 

// Topología MSPA agrupada
var aplicarMSPA_Agrupado = function(img) {
  var toList = ee.List.repeat(1, clasesNaturales.length);
  var binary = img.remap(clasesNaturales, toList, 0).rename('binary');
  
  var core = binary.focal_min({radius: edgeWidthPixels, kernelType: 'square', units: 'pixels'});
  var coreDilatado = core.focal_max({radius: edgeWidthPixels, kernelType: 'square', units: 'pixels'});
  var edge = coreDilatado.and(binary).and(core.not());
  
  var mspaClass = ee.Image(0)
    .where(binary.eq(1), 3) // Conectores e Islotes
    .where(edge.eq(1), 2)   // Borde / Transición
    .where(core.eq(1), 1)   // Core / Núcleo
    .where(binary.eq(0), 0) // Matriz de uso del suelo (no natural)
    .rename('mspa_class')
    .set('year', img.get('year'))
    .set('system:time_start', img.get('system:time_start'));
    
  return mspaClass.toByte();
};

var coleccionMSPA_Completa = coleccionOriginales.map(aplicarMSPA_Agrupado);

// Función para calcular el área (en hectáreas) de cada clase MSPA
var calcularAreas = function(img) {
  var areaImage = ee.Image.pixelArea().divide(10000); 
  var coreArea = areaImage.updateMask(img.eq(1)).rename('Core');
  var edgeArea = areaImage.updateMask(img.eq(2)).rename('Edge');
  var connArea = areaImage.updateMask(img.eq(3)).rename('Conectores');
  var imgAreas = ee.Image([coreArea, edgeArea, connArea]);
  
  var stats = imgAreas.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: roi.geometry(),
    scale: 30, 
    maxPixels: 1e13
  });
    
  return ee.Feature(null, {
    'year': ee.Number.parse(img.get('year')),
    'Core_ha': stats.get('Core'),
    'Edge_ha': stats.get('Edge'),
    'Conectores_ha': stats.get('Conectores')
  });
};

var serieTemporalAreas = ee.FeatureCollection(coleccionMSPA_Completa.map(calcularAreas));

// Generación del gráfico de la serie temporal
var chartAreas = ui.Chart.feature.byFeature({
  features: serieTemporalAreas,
  xProperty: 'year',
  yProperties: ['Core_ha', 'Edge_ha', 'Conectores_ha']
})
.setOptions({
  title: 'Dinámica de fragmentación estructural del Área de Estudio (1985-2024)',
  vAxis: {title: 'Área (Hectáreas)'},
  hAxis: {title: 'Año', format: '####'},
  lineWidth: 2,
  colors: ['#00441b', '#a1d99b', '#ff7f00'] 
});
print('Análisis de Áreas por Año:', chartAreas);

// --------------------------------------------------------
// 4 - EXTRACCIÓN DE NODOS DE HÁBITAT NATURAL Y COORDENADAS
// --------------------------------------------------------

var extraerNodosAnuales = function(img) {
  var year = ee.Number.parse(img.get('year'));
  
  var soloCore = img.select('mspa_class').eq(1);
  var coreMasked = soloCore.updateMask(soloCore);
  
  // Coordenadas de latitud y longitud 
  var imgConCoords = coreMasked.addBands(ee.Image.pixelLonLat());
  
  // Vectorizar los parches
  var parchesVector = coreMasked.reduceToVectors({
    geometry: roi.geometry(),
    crs: img.projection(),
    scale: 30, 
    geometryType: 'polygon',
    eightConnected: true, 
    labelProperty: 'mspa_class',
    maxPixels: 1e13
  });
  
  // Extraer la coordenada del primer píxel del interior a cada polígono
  var parchesConPunto = imgConCoords.reduceRegions({
    collection: parchesVector,
    reducer: ee.Reducer.first(),
    scale: 30,
    crs: img.projection()
  });
  
  var nodosConAtributos = parchesConPunto.map(function(feature) {
    var geom = feature.geometry();
    var areaHa = geom.area(1).divide(10000);
    
    // Tomamar la longitud y latitud del píxel interior
    var lon = feature.get('longitude');
    var lat = feature.get('latitude');
    
    return ee.Feature(null, {
      'id_nodo': feature.id(), 
      'area_ha': areaHa,
      'x_lon': lon,
      'y_lat': lat,
      'year': year
    });
  });
  
  return nodosConAtributos;
};

// Consolidar la tabla maestra uniendo los 40 años de datos con .flatten()
var tablaNodosMaestra = coleccionMSPA_Completa.map(extraerNodosAnuales).flatten();

// Exportar a Drive
Export.table.toDrive({
  collection: tablaNodosMaestra,
  description: 'Nodos_IPC_AreaDeEstudio_1985_2024',
  folder: 'Objetivo3', 
  fileFormat: 'CSV'
});

// Visualización en el mapa

// Función específica para extraer la geometría vectorial de los nodos para visualización al vuelo
var obtenerNodosParaVis = function(img) {
  var soloCore = img.select('mspa_class').eq(1);
  var coreMasked = soloCore.updateMask(soloCore);
  
  var parchesVector = coreMasked.reduceToVectors({
    geometry: roi.geometry(),
    crs: img.projection(),
    scale: 30,
    geometryType: 'polygon',
    eightConnected: true,
    labelProperty: 'mspa_class',
    maxPixels: 1e13
  });
  return parchesVector;
};

var visTopologia = {
  min: 0, max: 3, 
  palette: ['#e0e0e0', '#00441b', '#a1d99b', '#ff7f00'] // Matriz, Core, Edge, Conectores
};   

// Centrar el mapa y agregar el límite
Map.centerObject(roi, 11).setOptions("SATELLITE");

// Cargar las imágenes iniciales (1985)
var imgMapBiomasInicial = ee.Image(coleccionOriginales.filter(ee.Filter.eq('year', 1985)).first());
var imgMSPAInicial = ee.Image(coleccionMSPA_Completa.filter(ee.Filter.eq('year', 1985)).first());
var nodosIniciales = obtenerNodosParaVis(imgMSPAInicial);

// Definir estilo visual de los nodos (borde rojo semitransparente)
var estiloNodos = {
  color: 'FF0000', 
  fillColor: 'FF000044', 
  width: 2
};

// Configurar el orden inicial de las capas (Paso del slider a la posición 0, 1, 2, 3)
Map.layers().set(0, ui.Map.Layer(imgMapBiomasInicial, visClass, 'MapBiomas Área de Estudio 1985', false));
Map.layers().set(1, ui.Map.Layer(imgMSPAInicial, visTopologia, 'Topología MSPA Área de Estudio 1985', true));
Map.layers().set(2, ui.Map.Layer(nodosIniciales.style(estiloNodos), {}, 'Nodos Vectoriales Área de Estudio 1985', true));
Map.layers().set(3, ui.Map.Layer(roi.style({color: 'black', fillColor: '00000000', width: 2}), {}, 'Límite Área de Estudio', true));

// Crear el Slider de tiempo unificado
var slider = ui.Slider({
  min: 1985,
  max: 2024,
  value: 1985,
  step: 1,
  style: {width: '400px', position: 'bottom-center'},
  onChange: function(year) {
    
    // Filtrar las colecciones por el año seleccionado
    var imgMB_Actualizada = ee.Image(coleccionOriginales.filter(ee.Filter.eq('year', year)).first());
    var imgMSPA_Actualizada = ee.Image(coleccionMSPA_Completa.filter(ee.Filter.eq('year', year)).first());
    var nodosVis_Actualizado = obtenerNodosParaVis(imgMSPA_Actualizada);
    
    // Obtención del estado de visibilidad actual de las capas para no alterarlo al cambiar de año
    var showMB = Map.layers().get(0).getShown();
    var showMSPA = Map.layers().get(1).getShown();
    var showNodos = Map.layers().get(2).getShown();
    
    // Reemplazar las capas dinámicamente
    Map.layers().set(0, ui.Map.Layer(imgMB_Actualizada, visClass, 'MapBiomas Área de Estudio ' + year, showMB));
    Map.layers().set(1, ui.Map.Layer(imgMSPA_Actualizada, visTopologia, 'Topología MSPA Área de Estudio ' + year, showMSPA));
    Map.layers().set(2, ui.Map.Layer(nodosVis_Actualizado.style(estiloNodos), {}, 'Nodos Vectoriales Área de Estudio ' + year, showNodos));
  }
});

// Agregar el Slider
Map.add(slider);

// -----------------------------------
// 5- EXPORTACIÓN DE RASTERS PARA QGIS
// -----------------------------------

// Aislar las imágenes de la cobertura de MapBiomas
var imgMB_1985 = ee.Image(coleccionOriginales.filter(ee.Filter.eq('year', 1985)).first());
var imgMB_2024 = ee.Image(coleccionOriginales.filter(ee.Filter.eq('year', 2024)).first());

// Aislar las imágenes topológicas del MSPA
var imgMSPA_1985 = ee.Image(coleccionMSPA_Completa.filter(ee.Filter.eq('year', 1985)).first());
var imgMSPA_2024 = ee.Image(coleccionMSPA_Completa.filter(ee.Filter.eq('year', 2024)).first());

// Exportat a Drive en .TIF

// MapBiomas 1985
Export.image.toDrive({
  image: imgMB_1985.toByte(), 
  description: 'MapBiomas_AreaDeEstudio_1985',
  folder: 'Objetivo3',
  scale: 30,
  region: roi.geometry(),
  maxPixels: 1e13
});

// MapBiomas 2024
Export.image.toDrive({
  image: imgMB_2024.toByte(),
  description: 'MapBiomas_AreaDeEstudio_2024',
  folder: 'Objetivo3',
  scale: 30,
  region: roi.geometry(),
  maxPixels: 1e13
});

// MSPA 1985
Export.image.toDrive({
  image: imgMSPA_1985.toByte(),
  description: 'MSPA_AreaDeEstudio_1985',
  folder: 'Objetivo3',
  scale: 30,
  region: roi.geometry(),
  maxPixels: 1e13
});

// MSPA 2024
Export.image.toDrive({
  image: imgMSPA_2024.toByte(),
  description: 'MSPA_AreaDeEstudio_2024',
  folder: 'Objetivo3',
  scale: 30,
  region: roi.geometry(),
  maxPixels: 1e13
});

// ---------------------------------------------------
// PASO 6: CÁLCULO DE IPC, JOIN ESPACIAL Y EXPORTACIÓN
// ---------------------------------------------------

// Cargar la tabla generada en el script de R 
var tablaImportancia = ee.FeatureCollection("users/joaquinasad/Nodos_Importancia_2024_v2");

// Converción de la tabla CSV en Puntos espaciales reales usando x_lon e y_lat
var puntosIPC = tablaImportancia.map(function(feature) {
  var lon = ee.Number.parse(feature.get('x_lon'));
  var lat = ee.Number.parse(feature.get('y_lat'));
  return ee.Feature(ee.Geometry.Point([lon, lat]), feature.toDictionary());
});

// Extracción y vectorización de los parches Core del año 2024 al vuelo
var imgMSPA_2024 = ee.Image(coleccionMSPA_Completa.filter(ee.Filter.eq('year', 2024)).first());
var soloCore2024 = imgMSPA_2024.select('mspa_class').eq(1);
var coreMasked2024 = soloCore2024.updateMask(soloCore2024);

var nodosVector2024 = coreMasked2024.reduceToVectors({
  geometry: roi.geometry(),
  crs: imgMSPA_2024.projection(),
  scale: 30,
  geometryType: 'polygon',
  eightConnected: true,
  labelProperty: 'mspa_class',
  maxPixels: 1e13
});

// JOIN ESPACIAL (Intersección geométrica)
var filtroEspacial = ee.Filter.intersects({
  leftField: '.geo',
  rightField: '.geo'
});

var join = ee.Join.saveFirst('valores_ipc');
var nodosConIPC = join.apply(nodosVector2024, puntosIPC, filtroEspacial);

// Formatear la FeatureCollection para que 'Importancia_Relativa' sea una propiedad directa
var nodosFinales = nodosConIPC.map(function(feature) {
  var datosIPC = ee.Feature(feature.get('valores_ipc'));
  var importancia = ee.Number.parse(datosIPC.get('Importancia_Relativa'));
  return feature.set('Importancia_Relativa', importancia);
});

// Gradiente de colores para visualización
var visPriorizacion = {
  palette: ['#ffffcc', '#fd8d3c', '#e31a1c', '#800026'],
  min: 0,
  max: 1
};
var imgPriorizacion = ee.Image().float().paint({
  featureCollection: nodosFinales,
  color: 'Importancia_Relativa'
});

// Visualización
Map.centerObject(roi, 11);
Map.addLayer(imgPriorizacion, visPriorizacion, 'Priorización de Nodos (IPC) - 2024', true);

// Leyenda
var leyenda = ui.Panel({
  style: {position: 'bottom-right', padding: '8px 15px'}
});

var tituloLeyenda = ui.Label({
  value: 'Importancia Relativa (IPC)',
  style: {fontWeight: 'bold', fontSize: '14px', margin: '0 0 4px 0'}
});
leyenda.add(tituloLeyenda);

var palette = visPriorizacion.palette;
var nombres = ['Baja (0 - 0.25)', 'Media-Baja (0.25 - 0.50)', 'Media-Alta (0.50 - 0.75)', 'Alta (0.75 - 1)'];

for (var i = 0; i < palette.length; i++) {
  var colorBox = ui.Label({
    style: {backgroundColor: palette[i], padding: '8px', margin: '0 4px 4px 0'}
  });
  var descLabel = ui.Label({
    value: nombres[i],
    style: {margin: '0 0 4px 6px'}
  });
  var panelItem = ui.Panel({
    widgets: [colorBox, descLabel],
    layout: ui.Panel.Layout.Flow('horizontal')
  });
  leyenda.add(panelItem);
}
Map.add(leyenda);

// Exportar la capa del Índice de Probabilidad de Conectividad a Drive
Export.image.toDrive({
  image: imgPriorizacion.toFloat(),
  description: 'Priorizacion_IPC_2024_v2', 
  folder: 'Objetivo3',
  scale: 30,
  region: roi.geometry(),
  maxPixels: 1e13
});
